# encoding: utf-8

import numpy as np
import scipy.optimize as optimize

from roboscheduler.scheduler import Observer
from .cadence import Cadence


class Scheduler(object):
    """Top level scheduler for SDSS V plate program

    """

    def __init__(self, platePath=None, airmass_limit=2., dark_limit=0.35):
        self.platePath = platePath
        self._plates = None  # may not need this
        self._ras = None
        self._decs = None
        self._priorities = None
        self._lunations = None
        self._cadences = None
        self._programs = None
        self._lastObs = None
        self._carts = None

        # Mike did a great job setting up bright/dark/observability rules
        # we'll just borrow that
        self.Observer = Observer(observatory="apo", bright_twilight=-12)
        self.airmass_limit = airmass_limit
        self.dark_limit = dark_limit


    def getObs(plate):
        """place holder.

        at some point we need to query obs, I guess with plates,
        similar to how AS currently does
        """
        return 0

    def pullPlates(self):
        """Read in plates for scheduling

        For now this is a fits file, needs to be DB query soon
        """

        self._plates = fitsio.read(self.platePath)

        self._plateIDs = self._plates["id"]
        self._ras = self._plates["ra"]
        self._decs = self._plates["dec"]
        self._priorities = self._plates["priority"]
        self._skybrightness = self._plates["skybrightness"]
        self._programs = self._plates["program"]
        self._cadences = [Cadence(d) for d in self._plates["delta"]]
        self._lastObs = [self.getObs(i) for i in self._plateIDs]


    @property
    def plates(self):
        if self._plates is None:
            self.pullPlates()
        return self._plates


    @property
    def ras(self):
        if self._ras is None:
            self.pullPlates()
        return self._ras


    @property
    def decs(self):
        if self._decs is None:
            self.pullPlates()
        return self._decs


    @property
    def skybrightness(self):
        if self._skybrightness is None:
            self.pullPlates()
        return self._skybrightness


    @property
    def priorities(self):
        if self._priorities is None:
            self.pullPlates()
        return self._priorities


    @property
    def cadences(self):
        if self._cadences is None:
            self.pullPlates()
        return self._cadences


    def pullCarts(self):
        """Get cart info

        For now hardcode, but this will need to get carts from DB soon
        """

        self._carts = {
            1: "BOTH",
            2: "BOTH",
            3: "BOTH",
            4: "BOTH",
            5: "BOSS",
            6: "BOTH",
            7: "BOTH",
            8: "BOTH",
            9: "BOTH",
        }

    @property
    def carts(self):
        self.pullCarts()
        return self._carts
    
    def observable(self, mjd=None, check_cadence=True):
        """Return array of plates observable

        Parameters:
        ----------

        mjd : np.float64
            current MJD

        check_cadence: Boolean
            is it safe to ignore cadence requirements? Usually not,
            but if we already failed a check maybe we should.
        """

        (alt, az) = self.Observer.radec2altaz(mjd=mjd, ra=self.ras,
                                     dec=self.decs)
        airmass = self.Observer.alt2airmass(alt)
        skybrightness = self.Observer.skybrightness(mjd)

        observable = (alt > 0.) & (airmass < self.airmass_limit)\
                     & (skybrightness < self.skybrightness)

        if(check_cadence):
            for i, o in enumerate(observable):
                if o:
                    self.cadences[i](self.lastObs[i])

        iobservable = np.where(observable)[0]

        return iobservable

    def prioritize(self, iobservable):
        """prioritize observable plates
        """


    def _bright_dark_function(self, mjd=None, switch=None):
        if switch is None:
            switch = self.dark_limit

        return self.Observer.skybrightness(mjd) - switch

    def makeSlots(self, mjd):
        """Determine observable slots based on plates and time
        """
        night_start = self.Observer.evening_twilight(mjd)
        night_end = self.Observer.morning_twilight(mjd)
        nightLength = night_end - night_start

        bright_start = bool(self.Observer.skybrightness(night_start + 1 / 24) >= 0.35)
        bright_end = bool(self.Observer.skybrightness(night_end - 1 / 24) >= 0.35)
        dark_start = bool(self.Observer.skybrightness(night_start + 1 / 24) < 0.35)
        dark_end = bool(self.Observer.skybrightness(night_end - 1 / 24) < 0.35)

        short_slot = (30 + 20) / 60 / 24  # 30 min GG size
        dark_slot = (67 + 20) / 60 / 24

        split_night = False
        if bright_start and bright_end:
            bright_slots = int(nightLength // short_slot)
            rm_slots = 0
            dark_slots = 0
        elif dark_start and dark_end:
            # we WILL have an RM plate, assume it gets no overhead
            # even if it isn't first or last it's easier this way
            remainder = nightLength - 2/24
            rm_slots = 1
            dark_slots = int(remainder // dark_slot)
            extra = remainder % (dark_slots * dark_slot)
            if extra >= 67 / 60 / 24:
                # we can ignore a second overhead on a full darknight
                dark_slots += 1
            elif extra >= 30 / 60 / 24:
                # ignore the overhead and add a GG plate
                bright_slots = 1
            else:
                bright_slots = 0
        elif dark_start and bright_end:
            split_night = True
            split = optimize.brenth(self._bright_dark_function,
                              night_start + 1 / 24, night_end - 1 / 24,
                              args=self.dark_limit)
            bright_time = night_end - split
            dark_time = split - night_start

        elif bright_start and dark_end:
            split_night = True
            split = optimize.brenth(self._bright_dark_function,
                              night_start + 1 / 24, night_end - 1 / 24,
                              args=self.dark_limit)
            bright_time = split - night_start 
            dark_time = night_end - split
        else:
            raise Exception("You broke boolean algebra!")

        if split_night:
            bright_slots = int(bright_time // short_slot)
            if dark_time > 2 / 24:
                # again assume the free overhead is for RM exposure
                rm_slots = 1
                dark_time = dark_time - 2 / 24
            else:
                rm_slots = 0
            dark_slots = int(dark_time // dark_slot)

            extra = dark_time % (dark_slots * dark_slot)
            if extra >= 67 / 60 / 24 and rm_slots == 0:
                # we can give this slot the free overhead
                dark_slots += 1
            elif extra >= 30 / 60 / 24:
                # ignore the overhead and add a GG plate
                bright_slots += 1

        slots = np.sum([bright_slots, dark_slots, rm_slots])
        if slots > len(self.carts):
            used = dark_slots + rm_slots
            avail = len(self.carts) - used
            n = bright_slots - avail
            if n <= 0:
                bright_slots_short = bright_slots
                bright_slots_long = 0
            else:
                bright_slots_short = bright_slots - 2*n
                bright_slots_long = n
            if bright_slots_short < 0:
                # this came up briefly in testing
                bright_slots_short = 0
            if bright_slots_long > avail:
                # this also came up in testing with
                # small number of carts, shouldn't happen in practice but...
                bright_slots_long = avail
            assert bright_slots_short + bright_slots_long + dark_slots + rm_slots\
                    <= len(self.carts), "cart assignment made up extra carts!"
        else:
            bright_slots_short = bright_slots
            bright_slots_long = 0

        bright_lengths = [30 + 20 for i in range(bright_slots_short)]
        bright_lengths.extend([67 + 20 for i in range(bright_slots_long)])
        dark_lengths = [67 + 20 for i in range(dark_slots)]
        rm_lengths = [120 for i in range(rm_slots)]

        all_lengths = np.sum(bright_lengths) + np.sum(dark_lengths) + np.sum(rm_lengths)

        waste = nightLength - all_lengths / 60 /24

        return bright_lengths, dark_lengths, rm_lengths, waste
