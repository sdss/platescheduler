# encoding: utf-8

import numpy as np
import scipy.optimize as optimize
import fitsio

from roboscheduler.scheduler import Observer
from .cadence import assignCadence


class Scheduler(object):
    """Top level scheduler for SDSS V plate program

    """

    def __init__(self, platePath=None, airmass_limit=2., dark_limit=0.35):
        self.platePath = platePath
        self._plates = None
        self._plateIDtoField = None
        self._cadences = {}
        self._obs_query = []  # stash raw result probably
        self.obs_hist = {}
        self._carts = None

        self.bright_cadences = ["YSO", "RV6", "RV12"]

        # Mike did a great job setting up bright/dark/observability rules
        # we'll just borrow that
        self.Observer = Observer(observatory="apo", bright_twilight=-12)
        self.airmass_limit = airmass_limit
        self.dark_limit = dark_limit

        self.moon_threshold = 15.
        self._last_full_moon = None


    def getObs(self):
        """place holder.

        at some point we need to query obs, I guess with plates,
        similar to how AS currently does
        """
        return self._obs_query

    def pullPlates(self):
        """Read in plates for scheduling

        For now this is a fits file, needs to be DB query soon
        """

        self._plates = fitsio.read(self.platePath)
        self._plateIDtoField = dict()

        for p in self._plates:
            if not p["CADENCE"] in self._cadences:
                self._cadences[p["CADENCE"]] = assignCadence(p["CADENCE"])
            self._plateIDtoField[p["PLATE_ID"]] = p["FIELD"]
            self.obs_hist[p["FIELD"]] = []


    @property
    def plates(self):
        if self._plates is None:
            self.pullPlates()
        return self._plates


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
        if self._carts is None:
            self.pullCarts()
        return self._carts


    def _inverseMoon(self, mjd):
        return -1 * self.Observer.moon_illumination(mjd)


    def last_full_moon(self, mjd):
        d_moon_dt = optimize.approx_fprime(mjd, self.Observer.moon_illumination, 0.5)

        if d_moon_dt > 0:
            # waxing, go way back
            res = optimize.minimize(self._inverseMoon, mjd-15, bounds=[[mjd-40, mjd+0.1]])
        else:
            # waning
            res = optimize.minimize(self._inverseMoon, mjd, bounds=[[mjd-27, mjd+1]])
        self._last_full_moon = float(res.x)
        return self._last_full_moon


    def observable(self, plates, mjd=None, check_cadence=True, duration=60.):
        """Return array of plates observable

        Parameters:
        ----------

        mjd : np.float64
            current MJD

        check_cadence: Boolean
            is it safe to ignore cadence requirements? Usually not,
            but if we already failed a check maybe we should.

        duration: float, int
            duration in minutes; will be converted to decimal days
        """

        window_lst_start = self.Observer.lst(mjd) / 15.
        window_lst_end = self.Observer.lst(mjd + duration / 60 / 24) / 15.

        in_window = [True for p in plates]

        # debug = len(plates) < 4 and int(mjd) in [59146, 59147]
        debug = False

        if debug:
            print(window_lst_start, window_lst_end)

        for i, p in enumerate(plates):
            if debug:
                print(p["RA"], p["HA_MIN"], p["HA_MAX"])
            start = (p["RA"] + p["HA_MIN"]) / 15.
            if start < 0:
                start += 24
            elif start > 24:
                start -= 24
            end = (p["RA"] + p["HA_MAX"]) / 15.
            if end < 0:
                end += 24
            elif end > 24:
                end -= 24
            start_diff = float(lstDiff(start, window_lst_start))
            end_diff = float(lstDiff(end, window_lst_end))

            if debug:
                print(f"{start:.2f} {start_diff:.2f} {end:.2f} {end_diff:.2f}")

            if start_diff < 0:
                # start is less than window start
                in_window[i] = False
            elif end_diff > 0:
                in_window[i] = False

        moonra, moondec = self.Observer.moon_radec(mjd)

        moon_dist = np.power((plates["RA"] - moonra)*np.cos(plates["DEC"]*np.pi), 2)\
                  + np.power((plates["DEC"] - moondec), 2)

        # (alt, az) = self.Observer.radec2altaz(mjd=mjd, ra=plates["RA"],
        #                              dec=plates["DEC"])
        # airmass = self.Observer.alt2airmass(alt)
        # check skybrightness half way through so we aren't recording twilight
        # otherwise we'll lose dark programs
        skybrightness = self.Observer.skybrightness(mjd + duration / 60 / 24 / 2)

        observable = (moon_dist > self.moon_threshold) & in_window\
                     & (skybrightness <= plates["SKYBRIGHTNESS"])

        if debug:
            print("B  ", plates["PLATE_ID"], observable)

        if(check_cadence):
            for i, p in enumerate(plates):
                if observable[i] and p["PRIORITY"] < 10:
                    field = p["FIELD"]
                    cad = p["CADENCE"]
                    observable[i] = self.cadences[cad](mjd, hist=self.obs_hist[field], last_full_moon=self._last_full_moon)

        if debug:
            print("A  ", plates["PLATE_ID"], observable)

        return plates[observable]


    def prioritize(self, plates):
        """prioritize observable plates
        """
        base = plates["PRIORITY"] * 100.0

        # super boost pri 10 plates
        pri_ten = [1e4 if p > 950 else 0 for p in base]
        base = base + np.array(pri_ten)

        # gaussian de-weight for easier plates
        dec = base - 50.0 * np.power(np.exp(-1 * (plates["DEC"] - 33)), 2) / (2 * 20**2)

        # TO-DO: possibly add weight by some delta-times, passed from cadence check?

        return dec


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
        night_sched = {"start": night_start,
                       "end": night_end}

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
            night_sched["bright_start"] = night_start
            night_sched["bright_end"] = night_end
            night_sched["dark_start"] = 0
            night_sched["dark_end"] = 0
        elif dark_start and dark_end:
            # we WILL have an RM plate, assume it gets no overhead
            # even if it isn't first or last it's easier this way
            night_sched["bright_start"] = 0
            night_sched["bright_end"] = 0
            night_sched["dark_start"] = night_start
            night_sched["dark_end"] = night_end
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
                night_sched["dark_end"] = night_end - 30 / 60 / 24
                night_sched["bright_start"] = night_end - 30 / 60 / 24
                night_sched["bright_end"] = night_end
            else:
                bright_slots = 0
        elif dark_start and bright_end:
            split_night = True
            split = optimize.brenth(self._bright_dark_function,
                              night_start + 1 / 24, night_end - 1 / 24,
                              args=self.dark_limit)
            bright_time = night_end - split
            dark_time = split - night_start
            night_sched["bright_start"] = split
            night_sched["bright_end"] = night_end
            night_sched["dark_start"] = night_start
            night_sched["dark_end"] = split

        elif bright_start and dark_end:
            split_night = True
            split = optimize.brenth(self._bright_dark_function,
                              night_start + 1 / 24, night_end - 1 / 24,
                              args=self.dark_limit)
            bright_time = split - night_start 
            dark_time = night_end - split
            night_sched["bright_start"] = night_start
            night_sched["bright_end"] = split
            night_sched["dark_start"] = split
            night_sched["dark_end"] = night_end
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
                print("THIS SHOULDN'T HAPPEN")
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

            bright_total = (night_sched["bright_end"] - night_sched["bright_start"]) * 24 * 60
            bright_waste = bright_total - 50 * bright_slots_short - 87 * bright_slots_long
            # print(f"waste:        {mjd} {bright_slots_short}, {bright_slots_long}, {bright_waste:.1f}")
            while bright_waste > 37:
                # print(f"re-allocating {mjd} {bright_slots_short}, {bright_slots_long}, {bright_waste:.1f}")
                bright_slots_short -= 1
                bright_slots_long += 1
                bright_waste = bright_total - 50 * bright_slots_short - 87 * bright_slots_long
                # print(f"re-allocating {mjd} {bright_slots_short}, {bright_slots_long}, {bright_waste:.1f}")

            assert bright_slots_short + bright_slots_long + dark_slots + rm_slots\
                    <= len(self.carts), "cart assignment made up extra carts!"
        else:
            bright_slots_short = bright_slots
            bright_slots_long = 0

        gg_len = [30 + 20 for i in range(bright_slots_short)]
        long_bright =[67 + 20 for i in range(bright_slots_long)]
        dark_lengths = [67 + 20 for i in range(dark_slots)]
        rm_lengths = [120 for i in range(rm_slots)]

        all_lengths = np.sum(gg_len) + np.sum(long_bright) + np.sum(dark_lengths) + np.sum(rm_lengths)

        waste = nightLength - all_lengths / 60 /24

        return night_sched, gg_len, long_bright, dark_lengths, rm_lengths, waste


    def scheduleMjd(self, mjd):
        """Run the scheduling
        """
        night_sched, gg_len, long_bright, dark_lengths, rm_lengths, waste = \
                                                          self.makeSlots(mjd)

        bright_starts = list()
        dark_starts = list()
        tonight_ids = list()
        if night_sched["bright_start"] > 0:
            w_gg = np.where([c == "GG" for c in self.plates["CADENCE"]])
            w_bright_long = np.where([c in self.bright_cadences for c in \
                                      self.plates["CADENCE"]])
            gg_sched = 0

            now = night_sched["bright_start"]

            while len(bright_starts) < (len(gg_len) + len(long_bright)):
                if now + 15 / 60 / 24 >= night_sched["bright_end"]:
                    break

                if gg_sched < len(gg_len):
                    gg_obs = self.observable(self.plates[w_gg], mjd=now,
                                             check_cadence=True, duration=30.)
                else:
                    gg_obs = []
                gg_this = False
                if len(gg_obs) > 0:
                    slot_priorities = self.prioritize(gg_obs)

                    sorted_priorities = np.argsort(slot_priorities)[::-1]
                    i = 0
                    if len(sorted_priorities) > 0:
                        while gg_obs[sorted_priorities[i]]["PLATE_ID"] in tonight_ids:
                            i += 1
                            if i >= len(sorted_priorities) - 1:
                                # g-e for case of len(sorted_priorities) == 1
                                i = -1
                                break
                        if i != -1:
                            tonight_ids.append(gg_obs[sorted_priorities[i]]["PLATE_ID"])

                            bright_starts.append({"start": now,
                                                  "length": 50,
                                                  "plateid": gg_obs[sorted_priorities[i]]["PLATE_ID"]
                                })
                            now += 50 / 60 / 24
                            gg_sched += 1
                            gg_this = True
                if not gg_this:
                    # or we can always do a long plate if needed
                    long_obs = self.observable(self.plates[w_bright_long], mjd=now,
                                         check_cadence=True, duration=30.)
                    slot_priorities = self.prioritize(long_obs)
                    sorted_priorities = np.argsort(slot_priorities)[::-1]
                    i = 0
                    if len(sorted_priorities) > 0:
                        while long_obs[sorted_priorities[i]]["PLATE_ID"] in tonight_ids:
                            i += 1
                            if i >= len(sorted_priorities) - 1:
                                # g-e for case of len(sorted_priorities) == 1
                                i = -1
                                break
                        if i == -1:
                            bright_starts.append({"start": now,
                                                  "length": 50,
                                                  "plateid": None
                            })
                            now += 50 / 60 / 24
                        else:
                            tonight_ids.append(long_obs[sorted_priorities[i]]["PLATE_ID"])

                            bright_starts.append({"start": now,
                                                  "length": 87,
                                                  "plateid": long_obs[sorted_priorities[i]]["PLATE_ID"]
                                })
                            now += 87 / 60 / 24
                    else:
                        # everything else failed
                        bright_starts.append({"start": now,
                                              "length": 50,
                                              "plateid": None
                        })
                        now += 50 / 60 / 24

                # strip overhead appropriately
                # except petunia doesn't care about overhead... leave it alone for now
                # if night_sched["bright_start"] == night_sched["start"]:
                #     bright_starts[0]["length"] -= 20
                # if night_sched["bright_end"] == night_sched["end"]:
                #     bright_starts[-1]["length"] -= 20

        if night_sched["dark_start"] > 0:
            self.last_full_moon(mjd)
            w_rm = np.where([c == "RM" for c in self.plates["CADENCE"]])
            w_aqmes = np.where(["AQMES" in c for c in self.plates["CADENCE"]])

            now = night_sched["dark_start"]

            rm_sched = 0

            while len(dark_starts) < (len(rm_lengths) + len(dark_lengths)):
                if now + 15 / 60 / 24 >= night_sched["dark_end"]:
                    break

                if rm_sched < len(rm_lengths):
                    rm_obs = self.observable(self.plates[w_rm], mjd=now,
                                             check_cadence=True, duration=60.)
                else:
                    rm_obs = []
                rm_this = False
                if len(rm_obs) > 0:
                    slot_priorities = self.prioritize(rm_obs)

                    sorted_priorities = np.argsort(slot_priorities)[::-1]
                    i = 0
                    if len(sorted_priorities) > 0:
                        while rm_obs[sorted_priorities[i]]["PLATE_ID"] in tonight_ids:
                            i += 1
                            if i >= len(sorted_priorities) - 1:
                                # g-e for case of len(sorted_priorities) == 1
                                i = -1
                                break
                        if i != -1:
                            tonight_ids.append(rm_obs[sorted_priorities[i]]["PLATE_ID"])

                            dark_starts.append({"start": now,
                                                  "length": 120,
                                                  "plateid": rm_obs[sorted_priorities[i]]["PLATE_ID"]
                                })
                            now += 120 / 60 / 24
                            rm_sched += 1
                            rm_this = True
                if not rm_this:
                    aqmes_obs = self.observable(self.plates[w_aqmes], mjd=now,
                                         check_cadence=True, duration=30.)
                    slot_priorities = self.prioritize(aqmes_obs)
                    sorted_priorities = np.argsort(slot_priorities)[::-1]
                    i = 0
                    if len(sorted_priorities) > 0:
                        while aqmes_obs[sorted_priorities[i]]["PLATE_ID"] in tonight_ids:
                            i += 1
                            if i >= len(sorted_priorities) - 1:
                                # g-e for case of len(sorted_priorities) == 1
                                i = -1
                                break
                        if i == -1:
                            dark_starts.append({"start": now,
                                                  "length": 87,
                                                  "plateid": None
                            })
                            now += 87 / 60 / 24
                        else:
                            tonight_ids.append(aqmes_obs[sorted_priorities[i]]["PLATE_ID"])

                            dark_starts.append({"start": now,
                                                  "length": 87,
                                                  "plateid": aqmes_obs[sorted_priorities[i]]["PLATE_ID"]
                                })
                            now += 87 / 60 / 24
                    else:
                        # everything else failed
                        dark_starts.append({"start": now,
                                              "length": 87,
                                              "plateid": None
                        })
                        now += 87 / 60 / 24

        return bright_starts, dark_starts, night_sched


def isWithinSlot(starts, mjd, delta):

    w_within = np.where(np.abs(np.array(starts) - mjd) < delta / 60 / 24)

    if len(w_within[0]) == 0:
        return None

    return w_within[0][0]


def lstDiff(a, b):
    """Intelligently find signed difference in 2 lsts, b-a

    Assuming a clockwise coordinate system that wraps at 24=0

    A value ahead clockwise is "larger" for subtraction purposes,
    even if it is on the other side of zero.

    Parameters:
    -----------

    a : np.float32, float
        first lst in hours
    b : np.float32, float
        second lst in hours

    Returns:
    --------

    diff: float
        b-a, wrapped correctly

    Comments:
    ---------

    """

    if a < b:
        # if a is closer going through zero, wrap is used
        # meaning a is clockwise of a, technically ahead of it
        # so its direction is negative
        simple = b - a
        wrap = a + 24 - b

        if wrap < simple:
            return -wrap
        else:
            return simple
    else:
        # a is greater, but direction is tricky
        simple = b - a  # always negative now
        wrap = b + 24 - a  # always positive

        # if wrap < abs(simple), it is closer to go through zero again
        # in this case though that means b is clockwise a
        # so wrap is appropriately positive
        # otherwise simple is appropriately negative

        if wrap < abs(simple):
            return wrap
        else:
            return simple
