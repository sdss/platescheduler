# encoding: utf-8

import numpy as np
import scipy.optimize as optimize

from roboscheduler.scheduler import Observer
from .cadence import assignCadence


class Scheduler(object):
    """Top level scheduler for SDSS V plate program

    """

    def __init__(self, platePath=None, airmass_limit=2., dark_limit=0.35):
        self.platePath = platePath
        self._plates = None  # may not need this
        # self._ras = None
        # self._decs = None
        # self._priorities = None
        # self._lunations = None
        self._cadences = None
        # self._programs = None
        self._obs_query = None  # stash raw result probably
        self._obs_hist = None
        self._carts = None

        self.bright_cadences = ["YSO", "RV6", "RV12"]

        # Mike did a great job setting up bright/dark/observability rules
        # we'll just borrow that
        self.Observer = Observer(observatory="apo", bright_twilight=-12)
        self.airmass_limit = airmass_limit
        self.dark_limit = dark_limit


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

        for p in self._plates:
            if not p["CADENCE"] in self._cadences:
                self._cadences[p["CADENCE"]] = assignCadence(p["CADENCE"])


    @property    
    def obs_hist(self):
        # this is temporary
        # need to hook into DB and establish AB pairs, etc probably?
        if self._obs_hist is None:
            self._obs_hist = dict()
            hist = self.getObs()
            for h in hist:
                if h["FIELD"] in self._obs_hist:
                    self._obs_hist[h["FIELD"]].append(h["MJD"])
                else:
                    self._obs_hist[h["FIELD"]] = [h["MJD"]]


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


    def observable(self, plates, mjd=None, check_cadence=True):
        """Return array of plates observable

        Parameters:
        ----------

        mjd : np.float64
            current MJD

        check_cadence: Boolean
            is it safe to ignore cadence requirements? Usually not,
            but if we already failed a check maybe we should.
        """

        (alt, az) = self.Observer.radec2altaz(mjd=mjd, ra=plates["RA"],
                                     dec=plates["DEC"])
        airmass = self.Observer.alt2airmass(alt)
        skybrightness = self.Observer.skybrightness(mjd)

        observable = (alt > 0.) & (airmass < self.airmass_limit)\
                     & (skybrightness < plates["SKYBRIGHTNESS"])

        if(check_cadence):
            for p, o in zip(plates, observable):
                if o and p["PRIORITY"] < 10:
                    field = p["FIELD"]
                    cad = p["CADENCE"]
                    o = self.cadences[cad](self.obs_hist[field])


        return plates[observable]

    def prioritize(self, plates):
        """prioritize observable plates
        """
        base = plates["PRIORITY"] * 100.0

        # super boost pri 10 plates
        pri_ten = [1e4 if p > 950 else 0 for p in base]
        base = base + np.array(pri_ten)

        # gaussian de-weight for easier plates
        dec = base - 50.0 * np.power(np.exp(-(plates["DEC"] - 33), 2) / (2 * (20)**2))

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
        if night_sched["bright_start"] > 0:
            if len(long_bright) > 0:
                # find highest priority long slots and assign them
                w_bright_long = np.where([c in self.bright_cadences for c in \
                                 self.plates["CADENCE"]])
                
                long_slots = np.arange(night_sched["bright_start"], night_sched["bright_end"], 50 / 60 / 24)
                priorities = list()
                plate_ids = list()
                for l in long_slots:
                    b_long_obs = observable(self.plates[w_bright_long], mjd=l, check_cadence=True)
                    slot_priorities = self.prioritize(b_long_obs)
                    highest = np.argmax(slot_priorities)
                    priorities.append(slot_priorities[highest])
                    plate_ids.append(plates[highest]["PLATE_ID"])

                # sort fields into descending priority order
                # track index of each, so we can use plate_ids and long_slots
                priority_order = np.argsort(priorities)[::-1]

                long_starts = list()
                long_plates = list()
                i = 0  # for iterating in while loop
                while len(long_starts) < len(long_bright):
                    idx = priority_order[i]
                    long_starts.append(long_slots[idx])
                    long_plates.append(plate_ids[idx])
                    i += 1
                
                begin_night = isWithinSlot(long_starts, night_sched["bright_start"], 20 / 60 / 24)

                now = night_sched["bright_start"]

                if begin_night is not None:
                    bright_starts.append({"start": night_sched["bright_start"],
                                          "length": 87,
                                          "plateid": long_plates[begin_night]
                        })
                    now += 87 / 60 / 24

                while len(bright_starts) < (gg_len + long_bright):
                    long_check = isWithinSlot(long_starts, now, 50)

                    if long_check is not None:
                        bright_starts.append({"start": now,
                                              "length": 87,
                                              "plateid": long_plates[long_check]
                            })
                        now += 87 / 60 / 24
                    else:
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

            else:
                # all short, that's easier
                now = night_sched["bright_start"]
                while len(bright_starts) < (gg_len):
                    bright_starts.append({"start": now,
                                          "length": 50,
                                          "plateid": None
                        })
                    now = night_sched["start"] + 50 / 60 / 24

                # strip overhead appropriately
                # except petunia doesn't care about overhead... leave it alone for now
                # if night_sched["bright_start"] == night_sched["start"]:
                #     bright_starts[0]["length"] -= 20
                # if night_sched["bright_end"] == night_sched["end"]:
                #     bright_starts[-1]["length"] -= 20

            w_gg = np.where([c == "GG" for c in self.plates["CADENCE"]])

            for b in bright_starts:
                # now fill in GG plates as needed
                if b["plateid"] is None:
                    gg_obs = observable(self.plates[w_gg], mjd=b["start"], check_cadence=True)
                    slot_priorities = self.prioritize(gg_obs)
                    highest = np.argmax(slot_priorities)
                    b["plateid"] = gg_obs[highest]["PLATE_ID"]

        if night_sched["dark_start"] > 0:
            if len(rm_lengths) > 0:
                w_rm = np.where([c == "RM" for c in self.plates["CADENCE"]])
                
                rm_slots = np.arange(night_sched["dark_start"], night_sched["dark_end"], 80 / 60 / 24)
                priorities = list()
                plate_ids = list()
                for l in rm_slots:
                    b_long_obs = observable(self.plates[w_rm], mjd=l, check_cadence=True)
                    slot_priorities = self.prioritize(b_long_obs)
                    highest = np.argmax(slot_priorities)
                    priorities.append(slot_priorities[highest])
                    plate_ids.append(plates[highest]["PLATE_ID"])

                # sort fields into descending priority order
                # track index of each, so we can use plate_ids and long_slots
                priority_order = np.argsort(priorities)[::-1]

                rm_starts = list()
                rm_plates = list()
                i = 0  # for iterating in while loop
                while len(rm_starts) < len(rm_lengths):
                    idx = priority_order[i]
                    rm_starts.append(rm_slots[idx])
                    rm_plates.append(plate_ids[idx])
                    i += 1
                
                begin_night = isWithinSlot(rm_starts, night_sched["dark_start"], 60 / 60 / 24)

                now = night_sched["dark_start"]

                if begin_night is not None:
                    dark_starts.append({"start": night_sched["dark_start"],
                                          "length": 120,
                                          "plateid": rm_plates[begin_night]
                        })
                    now += 120 / 60 / 24

                while len(dark_starts) < (dark_lengths + rm_lengths):
                    long_check = isWithinSlot(rm_starts, now, 50)

                    if long_check is not None:
                        dark_starts.append({"start": now,
                                              "length": 120,
                                              "plateid": rm_plates[long_check]
                            })
                        now += 120 / 60 / 24
                    else:
                        dark_starts.append({"start": now,
                                              "length": 80,
                                              "plateid": None
                            })
                        now += 80 / 60 / 24
                
                # strip overhead appropriately
                # except petunia doesn't care about overhead... leave it alone for now

                # if night_sched["dark_start"] == night_sched["start"]:
                #     bright_starts[0]["length"] -= 20
                # if night_sched["dark_end"] == night_sched["end"]:
                #     bright_starts[-1]["length"] -= 20

            else:
                now = night_sched["dark_start"]
                while len(dark_starts) < (dark_lengths):
                    bright_starts.append({"start": now,
                                          "length": 80,
                                          "plateid": None
                        })
                    now = night_sched["start"] + 80 / 60 / 24

                # strip overhead appropriately
                # except petunia doesn't care about overhead... leave it alone for now

                # if night_sched["dark_start"] == night_sched["start"]:
                #     bright_starts[0]["length"] -= 20
                # if night_sched["dark_end"] == night_sched["end"]:
                #     bright_starts[-1]["length"] -= 20

            w_aqmes = np.where(["AQMES" in c for c in self.plates["CADENCE"]])

            for b in dark_starts:
                if b["plateid"] is None:
                    aqmes_plates = observable(self.plates[w_aqmes], mjd=b["start"], check_cadence=True)
                    slot_priorities = self.prioritize(aqmes_plates)
                    highest = np.argmax(slot_priorities)
                    b["plateid"] = aqmes_plates[highest]["PLATE_ID"]

    return bright_starts, dark_starts

    # don't forget to match plates at different hour angles


def isWithinSlot(starts, mjd, delta):

    w_within = np.where(np.abs(np.array(starts) - mjd) < delta)

    if len(w_within[0]) == 0:
        return None

    return w_within[0][0]
