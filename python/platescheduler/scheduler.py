# encoding: utf-8

from __future__ import print_function, division, absolute_import
from time import time
import os

import numpy as np
import scipy.optimize as optimize
# import fitsio
from petunia.plateDBtools.database.apo.platedb import ModelClasses as pdb
from petunia.plateDBtools.database.apo.apogeeqldb import ModelClasses as qldb

import sqlalchemy
from sqlalchemy import or_

import yaml

from .roboscheduler import Observer
from .cadence import assignCadence


def get_plates(session):
    '''DESCRIPTION: Reads in APOGEE-II plate information from platedb
    INPUT: None
    OUTPUT: apg -- list of objects with all APOGEE-II plate information'''
    start_time = time()

    from petunia.plateDBtools.database.connections.MyLocalConnection import db

    # session = db.Session()

    try:
        acceptedStatus = session.query(pdb.PlateStatus).filter(pdb.PlateStatus.label == "Accepted").one()
    except sqlalchemy.orm.exc.NoResultFound:
        raise Exception("Could not find 'Accepted' status in plate_status table")

    # ##################################
    # check up on this survey lable for live db
    # ##################################

    try:
        mwm = session.query(pdb.Survey).filter(pdb.Survey.label == "MWM").one()
    except sqlalchemy.orm.exc.NoResultFound:
        raise Exception("Could not find 'mwm' survey in survey table")

    try:
        bhm = session.query(pdb.Survey).filter(pdb.Survey.label == "BHM").one()
    except sqlalchemy.orm.exc.NoResultFound:
        raise Exception("Could not find 'bhm' survey in survey table")

    try:
        plateLoc = session.query(pdb.PlateLocation).filter(pdb.PlateLocation.label == "UW").one()
    except sqlalchemy.orm.exc.NoResultFound:
        raise Exception("Could not find 'APO' location in plate_location table")

    try:
        plateLoc2 = session.query(pdb.PlateLocation).filter(pdb.PlateLocation.label == "APO").one()
    except sqlalchemy.orm.exc.NoResultFound:
        raise Exception("Could not find 'APO' location in plate_location table")

    # Pull all relevant plate information for APOGEE plates
    protoList = list()
    with session.begin():
        protoList = session.query(pdb.Plate.plate_id)\
                .join(pdb.PlateToSurvey, pdb.Survey)\
                .join(pdb.PlateLocation)\
                .join(pdb.PlateToPlateStatus, pdb.PlateStatus)\
                .filter(or_(pdb.Survey.pk == bhm.pk, pdb.Survey.pk == mwm.pk))\
                .filter(or_(pdb.Plate.location == plateLoc, pdb.Plate.location == plateLoc2)).all()
                # .filter(pdb.Plate.location == plateLoc).all()
                # .filter(pdb.PlateStatus.pk == acceptedStatus.pk).all()
        locIDS = session.query(pdb.Plate.location_id)\
               .filter(or_(pdb.Survey.pk == bhm.pk, pdb.Survey.pk == mwm.pk))\
               .filter(pdb.Plate.plate_id.in_(protoList)).all()
        plate_query = session.query(pdb.Plate)\
               .join(pdb.PlateToSurvey, pdb.Survey)\
               .filter(or_(pdb.Survey.pk == bhm.pk, pdb.Survey.pk == mwm.pk))\
               .filter(pdb.Plate.location_id.in_(locIDS)).all()

    q1Time = time()
    print('[SQL]: plate query completed in {} s'.format(q1Time-start_time))

    exposedPlates = [p.plate_id for p in plate_query]
    exposures = session.query(sqlalchemy.func.floor(pdb.Exposure.start_time/86400+.3),\
                 pdb.Plate.plate_id, qldb.Quickred.snr_standard, qldb.Reduction.snr)\
                .join(pdb.ExposureFlavor).join(pdb.Observation).join(pdb.PlatePointing)\
                .join(pdb.Plate).outerjoin(qldb.Quickred).outerjoin(qldb.Reduction)\
                .filter(pdb.ExposureFlavor.label == 'Object')\
                .filter(pdb.Plate.plate_id.in_(exposedPlates)).all()

    plate_exps = {int(p): [] for p in exposedPlates}

    for e in exposures:
        mjd = int(e[0])  # sqlalchemy.func doesn't give an attribute
        if e.snr > 10:
            plate_exps[int(e.plate_id)].append(mjd)
        elif not e.snr and e.snr_standard > 10:
            plate_exps[int(e.plate_id)].append(mjd)

    field_exp_hist = dict()

    PLATE_ID = list()
    FIELD = list()
    RA = list()
    DEC = list()
    HA = list()
    HA_MIN = list()
    HA_MAX = list()
    SKYBRIGHTNESS = list()
    CADENCE = list()
    PRIORITY = list()

    for p in plate_query:
        if "cadence" not in p.design.designDictionary:
            print("skipping: ", p.plate_id)
            continue

        survey_mode = p.currentSurveyMode.definition_label
        field = str(p.name)

        # print(p.plate_id, survey_mode)

        PLATE_ID.append(int(p.plate_id))
        FIELD.append(field)
        if not field in field_exp_hist:
            field_exp_hist[field] = list()
        RA.append(float(p.firstPointing.center_ra))
        DEC.append(float(p.firstPointing.center_dec))
        HA.append(float(p.firstPointing.platePointing(p.plate_id).hour_angle))
        HA_MIN.append(float(p.firstPointing.platePointing(p.plate_id).ha_observable_min) - 7.5)
        HA_MAX.append(float(p.firstPointing.platePointing(p.plate_id).ha_observable_max) + 7.5)
        if survey_mode.lower() == "mwmlead":
            SKYBRIGHTNESS.append(1)
        else:
            SKYBRIGHTNESS.append(0.35)
        CADENCE.append(p.design.designDictionary["cadence"])
        PRIORITY.append(float(p.firstPointing.platePointing(p.plate_id).priority))

        if survey_mode.lower() == "mwmlead":
            plate_mjds = np.array(plate_exps[PLATE_ID[-1]])
            mjds = np.unique(plate_mjds)
            for m in mjds:
                day = np.where(plate_mjds)
                if len(day[0]) > 1:
                    # assume 2 exps count for a AB pair
                    field_exp_hist[field].append(m)
        else:
            for plug in p.pluggings:
                for s in plug.scienceExposures():
                    if "good" in s.status.label.lower():
                        field_exp_hist[field].append(s.mjd)

    # hist = dict()
    # for p, f in zip(PLATE_ID, FIELD):
    #     hist[p] = field_exp_hist[f]

    plate_types = [("PLATE_ID", np.int32),
                   ("FIELD", np.dtype("a20")),
                   ("RA", np.float64),
                   ("DEC", np.float64),
                   ("HA", np.float64),
                   ("HA_MIN", np.float64),
                   ("HA_MAX", np.float64),
                   ("SKYBRIGHTNESS", np.float64),
                   ("CADENCE", np.dtype("a20")),
                   ("PRIORITY", np.int32)]


    plates = np.zeros(len(PLATE_ID), dtype=plate_types)

    plates["PLATE_ID"] = np.array(PLATE_ID)
    plates["FIELD"] = np.array(FIELD)
    plates["RA"] = np.array(RA)
    plates["DEC"] = np.array(DEC)
    plates["HA"] = np.array(HA)
    plates["HA_MIN"] = np.array(HA_MIN)
    plates["HA_MAX"] = np.array(HA_MAX)
    plates["SKYBRIGHTNESS"] = np.array(SKYBRIGHTNESS)
    plates["CADENCE"] = np.array(CADENCE)
    plates["PRIORITY"] = np.array(PRIORITY)

    buildTime = time()
    print('[PY]: array built in {} s'.format(buildTime - q1Time))

    return plates, field_exp_hist


class Scheduler(object):
    """Top level scheduler for SDSS V plate program

    """

    def __init__(self, session=None, airmass_limit=2., dark_limit=0.35):
        if session is None:
            self.session = db.Session()
        else:
            self.session = session
        self._plates = None
        self._plateIDtoField = None
        self._cadences = {}
        self.obs_hist = None
        self._carts = None

        self.bright_cadences = ["YSO", "RV6", "RV12"]

        # Mike did a great job setting up bright/dark/observability rules
        # we'll just borrow that
        self.Observer = Observer(observatory="apo", bright_twilight=-12)
        self.airmass_limit = airmass_limit
        self.dark_limit = dark_limit

        self.moon_threshold = 15.
        self._last_full_moon = None

        # exp  times
        self.aqm_time = 60
        self.rm_time = 120
        self.apogee_time = 67
        self.gg_time = 33
        self.overhead = 20



    def pullPlates(self):
        """Read in plates for scheduling

        For now this is a fits file, needs to be DB query soon
        """

        # self._plates = fitsio.read(self.platePath)
        self._plateIDtoField = dict()
        self._plates, self.obs_hist = get_plates(self.session)

        for p in self._plates:
            if not p["CADENCE"] in self._cadences:
                self._cadences[p["CADENCE"]] = assignCadence(p["CADENCE"])
            self._plateIDtoField[p["PLATE_ID"]] = p["FIELD"]
            # print("{pid:6d} {field:10s} {bright:5.2f}".format(pid=p["PLATE_ID"], field=p["FIELD"], bright=p["SKYBRIGHTNESS"]))


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

        # self._carts = {
        #     1: {"type": "BOTH", "plate": 0},
        #     2: {"type": "BOTH", "plate": 0},
        #     3: {"type": "BOTH", "plate": 0},
        #     4: {"type": "BOTH", "plate": 0},
        #     5: {"type": "BOSS", "plate": 0},
        #     6: {"type": "BOTH", "plate": 0},
        #     7: {"type": "BOTH", "plate": 0},
        #     8: {"type": "BOTH", "plate": 0},
        #     9: {"type": "BOTH", "plate": 0},
        # }

        self._carts = yaml.load(open(os.path.expanduser("~/.cart_status.yml")))

        currentplug = self.session.query(pdb.Cartridge.number, pdb.Plate.plate_id)\
                                .join(pdb.Plugging).join(pdb.Plate).join(pdb.ActivePlugging)\
                                .order_by(pdb.Cartridge.number).all()

        for cart, plate in currentplug:
            assert cart in self._carts, "CART {} UNACOUNTED FOR. Update ~/.cart_status.yml".format(cart)
            self._carts[cart]["plate"] = plate


    @property
    def carts(self):
        if self._carts is None:
            self.pullCarts()
        return self._carts


    def assignCarts(self, bright_starts, dark_starts):
        """
        """

        plugged = {v["plate"]: k for k, v in self.carts.items()}

        available = {k: v["type"] for k, v in self.carts.items()} # definitely want a copy

        for s in bright_starts:
            if s["cart"] is None:
                if s["plate"] in plugged:
                    s["cart"] = plugged[s["plate"]]
                    del available[plugged[s["plate"]]]

        for s in dark_starts:
            if s["cart"] is None:
                if s["plate"] in plugged:
                    s["cart"] = plugged[s["plate"]]
                    del available[plugged[s["plate"]]]


        # THIS MAY NEED TO CHANGE!
        # for testing I set one cart BOSS only
        # if there are APOGEE only carts and no BOSS only, this order reverses
        # if there are both APOGEE only and BOSS only, we need more loops
        # (assign bright apogee only, then dark boss only, then triage "both")

        for s in bright_starts:
            if s["cart"] is None:
                for k, v in available.items():
                    if v in ["BOTH", "APOGEE"]:
                        s["cart"] = k
                        del available[k]
                        break

        for s in dark_starts:
            if s["cart"] is None:
                for k, v in available.items():
                    print(v)
                    if v in ["BOTH", "BOSS"]:
                        print(k, s["plate"])
                        s["cart"] = k
                        del available[k]
                        break

        # print(self.carts)
        # print(bright_starts + dark_starts)

        for s in bright_starts + dark_starts:
            if s["cart"] is None:
                raise Exception("no cart assigned for {}".format(s["plate"]))


    def _inverseMoon(self, mjd):
        return -1 * self.Observer.moon_illumination(mjd)


    def last_full_moon(self, mjd):
        mjd = np.array([mjd])
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

        if len(plates) == 0:
            return []

        window_lst_start = self.Observer.lst(mjd) / 15.
        window_lst_end = self.Observer.lst(mjd + duration / 60 / 24) / 15.

        in_window = [True for p in plates]

        for i, p in enumerate(plates):
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

        if(check_cadence):
            for i, p in enumerate(plates):
                if observable[i] and p["PRIORITY"] < 10:
                    field = p["FIELD"]
                    cad = p["CADENCE"]
                    observable[i] = self.cadences[cad](mjd, hist=self.obs_hist[field], last_full_moon=self._last_full_moon)

        return plates[observable]


    def prioritize(self, plates):
        """prioritize observable plates
        """
        if len(plates) == 0:
            return np.array([])
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

        short_slot = (self.gg_time + self.overhead) / 60 / 24  # 30 min GG size
        dark_slot = (self.aqm_time + self.overhead) / 60 / 24

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
            if extra >= self.apogee_time / 60 / 24:
                # we can ignore a second overhead on a full darknight
                dark_slots += 1
                bright_slots = 0
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
            if dark_time > self.rm_time / 60 / 24:
                # again assume the free overhead is for RM exposure
                rm_slots = 1
                dark_time = dark_time - self.rm_time / 60 / 24
            else:
                rm_slots = 0
            dark_slots = int(dark_time // dark_slot)

            extra = dark_time % (dark_slots * dark_slot)
            if extra >= self.apogee_time / 60 / 24 and rm_slots == 0:
                # we can give this slot the free overhead
                dark_slots += 1
            elif extra >= self.gg_time / 60 / 24:
                # ignore the overhead and add a GG plate
                bright_slots += 1

        slots = np.sum([bright_slots, dark_slots, rm_slots])
        bright_carts = [v["type"] for k, v in self.carts.items() if v["type"] in ["BOTH", "APOGEE"]]
        if slots > len(self.carts):
            used = dark_slots + rm_slots
            avail = min(len(self.carts) - used, len(bright_carts))
            print(mjd, avail, len(self.carts) - used, len(bright_carts))
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
            bright_waste = bright_total - (self.gg_time + self.overhead) * bright_slots_short\
                                        - (self.apogee_time + self.overhead) * bright_slots_long
            # print(f"waste:        {mjd} {bright_slots_short}, {bright_slots_long}, {bright_waste:.1f}")
            while bright_waste > 37:
                # print(f"re-allocating {mjd} {bright_slots_short}, {bright_slots_long}, {bright_waste:.1f}")
                bright_slots_short -= 1
                bright_slots_long += 1
                bright_waste = bright_total - (self.gg_time + self.overhead) * bright_slots_short\
                                            - (self.apogee_time + self.overhead) * bright_slots_long
                # print(f"re-allocating {mjd} {bright_slots_short}, {bright_slots_long}, {bright_waste:.1f}")

            assert bright_slots_short + bright_slots_long + dark_slots + rm_slots\
                    <= len(self.carts), "cart assignment made up extra carts!"
        else:
            bright_slots_short = bright_slots
            bright_slots_long = 0

        gg_len = [self.gg_time + self.overhead for i in range(bright_slots_short)]
        long_bright =[self.apogee_time + self.overhead for i in range(bright_slots_long)]
        dark_lengths = [self.aqm_time + self.overhead for i in range(dark_slots)]
        rm_lengths = [self.rm_time for i in range(rm_slots)]

        all_lengths = np.sum(gg_len) + np.sum(long_bright) + np.sum(dark_lengths) + np.sum(rm_lengths)

        waste = nightLength - all_lengths / 60 /24

        return night_sched, gg_len, long_bright, dark_lengths, rm_lengths, waste


    def scheduleMjd(self, mjd):
        """Run the scheduling
        """

        self.pullCarts()

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

                            bright_starts.append({"obsmjd": now,
                                                  "exposure_length": self.gg_time,
                                                  "plate": gg_obs[sorted_priorities[i]]["PLATE_ID"],
                                                  "cart": None
                                })
                            now += (self.gg_time + self.overhead) / 60 / 24
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
                            bright_starts.append({"obsmjd": now,
                                                  "exposure_length": self.gg_time,
                                                  "plate": -1,
                                                  "cart": -1
                            })
                            now += (self.gg_time + self.overhead) / 60 / 24
                        else:
                            tonight_ids.append(long_obs[sorted_priorities[i]]["PLATE_ID"])

                            bright_starts.append({"obsmjd": now,
                                                  "exposure_length": self.apogee_time,
                                                  "plate": long_obs[sorted_priorities[i]]["PLATE_ID"],
                                                  "cart": None
                                })
                            now += (self.apogee_time + self.overhead) / 60 / 24
                    else:
                        # everything else failed
                        bright_starts.append({"obsmjd": now,
                                              "exposure_length": self.gg_time,
                                              "plate": -1,
                                              "cart": -1
                        })
                        now += (self.gg_time + self.overhead) / 60 / 24

                # strip overhead appropriately
                # except petunia doesn't care about overhead... leave it alone for now
                # if night_sched["bright_start"] == night_sched["obsmjd"]:
                #     bright_starts[0]["exposure_length"] -= self.overhead
                # if night_sched["bright_end"] == night_sched["end"]:
                #     bright_starts[-1]["exposure_length"] -= self.overhead

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
                            this_plate = rm_obs[sorted_priorities[i]]

                            tonight_ids.append(this_plate["PLATE_ID"])

                            alt, az = self.Observer.radec2altaz(mjd=now,
                                                            ra=this_plate["RA"],
                                                            dec=this_plate["DEC"])
                            exp_time = darkExpTime(alt, default=self.rm_time)

                            dark_starts.append({"obsmjd": now,
                                                "exposure_length": exp_time,
                                                "plate": this_plate["PLATE_ID"],
                                                "cart": None
                                })
                            now += exp_time / 60 / 24
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
                            dark_starts.append({"obsmjd": now,
                                                "exposure_length": self.aqm_time,
                                                "plate": -1,
                                                "cart": -1
                            })
                            now += (self.aqm_time + self.overhead) / 60 / 24
                        else:
                            this_plate = aqmes_obs[sorted_priorities[i]]

                            tonight_ids.append(this_plate["PLATE_ID"])

                            alt, az = self.Observer.radec2altaz(mjd=now,
                                                            ra=this_plate["RA"],
                                                            dec=this_plate["DEC"])
                            exp_time = darkExpTime(alt, default=self.aqm_time)

                            dark_starts.append({"obsmjd": now,
                                                "exposure_length": exp_time,
                                                "plate": this_plate["PLATE_ID"],
                                                "cart": None
                                })
                            now += (exp_time + self.overhead) / 60 / 24
                    else:
                        # everything else failed
                        dark_starts.append({"obsmjd": now,
                                            "exposure_length": self.aqm_time,
                                            "plate": -1,
                                            "cart": -1
                        })
                        now += (self.aqm_time + self.overhead) / 60 / 24

        self.assignCarts(bright_starts, dark_starts)

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


def darkExpTime(alt, default=67):
    """Compute the expected exposure time for a dark plate
       at a given altitude

       alt: float, altitude of the plate

       default: float/int, normal visit length
    """

    airmass = 1. / np.sin(np.pi / 180. * float(alt))

    # compute additional time
    expected_time = default * airmass

    # extra_exp = expected_time / 15.

    return expected_time
