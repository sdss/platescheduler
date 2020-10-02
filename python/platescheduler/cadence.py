# encoding: utf-8

import numpy as np


def rmCadence(mjd, hist=[], **kwargs):
    """cadence requirements for RM

    Request: 10 epochs per month, a dark run is entirely between full moons

    For now: if it hasn't been 10, do it! Could add more complicated logic, e.g.
    wait a day if it's near beginning of run, cram at the end

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """
    last_full_moon = kwargs.get("last_full_moon", 0)

    if len(hist) == 0:
        return True

    deltas = np.array(hist) - last_full_moon
    this_month = np.where(deltas > 0)

    deltas = mjd - np.array(hist)

    return len(this_month[0]) < 5 and np.min(deltas) > 3


def AQMESmedium(mjd, hist=[], **kwargs):
    """cadence requirements for AQMES-Medium

    Request: 6 epochs, 1/month

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """
    last_full_moon = kwargs.get("last_full_moon", 0)

    if len(hist) == 0:
        return True

    if len(hist) > 6:
        return False

    deltas = np.array(hist) - last_full_moon
    this_month = np.where(deltas > 0)
    
    # make sure it's been at least since last run
    return len(this_month[0]) == 0


def single(mjd, hist=[], **kwargs):
    """cadence requirements for single-epoch

    Request: single epoch

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """

    return len(hist) == 0


def yso(mjd, hist=[], **kwargs):
    """cadence requirements for YSO

    Request: 3 epochs, ~3 days, then ~20 days

    mjd: float, or int should be ok

    hist: list, list of previous MJDs
    """

    if len(hist) == 0:
        return True
    elif len(hist) == 1:
        return mjd - hist[0] > 2
    elif len(hist) > 1:
        # this attempts to allow for more exposures than expected
        return mjd - hist[-1] > 19
    else:
        return False


def rv6(mjd, hist=[], **kwargs):
    """cadence requirements for 6-visit RV plates

    Request: 6 total, ~3 per month, ideally within 1 week

    mjd: float, or int should be ok

    hist: list, list of previous MJDs
    """

    if len(hist) == 0:
        return True

    if len(hist) > 6:
        return False

    deltas = mjd - np.array(hist)
    this_month = deltas[np.where(deltas < 15)]

    if len(deltas) > 2:
        return False

    # would allow for observations more than a week after previous
    return np.min(deltas) > 2


def rv12(mjd, hist=[], **kwargs):
    """cadence requirements for 12-visit RV plates

    Request: 12 total, ~3 per month, ideally within 1 week

    mjd: float, or int should be ok

    hist: list, list of previous MJDs
    """

    if len(hist) == 0:
        return True

    if len(hist) > 12:
        return False

    deltas = mjd - np.array(hist)
    this_month = deltas[np.where(deltas < 15)]

    if len(deltas) > 3:
        return False

    # would allow for observations more than a week after previous
    return np.min(deltas) > 1


def tess(mjd, hist=[], **kwargs):
    """cadence requirements for tess

    Request: 2 observations?

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """

    if len(hist) == 0:
        return True

    if len(hist) > 1:
        return False

    return mjd - hist[0] > 1


def assignCadence(name):
    """for now lets just have a dictrionary serve as a 
    look-up table
    """

    lookUpCadence = {"RM": rmCadence,
                     "AQMES-Medium": AQMESmedium,
                     "RV6": rv6,
                     "RV12": rv12,
                     "YSO": yso,
                     "TESS": tess}

    # match cadence, if missing use single
    # e.g. aqmes-wide, gg use single
    return lookUpCadence.get(name, single)