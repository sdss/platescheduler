# encoding: utf-8

import numpy as np


def rmCadence(mjd, hist=[]):
    """cadence requirements for RM

    Request: 10 epochs per month, month being dark run
    so check the last 15 days. This may bleed over too
    much from previous dark run? It shouldn't though,
    not with a bright run in-between

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """

    if len(hist) == 0:
        return True

    deltas = mjd - np.array(hist)
    this_month = np.where(deltas < 15)

    return len(this_month[0]) < 10


def AQMESmedium(mjd, hist=[]):
    """cadence requirements for AQMES-Medium

    Request: 6 epochs, 1/month, again assume 15-day run

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """

    if len(hist) == 0:
        return True

    if len(hist) > 6:
        return False

    delta = mjd - np.max(hist)
    
    # make sure it's been at least since last run
    return delta > 15


def single(mjd, hist=[]):
    """cadence requirements for single-epoch

    Request: single epoch

    mjd: float or int should be ok

    hist: list, list of previous MJDs
    """

    return len(hist) == 0


def yso(mjd, hist=[]):
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


def rv6(mjd, hist=[]):
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


def rv12(mjd, hist=[]):
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

    if len(deltas) > 2:
        return False

    # would allow for observations more than a week after previous
    return np.min(deltas) > 2


def tess(mjd, hist=[]):
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