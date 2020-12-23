#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import argparse
import os
import sys

import yaml
import numpy as np
import matplotlib.pyplot as plt
import fitsio
from astropy.time import Time

from platescheduler.scheduler import Scheduler
from platescheduler.weather import Weather


def plotWindow(plate_dict, sched, slots, ax):

    plates = sched.plates
    for i, s in enumerate(slots):
        ax.set_title(int(s["obsmjd"]))
        lst_s = sched.Observer.lst(s["obsmjd"])
        lst_e = sched.Observer.lst(s["obsmjd"] + s["obsmjd"] / 60 / 24)
        ax.scatter([lst_s, lst_e], [i, i], c="r", marker="|")

        if s["plateid"] is None:
            continue
        p = plate_dict[s["plateid"]]
        ax.plot([p["RA"]+p["HA_MIN"], p["RA"]+p["HA_MAX"]], [i, i], linewidth=2, c="k")
        ax.scatter(p["RA"], i, marker="^", c="b")
        ax.axvline(360)

        # print(i, s["plateid"],p["RA"]+p["HA_MIN"], p["RA"]+p["HA_MAX"], lst_s, lst_e, s["obsmjd"], p["CADENCE"])
    return


def monthSim(sched, weather, mjd_start=59135, mjd_end=59165):
    plate_to_cadence = {p["PLATE_ID"]: p["CADENCE"] for p in sched.plates}
    plate_to_cadence[-1] = "FAIL_CAD"

    # it's just easier
    plate_dict = {p["PLATE_ID"]: p for p in sched.plates}

    mjds = np.arange(mjd_start, mjd_end, 1)

    night_scheds = list()

    weather_lost = list()

    no_plate = {"bright": [], "dark": []}

    full_hist = []

    last_engineering = 0
    for m in mjds:
        if sched.Observer.moon_illumination(m) > 0.99:
            if m - last_engineering > 2:
                last_engineering = m
                print("SKIPPING {m} for eng".format(m=m))
                continue

        bright_starts, dark_starts, night_sched = sched.scheduleMjd(m)

        night_scheds.append(night_sched)

    #     plt.figure(figsize=(6,6))
    #     ax1 = plt.subplot(211)
    #     ax2 = plt.subplot(212)

    #     plotWindow(plate_dict, sched, dark_starts, ax1)
    #     plotWindow(plate_dict, sched, bright_starts, ax2)

    #     plt.show()

        for b in bright_starts:
            # try:
            #     print(b["plate"], sched._plateIDtoField[b["plate"]], plate_to_cadence[b["plate"]])
            # except KeyError:
            #     pass
            if weather.clear(b["obsmjd"], returnNext=False):
                # print("WEATHER", b["obsmjd"])
                weather_lost.append(b)
                continue
            if b["plate"] is None or b["plate"] == -1:
                # print(f"NO bright PLATE {float(sched.Observer.lst(b['start']))/15:5.1f} {b['start']:9.3f}")
                lst = float(sched.Observer.lst(b['obsmjd']))/15
                no_plate["bright"].append(lst)
                field = "{lst:4.2f}".format(lst=lst)
                pid = -1
                cadence = "BAD"
            else:
                field = sched._plateIDtoField[b["plate"]]
                pid = b["plate"]
                cadence = plate_to_cadence[b["plate"]]

            full_hist.append({"plate": pid,
                              "field": field,
                              "cadence": cadence,
                              "mjd": b["obsmjd"],
                              "len": b["exposure_length"]})
            if pid > 0:
                sched.obs_hist[field].append(b["obsmjd"])
                plate_dict[pid]["SN_TOT"] += 1601

        for b in dark_starts:
            if weather.clear(b["obsmjd"], returnNext=False):
                # print("WEATHER", b["obsmjd"])
                weather_lost.append(b)
                continue
            if b["plate"] is None or b["plate"] == -1:
                print("NO  dark  PLATE {0:5.1f} {1:9.3f}".format(float(sched.Observer.lst(b['obsmjd']))/15, b['obsmjd']))
                lst = float(sched.Observer.lst(b['obsmjd']))/15
                no_plate["dark"].append(lst)
                field = "{lst:4.2f}".format(lst=lst)
                pid = -1
                cadence = "BAD"
            else:
                field = sched._plateIDtoField[b["plate"]]
                pid = b["plate"]
                cadence = plate_to_cadence[b["plate"]]

            full_hist.append({"plate": pid,
                              "field": field,
                              "cadence": cadence,
                              "mjd": b["obsmjd"],
                              "len": b["exposure_length"]})
            if pid > 0:
                sched.obs_hist[field].append(b["obsmjd"])

    return night_scheds, weather_lost, no_plate, full_hist


def brightField(f):
    if "GG" in f:
        return True
    elif "RV" in f:
        return True
    elif "YSO" in f:
        return True
    else:
        return False


def darkField(f):
    if "AQ" in f:
        return True
    else:
        return False


def plotMonth(sched, night_scheds, weather_lost, no_plate, full_hist, name="oops", seed=1):
    field_to_ra = {p["FIELD"]: p["RA"] + p["HA"] for p in sched.plates}

    not_used = [f for f, o in sched.obs_hist.items() if len(o) == 0]

    # for n in not_used:
    #     print(f"{n:15s}: {field_to_ra[n]/15:4.1f}")

    bright_missed = [field_to_ra[f]/15 for f in not_used if brightField(f)]
    dark_missed = [field_to_ra[f]/15 for f in not_used if darkField(f)]
    plt.figure(figsize=(7, 5))

    bins = np.arange(0, 24, 1)
    plt.hist(no_plate["dark"], bins=bins, color="b", alpha=0.5)
    plt.hist(no_plate["bright"], bins=bins, color="r", alpha=0.5)
    plt.hist(bright_missed, bins=bins, facecolor="None", edgecolor="r", linewidth=1, linestyle="--")
    plt.hist(dark_missed, bins=bins, facecolor="None", edgecolor="b", linewidth=1, linestyle="--")
    plt.savefig("lst_probs_{seed}.pdf".format(seed=seed))

    plt.figure(figsize=(10, 20))
    plt.gca().invert_yaxis()
    for s in night_scheds:
        if int(s["bright_start"]) > 0:
            plt.scatter(s["bright_start"]%1, int(s["bright_start"]), marker="|", c="r", alpha=0.7, s=60)
            plt.scatter(s["bright_end"]%1, int(s["bright_end"]), marker="|", c="r", alpha=0.7, s=60)
        if int(s["dark_start"]) > 0:
            plt.scatter(s["dark_start"]%1, int(s["dark_start"]), marker="|", c="b", alpha=0.7, s=60)
            plt.scatter(s["dark_end"]%1, int(s["dark_end"]), marker="|", c="b", alpha=0.7, s=60)

    for o in full_hist:
        cad = o["cadence"]
        if cad in ["GG", "YSO", "RV6", "RV12"]:
            # print(o)
            color = "r"
        elif cad == "RM":
            color = "grey"
        elif cad == "BAD":
            color = "k"
        else:
            color = "b"
        m = o["mjd"]
        plt.text(m%1 - 0.005, int(m)-0.1, o["field"][:10].replace("_", " "), fontsize=9)
        plt.plot([m%1, m%1 + o["len"]/60/24], [int(m), int(m)], linewidth=4, c=color)
    #     print(f"{m:.3f}", o["len"], o["field"], o["plate"])

    for w in weather_lost:
        m = w["obsmjd"]
    #     print(m)
        plt.text(m%1 - 0.005, int(m)-0.05, "weather", fontsize=9)
        plt.plot([m%1, m%1 + o["len"]/60/24], [int(m), int(m)], linewidth=1, c="k", linestyle="--")

    plt.savefig("{name}_plate_proj_{seed}.pdf".format(name=name, seed=seed))


def wrapForMonths(name="", platePath=None, mjd_start=59135, mjd_end=59165,
                  histFile=None, seed=1):
    sched = Scheduler(session=1, platePath=platePath)

    # force plates to read in
    sched.plates

    sched._plateIDtoField[-1] = "FAIL"

    if histFile is not None:
        sched.obs_hist = yaml.load(open(histFile),
                                   Loader=yaml.Loader)

    weather = Weather(mjd_start=mjd_start,
                      mjd_end=mjd_end,
                      seed=seed, fclear=0.5)

    months = (mjd_end-mjd_start) // 30
    last_month = (mjd_end-mjd_start) - months * 30

    mjd_now = mjd_start

    for m in range(months):
        name = "month_{}_".format(m) + name
        args = monthSim(sched, weather, mjd_start=mjd_now, mjd_end=mjd_now+30)
        plotMonth(sched, *args, name=name,  seed=i)
        mjd_now += 30

    if mjd_end > mjd_now:
        name = "month_{}".format(months) + name
        args = monthSim(sched, weather, mjd_start=mjd_now, mjd_end=mjd_end+1)
        plotMonth(sched, *args, name=name,  seed=i)

    with open("{}_sim_{}.dat".format(name, seed), "w") as outfile:
        for f, m in sched.obs_hist.items():
            print("{f:18s} | {num:3d} | {mjds}".format(f=f, num=len(m), mjds=", ".join([str(int(i)) for i in m])), file=outfile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Simulate a SDSS V plate observing.')

    parser.add_argument("-n", "--name", dest="name", type=str,
                        required=False, help="name this sim",
                        default=None)
    parser.add_argument("-f", "--file", dest="file", type=str,
                        required=False, help="plate input fits file",
                        default=None)
    parser.add_argument("-p", "--prev", dest="prev", type=str,
                        required=False, help="file containing obs hist",
                        default=None)
    parser.add_argument("-i", "--iter", dest="iter", type=int,
                        required=False, help="number of sims to run",
                        default=1)
    parser.add_argument("-s", "--start", dest="start", type=int,
                        required=False, help="mjd to start on",
                        default=None)
    parser.add_argument("-e", "--end", dest="end", type=int,
                        required=False, help="mjd to end on",
                        default=None)

    args = parser.parse_args()
    inputFile = args.file
    histFile = args.prev
    nIter = args.iter
    name = args.name
    start = args.start
    end = args.end

    if inputFile is None:
        inputFile = os.getenv('FIVEPLATE_FITS_INPUT')
        print("Input", inputFile)

    assert inputFile is not None, "Input must be specified"

    if histFile is None:
        histFile = os.getenv('FIVEPLATE_HISTORY')

    if name is None:
        name = "forgot-to-name"

    if start is None:
        start = int(Time(Time.now(), format="mjd").value)

    if end is None:
        end = start + 30

    assert end > start, "start must be before end!"

    for i in range(nIter):
        wrapForMonths(name=name, platePath=inputFile, mjd_start=start, mjd_end=end,
                      histFile=histFile, seed=i)
