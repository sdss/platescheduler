#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import fitsio

from platescheduler.scheduler import Scheduler
from observesim.weather import Weather


def plotWindow(plate_dict, sched, slots, ax):
    
    plates = sched.plates
    for i, s in enumerate(slots):
        ax.set_title(int(s["obsmjd"]))
        lst_s = sched.Observer.lst(s["obsmjd"])
        lst_e = sched.Observer.lst(s["obsmjd"] + s["obsmjd"] / 60 / 24)
        ax.scatter([lst_s, lst_e], [i,i], c="r", marker="|")
        
        if s["plateid"] is None:
            continue
        p = plate_dict[s["plateid"]]
        ax.plot([p["RA"]+p["HA_MIN"], p["RA"]+p["HA_MAX"]], [i,i], linewidth=2, c="k")
        ax.scatter(p["RA"], i, marker="^", c="b")
        ax.axvline(360)
        
#         print(i, s["plateid"],p["RA"]+p["HA_MIN"], p["RA"]+p["HA_MAX"], lst_s, lst_e, s["obsmjd"], p["CADENCE"])
    return


def monthSim(seed=1):
    sched = Scheduler(session=1)
    plate_to_cadence = {p["PLATE_ID"]: p["CADENCE"] for p in sched.plates}
    plate_to_cadence[-1] = "FAIL_CAD"

    # it's just easier
    plate_dict = {p["PLATE_ID"]: p for p in sched.plates}

    sched._plateIDtoField[-1] = "FAIL"

    # first day!
    start = 59140
    # take out a full moon night somewhere

    mjds = np.arange(start, start+30, 1)

    weather = Weather(mjd_start=start,
                      mjd_end=mjds[-1],
                      seed=seed, fclear=0.5)

    night_scheds = list()

    weather_lost = list()
    
    no_plate = {"bright":[], "dark":[]}

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
    #             print("WEATHER", b["obsmjd"])
                weather_lost.append(b)
                continue
            if b["plate"] is None or b["plate"] == -1:
#                 print(f"NO bright PLATE {float(sched.Observer.lst(b['start']))/15:5.1f} {b['start']:9.3f}")
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

        for b in dark_starts:
            if weather.clear(b["obsmjd"], returnNext=False):
    #             print("WEATHER", b["obsmjd"])
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

    return sched, night_scheds, weather_lost, no_plate, full_hist

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

def plotMonth(sched, night_scheds, weather_lost, no_plate, full_hist, seed=1):
    field_to_ra = {p["FIELD"]: p["RA"] + p["HA"] for p in sched.plates}
    
    not_used = [f for f, o in sched.obs_hist.items() if len(o) == 0]
    
    # for n in not_used:
    #     print(f"{n:15s}: {field_to_ra[n]/15:4.1f}")
    
    bright_missed = [field_to_ra[f]/15 for f in not_used if brightField(f)]
    dark_missed = [field_to_ra[f]/15 for f in not_used if darkField(f)]
    plt.figure(figsize=(7,5))
    
    bins = np.arange(0, 24, 1)
    plt.hist(no_plate["dark"], bins=bins, color="b", alpha=0.5)
    plt.hist(no_plate["bright"], bins=bins, color="r", alpha=0.5)
    plt.hist(bright_missed, bins=bins, facecolor="None", edgecolor="r", linewidth=1, linestyle="--")
    plt.hist(dark_missed, bins=bins, facecolor="None", edgecolor="b", linewidth=1, linestyle="--")
    plt.savefig("lst_probs_{seed}.pdf".format(seed=seed))
    
    plt.figure(figsize=(10,20))
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
            print(o)
            color = "r"
        elif cad == "RM":
            color = "grey"
        elif cad == "BAD":
            color = "k"
        else:
            color = "b"
        m = o["mjd"]
        plt.text(m%1 - 0.005, int(m)-0.1, o["field"][:10].replace("_", " "), fontsize=9)
        plt.plot([m%1, m%1 + o["len"]/60/24 -20/60/24], [int(m), int(m)], linewidth=4, c=color)
    #     print(f"{m:.3f}", o["len"], o["field"], o["plate"])

    for w in weather_lost:
        m = w["obsmjd"]
    #     print(m)
        plt.text(m%1 - 0.005, int(m)-0.05, "weather", fontsize=9)
        plt.plot([m%1, m%1 + o["len"]/60/24], [int(m), int(m)], linewidth=1, c="k", linestyle="--")

    plt.savefig("nov_plate_proj_{seed}.pdf".format(seed=seed))

    fields = defaultdict(list)
    for o in full_hist:
        fields[o["field"]].append(int(o["mjd"]))
    with open("nov_sim_{}.dat".format(seed), "w") as outfile:
        for f, m in fields.items():
            if f == "EMPTY":
                continue
            print("{f:18s} | {num:3d} | {mjds}".format(f=f, num=len(m), mjds=", ".join([str(i) for i in m])), file=outfile)


if __name__ == "__main__":
    for i in range(5):
        args = monthSim(seed=i)
        plotMonth(*args, seed=i)

    seed = 0
    args = monthSim(seed=seed)
    plotMonth(*args, seed=seed)


