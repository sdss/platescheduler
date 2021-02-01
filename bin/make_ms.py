#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy.time import Time

from platescheduler.scheduler import Scheduler

sched = Scheduler(session=1, platePath="five_plates.fits")
# bright_starts, dark_starts, night_sched = sched.scheduleMjd(59185)
# print(night_sched)


# # startTime = Time(59185, format="mjd").datetime
# startTime = Time(59185.048022177776, format="mjd").datetime
# endTime = Time(59185.525353942721, format="mjd").datetime

# # startTime.strftime("%a %d %b %Y")

# print(startTime.time(), startTime.hour, startTime.minute)

def lst_to_hhmm(lst):
    hour = int(lst // 1)
    minute = int((lst - hour) * 60 // 1)
    return f"{hour:02d}:{minute:02d}"


# print("16.5 -> ", lst_to_hhmm(16.5))
# print("6.1  -> ", lst_to_hhmm(6.1))

last_engineering = 0


def night_string(m):
    global last_engineering
    m = int(m)
    engineering = "---- "
    if sched.Observer.moon_illumination(m) > 0.98:
        engineering = "Day  "
        if m - last_engineering > 2:
            last_engineering = m
            engineering = "Night"

    bright_starts, dark_starts, night_sched = sched.scheduleMjd(m)

    startTime = Time(night_sched["start"], format="mjd").datetime
    endTime = Time(night_sched["end"], format="mjd").datetime

    bright_start = Time(night_sched["bright_start"], format="mjd").datetime
    bright_end = Time(night_sched["bright_end"], format="mjd").datetime

    dark_start = Time(night_sched["dark_start"], format="mjd").datetime
    dark_end = Time(night_sched["dark_end"], format="mjd").datetime

    lst_start = float(sched.Observer.lst(night_sched["start"]))/15
    lst_end = float(sched.Observer.lst(night_sched["end"]))/15

    lst_start = lst_to_hhmm(lst_start)
    lst_end = lst_to_hhmm(lst_end)

    aTime = Time(m-1, format="mjd").datetime

    # day dd MON yyyy
    day_format = aTime.strftime("%a %d %b %Y")
    moon = sched.Observer.moon_illumination(np.median([night_sched["start"], night_sched["end"]]))

    form_str = f"{day_format:18s} {m:5d}  {float(moon):.2f}  " +               f"{int(startTime.hour):02d}:{int(startTime.minute):02d}  " +               f"{int(endTime.hour):02d}:{int(endTime.minute):02d}   " +               f"{lst_start}  {lst_end}   {engineering}   XX  " 

    dark_time = (night_sched["dark_end"] - night_sched["dark_start"]) * 24
    bright_time = (night_sched["bright_end"] - night_sched["bright_start"]) * 24

    if bright_time and not dark_time:
        time_sum = f"MWM ({bright_time:4.1f} h)"
    elif dark_time and not bright_time:
        time_sum = f"MWM ({dark_time:4.1f} h)"
    elif night_sched["bright_start"] < night_sched["dark_start"]:
        time_sum = f"MWM until {bright_end.hour:02d}:{bright_end.minute:02d} ({bright_time:4.1f} h), " +                   f"then BHM ({dark_time:4.1f} h)"
    elif night_sched["dark_start"] < night_sched["bright_start"]:
        time_sum = f"BHM until {dark_end.hour:02d}:{dark_end.minute:02d} ({dark_time:4.1f} h), " +                   f"then MWM ({bright_time:4.1f} h)"
    else:
        assert False, "you broke something"

    form_str += time_sum
    return form_str


start_mjd = 59185

with open(f"ms_{start_mjd:5d}.dat", "w") as outfile:
    mjds = np.arange(start_mjd, start_mjd + 90, 1)
    for m in mjds:
        print(night_string(m), file=outfile)
