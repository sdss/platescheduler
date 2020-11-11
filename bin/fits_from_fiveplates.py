#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys

import numpy as np
import fitsio


def windowWidth(dec, designHA):
    # dec and designHA in degrees!
    # compute plates obs window in hours
    obs_win = 6 - ( np.abs(dec - 40) + np.abs(designHA))/15  # *1h
    return obs_win*15


def parseLine(line):
    els = line.split()
    field = els[0].strip()
    ra = float(els[1])
    dec = float(els[2])
    ha = float(els[5])
    cadence = els[6].strip()
    priority = int(els[7][:1])//2
    return f"{field}, {ra}, {dec}, {ha}, {cadence}, {priority}"


def convert(platerun_root, outPath):
    five_dirs = os.listdir(platerun_root)
    use_files = list()
    for d in five_dirs:
        test_file = platerun_root + d + "/" + d + ".summary"
        if os.path.isfile(test_file):
            use_files.append(test_file)

    full_csv = ["FIELD, RA, DEC, HA, CADENCE, PRIORITY"]

    for f in use_files:
        with open(f, "r") as sumFile:
            for l in sumFile.readlines():
                if l[0:5] == "Field"\
                   or len(l.strip()) == 0\
                   or l[0] == "#":
                    continue
                full_csv.append(parseLine(l))

    plate_array = np.genfromtxt(full_csv, delimiter=",", dtype=None, names=True)

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

    plates = np.zeros(len(plate_array), dtype=plate_types)
    plates["PLATE_ID"] = 2e4 + np.arange(len(plate_array))
    plates["FIELD"] = np.array([c.decode().strip() for c in plate_array["FIELD"]], dtype=np.dtype("a20"))
    plates["RA"] = plate_array["RA"]
    plates["DEC"] = plate_array["DEC"]
    plates["HA"] = plate_array["HA"]
    obs_win = windowWidth(plate_array["DEC"], plate_array["HA"])/2
    plates["HA_MIN"] = plate_array["HA"] - obs_win
    plates["HA_MAX"] = plate_array["HA"] + obs_win
    plates["SKYBRIGHTNESS"] = np.array([1 if c.decode().strip() in ["RV6", "RV12", "YSO", "GG"] else 0.35 for c in plate_array["CADENCE"]])
    plates["CADENCE"] = np.array([c.decode().strip() for c in plate_array["CADENCE"]], dtype=np.dtype("a20"))
    plates["PRIORITY"] = plate_array["PRIORITY"]

    fitsio.write(outPath, plates, clobber=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Create fiveplates plate fits input')

    parser.add_argument("-r", "--root", dest="root", type=str,
                        required=True, help="platescheduler root")
    parser.add_argument("-o", "--out", dest="out", type=str,
                        required=False, help="output path",
                        default=None)

    args = parser.parse_args()
    root = args.root
    outPath = args.out

    if outPath is None:
        outPath = "five_plates.fits"

    convert(root, outPath)
