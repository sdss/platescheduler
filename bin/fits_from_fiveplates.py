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
    return f"{field:20s}, {ra:11.6f}, {dec:+12.6f}, {ha:+5.1f}, {cadence:13s}, {priority:2d}"


def mkCsv(platerun_root):
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

    return full_csv


def mkFits(csv, outPath):
    """Creates and saves a fits file to output
       using csv.

       csv can be a path to a csv file or the
       csv text created by mkCsv
    """

    plate_array = np.genfromtxt(csv, delimiter=",", dtype=None, names=True)

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


def writeFullCsv(platerun_root, csvPath):
    csv = mkCsv(platerun_root)

    with open(csvPath, "w") as outFile:
        for l in csv:
            print(l, file=outFile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description="""Create fiveplates plate fits input or intermediate csv
                       files.""")

    parser.add_argument("-r", "--root", dest="root", type=str,
                        required=False, help="platescheduler root",
                        default=None)
    parser.add_argument("-c", "--csv", dest="csv", type=str,
                        required=False, help="csv file to use for input",
                        default=None)
    parser.add_argument("-o", "--out", dest="out", type=str,
                        required=False, help="output path, WITHOUT file extension",
                        default=None)

    args = parser.parse_args()
    root = args.root
    outPath = args.out
    csv = args.csv

    if outPath is None:
        outPath = "five_plates"

    fitsPath = outPath + ".fits"

    if root is not None:
        csvPath = outPath + ".csv"
        writeFullCsv(root, csvPath)
        mkFits(csvPath, fitsPath)
    elif csv is not None:
        mkFits(csv, fitsPath)
    else:
        print(" ERROR: incorrect usage. \n",
              "You must specify five_plates root with -r \n",
              "or an input csv file with -c")
