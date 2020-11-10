# platescheduler

![Versions](https://img.shields.io/badge/python->3.7-blue)
[![Documentation Status](https://readthedocs.org/projects/sdss-platescheduler/badge/?version=latest)](https://sdss-platescheduler.readthedocs.io/en/latest/?badge=latest)
[![Travis (.org)](https://img.shields.io/travis/sdss/platescheduler)](https://travis-ci.org/sdss/platescheduler)
[![codecov](https://codecov.io/gh/sdss/platescheduler/branch/master/graph/badge.svg)](https://codecov.io/gh/sdss/platescheduler)

autoscheduler for sdss5 plate program, now with simulation tools!

## Install

To use the simulations, checkout a copy of the repo and switch to the testFromFits branch:

    git clone https://github.com/sdss/platescheduler.git 
    git checkout testFromFits

To install locally, it is recommended to create a virtual environment to contain this project. From within the environment, it is recommended you install in editable mode so you can `git pull` my bug fixes:

    pip install -e .

Alternatively, a moduleFile setup would work if all dependencies specified in setup.cfg are satisfied. Just make sure to add `bin/` to `$PATH` in the module file.

## Usage

To run a simulation, a five_plates fits input file containing all plates to include must be generated. To make this file, a `fits_from_fiveplates.py` executable is added to your path during pip installation. It takes one required argument, -r or --root: the path to five_plates/plateruns, and one optional argument -o or --out: the file name or path to save the output fits file. By default a `five_plats.fits` file is created in the current directory. 

Example usage:

    fits_from_fiveplates.py -r /home/john/software/five_plates/plateruns 

Will create a `five_plates.fits` in the current directory. It may be useful to create an empty sim directory and work from that directory, so all sim files are created in a convenient location.

`plate_sim.py` was also added to your path. This script uses the latest platescheduler logic and a `five_plates.fits` file to simulate 30 days of plate observations. 

At minimum a `five_plates.fits` file is needed (though it can be named anything). The file can be specified with a -f or --file argument or be set in a FIVEPLATE_FITS_INPUT environment variable (useful if it is deep in the five_plates product). 

The simulation defaults to MJD 59135, which is probably not useful. Specify a start MJD with -m or --mjd arguments. 

A name for this sim can be specified with a -n or --name argument. The number of iterations can be specified with a -i or --iter arguments; N simulations will be run with different random seeds for weather. 

Example usage:
    
    export FIVEPLATE_FITS_INPUT=/home/john/software/five_plates/five_plates.fits

    plate_sim.py -n test -p five_plates_hist.yaml -m 59160 -i 2

OR if working in dedicated directory with a `five_plates.fits` file just created:

    plate_sim.py -n test -f five_plates.fits -p five_plates_hist.yaml -m 59160 -i 2

