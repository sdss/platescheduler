# platescheduler

![Versions](https://img.shields.io/badge/python->3.7-blue)
[![Documentation Status](https://readthedocs.org/projects/sdss-platescheduler/badge/?version=latest)](https://sdss-platescheduler.readthedocs.io/en/latest/?badge=latest)
[![Travis (.org)](https://img.shields.io/travis/sdss/platescheduler)](https://travis-ci.org/sdss/platescheduler)
[![codecov](https://codecov.io/gh/sdss/platescheduler/branch/master/graph/badge.svg)](https://codecov.io/gh/sdss/platescheduler)

autoscheduler for sdss5 plate program, now with master schedule tools!

## Install

To use the simulations, checkout a copy of the repo and switch to the makeMS branch:

    git clone https://github.com/sdss/platescheduler.git 
    git checkout makeMS

To install locally, it is recommended to create a virtual environment to contain this project. From within the environment, it is recommended you install in editable mode so you can `git pull` my bug fixes:

    pip install -e .

Alternatively, a moduleFile setup would work if all dependencies specified in setup.cfg are satisfied. Just make sure to add `bin/` to `$PATH` in the module file.

## Usage

To make a master scheduler, a `make_ms` script is supplied in the `bin` directory, available after a pip install.

By default, a schedule for the next 30 days will be created in the current directory. 4 input parameters are provided to tune this output:

`-s`, `--start`: mjd to start on
`-e`, `--end`: mjd to end on
`-o`, `--out`: path to output directory
`-n`, `--name`: name of the output file (not include extension)

Example usage:

    make_ms -o ~/Downloads/ -s 59200 -e 59300 -n testMS 

Will create a `testMs.dat` file in ~/Downloads covering the 100 days from 59200 to 59300.
