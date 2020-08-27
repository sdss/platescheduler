#!/usr/bin/env python
# encoding: utf-8
#
# plate_plans_db.py
#
# Originally created by Demitri Muna in 2004
# Rewritten by José Sánchez-Gallego on 14 Jun 2017.


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import hashlib
import os
import six
import warnings

import peewee

from astropy import table

from sdssdb.peewee.operationsdb import platedb
from sdssdb.peewee.operationsdb import database

import numpy as np
import fitsio


@platedb.database.atomic()
def _load_design(design_id, definition, overwrite=True):
    """Loads a design into the DB."""

    design_dbo, created = platedb.Design.get_or_create(pk=design_id)

    if not created:
        print('found design {0} in the DB.'.format(design_id))
        if overwrite:
            warnings.warn('overwriting design for design_id={0}.'.format(design_id),
                          UserWarning)
        else:
            return
    else:
        print('creating new Design for design_id={0}.'.format(design_id))

    # definition = definition_from_id(design_id)

    # Delete old values (if present; easier than syncing).
    for xx in design_dbo.values:
        xx.delete_instance()

    for key in definition.keys():

        design_value_dbo = platedb.DesignValue()
        design_value_dbo.design = design_dbo
        design_value_dbo.value = definition[key]

        design_field_dbo, created = platedb.DesignField.get_or_create(label=key)
        if created:
            print('created new row in DesignField with value label={0!r}'.format(key))

        design_value_dbo.field = design_field_dbo

        design_value_dbo.save()

    pointing_dbo, created = platedb.Pointing.get_or_create(design_pk=design_dbo.pk,
                                                               pointing_no=1)
    pointing_dbo.design = design_dbo
    pointing_dbo.center_ra = definition['racen']
    pointing_dbo.center_dec = definition['deccen']
    pointing_dbo.pointing_no = 1

    pointing_dbo.save()

    print("!! pointing!", pointing_dbo)

    # Handle inputs (PlateInput)

    # Delete any old inputs if present.
    # for xx in design_dbo.inputs:
    #     xx.delete_instance()

    # priority_list = 'priority' in definition and list(map(int, definition['priority'].split()))

    # for key in definition:
    #     if not key.startswith('plateinput'):
    #         continue
    #     input_number = int(key.strip('plateinput'))
    #     priority = priority_list[input_number - 1] if priority_list else None
    #     filepath = definition[key]

    #     plate_input_dbo = platedb.PlateInput()
    #     plate_input_dbo.design = design_dbo
    #     plate_input_dbo.input_number = input_number
    #     plate_input_dbo.priority = priority
    #     plate_input_dbo.filepath = filepath

    #     plate_input_full_path = os.path.join(os.environ['PLATELIST_DIR'], 'inputs', filepath)

    #     if not os.path.exists(plate_input_full_path):
    #         warnings.warn('cannot find plateInput {0}. '
    #                       'MD5 check sum will be null.'.format(filepath), UserWarning)
    #     else:
    #         plate_input_dbo.md5_checksum = hashlib.md5(
    #             open(plate_input_full_path).read().encode('utf-8')).hexdigest()

    #     print('added plateInput file {0}'.format(filepath))
    #     plate_input_dbo.save()

    # # Create Pointings
    # no_pointings = int(definition['npointings'])

    # # If the design already has pointings it must be because we are loading a
    # # plate from a design already loaded. In that case we check that the number
    # # of pointings loaded matches.
    # if len(design_dbo.pointings) > 0 and len(design_dbo.pointings) != no_pointings:
    #     # If the number of pointins disagree but they do not have plate_pointings
    #     # associated, we can remove the pointings and start from scratch.
    #     no_plate_pointings = np.sum([len(pointing.plate_pointings)
    #                                  for pointing in design_dbo.pointings])

    #     if no_plate_pointings > 0:
    #         raise RuntimeError('design_id {0} has pointings with '
    #                            'already created plate_pointings. '
    #                            'This requires manual intervention.'
    #                            .format(design_id))
    #     else:
    #         for pointing_dbo in design_dbo.pointings:
    #             pointing_dbo.delete_instance()

    # for pno in range(1, no_pointings + 1):
    #     pointing_dbo, created = platedb.Pointing.get_or_create(design_pk=design_dbo.pk,
    #                                                            pointing_no=pno)
    #     pointing_dbo.design = design_dbo
    #     pointing_dbo.center_ra = definition['racen'].split()[pno - 1]
    #     pointing_dbo.center_dec = definition['deccen'].split()[pno - 1]
    #     pointing_dbo.pointing_no = pno

    #     pointing_dbo.save()

    #     if created:
    #         print('created pointing #{0} for design {1}'.format(pno, design_id))
    #     else:
    #         print('found pointing #{0} for design {1} in DB'.format(pno, design_id))


@platedb.database.atomic()
def _load_plate(plate_id, plateplans_line, overwrite=True):
    """Loads a plate and plate_pointing info to the DB."""

    # Does this plate exist in the database?
    plate_dbo, created = platedb.Plate.get_or_create(plate_id=plate_id)

    if not created:
        print('found plate {0} in the DB.'.format(plate_id))
        if overwrite:
            warnings.warn('overwriting plate for plate_id={0}.'.format(plate_id),
                          UserWarning)
        else:
            return
    else:
        print('creating new Plate for plate_id={0}.'.format(plate_id))

    plate_dbo.location_id = plateplans_line['locationid']
    plate_dbo.temperature = 15.0
    plate_dbo.epoch = round(2020.75, 6)
    plate_dbo.center_ra = round(plateplans_line['raCen'], 6)
    plate_dbo.center_dec = round(plateplans_line['decCen'], 6)
    # plate_dbo.rerun = plateplans_line['rerun']
    # plate_dbo.chunk = plateplans_line['chunk']
    plate_dbo.chunk = "2020.8.a.mwm-bhm-fake"

    if plateplans_line['name'] != "''" and len(plateplans_line['name']) > 0:
        plate_dbo.name = plateplans_line['name']

    # plate_dbo.comment = plateplans_line['comments']
    plate_dbo.comment = "999"

    # Tile info
    # tileid = plateplans_line['tileid']
    # if tileid > -1:
    #     plate_dbo.tile_id = tileid
    #     tile_dbo, created = platedb.Tile.get_or_create(id=tileid)
    #     if not created:
    #         print('found tile {0} in the DB'.format(tileid))
    #     else:
    #         print('created new tile with id={0}'.format(tileid))

    #     tile_dbo.save()

    #     plate_dbo.tile = tile_dbo

    # plate_dbo.epoch = round(plateplans_line['epoch'], 6)
    # plate_dbo.center_ra = round(plateplans_line['raCen'], 6)
    # plate_dbo.center_dec = round(plateplans_line['decCen'], 6)

    # Assigns the platerun
    try:
        platerun_dbo = platedb.PlateRun.get(label="2020.8.a.mwm-bhm-fake")
    except peewee.DoesNotExist:
        raise ValueError('cannot found a PlateRun row for plate {0}. '
                         'The design should already be in the DB.'.format(plate_id))

    plate_dbo.plate_run = platerun_dbo

    # Sets the plate status to design.
    design_status = platedb.PlateStatus.get(label='Design')
    plate_dbo.statuses.clear()  # First remove statuses
    plate_dbo.statuses.add([design_status])

    # Handle survey relationship
    plate_dbo.surveys.clear()  # First remove surveys
    for survey in six.u(plateplans_line['survey']).split('-'):
        plate_dbo.surveys.add([platedb.Survey.get(plateplan_name=survey)])

    # Ensure "design" foreign key constraint is met (lookup design from db).
    design_id = plateplans_line['designid']
    try:
        # Look for existing design in the database.
        design_dbo = platedb.Design.get(pk=design_id)
        print('found design {0} for plate {1}.'.format(design_id, plate_id))
    except peewee.DoesNotExist:
        raise ValueError('cannot found a Design for plate {0}. '
                         'The design should already be in the DB.'.format(plate_id))

    plate_dbo.design = design_dbo

    # The default survey mode key needs to also be written to the plate table
    defaultsurveymode_dbo = platedb.DesignValue.select().join(
        platedb.DesignField).where((platedb.DesignValue.design_pk == design_dbo.pk) &
                                   (platedb.DesignField.label == 'defaultsurveymode'))

    if len(defaultsurveymode_dbo) == 0:
        warnings.warn('cannot find defaultsurveymode for '
                      'design {0} for plate {1}. '
                      'Not setting current_survey_mode'.format(design_dbo.pk, plate_id))
    else:
        defaultsurveymode = defaultsurveymode_dbo[0].value

        survey_mode_pk = platedb.SurveyMode.select(platedb.SurveyMode.pk).where(
            platedb.SurveyMode.definition_label ** defaultsurveymode).scalar()

        if not survey_mode_pk:
            raise RuntimeError('The database is missing an entry in \'survey_mode\' '
                               'for the entry {0!r}.'.format(defaultsurveymode))

        plate_dbo.current_survey_mode_pk = survey_mode_pk

    plate_dbo.save()

    # PlatePointings
    # The actual instance of a telescope pointing - the parameters of Pointing
    # plus an actual plate and hour angle.

    for pointing_dbo in plate_dbo.design.pointings:

        plate_pointing_dbo, created = platedb.PlatePointing.get_or_create(
            pointing_pk=pointing_dbo.pk, plate_pk=plate_dbo.pk,
            defaults={'pointing_name': 'A'})

        if not created:
            print('found plate_pointing for plate_id={0} in DB.'.format(plate_id))
            return

        pno = pointing_dbo.pointing_no

        plate_pointing_dbo.design = design_dbo
        plate_pointing_dbo.pointing = pointing_dbo
        plate_pointing_dbo.plate = plate_dbo
        plate_pointing_dbo.hour_angle = plateplans_line['ha']
        plate_pointing_dbo.ha_observable_min = plateplans_line['ha_min']
        plate_pointing_dbo.ha_observable_max = plateplans_line['ha_max']

        # pointing_name = platedb.DesignValue.select().join(
        #     platedb.DesignField).where((platedb.DesignValue.design_pk == design_dbo.pk) &
        #                                (platedb.DesignField.label == 'pointing_name'))

        # if len(pointing_name) == 0:
        #     raise ValueError('cannot find pointing_name for '
        #                      'design {0} for plate {1}'.format(design_dbo.pk, plate_id))

        # plate_pointing_dbo.pointing_name = pointing_name[0].value.split()[pno - 1]

        # Sets the priority to 5 for MaNGA and APOGEE, 2 for eBOSS
        survey = six.u(plateplans_line['survey'])
        if 'manga' in survey or 'apogee' in survey:
            plate_pointing_dbo.priority = 5
        else:
            plate_pointing_dbo.priority = 5

        plate_pointing_dbo.save()
        print('created plate_pointing for plate_id={}.'
                  .format(plate_id))


def plate_plans_db(plate_id, design_id, plate_line, design_defs, overwrite=True):
    platerun = "2020.8.a.mwm-bhm-fake"

    # Checks the connection
    conn_status = platedb.database.connected
    if conn_status:
        print('database connection is open.')
    else:
        raise RuntimeError('cannot connect to the database. Review you connection settings.')

    pr, created = platedb.PlateRun.get_or_create(label=platerun, year=2020)

    if not created:
        print('platerun {0} already is already in the DB.'.format(platerun))
    else:
        print('added platerun {0} to the plate_run table.'.format(platerun))

    # design_ids = np.unique(run_lines['designid'])

    # design_ids = [9100]

    # for design_id in design_ids:

    
    print('loading design_id={0}'.format(design_id))
    _load_design(design_id, design_defs, overwrite=overwrite)

    # if load_addenda:
    #     log.important('loading plateDefinitionAddendas ...')
    #     plate_addenda_db(design_ids, design_mode=True, log=log)

    # plate_ids = np.unique(run_lines['plateid'])

    # for plate_id in plate_ids:

    #     plate_line = run_lines[run_lines['plateid'] == plate_id][0]

    print('loading plate_id={0}'.format(plate_id))
    _load_plate(plate_id, plate_line, overwrite=overwrite)

    # log.important('populating observing ranges for {0} ... '.format(platerun))
    # populate_obs_range(plate_ids, log=log)


def fitsToDb(line, designid):
    plateplans_line = {}

    if line["CADENCE"] in ["YSO", "RV6", "RV12", "GG"]:
        survey = "mwm"
    else:
        survey = "bhm"

    plateplans_line['survey'] = survey
    plateplans_line['ha'] = line["HA"]
    plateplans_line['ha_min'] = line["HA_MIN"]
    plateplans_line['ha_max'] = line["HA_MAX"]
    plateplans_line['designid'] = designid
    plateplans_line['locationid'] = designid
    plateplans_line['raCen'] = line["RA"]
    plateplans_line['decCen'] = line["DEC"]
    plateplans_line['name'] = line["FIELD"]

    design_defs = {"cadence": line["CADENCE"],
                   "racen": line["RA"],
                   "deccen": line["DEC"]}


    plate_plans_db(line["PLATE_ID"], designid, plateplans_line, design_defs)


if __name__ == "__main__":
    database.connect_from_parameters(dbname="apodb", user="sdssdb_admin", 
                                     host="localhost", port="5500")

    plates = fitsio.read("/home/john/Downloads/first_plates.fits")

    for i, p in enumerate(plates):
        fitsToDb(p, 91000+i)
