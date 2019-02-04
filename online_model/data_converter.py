"""
Module online_model.data_converter
---------------------------------------
Functions to convert the online model output files.
"""
import re
import os
import pytz
from datetime import datetime
import numpy as np

import online_model.constants as const
from tfs_files import tfs_pandas as tfs
from collections import OrderedDict
from utils.contexts import suppress_exception

# Globals for finding data in output file ######################################
PRE_BEAMPROCESS_STRING = "beamprocess"  # word before beamprocess name (line 0)
PRE_KNOB_STRING = "knob"                # word before knobname (line 1)
PRE_TIME_STRING = "time"                # word before time (line 1)
TRIM_VALUE_NAME = r"bbmext\.knob\d*"    # matcher for trim value name
OPTIC_NAME = "OPTIC_NAME"               # ID for the optic name
FILLNUMBER = "Fillnumber"               # Header in orbit file

BEAMPROCESS_HEADER = const.get_beamprocess_header()
TRIM_HEADER = const.get_trim_header()
TRIMHINT_HEADER = const.get_trimhint_header()
TIME_HEADER = const.get_time_header()
KNOB_HEADER = const.get_knob_header()


NAME_COLUMN = const.get_name_column(with_type=True)
CIRCUIT_COLUMN = const.get_circuit_column(with_type=True)
VALUE_COLUMN = const.get_value_column(with_type=True)
DELTA_COLUMN = const.get_delta_column(with_type=True)
TRIM_COLUMN = const.get_trim_column(with_type=True)
COLUMNS_K = OrderedDict([NAME_COLUMN, CIRCUIT_COLUMN, VALUE_COLUMN, DELTA_COLUMN])
RE_K_LINE = r"\s*(?P<NAME>\S+)\s*=\s*(?P<VALUE>[^;]+)[^:]+:\s*(?P<CIRCUIT>[^/]+/[^.]+)\.[^\d+\-]+(?P<DELTA>[^.]+\.[^.]+)\."
RE_TRIM = r"\s*(?P<trimname>" + TRIM_VALUE_NAME + r")\s*=\s*(?P<trimvalue>[^;]+)\s*;"


def knobs_k_to_tfs(knob_file=const.get_extractor_knobs_filename(), trim=None):
    """ Convert from knobs.madx output file to a proper tfs file.

    Assumes no one changes the output order. See global variables is this file.

    Args:
        knob_file: path to online model extractor K-Type output file.
        trim: trim value if available.

    Returns:
        List of TfsDataFrames, one for each knob.
    """
    with open(knob_file, "r") as f:
        content = f.readlines()

    ls = content[0].split()
    beamprocess = ls[ls.index(PRE_BEAMPROCESS_STRING) + 1]
    del content[0]

    df_list = []
    while len(content) > 0:
        df = tfs.TfsDataFrame(columns=COLUMNS_K.keys())

        # headers
        df.headers[BEAMPROCESS_HEADER] = beamprocess

        ls = content[0].split()
        df.headers[KNOB_HEADER] = ls[ls.index(PRE_KNOB_STRING) + 1]

        time_idx = ls.index(PRE_TIME_STRING)
        df.headers[TIME_HEADER] = convert_local_time_to_utc(
            " ".join(ls[time_idx + 1:time_idx+3])[:-1]
        )

        if trim is not None:
            df.headers[TRIM_HEADER] = trim
            df.headers[TRIMHINT_HEADER] = "VALUE = BEFORE_VALUE + DELTA * TRIM"

        # get data
        pattern = re.compile(RE_K_LINE)
        idx = 0
        for idx, line in enumerate(content[1:]):
            if line.startswith("!"):
                break  # next knob, keep last line

            match = re.match(pattern, line)
            for col in COLUMNS_K.keys():
                    df.loc[idx, col] = match.group(col)
        else:
            idx += 1  # no further knob, delete last line

        for col in COLUMNS_K.keys():
            with suppress_exception(KeyError):
                df[col] = df[col].astype(COLUMNS_K[col])

        df_list.append(df)
        del content[:idx+1]
    return df_list


def knobs_kvalues_to_tfs(knob_file=const.get_extractor_knobs_filename()):
    """ Convert data from KNOBVALUES extraction to tfs file.

    Args:
        knob_file: path to online model extractor KNOBVALUE-Type output file.

    Returns:
        TfsDataFrame with trimvalues and additional time an beamprocess in headers.
    """
    with open(knob_file, "r") as f:
        content = f.readlines()

    df = tfs.TfsDataFrame(columns=[NAME_COLUMN[0], TRIM_COLUMN[0]])

    ls = content[0].split()
    df.headers[BEAMPROCESS_HEADER] = ls[ls.index(PRE_BEAMPROCESS_STRING) + 1]
    del content[0]

    while len(content) > 0:

        # get time
        ls = content[0].split()
        if TIME_HEADER not in df.headers:
            time_idx = ls.index(PRE_TIME_STRING)
            df.headers[TIME_HEADER] = convert_local_time_to_utc(
                " ".join(ls[time_idx + 1:time_idx+3])[:-1]
            )

        # get knob name
        knob = ls[ls.index(PRE_KNOB_STRING) + 1]

        if len(content) == 1:
            raise RuntimeError("Error in knob '{:s}':".format(knob) +
                               "No data found. See log for details.")

        # get trim
        match = re.match(RE_TRIM, content[1])
        if not match:
            raise RuntimeError("Error in knob '{:s}':".format(knob) +
                               "No trim found. See log for details.")
        trim = np.float64(match.group("trimvalue"))

        # add to dataframe
        df.loc[len(df), :] = (knob, trim)  # do not use append, it deletes headers

        for col, typ in (NAME_COLUMN, TRIM_COLUMN):
            df[col] = df[col].astype(typ)

        # clear lines
        idx = 0
        for idx, line in enumerate(content[1:]):
            if line.startswith("!"):
                break  # next knob, keep last line
        else:
            idx += 1  # no further knob, delete last line
        del content[:idx+1]
    return df


def get_optics(optics_file):
    """ Extracts Optics from optics file (-e switch in online model extractor). """
    with open(optics_file, "r") as f:
        content = f.readlines()
    for line in content:
        ls = line.split()
        if ls[1] == OPTIC_NAME:
            return ls[2]
    raise RuntimeError("No optics found in '{:s}'".format(optics_file))


def get_fill_from_orbitfile(orbit_file):
    """ Extracts the fill number from orbit file (-oe switch in online model extractor).

    Args:
        orbit_file: Path to orbit file.
    """
    return int(tfs.read_tfs(orbit_file).headers[FILLNUMBER])


# Time Functions ###############################################################


def get_utc_now():
    """ Returns current utc time as timeformat """
    return datetime.now(tz=pytz.utc).strftime(const.get_time_format())


def convert_utc_time_to_local(time):
    """ Converts UTC time in isoformat to local time """
    local = const.get_experiment_timezone()
    time = pytz.utc.localize(datetime.strptime(time, const.get_time_format()))
    return time.astimezone(local).strftime(const.get_time_format())


def convert_local_time_to_utc(time):
    """ Converts local time in isoformat to UTC time """
    local = const.get_experiment_timezone()
    time = local.localize(datetime.strptime(time, const.get_time_format()))
    return time.astimezone(pytz.utc).strftime(const.get_time_format())


def post_trim_extract(change_file):
    """ Extract trim value and move file to better filename """
    new_filename = os.path.join(os.path.dirname(change_file), const.get_change_filename_out())
    os.rename(change_file, new_filename)
    return _get_trim_from_changefile(new_filename), new_filename


# Private Functions ############################################################


def _get_trim_from_changefile(change_file):
    """ Extract the trim value from the change_file """
    pattern = re.compile(RE_TRIM)
    with open(change_file, "r") as f:
        for line in f.readlines():
            match = re.match(pattern, line)
            if match:
                return np.float64(match.group("trimvalue"))
    raise StandardError("Trim Value could not be extracted from '{:s}'".format(change_file))


