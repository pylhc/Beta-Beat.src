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
RE_LINE = r"\s*(?P<NAME>\S+)\s*=\s*(?P<VALUE>[^;]+)[^:]+:\s*(?P<CIRCUIT>[^/]+/[^.]+)\.[^\d+\-]+(?P<DELTA>[^.]+\.[^.]+)\."
COLUMNS = OrderedDict([("NAME", str), ("CIRCUIT", str),
                       ("VALUE", np.float64), ("DELTA", np.float64)])

RE_TRIM = r"\s*(?P<trimname>" + TRIM_VALUE_NAME + r")\s*=\s*(?P<trimvalue>[^;]+)\s*;"


def knob_k_to_tfs(knob_file=const.get_default_knob_filename(), trim=None):
    """ Convert from knob.madx output file to a proper tfs file.

    Assumes no one changes the output order. See global variables is this file.

    Args:
        knob_file: path to online model extractor K-Type output file.
        trim: trim value if available.
    """
    with open(knob_file, "r") as f:
        content = f.readlines()

    df = tfs.TfsDataFrame(index=range(len(content[2:])), columns=COLUMNS.keys())

    # get headers
    ls = content[0].split()
    df.headers["BEAMPROCESS"] = ls[ls.index(PRE_BEAMPROCESS_STRING) + 1]

    ls = content[1].split()
    df.headers["KNOB"] = ls[ls.index(PRE_KNOB_STRING) + 1]

    time_idx = ls.index(PRE_TIME_STRING)
    df.headers["UTCTIME"] = convert_local_time_to_utc(
        " ".join(ls[time_idx + 1:time_idx+3])[:-1]
    )

    if trim is not None:
        df.headers["TRIM"] = trim
        df.headers["TRIMHINT"] = "VALUE = BEFORE_VALUE + DELTA * TRIM"

    # get data
    pattern = re.compile(RE_LINE)
    for idx, line in enumerate(content[2:]):
        match = re.match(pattern, line)
        for col in COLUMNS.keys():
                df.loc[idx, col] = match.group(col)

    for col in COLUMNS.keys():
        with suppress_exception(KeyError):
            df[col] = df[col].astype(COLUMNS[col])

    return df


def convert_utc_time_to_local(time):
    """ Converts UTC time in isoformat to local time """
    utc = pytz.timezone("UTC")
    local = const.get_experiment_timezone()
    time = utc.localize(datetime.strptime(time, const.get_time_format()))
    return time.astimezone(local).strftime(const.get_time_format())


def convert_local_time_to_utc(time):
    """ Converts local time in isoformat to UTC time """
    utc = pytz.timezone("UTC")
    local = const.get_experiment_timezone()
    time = local.localize(datetime.strptime(time, const.get_time_format()))
    return time.astimezone(utc).strftime(const.get_time_format())


def post_trim_extract(change_file):
    """ Extract trim value and move file to better filename """
    new_filename = os.path.join(os.path.dirname(change_file), const.get_default_change_filename())
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
