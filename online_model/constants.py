"""
Module online_model.constants
---------------------------------------
Functions to retrieve constants in online model extractor.
"""
import os
import pytz
import numpy as np


def get_om_extractor():
    """ AFS-path to online model extractor """
    return os.path.join("/", "afs", "cern.ch", "eng", "lhc_online_model", "dev", "bin",
                        "lhc-model-extractor.sh")


def get_extractor_output_filename():
    """ Default output filename of extractor """
    return "modifiers.madx"


def get_extractor_knobs_filename():
    """ Default knob filename (output of extractor) """
    return "knobs.madx"


def get_extractor_log_filename():
    """ Default log filename (output of extractor) """
    return "extractor.log"


def get_default_knob_summary_filename(idx):
    """ Default knob summary filename (output of extractor) """
    return "knobs.madxsummary{:d}".format(idx)


def get_default_orbit_filename(beam):
    """ Default orbit filename (output of extractor """
    return "LHCB{:d}.orbit1.tfs".format(beam)


def get_experiment_timezone():
    """ Get time zone for measurement data. """
    return pytz.timezone("Europe/Zurich")


def get_time_format():
    """ Default time format """
    return "%Y-%m-%d %H:%M:%S.%f"


def get_ktype(ktype):
    """ Maps from more understandable names to modelextractor names. """
    return {
        "k": "K",
        "trim": "KNOBVALUE",
    }[ktype]


def get_change_filename_out():
    """ Default changeparameters filename """
    return "changeparameters_fromknobs.madx"


def get_knob_tfs_filename(knob):
    """ Default changeparameters filename """
    return "{:s}.tfs".format(knob.replace("/", ".").replace("-", "_").lower())


def get_overview_filename():
    """ Default changeparameters filename """
    return "overview.tfs"


def get_optics_header():
    return "OPTICS"


def get_beamprocess_header():
    return "BEAMPROCESS"


def get_trim_header():
    return "TRIM"


def get_trimhint_header():
    return "TRIMHINT"


def get_time_header():
    return "UTCTIME"


def get_knob_header():
    return" KNOB"


def get_fill_header():
    return "FILL"


def get_name_column(with_type=False):
    name = "NAME"
    if with_type:
        return (name, str)
    return


def get_circuit_column(with_type=False):
    name = "CIRCUIT"
    if with_type:
        return (name, str)
    return


def get_value_column(with_type=False):
    name = "VALUE"
    if with_type:
        return (name, np.float64)
    return


def get_delta_column(with_type=False):
    name = "DELTA"
    if with_type:
        return (name, np.float64)
    return


def get_trim_column(with_type=False):
    name = "TRIM"
    if with_type:
        return (name, np.float64)
    return


