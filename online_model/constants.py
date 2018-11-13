"""
Module online_model.constants
---------------------------------------
Functions to retrieve constants in online model extractor.
"""
import os
import pytz


def get_om_extractor():
    """ AFS-path to online model extractor """
    return os.path.join("/", "afs", "cern.ch", "eng", "lhc_online_model", "dev", "bin",
                        "lhc-model-extractor.dev_knobs.sh")


def get_default_knob_filename():
    """ Default knob filename (output of extractor) """
    return "knobs.madx"


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


def get_default_change_filename():
    """ Default changeparameters filename """
    return "changeparameters_fromknob.madx"