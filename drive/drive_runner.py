from __future__ import print_function
import __init__  # noqa
import os
import sys
from utils import iotools
import subprocess


DRIVING_TERMS_FILENAME = "DrivingTerms"
DRIVING_TERMS_TEMPL = """
%(SDDS_PATH)s %(START_TURN)i %(END_TURN)i
"""

DRIVE_INP_FILENAME = "Drive.inp"
DRIVE_INP_TEMPL = """
KICK=1
CASE(1[H], 0[V])=1
KPER(KICK PERCE.)=0.5
TUNE X=%(TUNE_X)s
TUNE Y=%(TUNE_Y)s
PICKUP START=0
PICKUP END=538
ISTUN=%(ISTUN)s
LABEL RUN (1[yes])=0
WINDOWa1=0.04
WINDOWa2=0.1
WINDOWb1=0.4
WINDOWb2=0.45
%(NAT_TUNE_TEXT)s
"""
NAT_TUNE_TMPL = """
NATURAL X=%(NAT_TUNE_X)s
NATURAL Y=%(NAT_TUNE_Y)s
"""


def run_drive(sdds_file_path, start_turn, end_turn, tune_x, tune_y,
              nat_tune_x=None, nat_tune_y=None, clean_up=False,
              stdout=None, stderr=None, tune_window=0.01):
    output_path = os.path.abspath(os.path.dirname(sdds_file_path))
    driving_terms_path = _generate_driving_terms(sdds_file_path, output_path,
                                                 start_turn, end_turn)
    drive_inp_path = _generate_drive_inp(output_path, tune_x, tune_y,
                                         nat_tune_x, nat_tune_y, tune_window)
    _run_drive_exe(output_path, stdout, stderr)
    if clean_up:
        iotools.delete_item(driving_terms_path)
        iotools.delete_item(drive_inp_path)


def _generate_driving_terms(sdds_file_path, output_path, start_turn, end_turn):
    driving_terms_path = os.path.join(output_path, DRIVING_TERMS_FILENAME)
    driving_terms_text = DRIVING_TERMS_TEMPL % {"SDDS_PATH": sdds_file_path,
                                                "START_TURN": start_turn,
                                                "END_TURN": end_turn}
    print("Writing DrivingTerms file into: ", driving_terms_path)
    iotools.write_string_into_new_file(driving_terms_path, driving_terms_text)
    return driving_terms_path


def _generate_drive_inp(output_path, tune_x, tune_y,
                        nat_tune_x, nat_tune_y, tune_window):
    drive_inp_path = os.path.join(output_path, DRIVE_INP_FILENAME)
    if nat_tune_x is not None and nat_tune_y is not None:
        nat_tune_text = NAT_TUNE_TMPL % {"NAT_TUNE_X": str(nat_tune_x),
                                         "NAT_TUNE_Y": str(nat_tune_y)}
    else:
        nat_tune_text = ""
    drive_inp_text = DRIVE_INP_TEMPL % {"TUNE_X": str(tune_x),
                                        "TUNE_Y": str(tune_y),
                                        "ISTUN": str(tune_window),
                                        "NAT_TUNE_TEXT": nat_tune_text}
    print("Writing Drive.inp file into: ", drive_inp_path)
    iotools.write_string_into_new_file(drive_inp_path, drive_inp_text)
    return drive_inp_path


def _run_drive_exe(output_path, stdout=None, stderr=None):
    drive_exe_path = _get_so_independent_drive_exe_path()
    process = subprocess.Popen([drive_exe_path, output_path],
                               stdin=subprocess.PIPE,
                               stdout=stdout, stderr=stderr)
    process.communicate()
    return_code = process.wait()
    return return_code


def _get_so_independent_drive_exe_path():
    this_dir_path = os.path.abspath(os.path.dirname(__file__))
    if sys.platform.startswith("linux"):
        return os.path.join(this_dir_path, "Drive_God_lin")
    elif sys.platform == "win32":  # Windows...
        return os.path.join(this_dir_path, "Drive_God_lin_win.exe")
    else:
        raise OSError("No drive version for given os: " + sys.platform)


if __name__ == '__main__':
    print >> sys.stderr, "Standalone usage not implemented yet."
    sys.exit(-1)
