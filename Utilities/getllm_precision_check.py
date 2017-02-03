from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../")
))
from madx import madx_wrapper
from drive import drive_runner
from GetLLM import GetLLM
from Python_Classes4MAD import metaclass
import tempfile
import iotools
import ADDbpmerror
from contextlib import contextmanager

HOR, VER = 0, 1

MADX_SCRIPT = """
title, "Tracking test for GetLLM";

!@require lhc_runII_2016
!@require tracking

option, -echo;
exec, full_lhc_def("%(MODIFIERS)s", %(BEAM)s);
option, echo;

! exec, high_beta_matcher();

exec, match_tunes(64.%(IQX)s, 59.%(IQY)s, %(BEAM)s);

exec, do_twiss_monitors(
    LHCB%(BEAM)i,
    "%(TWISS)s",
    0.0
);
exec, do_twiss_elements(
    LHCB%(BEAM)i,
    "%(TWISS_ELEMENTS)s",
    0.0
);

exec, do_track_single_particle(
    %(KICK_X)s, %(KICK_Y)s,
    %(NUM_TURNS)s, "%(TRACK_PATH)s"
);

stop;

exec, twiss_ac_dipole(
    %(QX)s, %(QY)s, %(DQX)s, %(DQX)s,
    %(BEAM)s, "%(TWISS_AC)s", 0.0
);
"""

OPTICS = "40cm"

MODIFIERS = {
    "injection": "call file=\"runII_2016/opt_inj.madx\";",
    "40cm": "call file=\"runII_2016/opt_400_10000_400_3000.madx\";",
}

ERR_DEF_FILES = {
    "injection": "0450GeV",
    "40cm": "6500GeV",
}
ERR_DEF_PATH = os.path.join(
    "/afs", "cern.ch", "work", "o", "omc",
    "Error_definition_files", ERR_DEF_FILES[OPTICS]
)

BEAM = 1
NUM_TURNS = 2200
QX = 0.28
QY = 0.31
DQX = 0.27
DQY = 0.32
KICK_X = 0.00035
KICK_Y = 0.0002


def print_getllm_precision():
    output_dir = tempfile.mkdtemp()
    try:
        _run_tracking_model(output_dir)
        _do_analysis(output_dir)
        _comprare_results(output_dir)
    finally:
        pass
        # _clean_up_files(output_dir)


def _run_tracking_model(directory):
    print("Creating model and tracking...")
    madx_script = _get_madx_script(BEAM, directory)
    madx_wrapper.resolve_and_run_string(madx_script, log_file=os.devnull)
    track_path = _get_track_path(directory, one=True)
    tbt_path = _get_tbt_path(directory)
    with silence():
        ADDbpmerror.convert_files(infile=track_path, outfile=tbt_path)


def _get_madx_script(beam, directory):
    modifiers_file_path = _get_modifiers_file(OPTICS, directory)
    twiss_path = _get_twiss_path(directory)
    twiss_elem_path = _get_twiss_elem_path(directory)
    twiss_ac_path = _get_twiss_ac_path(directory)
    track_path = _get_track_path(directory)
    madx_script = MADX_SCRIPT % {
        "MODIFIERS": modifiers_file_path,
        "BEAM": beam,
        "IQX": str(QX).replace("0.", ""),
        "IQY": str(QY).replace("0.", ""),
        "QX": QX,
        "QY": QY,
        "DQX": DQX,
        "DQY": DQY,
        "TWISS": twiss_path,
        "TWISS_ELEMENTS": twiss_elem_path,
        "TWISS_AC": twiss_ac_path,
        "NUM_TURNS": NUM_TURNS,
        "TRACK_PATH": track_path,
        "KICK_X": KICK_X,
        "KICK_Y": KICK_Y,
    }
    return madx_script


def _get_twiss_path(directory):
    return os.path.join(directory, "twiss.dat")


def _get_twiss_elem_path(directory):
    return os.path.join(directory, "twiss_elements.dat")


# TODO: This into madx template
def _get_twiss_centre_path(directory):
    return os.path.join(directory, "twiss_elements_centre.dat")


def _get_twiss_ac_path(directory):
    return os.path.join(directory, "twiss_ac.dat")


def _get_track_path(directory, one=False):
    if not one:
        return os.path.join(directory, "track")
    return os.path.join(directory, "trackone")


def _get_tbt_path(directory):
    return os.path.join(directory, "ALLBPMs")


def _do_analysis(directory):
    print("Performing analysis...")
    tbt_path = _get_tbt_path(directory)
    print("    -> Running drive...")
    with silence():
        drive_runner.run_drive(tbt_path, 0, NUM_TURNS,
                               DQX, DQY, QX, QY,
                               stdout=open(os.devnull, "w"))
    twiss_path = os.path.join(directory, "twiss.dat")
    err_def_path = _copy_error_def_file(directory)
    print("    -> Running GetLLM...")
    with silence():
        GetLLM.main(directory, tbt_path, twiss_path, bpmu="mm",
                    errordefspath=err_def_path)


def _copy_error_def_file(directory):
    new_err_def_path = os.path.join(directory, "error_deff.txt")
    iotools.copy_item(ERR_DEF_PATH, new_err_def_path)
    return new_err_def_path


def _get_modifiers_file(optics, directory):
    modifiers_file_path = os.path.join(directory, "modifiers.madx")
    with open(modifiers_file_path, "w") as modifiers_file:
        modifiers_file.write(MODIFIERS[optics])
    return modifiers_file_path


def _comprare_results(directory):
    _print_free_results(directory)
    if os.path.isfile(_get_twiss_ac_path(directory)):
        _print_ac_results(directory)


def _print_free_results(directory):
    print("+++++++ Results for free files +++++++\n")
    model_data = metaclass.twiss(_get_twiss_path(directory))
    betax_data = _get_beta_data(directory, HOR)
    betay_data = _get_beta_data(directory, VER)
    _compare_betas(model_data, betax_data, betay_data)
    print("++++++++++++++++++++++++++++++++++++++\n")


def _print_ac_results(directory):
    print("+++++++ Results for AC free files +++++++\n")
    model_data = metaclass.twiss(_get_twiss_ac_path(directory))
    betax_data = _get_beta_data(directory, HOR, free=False)
    betay_data = _get_beta_data(directory, VER, free=False)
    _compare_betas(model_data, betax_data, betay_data)
    print("++++++++++++++++++++++++++++++++++++++\n")


PLANE_SUFFIX = {HOR: "x", VER: "y"}


def _get_beta_data(directory, plane, free=True):
    suffix = PLANE_SUFFIX[plane]
    getbeta = os.path.join(directory, "getbeta" + suffix + ".out")
    getbetafree = os.path.join(directory, "getbeta" + suffix + "_free.out")
    if not free:
        return metaclass.twiss(getbeta)
    return _get_twiss_for_one_of(getbetafree, getbeta)


def _get_twiss_for_one_of(*paths):
    for path in paths:
        if os.path.isfile(path):
            return metaclass.twiss(path)
    raise IOError("None of the files exist:\n\t" + "\n\t".join(path))


def _compare_betas(model, betax, betay):
    print("-> Beta-beating:")
    x_model_beta_values = []
    y_model_beta_values = []
    x_beta_values = []
    y_beta_values = []
    for bpm_name in model.NAME:
        if bpm_name not in betax.indx or bpm_name not in betay.indx:
            continue
        x_indx = betax.indx[bpm_name]
        x_beta_values.append((betax.BETX[x_indx], betax.ERRBETX[x_indx]))
        y_indx = betay.indx[bpm_name]
        y_beta_values.append((betay.BETY[y_indx], betay.ERRBETY[y_indx]))
        mdl_indx = model.indx[bpm_name]
        x_model_beta_values.append(model.BETX[mdl_indx])
        y_model_beta_values.append(model.BETY[mdl_indx])
    x_beta_values, x_beta_errs = zip(*x_beta_values)
    y_beta_values, x_beta_errs = zip(*y_beta_values)
    x_beta_values = np.array(x_beta_values)
    y_beta_values = np.array(y_beta_values)
    x_model_beta_values = np.array(x_model_beta_values)
    y_model_beta_values = np.array(y_model_beta_values)
    x_beta_beat = ((x_beta_values - x_model_beta_values) /
                   x_model_beta_values)
    y_beta_beat = ((y_beta_values - y_model_beta_values) /
                   y_model_beta_values)
    print("    Horizontal:")
    print("    -RMS beating:", np.std(x_beta_beat))
    print("    -Peak beating:", x_beta_beat[np.argmax(np.abs(x_beta_beat))])
    print("    Vertical:")
    print("    -RMS beating:", np.std(y_beta_beat))
    print("    -Peak beating:", y_beta_beat[np.argmax(np.abs(y_beta_beat))])
    print("")


def _beating(model, meas):
    return (meas - model) / model


def _clean_up_files(ouput_dir):
    iotools.delete_item(ouput_dir)


@contextmanager
def silence():
    stdout = sys.stdout
    stderr = sys.stderr
    devnull = open(os.devnull, "w")
    sys.stdout = devnull
    sys.stderr = devnull
    try:
        yield
    finally:
        sys.stdout = stdout
        sys.stderr = stderr
        devnull.close()


if __name__ == "__main__":
    print_getllm_precision()
