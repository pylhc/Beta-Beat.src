from __future__ import print_function
import sys
import os
import shutil
import numpy as np
import argparse
import time

from os.path import abspath, join, dirname, pardir
new_path = abspath(join(dirname(abspath(__file__)), pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

import madx_wrapper
from drive import drive_runner
from GetLLM import GetLLM
from Python_Classes4MAD import metaclass
from utils import iotools, ADDbpmerror, logging_tools
from utils.contexts import silence
import hole_in_one
from harmonic_analysis.io_handlers import input_handler as hio_input_handler
from model import manager

LOGGER = logging_tools.get_logger(__name__)

HOR, VER = 0, 1
PLANE_SUFFIX = {HOR: "x", VER: "y"}

THIS_DIR = os.path.dirname(__file__)
FILES_PATH = os.path.abspath(os.path.join(THIS_DIR,
                                          "getllm_precision_check"))
MADX_PATH = os.path.abspath(os.path.join(THIS_DIR, "..",
                                         "binaries", "madx_dev"))

FREE_BB_MAX = 0.25e-2  #  maximal bete beating error for free kick 0.25 %
FREE_BB_PEAK_MAX = 0.5e-2  #  maximal bete beating peak error for free kick 0.5 %
OK_FORMAT = "\33[38;2;0;255;0m"
NOTOK_FORMAT = "\33[38;2;255;0;0m"

MADX_SCRIPT = """
title, "Tracking test for GetLLM";

!@require lhc_runII_2016
!@require tracking

load_main_sequence(): macro = {
    call, file = "%(FILES_PATH)s/lhc_as-built.seq";
};

load_beam4_and_slice(): macro = {
    if (%(BEAM)i == 2){
        call, file = "%(FILES_PATH)s/lhcb4_as-built.seq";
    }
    call, file="%(FILES_PATH)s/slice.madx";
}


option, -echo;

exec, define_lhc_links();
exec, load_main_sequence();

beam, sequence=LHCB1, particle=proton, energy=6500,
    kbunch=1, npart=1.15E11, bv=1;
beam, sequence=LHCB2, particle=proton, energy=6500,
    kbunch=1, npart=1.15E11, bv=-1;

call, file = "%(MODIFIERS)s";
exec, cycle_sequences();
exec, set_default_crossing_scheme();


use, sequence=LHCB%(BEAM)i;
option, echo;

if(%(DO_COUPLING)s == 1){
    exec, coupling_knob(1);
    b1_re_ip7_knob = b1_re_ip7_knob - 0.01;
    b1_im_ip7_knob = b1_im_ip7_knob - 0.002;
};


exec, match_tunes(64.%(IQX)s, 59.%(IQY)s, %(BEAM)i);


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

if(%(DO_ACD)s == 1){
    exec, twiss_ac_dipole(
        %(QX)s, %(QY)s,
        %(DQX)s, %(DQY)s,
        %(BEAM)i, "%(TWISS_AC_OR_ADT)s", 0.0
    );

    ! Uninstall AC-dipole matrix...
    seqedit, sequence=LHCB%(BEAM)i;
        remove, element=hacmap;
        remove, element=vacmap;
    endedit;

     exec, load_beam4_and_slice();

    ! ...and install as element for tracking
    exec, install_acd_as_element(
        %(QX)s, %(QY)s,
        %(DQX)s, %(DQY)s,
        %(BEAM)i,
        %(RAMP1)s, %(RAMP2)s, %(RAMP3)s, %(RAMP4)s
    );
}elseif(%(DO_ADT)s == 1){
    exec, twiss_adt(
        %(QX)s, %(QY)s,
        %(DQX)s, %(DQY)s,
        %(BEAM)i, "%(TWISS_AC_OR_ADT)s", 0.0
    );

    ! Uninstall ADT matrix...
    seqedit, sequence=LHCB%(BEAM)i;
        remove, element=hacmap;
        remove, element=vacmap;
    endedit;

     exec, load_beam4_and_slice();

    ! ...and install as element for tracking
    exec, install_adt_as_element(
        %(QX)s, %(QY)s,
        %(DQX)s, %(DQY)s,
        %(BEAM)i,
        %(RAMP1)s, %(RAMP2)s, %(RAMP3)s, %(RAMP4)s
    );
}else {
    exec, load_beam4_and_slice();
    use, sequence=lhcb%(BEAM)i;
}


exec, do_madx_track_single_particle(
    %(KICK_X)s, %(KICK_Y)s, 0.0, 0.0,
    %(NUM_TURNS)s, "%(TRACK_PATH)s"
);
"""


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--acd",
        help="Activate excitation with ACD.",
        dest="acd",
        action="store_true",
    )
    parser.add_argument(
        "--adt",
        help="Activate excitation with ADT.",
        dest="adt",
        action="store_true",
    )
    parser.add_argument(
        "--optics",
        help="Optics can be 'injection' or '40cm'",
        dest="optics",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--analyse",
        help="Analyse method can be 'sussix', 'harpy_bpm', 'harpy_svd' or 'harpy_fast'",
        dest="analyse",
        required=True,
        choices=("sussix", "harpy_bpm", "harpy_svd", "harpy_fast"),
        type=str,
    )
    parser.add_argument(
        "--coupling",
        help="Include coupling in the simulation",
        dest="coupling",
        action="store_true",
    )
    parser.add_argument(
        "--usetracking",
        help="Path to previously created tracking model",
        dest="usetracking",
        default=None,
        type=str
    )
    parser.add_argument(
        "--nodelete",
        help="Delete output folder",
        dest="nodelete",
        action="store_true",
    )
    parser.add_argument(
        "--addnoise",
        help="Add noise",
        dest="addnoise",
        action="store_true",
    )
    parser.add_argument(
        "--nosilent",
        help="Don't suppress output",
        dest="nosilent",
        action="store_true",
    )
    return parser.parse_args()


MODIFIERS = {
    "injection": "call file=\"" + os.path.join(
        FILES_PATH,
        "opt_inj.madx"
    ) + "\";",
    "40cm": "call file=\"" + os.path.join(
        FILES_PATH,
        "opt_400_10000_400_3000.madx"
    ) + "\";",
}

ERR_DEF_FILES = {
    "injection": "0450GeV",
    "40cm": "6500GeV",
}
ERR_DEF_PATH = os.path.join(
    FILES_PATH,
    "Error_definition_files"
)

BEAM = 1
QX = 0.28
QY = 0.31
DQX = 0.27
DQY = 0.32
TUNE_WINDOW = 0.005

KICK_X = 1e-6
KICK_Y = 1e-6

RAMP1 = 0  # Ramp up start turn
RAMP2 = 2001  # Ramp up end turn
RAMP3 = 8561  # Ramp down start turn
RAMP4 = 10500  # Ramp down end turn

NUM_TURNS = RAMP3


def print_getllm_precision(options):
    if options.usetracking:
        output_dir = os.path.abspath(options.usetracking)
    else:
        output_dir = os.path.join(THIS_DIR, "test_getllm_" + time.strftime("%Y_%m_%d_%Hh_%Mm_%Ss"))
    _create_folders(output_dir, options)
    
    if options.acd and options.adt:
        raise IOError(
            "ADT and AC-dipole are both set to 1. Please select only one"
        )

    _run_tracking_model(output_dir, options)
    _do_analysis(output_dir, options)
    _comprare_results(output_dir, options)
    _clean_up_files(output_dir, options)
    LOGGER.info("Test finished")


def _run_tracking_model(directory, options):
    tbt_path = _get_tbt_path(directory)

    if options.usetracking:
        if not os.path.exists(tbt_path):
            raise IOError("Tracking data '" + options.usetracking + "' not found!")
        LOGGER.info("Using previously created model and tracking: " + options.usetracking)
        return
    print("Creating model and tracking...")
    madx_script = _get_madx_script(BEAM, directory, options)
    madx_wrapper.resolve_and_run_string(
        madx_script,
        madx_path=MADX_PATH,
        output_file=os.path.join(directory, 'job.test.madx'),
        log_file=os.path.join(directory, 'madx_log.txt')
    )
    track_path = _get_track_path(directory, one=True)
    with silence():
        if options.addnoise:
            ADDbpmerror.convert_files(infile=track_path, outfile=tbt_path, 
                                      x_error=1e-4, y_error=1e-4)
        else:
            ADDbpmerror.convert_files(infile=track_path, outfile=tbt_path)
        headers = [
            "\n".join(["#SDDSASCIIFORMAT v1",
                       "#Beam: LHCB" + str(BEAM),
                       "#Created: " +
                           time.strftime("%Y-%m-%d at %H:%M:%S") + " By: ADDbpmerror",
                       "#Number of turns: " + str(NUM_TURNS),
                       "#Number of horizontal monitors: 522",
                       "#Number of vertical monitors: 522",
                       "#Acquisition date: " + time.strftime("%Y-%m-%d at %H:%M:%S"),
                       "#dpp: 0.0",
                       ])
        ]
        lines = []
        with open(tbt_path, "r") as tbt_data:
            for line in tbt_data:
                if line.startswith("#"):
                    headers.append(line)
                else:
                    lines.append(line)
        with open(tbt_path, "w") as tbt_data:
            for header in headers:
                tbt_data.write(header)
            for line in reversed(lines):
                tbt_data.write(line)


def _get_madx_script(beam, directory, options):
    modifiers_file_path = _get_modifiers_file(options.optics, directory)
    twiss_path = _get_twiss_path(directory)
    twiss_elem_path = _get_twiss_elem_path(directory)
    kick_x = 0
    kick_y = 0
    do_acd = 0
    do_adt = 0
    if options.acd:
        twiss_ac_or_adt_path = _get_twiss_ac_path(directory)
        do_acd = 1
    elif options.adt:
        twiss_ac_or_adt_path = _get_twiss_adt_path(directory)
        do_adt = 1
    else:
        twiss_ac_or_adt_path = _get_twiss_path(directory)
        kick_x = KICK_X
        kick_y = KICK_Y
    do_coupling = 0
    if options.coupling:
        do_coupling = 1
    track_path = _get_track_path(directory)
    madx_script = MADX_SCRIPT % {
        "FILES_PATH": FILES_PATH,
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
        "TWISS_AC_OR_ADT": twiss_ac_or_adt_path,
        "NUM_TURNS": NUM_TURNS,
        "TRACK_PATH": track_path,
        "DO_ACD": do_acd,
        "DO_ADT": do_adt,
        "DO_COUPLING": do_coupling,
        "KICK_X": kick_x,
        "KICK_Y": kick_y,
        "RAMP1": RAMP1,
        "RAMP2": RAMP2,
        "RAMP3": RAMP3,
        "RAMP4": RAMP4,
    }
    return madx_script


def _do_analysis(directory, options):
    LOGGER.info("Performing analysis...")

    if options.analyse == 'sussix':
        LOGGER.info("    -> Running drive...")
        with silence():
            _run_drive(directory, options)
    else:
        LOGGER.info("    -> Running hole_in_one...")
        with silence():
            _run_harpy(directory, options)

    tbt_path = _get_tbt_path(directory)
    err_def_path = _copy_error_def_file(directory, options)
    #shutil.copy(err_def_path)

    LOGGER.info("    -> Running GetLLM...")
    # with silence():
    accel = manager.get_accel_instance(accel="lhc", lhc_mode="lhc_runII_2016", beam=BEAM, model_dir=directory)
    GetLLM.main(accel, directory, directory, tbt_path,
                bpmu="mm")
    LOGGER.info("")


def _run_drive(directory, options):
    tbt_path = _get_tbt_path(directory)

    if not options.acd and not options.adt:
        drive_runner.run_drive(tbt_path, RAMP2, RAMP3,
                               QX, QY, QX, QY,
                               stdout=open(os.devnull, "w"),
                               tune_window=TUNE_WINDOW)
    else:
        drive_runner.run_drive(tbt_path, RAMP2, RAMP3,
                               DQX, DQY, QX, QY,
                               stdout=open(os.devnull, "w"),
                               tune_window=TUNE_WINDOW)


def _run_harpy(directory, options):
    file_path = _get_tbt_path(directory)
    twiss_path = _get_twiss_path(directory)

    hio_args = [
        "--file", file_path,
        "--model", twiss_path,
        "--outputdir", directory,
        "clean",
        "--noresync",
        "harpy",
        "--startturn", str(RAMP2),
        "--endturn", str(RAMP3),
        "--tunex", str(DQX) if any([options.acd, options.adt]) else str(QX),
        "--tuney", str(DQY) if any([options.acd, options.adt]) else str(QY),
        "--nattunex", str(QX),
        "--nattuney", str(QY),
        "--tolerance", str(TUNE_WINDOW),
        "--harpy_mode", options.analyse.replace("harpy_", ""),
    ]
    hole_in_one.run_all(*hio_input_handler.parse_args(hio_args))

    shutil.move(file_path + '.linx', file_path + '_linx')
    shutil.move(file_path + '.liny', file_path + '_liny')


def _get_twiss_path(directory):
    return os.path.join(directory, "twiss.dat")


def _get_twiss_elem_path(directory):
    return os.path.join(directory, "twiss_elements.dat")


# TODO: This into madx template
def _get_twiss_centre_path(directory):
    return os.path.join(directory, "twiss_elements_centre.dat")


def _get_twiss_ac_path(directory):
    return os.path.join(directory, "twiss_ac.dat")


def _get_twiss_adt_path(directory):
    return os.path.join(directory, "twiss_adt.dat")


def _get_track_path(directory, one=False):
    if not one:
        return os.path.join(directory, "track")
    return os.path.join(directory, "trackone")


def _get_tbt_path(directory):
    return os.path.join(directory, "ALLBPMs")


def _get_harmonic_path(directory, plane):
    return os.path.join(directory, "ALLBPMs_lin" + PLANE_SUFFIX[plane])


def _get_method_path(directory, options):
    return os.path.join(directory, "results_" + options.analyse)


def _get_llm_results_path(directory, options):
    return os.path.join(_get_method_path(directory, options), "getllm_" + options.analyse)


def _move_nonmadx_files(directory, options):
    dest_dir = _get_method_path(directory, options)
    getllm_dir = _get_llm_results_path(directory, options)

    # cleanup results folder in case there is one already
    shutil.rmtree(os.path.join(dest_dir, 'BPM'), ignore_errors=True)

    # move getLLM output to subfolder
    for f in os.listdir(directory):
        if (not os.path.isdir(os.path.join(directory, f))) and f.startswith('get') and f.endswith('.out'):
            shutil.move(os.path.join(directory, f), os.path.join(getllm_dir, f))
            LOGGER.warning("move {} to {}".format(os.path.join(directory, f), os.path.join(getllm_dir, f)))

    # move all except for madx files into results folder
    if options.analyse == "sussix":
        for f in os.listdir(directory):
            if any([f.startswith(prefix) for prefix in ["BPM", "ALLBPMs_", "Driv",
                                                        "results.txt", "sussix"]]):
                shutil.move(os.path.join(directory, f), os.path.join(dest_dir, f)) 
    else:
        for f in os.listdir(directory):
            if any([f.startswith(prefix) for prefix in ["ALLBPMs.", "ALLBPMs_", "results.txt"]]):
                shutil.move(os.path.join(directory, f), os.path.join(dest_dir, f))


def _copy_error_def_file(directory, options):
    new_err_def_path = os.path.join(directory, "error_deffs.txt")
    LOGGER.warning("copying errordefs file to '{}'".format(new_err_def_path))
    iotools.copy_item(os.path.join(
        ERR_DEF_PATH, ERR_DEF_FILES[options.optics]
    ), new_err_def_path)
    return new_err_def_path


def _create_folders(directory, options):
    if not os.path.exists(directory):
        os.mkdir(directory)

    output_method = _get_method_path(directory, options)
    if not os.path.exists(output_method):
        os.mkdir(output_method)

    output_getllm = _get_llm_results_path(directory, options)
    if not os.path.exists(output_getllm):
        os.mkdir(output_getllm)


def _get_modifiers_file(optics, directory):
    modifiers_file_path = os.path.join(directory, "modifiers.madx")
    with open(modifiers_file_path, "w") as modifiers_file:
        modifiers_file.write(MODIFIERS[optics])
    return modifiers_file_path


def _comprare_results(directory, options):
    _print_results(directory, free=not any([options.acd, options.adt]))
    output_method = _get_method_path(directory, options)
    orig_stdout = sys.stdout
    with open(os.path.join(output_method, 'results.txt'), "w+") as f:
        sys.stdout = f
        _print_results(directory, free=not any([options.acd, options.adt]))
    sys.stdout = orig_stdout


def _print_results(directory, free=True):
    if free:
        LOGGER.info("+++++++ Results for free files +++++++\n")
    else:
        LOGGER.info("+++++++ Results for excited files +++++++\n")
    meas_tunex = _get_tune_data(directory, HOR, free=free)
    meas_tuney = _get_tune_data(directory, VER, free=free)
    if free:
        _compare_tunes(meas_tunex, meas_tuney, QX, QY)
    else:
        meas_nat_tunex = _get_tune_data(directory, HOR)
        meas_nat_tuney = _get_tune_data(directory, VER)
        _compare_tunes(meas_tunex, meas_tuney, DQX, DQY)
        _compare_nat_tunes(meas_nat_tunex, meas_nat_tuney, QX, QY)
    betax_data = _get_beta_data(directory, HOR, free=free)
    betay_data = _get_beta_data(directory, VER, free=free)
    _compare_betas(betax_data, betay_data)
    phasex_data = _get_phase_data(directory, HOR, free=free)
    phasey_data = _get_phase_data(directory, VER, free=free)
    _compare_phases(phasex_data, phasey_data)
    ampbetax_data = _get_beta_amp_data(directory, HOR, free=free)
    ampbetay_data = _get_beta_amp_data(directory, VER, free=free)
    _compare_amp_betas(ampbetax_data, ampbetay_data)

    _compare_coupling(directory, free=free)
    LOGGER.info("++++++++++++++++++++++++++++++++++++++\n")


def _get_tune_data(directory, plane, free=True):
    har_data = metaclass.twiss(_get_harmonic_path(directory, plane))
    plane_index = {HOR: "1", VER: "2"}
    if not free:
        return (getattr(har_data, "Q" + plane_index[plane]),
                getattr(har_data, "Q" + plane_index[plane] + "RMS"))
    return (getattr(har_data, "NATQ" + plane_index[plane]),
            getattr(har_data, "NATQ" + plane_index[plane] + "RMS"))


def _get_beta_data(directory, plane, free=True):
    suffix = PLANE_SUFFIX[plane]
    getbeta = os.path.join(directory, "getbeta" + suffix + ".out")
    getbetafree = os.path.join(directory, "getbeta" + suffix + "_free.out")
    if not free:
        return metaclass.twiss(getbeta)
    return _get_twiss_for_one_of(getbetafree, getbeta)


def _get_beta_amp_data(directory, plane, free=True):
    suffix = PLANE_SUFFIX[plane]
    getampbeta = os.path.join(directory, "getampbeta" + suffix + ".out")
    getampbetafree = os.path.join(directory, "getampbeta" + suffix + "_free.out")
    if not free:
        return metaclass.twiss(getampbeta)
    return _get_twiss_for_one_of(getampbetafree, getampbeta)


def _get_coupling_data(directory, free=True):
    getcouple = os.path.join(directory, 'getcouple.out')
    getcouplefree = os.path.join(directory, 'getcouple_free.out')
    if not free:
        return metaclass.twiss(getcouple)
    return _get_twiss_for_one_of(getcouplefree, getcouple)


def _get_coupling_twiss(directory):
    ctwiss = metaclass.twiss(os.path.join(directory, 'twiss.dat'))
    return ctwiss


def _compare_coupling(directory, free=True):
    diff_between_re = []
    diff_between_im = []
    error = []
    ctwiss = metaclass.twiss(_get_twiss_path(directory))
    ctwiss.Cmatrix()
    cdata = _get_coupling_data(directory, free)

    for i in range(0, len(cdata.F1001R)):
        for j in range(0, len(ctwiss.F1001R)):
            if cdata.NAME[i] in ctwiss.NAME[j]:
                diff_between_re.append((cdata.F1001R[i] - ctwiss.F1001R[j]))
                diff_between_im.append((cdata.F1001I[i] - ctwiss.F1001I[j]))

    for i in range(0, len(diff_between_re)):
        error.append(np.sqrt(diff_between_re[i]**2 + diff_between_im[i]**2))
    LOGGER.info("    Average difference of f1001: " + str(np.mean(error)))


def _get_phase_data(directory, plane, free=True):
    suffix = PLANE_SUFFIX[plane]
    getphase = os.path.join(directory,
                            "getphasetot" + suffix + ".out")
    getphasefree = os.path.join(directory,
                                "getphasetot" + suffix + "_free.out")
    if not free:
        return metaclass.twiss(getphase)
    return _get_twiss_for_one_of(getphasefree, getphase)


def _get_twiss_for_one_of(*paths):
    for path in paths:
        if os.path.isfile(path):
            return metaclass.twiss(path)
    raise IOError("None of the files exist:\n\t" + "\n\t".join(paths))


def _compare_tunes(meas_tunex, meas_tuney, model_qx, model_qy):
    LOGGER.info("-> Tunes:")
    value_tune_x, rms_tune_x = meas_tunex
    value_tune_y, rms_tune_y = meas_tuney
    LOGGER.info("    Horizontal:")
    LOGGER.info("    -Measured tune:" + str(value_tune_x) +  "+-" + str(rms_tune_x))
    LOGGER.info("    -Design tune:" + str(model_qx))
    LOGGER.info("    Vertical:")
    LOGGER.info("    -Measured tune:" + str(value_tune_y) + "+-" + str(rms_tune_y))
    LOGGER.info("    -Design tune:" + str(model_qy))
    LOGGER.info("")


def _compare_nat_tunes(meas_tunex, meas_tuney, model_qx, model_qy):
    LOGGER.info("-> Natural Tunes:")
    value_tune_x, rms_tune_x = meas_tunex
    value_tune_y, rms_tune_y = meas_tuney
    LOGGER.info("    Horizontal:")
    LOGGER.info("    -Measured natural tune:" + str(value_tune_x) + "+-" + str(rms_tune_x))
    LOGGER.info("    -Design natural tune:" + str(model_qx))
    LOGGER.info("    Vertical:")
    LOGGER.info("    -Measured natural tune:" + str(value_tune_y) + "+-" + str(rms_tune_y))
    LOGGER.info("    -Design natural tune:" + str( model_qy))
    LOGGER.info("")

def _print_bb(rms_x_bb, peak_x_bb):
    if rms_x_bb > FREE_BB_MAX:
        form = NOTOK_FORMAT
    else:
        form = OK_FORMAT
    LOGGER.info("    - RMS beating: {}{:.7f}\33[0m".format(form, rms_x_bb))
    if abs(peak_x_bb) > FREE_BB_PEAK_MAX:
        form = NOTOK_FORMAT
    else:
        form = OK_FORMAT
    LOGGER.info("    -Peak beating: {}{:.7f}\33[0m".format(form, peak_x_bb))

def _compare_betas(betax, betay):
    LOGGER.info("-> Beta from phase:")
    x_beta_beat = ((betax.BETX - betax.BETXMDL) /
                   betax.BETXMDL)
    y_beta_beat = ((betay.BETY - betay.BETYMDL) /
                   betay.BETYMDL)

    rms_x_bb = np.std(x_beta_beat)
    peak_x_bb = x_beta_beat[np.argmax(np.abs(x_beta_beat))]
    rms_y_bb = np.std(y_beta_beat)
    peak_y_bb = y_beta_beat[np.argmax(np.abs(y_beta_beat))]

    LOGGER.info("    Horizontal:")
    _print_bb(rms_x_bb, peak_x_bb)
    LOGGER.info("    Vertical:")
    _print_bb(rms_y_bb, peak_y_bb)
    LOGGER.info("")


def _compare_amp_betas(ampbetax, ampbetay):
    LOGGER.info("-> Beta from amplitude:")
    x_beta_beat = ((ampbetax.BETX - ampbetax.BETXMDL) /
                   ampbetax.BETXMDL)
    y_beta_beat = ((ampbetay.BETY - ampbetay.BETYMDL) /
                   ampbetay.BETYMDL)
    LOGGER.info("    Horizontal:")
    LOGGER.info("    -RMS beating:" + str(np.std(x_beta_beat)))
    LOGGER.info("    -Peak beating:" + str(x_beta_beat[np.argmax(np.abs(x_beta_beat))]))
    LOGGER.info("    Vertical:")
    LOGGER.info("    -RMS beating:" + str(np.std(y_beta_beat)))
    LOGGER.info("    -Peak beating:" + str(y_beta_beat[np.argmax(np.abs(y_beta_beat))]))
    LOGGER.info("")


def _compare_phases(phasex, phasey):
    LOGGER.info("-> Phase:")
    x_phase_err = phasex.PHASEX - phasex.PHXMDL
    y_phase_err = phasey.PHASEY - phasey.PHYMDL
    LOGGER.info("    Horizontal:")
    LOGGER.info("    -RMS error:" + str(np.std(x_phase_err)))
    LOGGER.info("    -Peak error:" + str(x_phase_err[np.argmax(np.abs(x_phase_err))]))
    LOGGER.info("    Vertical:")
    LOGGER.info("    -RMS error:" + str(np.std(y_phase_err)))
    LOGGER.info("    -Peak error:" + str(y_phase_err[np.argmax(np.abs(y_phase_err))]))
    LOGGER.info("")


def _clean_up_files(directory, options):
    LOGGER.info('Cleaning up...')
    if options.nodelete:
        _move_nonmadx_files(directory, options)
        LOGGER.info("Output folder:" + directory)
    else:
        iotools.delete_item(directory)


if __name__ == "__main__":
    _options = _parse_args()
    print_getllm_precision(_options)
