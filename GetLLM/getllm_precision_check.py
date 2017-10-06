from __future__ import print_function
import sys
import os
import numpy as np
import argparse
import time

sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
))

from madx import madx_wrapper
from drive import drive_runner
from GetLLM import GetLLM
from Python_Classes4MAD import metaclass
from Utilities import iotools, ADDbpmerror
from Utilities.contexts import silence


HOR, VER = 0, 1
PLANE_SUFFIX = {HOR: "x", VER: "y"}

THIS_DIR = os.path.dirname(__file__)
FILES_PATH = os.path.abspath(os.path.join(THIS_DIR,
                                          "getllm_precision_check"))
MADX_PATH = os.path.abspath(os.path.join(THIS_DIR, "..",
                                         "binaries", "madx_dev"))

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
        "--coupling",
        help="Include coupling in the simulation",
        dest="coupling",
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
    "/afs", "cern.ch", "work", "o", "omc",
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

RAMP1 = 50  # Ramp up start turn
RAMP2 = 2050  # Ramp up end start turn
RAMP3 = 4250  # Ramp down start turn
RAMP4 = 6250  # Ramp down start turn

NUM_TURNS = RAMP3


def print_getllm_precision(options):
    output_dir = os.path.join(THIS_DIR,
                              "test_getllm" + time.strftime("%d_%m_%Y_%h_%s"))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    try:
        _run_tracking_model(output_dir, options)
        _do_analysis(output_dir, options)
        _comprare_results(output_dir)
    finally:
        _clean_up_files(output_dir)
        print("Test finished")
        print("Output directory:", os.path.abspath(output_dir))


def _run_tracking_model(directory, options):
    print("Creating model and tracking...")
    madx_script = _get_madx_script(BEAM, directory, options)
    with silence():
        madx_wrapper.resolve_and_run_string(
            madx_script,
            madx_path=MADX_PATH,
            output_file=os.path.join(directory, 'job.test.madx'),
            log_file=os.path.join(directory, 'madx_log.txt')
        )
    track_path = _get_track_path(directory, one=True)
    tbt_path = _get_tbt_path(directory)
    with silence():
        ADDbpmerror.convert_files(infile=track_path, outfile=tbt_path)
        headers = []
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
    if options.acd and options.adt:
        raise IOError(
            "ADT and AC-dipole are both set to 1. Please select only one"
        )
    elif options.acd:
        twiss_ac_or_adt_path = _get_twiss_ac_path(directory)
        do_acd = 1
    elif options.adt:
        twiss_ac_or_adt_path = _get_twiss_adt_path(directory)
        do_adt = 1
    if do_acd == 0 and do_adt == 0:
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


def _do_analysis(directory, options):
    print("Performing analysis...")
    tbt_path = _get_tbt_path(directory)
    print("    -> Running drive...")
    with silence():
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
    twiss_path = os.path.join(directory, "twiss.dat")
    err_def_path = _copy_error_def_file(directory, options)
    print("    -> Running GetLLM...")
    with silence():
        accel = "LHCB" + str(BEAM)
        GetLLM.main(directory, tbt_path, twiss_path, bpmu="mm", accel=accel,
                    errordefspath=err_def_path)


def _copy_error_def_file(directory, options):
    new_err_def_path = os.path.join(directory, "error_deff.txt")
    iotools.copy_item(os.path.join(
        ERR_DEF_PATH, ERR_DEF_FILES[options.optics]
    ), new_err_def_path)
    return new_err_def_path


def _get_modifiers_file(optics, directory):
    modifiers_file_path = os.path.join(directory, "modifiers.madx")
    with open(modifiers_file_path, "w") as modifiers_file:
        modifiers_file.write(MODIFIERS[optics])
    return modifiers_file_path


def _comprare_results(directory):
    _print_results(directory)
    if os.path.isfile(_get_twiss_ac_path(directory)):
        _print_results(directory, free=False)
    if os.path.isfile(_get_twiss_adt_path(directory)):
        _print_results(directory, free=False)


def _print_results(directory, free=True):
    if free:
        print("+++++++ Results for free files +++++++\n")
    else:
        print("+++++++ Results for excited files +++++++\n")
    meas_tunex = _get_tune_data(directory, HOR, free=free)
    meas_tuney = _get_tune_data(directory, VER, free=free)
    model_qx, model_qy = QX, QY
    if not free:
        model_qx, model_qy = DQX, DQY
    _compare_tunes(meas_tunex, meas_tuney, model_qx, model_qy)
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
    print("++++++++++++++++++++++++++++++++++++++\n")


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
    getcouplePathTwiss = os.path.join(directory, 'twiss.dat')
    ctwiss = metaclass.twiss(getcouplePathTwiss)
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
    print("    Average difference of f1001: ", np.mean(error))


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
    raise IOError("None of the files exist:\n\t" + "\n\t".join(path))


def _compare_tunes(meas_tunex, meas_tuney, model_qx, model_qy):
    print("-> Tunes:")
    value_tune_x, rms_tune_x = meas_tunex
    value_tune_y, rms_tune_y = meas_tuney
    print("    Horizontal:")
    print("    -Measured tune:", value_tune_x, "+-", rms_tune_x)
    print("    -Design tune:", model_qx)
    print("    Vertical:")
    print("    -Measured tune:", value_tune_y, "+-", rms_tune_y)
    print("    -Design tune:", model_qy)
    print("")


def _compare_betas(betax, betay):
    print("-> Beta from phase:")
    x_beta_beat = ((betax.BETX - betax.BETXMDL) /
                   betax.BETXMDL)
    y_beta_beat = ((betay.BETY - betay.BETYMDL) /
                   betay.BETYMDL)
    print("    Horizontal:")
    print("    -RMS beating:", np.std(x_beta_beat))
    print("    -Peak beating:", x_beta_beat[np.argmax(np.abs(x_beta_beat))])
    print("    Vertical:")
    print("    -RMS beating:", np.std(y_beta_beat))
    print("    -Peak beating:", y_beta_beat[np.argmax(np.abs(y_beta_beat))])
    print("")


def _compare_amp_betas(ampbetax, ampbetay):
    print("-> Beta from amplitude:")
    x_beta_beat = ((ampbetax.BETX - ampbetax.BETXMDL) /
                   ampbetax.BETXMDL)
    y_beta_beat = ((ampbetay.BETY - ampbetay.BETYMDL) /
                   ampbetay.BETYMDL)
    print("    Horizontal:")
    print("    -RMS beating:", np.std(x_beta_beat))
    print("    -Peak beating:", x_beta_beat[np.argmax(np.abs(x_beta_beat))])
    print("    Vertical:")
    print("    -RMS beating:", np.std(y_beta_beat))
    print("    -Peak beating:", y_beta_beat[np.argmax(np.abs(y_beta_beat))])
    print("")


def _compare_phases(phasex, phasey):
    print("-> Phase:")
    x_phase_err = phasex.PHASEX - phasex.PHXMDL
    y_phase_err = phasey.PHASEY - phasey.PHYMDL
    print("    Horizontal:")
    print("    -RMS error:", np.std(x_phase_err))
    print("    -Peak error:", x_phase_err[np.argmax(np.abs(x_phase_err))])
    print("    Vertical:")
    print("    -RMS error:", np.std(y_phase_err))
    print("    -Peak error:", y_phase_err[np.argmax(np.abs(y_phase_err))])
    print("")


def _clean_up_files(ouput_dir):
    print('Cleaning up...')
    iotools.delete_item(ouput_dir)


if __name__ == "__main__":
    _options = _parse_args()
    print_getllm_precision(_options)
