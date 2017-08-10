from __future__ import print_function
import sys
import os
import time
import logging

import numpy as np
from scipy.fftpack import fft as scipy_fft

sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    ".."
)))

import clean
from io_handlers import input_handler, output_handler
from harpy import harpy
from Utilities.twiss_to_tbt import generate
from Utilities import tfs_pandas
from Utilities.contexts import timeit 
from model import manager
from sdds_files import turn_by_turn_reader

LOGGER = logging.getLogger(__name__)
LOG_SUFFIX = ".log"


def run_all(main_input, clean_input, harpy_input):
    with timeit(lambda time: LOGGER.info("Total time for file: %s", time)):
        if (not main_input.write_raw and
                clean_input is None and harpy_input is None):
            LOGGER.error("No file has been choosen to be writen!")
            return
        _setup_file_log_handler(main_input)
        tbt_files = turn_by_turn_reader.read_tbt_file(main_input.file)
        for tbt_file in tbt_files:
            run_all_for_file(tbt_file, main_input, clean_input, harpy_input)


def run_all_for_file(tbt_file, main_input, clean_input, harpy_input):
    if main_input.write_raw:
        _write_raw_file(tbt_file, main_input)

    if clean_input is not None or harpy_input is not None:
        clean_writer = output_handler.CleanedAsciiWritter(main_input, tbt_file.date)
        for plane in ("x", "y"):
            bpm_names = np.array(getattr(tbt_file, "bpm_names_" + plane))
            bpm_data = getattr(tbt_file, "samples_matrix_" + plane)

            all_bad_bpms = []
            usv = None
            if clean_input is not None:
                with timeit(lambda time: LOGGER.debug("Time for filtering: %s", time)):
                    bpm_names, bpm_data, bad_bpms_clean = clean.clean(
                        bpm_names, bpm_data, clean_input, tbt_file.date,
                    )
                with timeit(lambda time: LOGGER.debug("Time for SVD clean: %s", time)):
                    bpm_names, bpm_data, bpm_res, bad_bpms_svd, usv = clean.svd_clean(
                        bpm_names, bpm_data, clean_input,
                    )
                all_bad_bpms.extend(bad_bpms_clean)
                all_bad_bpms.extend(bad_bpms_svd)
                setattr(clean_writer, "bpm_names_" + plane, bpm_names)
                setattr(clean_writer, "samples_matrix_" + plane, bpm_data)

            if plane == "x":
                computed_dpp = calc_dp_over_p(main_input, bpm_names, bpm_data)

            if harpy_input is not None:
                with timeit(lambda time: LOGGER.debug("Time for harmonic_analysis: %s", time)):
                    drive_results, bad_bpms_fft = harmonic_analysis(
                        bpm_names, bpm_data, usv,
                        plane, main_input, harpy_input,
                    )
                    all_bad_bpms.extend(bad_bpms_fft)
                    #TODO: Writing of harpy should be done in output_handler
                    drive_results.write_full_results()

            output_handler.write_bad_bpms(
                main_input.file,
                all_bad_bpms,
                main_input.outputdir, plane
            )

        if clean_input.write_clean:
            clean_writer.dpp = computed_dpp
            clean_writer.write()


def harmonic_analysis(bpm_names, bpm_data, usv,
                      plane, main_input, harpy_input):
    tunes = harpy_input.tunex, harpy_input.tuney, harpy_input.tunez
    if harpy_input.nattunex is None or harpy_input.nattuney is None:
        nattunes = None
    else:
        nattunes = harpy_input.nattunex, harpy_input.nattuney, harpy_input.nattunez
    allowed = _get_allowed_length(
        rang=[0, bpm_data.shape[1]]
    )[-1]
    bpm_data = bpm_data[:, :allowed]
    output_file = output_handler.get_outpath_with_suffix(
        main_input.file, main_input.outputdir, ".lin" + plane
    )
    if harpy_input.harpy_mode == "bpm":
        drivemat = harpy.init_from_matrix(
            bpm_names, bpm_data, tunes, plane.upper(),
            output_file, main_input.model, nattunes=nattunes,
            tolerance=harpy_input.tolerance,
            start_turn=0, end_turn=None, sequential=harpy_input.sequential,
        )
    elif harpy_input.harpy_mode == "svd" or harpy_input.harpy_mode == "fast":
        if usv is None:
            raise ValueError("Running harpy SVD mode but not svd clean was run."
                             " Set 'clean' flag to use SVD mode.")
        usv = (usv[0], usv[1], usv[2][:, :allowed])
        fast = harpy_input.harpy_mode == "fast"
        drivemat = harpy.init_from_svd(
            bpm_names, bpm_data, usv, tunes, plane.upper(),
            output_file, main_input.model, nattunes=nattunes,
            tolerance=harpy_input.tolerance,
            start_turn=0, end_turn=None, sequential=harpy_input.sequential, fast=fast
        )
    bpms_after_fft = []
    for i in range(len(drivemat.bpm_results)):
        bpms_after_fft.append(drivemat.bpm_results[i].name)
    bad_bpms_fft = []
    for bpm_name in bpm_names:
        if bpm_name not in bpms_after_fft:
            bad_bpms_fft.append(bpm_name + " Could not find the main resonance")
    
    return drivemat, bad_bpms_fft


def calc_dp_over_p(main_input, bpm_names, bpm_data):
    model_twiss = tfs_pandas.read_tfs(main_input.model)
    model_twiss.set_index("NAME", inplace=True)
    sequence = model_twiss.headers["SEQUENCE"].lower().replace("b1", "").replace("b2", "")
    accel_cls = manager.get_accel_class(sequence)
    arc_bpms_mask = accel_cls.get_arc_bpms_mask(bpm_names)
    arc_bpm_data = bpm_data[arc_bpms_mask]
    arc_bpm_names = bpm_names[arc_bpms_mask]
    dispersions = model_twiss.loc[arc_bpm_names, "DX"] * 1e3  # We need it in mm
    closed_orbits = np.mean(arc_bpm_data, axis=1)
    numer = np.sum(dispersions * closed_orbits)
    denom = np.sum(dispersions ** 2)
    if denom == 0.:
        raise ValueError("Cannot compute dpp probably no arc BPMs.")
    dp_over_p = numer / denom
    return dp_over_p


def _get_allowed_length(rang=[300, 10000], p2max=14, p3max=9, p5max=6):
    ind = np.indices((p2max, p3max, p5max))
    nums = (np.power(2, ind[0]) *
            np.power(3, ind[1]) *
            np.power(5, ind[2])).reshape(p2max * p3max * p5max)
    nums = nums[(nums > rang[0]) & (nums <= rang[1])]
    return np.sort(nums)


def _setup_file_log_handler(main_input):
    file_handler = logging.FileHandler(
        output_handler.get_outpath_with_suffix(
            main_input.file,
            main_input.outputdir,
            LOG_SUFFIX
        ),
        mode="w",
    )
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    if __name__ == "__main__":
        logging.getLogger("").addHandler(file_handler)
    else:
        LOGGER.addHandler(file_handler)


def _set_up_logger():
    main_logger = logging.getLogger("")
    main_logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    main_logger.addHandler(console_handler)


if __name__ == "__main__":
    _set_up_logger()
    run_all(*input_handler.parse_args())
