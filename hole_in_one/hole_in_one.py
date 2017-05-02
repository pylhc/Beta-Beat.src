from __future__ import print_function
import sys
import os
import numpy as np
import time
import logging
from datetime import datetime
import clean
from input_handler import parse_args
from harpy import harpy

sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    ".."
)))

print(sys.path)
from Utilities import tfs_pandas
from model import manager
from sdds_files import turn_by_turn_reader

RAW_SUFFIX = ".raw"
CLEAN_SUFFIX = ".clean"
LOG_SUFFIX = ".log"
LOGGER = logging.getLogger(__name__)


def run_all(main_input, clean_input, harpy_input):
    start_time = time.time()
    if (not main_input.write_raw and
            clean_input is None and harpy_input is None):
        LOGGER.error("No file has been choosen to be writen!")
        return
    if main_input.outputdir is None:
        outdir = os.path.dirname(main_input.file)
        main_input.outputdir = outdir
        LOGGER.info("Setting outdir to: " + outdir)
    _setup_file_log_handler(main_input)
    tbt_files = turn_by_turn_reader.read_tbt_file(main_input.file)
    for tbt_file in tbt_files:
        if main_input.write_raw:
            LOGGER.debug("Writing raw sdds")
            tbt_file.write_to_ascii(
                main_input.model,
                get_outpath_with_suffix(main_input.file,
                                        main_input.outputdir,
                                        RAW_SUFFIX)
            )

        if clean_input is not None:
            if not clean_input.write_clean and harpy_input is None:
                LOGGER.warn("Clean but set to run and not write_clean." +
                            "Will not produce output.")
            dict_to_write_clean_sdds = {}
            computed_dpp = None
            clean_input.noresync = (clean_input.noresync and
                                    not resync_from_date(tbt_file.date))
            for plane in ["x", "y"]:
                bpm_names = np.array(getattr(tbt_file, "bpm_names_" + plane))
                bpm_data = getattr(tbt_file, "samples_matrix_" + plane)
                bpm_names, bpm_data, bad_bpms = clean.clean(
                    bpm_names, bpm_data,
                    clean_input,
                )
                bpm_names, bpm_data, bpm_res, bad_bpms_svd = clean.svd_clean(
                    bpm_names, bpm_data,
                    clean_input,
                )

                if clean_input.write_clean:
                    if plane == "x":
                        computed_dpp = calc_dp_over_p(main_input, bpm_names, bpm_data)
                    dict_to_write_clean_sdds["BPMs_" + plane] = bpm_names
                    dict_to_write_clean_sdds["SAMPLES_" + plane] = bpm_data
                    LOGGER.debug("Writing clean sdds")
                bad_bpms_fft = []

                if harpy_input is not None:
                    allowed = _get_allowed_length(
                        rang=[0, bpm_data.shape[1]]
                    )[-1]
                    bpm_data = bpm_data[:, :allowed]
                    bad_bpms_fft = harmonic_analysis(
                        bpm_names, bpm_data, bpm_res,
                        plane, main_input, harpy_input,
                    )
                write_bad_bpms_into_file(
                    main_input.file,
                    bad_bpms + bad_bpms_svd + bad_bpms_fft,
                    main_input.outputdir, plane
                )
            if clean_input.write_clean:
                headers_dict = {}
                if computed_dpp is not None:
                    headers_dict["dpp"] = computed_dpp
                turn_by_turn_reader.write_ascii_file(
                    main_input.model,
                    get_outpath_with_suffix(main_input.file,
                                            main_input.outputdir,
                                            CLEAN_SUFFIX),
                    list(dict_to_write_clean_sdds["BPMs_x"]),
                    dict_to_write_clean_sdds["SAMPLES_x"],
                    list(dict_to_write_clean_sdds["BPMs_y"]),
                    dict_to_write_clean_sdds["SAMPLES_y"],
                    tbt_file.date,
                    headers_dict,
                )

    LOGGER.debug(">> Total time for file: {0}s".format(
        time.time() - start_time
    ))


def harmonic_analysis(bpm_names, bpm_data, bpm_res,
                      plane, main_input, harpy_input):
    time_start = time.time()
    tunes = harpy_input.tunex, harpy_input.tuney, harpy_input.tunez
    if harpy_input.nattunex is None or harpy_input.nattuney is None:
        nattunes = None
    else:
        nattunes = harpy_input.nattunex, harpy_input.nattuney
    output_file = get_outpath_with_suffix(main_input.file,
                                          main_input.outputdir,
                                          ".lin" + plane)
    drivemat = harpy.init_from_matrix(
        bpm_names, bpm_data, tunes, plane.upper(),
        output_file, main_input.model, nattunes=nattunes,
        tolerance=harpy_input.tolerance,
        start_turn=0, end_turn=None, sequential=harpy_input.sequential,
    )
    bpms_after_fft = []
    for i in range(len(drivemat.bpm_results)):
        drivemat.bpm_results[i].bpm_resolution = bpm_res[np.where(bpm_names == drivemat.bpm_results[i].name)[0]][0]
        bpms_after_fft.append(drivemat.bpm_results[i].name)
    bad_bpms_fft = []
    for bpm_name in bpm_names:
        if bpm_name not in bpms_after_fft:
            bad_bpms_fft.append(bpm_name + " Could not find the main resonance")
    drivemat.write_full_results()
    LOGGER.debug(">> Time for harmonic_analysis: {0}s".format(
        time.time() - time_start)
    )
    return bad_bpms_fft


def calc_dp_over_p(main_input, bpm_names, bpm_data):
    model_twiss = tfs_pandas.read_tfs(main_input.model)
    model_twiss.set_index("NAME", inplace=True)
    sequence = model_twiss.headers["SEQUENCE"].lower().replace("b1", "").replace("b2", "")
    AccelClass = manager.get_accel_class(sequence)
    arc_bpms_mask = AccelClass.get_arc_bpms_mask(bpm_names)
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


def resync_from_date(acqdate):
    LOGGER.debug("Will resynchronize BPMs")
    return acqdate > datetime(2016, 4, 1)


def write_bad_bpms_into_file(bin_path, bad_bpms_with_reasons, output_dir, plane):
    with open(get_outpath_with_suffix(bin_path, output_dir, ".bad_bpms_" + plane), 'w') as f:
        for line in bad_bpms_with_reasons:
            f.write(line + '\n')


def get_outpath_with_suffix(path, output_dir, suffix):
    return os.path.join(output_dir, os.path.basename(path) + suffix)


def _setup_file_log_handler(main_input):
    file_handler = logging.FileHandler(
        get_outpath_with_suffix(
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
    run_all(*parse_args())
