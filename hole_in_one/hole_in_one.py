from __future__ import print_function
import sys
import os
import logging
from collections import OrderedDict

import numpy as np
import pandas as pd
import clean
import svd_harpy
from io_handlers import input_handler, output_handler

sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    ".."
)))

from Utilities import tfs_pandas as tfs
from Utilities.contexts import timeit 
from model import manager
from sdds_files import turn_by_turn_reader


LOGGER = logging.getLogger(__name__)
LOG_SUFFIX = ".log"


def run_all(main_input, clean_input, harpy_input):
    with timeit(lambda spanned: LOGGER.info("Total time for file: %s", spanned)):
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
        output_handler.write_raw_file(tbt_file, main_input)

    if clean_input is not None or harpy_input is not None:
        clean_writer = output_handler.CleanedAsciiWritter(main_input, tbt_file.date)
        model_tfs = tfs.read_tfs(main_input.model).loc[:, ['NAME', 'S', 'DX']] # dispersion in meters
        for plane in ("x", "y"):
            bpm_names = np.array(getattr(tbt_file, "bpm_names_" + plane))
            bpm_data = getattr(tbt_file, "samples_matrix_" + plane)
            bpm_names, bpm_data, bpms_not_in_model = get_only_model_bpms(bpm_names, bpm_data, model_tfs)
            all_bad_bpms = []
            usv = None
            if clean_input is not None:
                with timeit(lambda spanned: LOGGER.debug("Time for filtering: %s", spanned)):
                    bpm_names, bpm_data, bad_bpms_clean = clean.clean(
                        bpm_names, bpm_data, clean_input, tbt_file.date,
                    )
                with timeit(lambda spanned: LOGGER.debug("Time for SVD clean: %s", spanned)):
                    bpm_names, bpm_data, bpm_res, bad_bpms_svd, usv = clean.svd_clean(
                        bpm_names, bpm_data, clean_input,
                    )
                all_bad_bpms.extend(bpms_not_in_model)
                all_bad_bpms.extend(bad_bpms_clean)
                all_bad_bpms.extend(bad_bpms_svd)
                setattr(clean_writer, "bpm_names_" + plane, bpm_names)
                setattr(clean_writer, "samples_matrix_" + plane, bpm_data)

            if plane == "x":
                computed_dpp = calc_dp_over_p(main_input, bpm_names, bpm_data)

            if harpy_input is not None:
                with timeit(lambda spanned: LOGGER.debug("Time for orbit_analysis: %s", spanned)):
                    lin_frame = get_orbit_data(bpm_names, bpm_data, bpm_res, model_tfs)
                    bpm_data = (bpm_data.T - np.mean(bpm_data, axis=1)).T
                with timeit(lambda spanned: LOGGER.debug("Time for harmonic_analysis: %s", spanned)):
                    lin_result, spectrum, bad_bpms_fft = harmonic_analysis(
                        bpm_names, bpm_data, usv,
                        plane, harpy_input, lin_frame, model_tfs,
                    )
                    lin_result = _rescale_amps_to_main_line(lin_result, plane)
                    all_bad_bpms.extend(bad_bpms_fft)
                    headers = _compute_headers(lin_result, plane)
                    output_handler.write_harpy_output(
                        main_input,
                        lin_result,
                        headers,
                        spectrum,
                        plane
                    )

            output_handler.write_bad_bpms(
                main_input.file,
                all_bad_bpms,
                main_input.outputdir, plane
            )

        if clean_input.write_clean:
            clean_writer.dpp = computed_dpp
            clean_writer.write()


def get_only_model_bpms(bpm_names, bpm_data, model):
    bpms = np.ones([len(bpm_names)], dtype=bool)
    no_model=set(bpm_names) - set(model.loc[:,'NAME'].values)
    bpms_not_in_model=[]
    for bpm in no_model:
        bpms_not_in_model.append(bpm + "Not found in model")
    for i in range(len(bpm_names)):
        if bpm_names[i] in no_model:
            bpms[i] = False
    return bpm_names[bpms], bpm_data[bpms,:], bpms_not_in_model


def get_orbit_data(bpm_names, bpm_data, bpm_res, model):
    di = {'NAME': bpm_names, 
          'PK2PK': np.max(bpm_data,axis=1)-np.min(bpm_data,axis=1),
          'CO': np.mean(bpm_data,axis=1),
          'CORMS': np.std(bpm_data,axis=1) / np.sqrt(bpm_data.shape[1]),
          'BPM_RES': bpm_res
         }
    return pd.merge(model, pd.DataFrame.from_dict(di), on='NAME', how='inner')


def harmonic_analysis(bpm_names, bpm_data, usv, plane, harpy_input, panda, model_tfs):
    if usv is None:
        if harpy_input.harpy_mode == "svd" or harpy_input.harpy_mode == "fast":
            raise ValueError("Running harpy SVD mode but not svd clean was run."
                             " Set 'clean' flag to use SVD mode.")
    else:
        allowed = _get_allowed_length(rang=[0, bpm_data.shape[1]])[-1]
        bpm_data = bpm_data[:, :allowed]
        usv = (usv[0], usv[1], usv[2][:, :allowed])
    lin_result, spectrum, bad_bpms_fft = svd_harpy.svd_harpy(
        bpm_names, bpm_data, usv, plane.upper(), harpy_input, panda, model_tfs) # TODO lin_file header?
    return lin_result, spectrum, bad_bpms_fft


def _rescale_amps_to_main_line(panda, plane):
    cols = [col for col in panda.columns.values if col.startswith('AMP')]
    cols.remove('AMP' + plane.upper())
    panda.loc[:, cols] = panda.loc[:, cols].div(
        panda.loc[:, 'AMP' + plane.upper()],
        axis="index"
    )
    return panda


def calc_dp_over_p(main_input, bpm_names, bpm_data):
    model_twiss = tfs.read_tfs(main_input.model)
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


def _compute_headers(panda, plane):
    plane_number = {"x": "1", "y": "2"}[plane]
    headers = OrderedDict()
    tunes = panda.loc[:, "TUNE" + plane.upper()]
    headers["Q" + plane_number] = np.mean(tunes)
    headers["Q" + plane_number + "RMS"] = np.std(tunes)
    try:
        nattunes = panda.loc[:, "NATTUNE" + plane.upper()]
        headers["NATQ" + plane_number] = np.mean(nattunes)
        headers["NATQ" + plane_number + "RMS"] = np.std(nattunes)
    except KeyError:
        pass  # No natural tunes
    try:
        ztunes = panda.loc[:, "TUNEZ" + plane.upper()]
        headers["Q3"] = np.mean(ztunes)
        headers["Q3RMS"] = np.std(ztunes)
    except KeyError:
        pass  # No tune z
    # TODO: DPPAMP
    return headers


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
