from __future__ import print_function
import sys
import os
import logging
from collections import OrderedDict

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import clean
import harpy
from io_handlers import input_handler, output_handler

from Utilities import tfs_pandas as tfs
from Utilities.contexts import timeit 
from model import manager
from sdds_files import turn_by_turn_reader


LOGGER = logging.getLogger(__name__)
LOG_SUFFIX = ".log"


def run_all(main_input, clean_input, harpy_input, to_log):
    with timeit(lambda spanned: LOGGER.info("Total time for file: %s", spanned)):
        if (not main_input.write_raw and
                clean_input is None and harpy_input is None):
            LOGGER.error("No file has been choosen to be writen!")
            return
        _setup_file_log_handler(main_input)
        LOGGER.debug(to_log)
        tbt_files = turn_by_turn_reader.read_tbt_file(main_input.file)
        for tbt_file in tbt_files:
            run_all_for_file(tbt_file, main_input, clean_input, harpy_input)


def run_all_for_file(tbt_file, main_input, clean_input, harpy_input):

    if main_input.write_raw:
        output_handler.write_raw_file(tbt_file, main_input)

    if clean_input is not None or harpy_input is not None:
        clean_writer = output_handler.CleanedAsciiWritter(main_input, tbt_file.date)
        model_tfs = tfs.read_tfs(main_input.model).loc[:, ('NAME', 'S', 'DX')] # dispersion in meters
        for plane in ("x", "y"):
            bpm_data = getattr(tbt_file, "samples_matrix_" + plane)
            if (clean_input.startturn > 0 or
                    clean_input.endturn < bpm_data.shape[1]):
                start = max(0, clean_input.startturn)
                end = min(clean_input.endturn, new_bpm_data.shape[1])
                bpm_data = bpm_data[:, start:end]
            bpm_data, bpms_not_in_model = _get_only_model_bpms(bpm_data, model_tfs)
            all_bad_bpms = []
            usv = None
            if clean_input is not None:
                with timeit(lambda spanned: LOGGER.debug("Time for filtering: %s", spanned)):
                    bpm_data, bad_bpms_clean = clean.clean(
                        bpm_data, clean_input, tbt_file.date,
                    )
                with timeit(lambda spanned: LOGGER.debug("Time for SVD clean: %s", spanned)):
                    bpm_data, bpm_res, bad_bpms_svd, usv = clean.svd_clean(
                        bpm_data, clean_input,
                    )
                all_bad_bpms.extend(bpms_not_in_model)
                all_bad_bpms.extend(bad_bpms_clean)
                all_bad_bpms.extend(bad_bpms_svd)
                setattr(clean_writer, "samples_matrix_" + plane, bpm_data)

            if plane == "x":
                computed_dpp = calc_dp_over_p(main_input, bpm_data)

            if harpy_input is not None:
                lin_frame = pd.DataFrame(index=bpm_data.index,
                                         data={"NAME": bpm_data.index})
                bpm_positions = model_tfs.set_index("NAME").loc[bpm_data.index, "S"]
                lin_frame["S"] = bpm_positions
                with timeit(lambda spanned: LOGGER.debug("Time for orbit_analysis: %s", spanned)):
                    lin_frame = _get_orbit_data(lin_frame, bpm_data, bpm_res, model_tfs)

                bad_bpms_fft = _do_harpy(main_input, harpy_input, bpm_data,
                                         model_tfs, usv, lin_frame, plane,
                                         computed_dpp)
                all_bad_bpms.extend(bad_bpms_fft)

            output_handler.write_bad_bpms(
                main_input.file,
                all_bad_bpms,
                main_input.outputdir, plane
            )

        if clean_input.write_clean:
            clean_writer.dpp = computed_dpp
            clean_writer.write()


def _get_only_model_bpms(bpm_data, model):
    model_indx = model.set_index("NAME").index
    bpm_data_in_model = bpm_data.loc[model_indx.intersection(bpm_data.index)]
    not_in_model = bpm_data.index.difference(model_indx)
    bpms_not_in_model = []
    for bpm in not_in_model:
        bpms_not_in_model.append("{} not found in model".format(bpm))
    return bpm_data_in_model, bpms_not_in_model


def _get_orbit_data(lin_frame, bpm_data, bpm_res, model):
    lin_frame['PK2PK'] = np.max(bpm_data, axis=1) - np.min(bpm_data, axis=1)
    lin_frame['CO'] = np.mean(bpm_data, axis=1)
    lin_frame['CORMS'] = np.std(bpm_data, axis=1) / np.sqrt(bpm_data.shape[1])
    lin_frame['BPM_RES'] = bpm_res
    return lin_frame


def _do_harpy(main_input, harpy_input, bpm_data, model_tfs, usv, lin_frame, plane, computed_dpp):
    with timeit(lambda spanned: LOGGER.debug("Time for harmonic_analysis: %s", spanned)):
        bpm_data = (bpm_data.T - np.mean(bpm_data, axis=1)).T
        lin_result, spectrum, bad_bpms_fft = harmonic_analysis(
            bpm_data, usv,
            plane, harpy_input, lin_frame, model_tfs,
        )
        lin_result = _rescale_amps_to_main_line(lin_result, plane)
        lin_result.sort_values('S', axis=0, ascending=True, inplace=True, kind='mergesort')
        headers = _compute_headers(lin_result, plane, computed_dpp)
        output_handler.write_harpy_output(
            main_input,
            lin_result,
            headers,
            spectrum,
            plane
        )
    return bad_bpms_fft


def harmonic_analysis(bpm_data, usv, plane, harpy_input, panda, model_tfs):
    if usv is None:
        if harpy_input.harpy_mode == "svd" or harpy_input.harpy_mode == "fast":
            raise ValueError(
                "Running harpy SVD mode but not svd clean was run."
                " Set 'clean' flag to use SVD mode."
            )
    else:
        allowed = _get_allowed_length(rang=[0, bpm_data.shape[1]])[-1]
        bpm_data = bpm_data.loc[:, :allowed]
        usv = (usv[0], usv[1], usv[2][:, :allowed])
    lin_result, spectrum, bad_bpms_fft = harpy.harpy(
        bpm_data, usv, plane.upper(), harpy_input, panda, model_tfs)
    return lin_result, spectrum, bad_bpms_fft


def _rescale_amps_to_main_line(panda, plane):
    cols = [col for col in panda.columns.values if col.startswith('AMP')]
    cols.remove('AMP' + plane.upper())
    panda.loc[:, cols] = panda.loc[:, cols].div(
        panda.loc[:, 'AMP' + plane.upper()],
        axis="index"
    )
    return panda


def calc_dp_over_p(main_input, bpm_data):
    model_twiss = tfs.read_tfs(main_input.model)
    model_twiss.set_index("NAME", inplace=True)
    sequence = model_twiss.headers["SEQUENCE"].lower().replace("b1", "").replace("b2", "")
    accel_cls = manager.get_accel_class(sequence)
    arc_bpms_mask = accel_cls.get_arc_bpms_mask(bpm_data.index)
    arc_bpm_data = bpm_data[arc_bpms_mask]
    # We need it in mm:
    dispersions = model_twiss.loc[arc_bpm_data.index, "DX"] * 1e3
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


def _compute_headers(panda, plane, computed_dpp):
    plane_number = {"x": "1", "y": "2"}[plane]
    headers = OrderedDict()
    tunes = panda.loc[:, "TUNE" + plane.upper()]
    headers["Q" + plane_number] = np.mean(tunes)
    headers["Q" + plane_number + "RMS"] = np.std(tunes)
    headers["DPP"] = computed_dpp
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
