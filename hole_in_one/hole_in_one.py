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

from utils import tfs_pandas as tfs
from utils.contexts import timeit
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
    tbt_file = _cut_tbt_file(tbt_file,
                             main_input.startturn,
                             main_input.endturn)
    file_date = tbt_file.date
       
    bpm_datas = {"x": tbt_file.samples_matrix_x,
                 "y": tbt_file.samples_matrix_y}

    if main_input.write_raw:
        output_handler.write_raw_file(tbt_file, main_input)
    elif clean_input is None and harpy_input is None:
        raise ValueError("No clean or harpy or --write_raw...")
    # Defaults in case we dont run clean
    usvs = {"x": None, "y": None}
    all_bad_bpms = {"x": [], "y": []}
    bpm_ress = {"x": None, "y": None}
    dpp = 0.0

    model_tfs = tfs.read_tfs(main_input.model).loc[:, ('NAME', 'S', 'DX')]

    if clean_input is not None:
        usvs, all_bad_bpms, bpm_ress, dpp = _do_clean(
            main_input, clean_input,
            bpm_datas, file_date, model_tfs
        )

    if harpy_input is not None:
        all_bad_bpms = _do_harpy(main_input, harpy_input, bpm_datas, usvs,
                                 model_tfs, bpm_ress, dpp, all_bad_bpms)

    for plane in ("x", "y"):
        output_handler.write_bad_bpms(
            main_input.file,
            all_bad_bpms[plane],
            main_input.outputdir,
            plane
        )


def _do_clean(main_input, clean_input, bpm_datas, file_date, model_tfs):
    usvs, all_bad_bpms, bpm_ress = {}, {}, {}
    clean_writer = output_handler.CleanedAsciiWritter(main_input, file_date)
    for plane in ("x", "y"):
        bpm_data = bpm_datas[plane]
        bpm_data, bpms_not_in_model = _get_only_model_bpms(bpm_data, model_tfs)
        if bpm_data.empty:
            raise AssertionError("Check BPMs names! None of the BPMs was found in the model!")    
    
        bad_bpms = []
        usv = None
        
        with timeit(lambda spanned: LOGGER.debug("Time for filtering: %s", spanned)):
            bpm_data, bad_bpms_clean = clean.clean(
                bpm_data, clean_input, file_date,
            )
        with timeit(lambda spanned: LOGGER.debug("Time for SVD clean: %s", spanned)):
            bpm_data, bpm_res, bad_bpms_svd, usv = clean.svd_clean(
                bpm_data, clean_input,
            )
        bpm_ress[plane] = bpm_res
        bad_bpms.extend(bpms_not_in_model)
        bad_bpms.extend(bad_bpms_clean)
        bad_bpms.extend(bad_bpms_svd)
        all_bad_bpms[plane] = bad_bpms
        setattr(clean_writer, "samples_matrix_" + plane, bpm_data)

        if plane == "x":
            dpp = _calc_dp_over_p(main_input, bpm_data)
        bpm_datas[plane] = bpm_data
        usvs[plane] = usv

    if clean_input.write_clean:
        clean_writer.dpp = dpp
        clean_writer.write()
    return usvs, all_bad_bpms, bpm_ress, dpp


def _do_harpy(main_input, harpy_input, bpm_datas, usvs, model_tfs, bpm_ress, dpp, all_bad_bpms):
    lin_frames = {}
    for plane in ("x", "y"):
        bpm_data, usv = bpm_datas[plane], usvs[plane]
        lin_frames[plane] = pd.DataFrame(
            index=bpm_data.index,
            data=OrderedDict([
                ("NAME", bpm_data.index),
                ("S", model_tfs.set_index("NAME").loc[bpm_data.index, "S"])
            ])
        )
        with timeit(lambda spanned: LOGGER.debug("Time for orbit_analysis: %s", spanned)):
            lin_frames[plane] = _get_orbit_data(lin_frames[plane], bpm_data, bpm_ress[plane])
        bpm_datas[plane], usvs[plane] = _prepare_data_for_harpy(bpm_data, usv)

    with timeit(lambda spanned: LOGGER.debug("Time for harmonic_analysis: %s", spanned)):
        harpy_iterator = harpy.harpy(
            harpy_input,
            bpm_datas["x"], usvs["x"],
            bpm_datas["y"], usvs["y"],
        )

    for plane in ("x", "y"):
        harpy_results, spectr, bad_bpms_summaries = harpy_iterator.next()
        lin_frame = lin_frames[plane]
        lin_frame = lin_frame.loc[harpy_results.index].join(harpy_results)
        if harpy_input.is_free_kick:
            bpm_data = bpm_datas[plane]
            lin_frame = _kick_phase_correction(bpm_data, lin_frame, plane)
        _sync_phase(bpm_data, lin_frame, plane)    
        lin_frame = _rescale_amps_to_main_line(lin_frame, plane)
        lin_frame = _add_resonances_noise(lin_frame, plane, bpm_ress[plane])
        lin_frame = lin_frame.sort_values('S', axis=0, ascending=True)
        headers = _compute_headers(lin_frame, plane, dpp)
        output_handler.write_harpy_output(
            main_input,
            lin_frame,
            headers,
            spectr,
            plane
        )
        all_bad_bpms[plane].extend(bad_bpms_summaries)
    return all_bad_bpms


def _prepare_data_for_harpy(bpm_data, usv):
    bpm_data = bpm_data.subtract(bpm_data.mean(axis=1), axis=0)
    allowed = _get_allowed_length(rang=[0, bpm_data.shape[1]])[-1]
    bpm_data = bpm_data.loc[:, :allowed - 1]
    if usv is not None:
        usv = (usv[0], usv[1], usv[2][:, :allowed - 1])
    return bpm_data, usv


def _cut_tbt_file(tbt_file, start_turn, end_turn):
    bpm_data_x = tbt_file.samples_matrix_x
    bpm_data_y = tbt_file.samples_matrix_y
    start = max(0, start_turn)
    end = min(end_turn, bpm_data_x.shape[1])
    num_turns = end - start
    tbt_file.samples_matrix_x = pd.DataFrame(
        index=bpm_data_x.index,
        data=bpm_data_x.iloc[:, start:end].values
    )
    tbt_file.samples_matrix_y = pd.DataFrame(
        index=bpm_data_y.index,
        data=bpm_data_y.iloc[:, start:end].values
    )
    tbt_file.num_turns = num_turns
    return tbt_file


def _get_only_model_bpms(bpm_data, model):
    model_indx = model.set_index("NAME").index
    bpm_data_in_model = bpm_data.loc[model_indx.intersection(bpm_data.index)]
    not_in_model = bpm_data.index.difference(model_indx)
    bpms_not_in_model = []
    for bpm in not_in_model:
        bpms_not_in_model.append("{} not found in model".format(bpm))
    return bpm_data_in_model, bpms_not_in_model


def _get_orbit_data(lin_frame, bpm_data, bpm_res):
    lin_frame['PK2PK'] = np.max(bpm_data, axis=1) - np.min(bpm_data, axis=1)
    lin_frame['CO'] = np.mean(bpm_data, axis=1)
    lin_frame['CORMS'] = np.std(bpm_data, axis=1) / np.sqrt(bpm_data.shape[1])
    if bpm_res is not None:
        lin_frame['BPM_RES'] = bpm_res.loc[lin_frame.index]
        # TODO: Magic number 10?:
        lin_frame['NOISE'] = bpm_res / np.sqrt(bpm_data.shape[1]) / 10.0
        lin_frame['AVG_NOISE'] = lin_frame['NOISE']
    else:
        lin_frame['NOISE'] = 0.0
    return lin_frame

def _sync_phase(bpm_data_orig, lin_frame, plane):
    """ Produces MUXSYNC and MUYSYNC column that is MUX/Y but 
        shifted such that for bpm at index 0 is always 0.
        It allows to compare phases of consecutive measurements
        and if some measurements stick out remove them from the data set.
        author: skowron  
        """
    uplane = plane.upper()
    phase = np.copy(lin_frame.loc[:, 'MU' + uplane].values)
     
    phase0 = phase[0]
    for i in range(len(phase)):
        p = phase[i]
        p = p - phase0
        tmp = p
        while p < -0.5:
            p = p + 1
        while p >  0.5:
            p = p - 1
        #print('sync %d  in %f shifted %f final %f' % (i, phase[i],tmp,p))    
        phase[i] = p
    lin_frame['MU' + uplane + 'SYNC'] = phase


def _kick_phase_correction(bpm_data_orig, lin_frame, plane):
    uplane = plane.upper()
    bpm_data = bpm_data_orig.loc[lin_frame.index,:]
    damp, dstd = _get_damping(bpm_data)
    LOGGER.debug("Damping factor X: {0:2.2e} +- {1:2.2e}".format(damp, dstd))
    int_range = np.arange(0.0, bpm_data.shape[1])
    amp = lin_frame.loc[:, 'PK2PK'].values / 2
    tune = lin_frame.loc[:, 'TUNE' + uplane].values * 2 * np.pi
    phase = lin_frame.loc[:, 'MU' + uplane].values * 2 * np.pi
    damp_range = damp * int_range
    phase_range = np.outer(tune, int_range) + np.outer(phase, np.ones(bpm_data.shape[1]))

    e1 = np.sum(np.exp(2 * damp_range) * np.sin(2 * phase_range), axis=1) * amp/2
    e2 = np.sum(bpm_data * np.exp(damp_range) * np.sin(phase_range), axis=1)
    e3 = np.sum(bpm_data * np.exp(damp_range) * np.cos(phase_range), axis=1)
    e4 = np.sum(np.exp(2 * damp_range) * np.cos(2 * phase_range), axis=1) * amp
    cor = (e1 - e2) / ((e3 - e4) * 2 * np.pi)
    lin_frame['MU' + uplane] = lin_frame.loc[:, 'MU' + uplane].values + cor
    return lin_frame


def _get_damping(bpm_data):
    coefs = np.polyfit(np.arange(bpm_data.shape[1]),
                       np.maximum.accumulate(np.log(np.abs(bpm_data[::-1]))).T,
                       1)
    return np.mean(coefs[0, :]), np.std(coefs[0, :])


def _rescale_amps_to_main_line(panda, plane):
    cols = [col for col in panda.columns.values if col.startswith('AMP')]
    cols.remove('AMP' + plane.upper())
    panda.loc[:, cols] = panda.loc[:, cols].div(
        panda.loc[:, 'AMP' + plane.upper()],
        axis="index"
    )
    return panda


def _add_resonances_noise(lin_frame, plane, bpm_res):
    uplane = plane.upper()
    cols = [col for col in lin_frame.columns.values if col.startswith('AMP')]
    cols.remove('AMP' + uplane)
    noise_scaled = lin_frame.loc[:, 'NOISE'] / lin_frame.loc[:, 'AMP' + uplane]
    lin_frame.loc[:, "NOISE_SCALED"] = noise_scaled
    if np.max(noise_scaled) == 0.0:  # Do not calculated errors when no noise was calculated (all zeros)
        return lin_frame
    lin_frame.loc[:, "ERR_AMP" + uplane] = lin_frame.loc[:, 'NOISE']
    lin_frame.loc[:, "ERR_MU" + uplane] = _get_spectral_phase_error(
        lin_frame.loc[:, "AMP" + uplane],
        lin_frame.loc[:, "NOISE"],
    )
    if "NATTUNE" + uplane in lin_frame.columns:
        lin_frame.loc[:, "ERR_NATAMP" + uplane] = lin_frame.loc[:, 'NOISE']
        lin_frame.loc[:, "ERR_NATMU" + uplane] = _get_spectral_phase_error(
            lin_frame.loc[:, "NATAMP" + uplane],
            lin_frame.loc[:, "NOISE"],
        )
    for col in cols:
        this_amp = lin_frame.loc[:, col]
        lin_frame.loc[:, "ERR_" + col] = noise_scaled * np.sqrt(1 + np.square(this_amp))
        err_phase_col = "ERR_" + col.replace("AMP", "PHASE")
        if col == "AMPZ":
            err_phase_col = "ERR_MUZ"
        lin_frame.loc[:, err_phase_col] = _get_spectral_phase_error(this_amp, noise_scaled)
    return lin_frame


def _get_spectral_phase_error(amplitude, noise):
    """
    When the error is too big (> 2*pi*0.25 more or less) the noise is not
    Gaussian anymore and it close to an uniform distribution, so we set the
    error to be 0.3 as it is the standard deviation of uniformly
    distributed phases. This is an approximation that should keep the error
    of the error below 20%.
    """
    error = noise / (np.where(amplitude > 0.0, amplitude, 1e-15) * 2 * np.pi)
    return np.where(error > 0.25, 0.3, error)


def _calc_dp_over_p(main_input, bpm_data):
    model_twiss = tfs.read_tfs(main_input.model)
    model_twiss.set_index("NAME", inplace=True)
    sequence = model_twiss.headers["SEQUENCE"].lower().replace("b1", "").replace("b2", "")
    if sequence != "lhc":
        return 0.0  # TODO: What do we do with other accels.
    accel_cls = manager.get_accel_class(accel=sequence)
    arc_bpms_mask = accel_cls.get_element_types_mask(bpm_data.index, types=["arc_bpm"])
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
