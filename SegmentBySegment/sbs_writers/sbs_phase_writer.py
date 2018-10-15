import math
import os

import numpy as np

from SegmentBySegment.sbs_writers import sbs_beta_writer
from SegmentBySegment.sbs_writers.sbs_beta_writer import intersect
from tfs_files import tfs_file_writer

FIRST_BPM_B1 = "BPMSW.1L2.B1"
FIRST_BPM_B2 = "BPMSW.1L8.B2"


def write_phase(element_name, measured_hor_phase, measured_ver_phase, measured_hor_beta, measured_ver_beta,
                    propagated_models, save_path):

    file_phase_x, file_phase_y = _get_phase_tfs_files(element_name, save_path)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected
    model_back_cor = propagated_models.corrected_back_propagation

    bpms_list_x = intersect([model_propagation, model_cor, model_back_propagation, model_back_cor, measured_hor_phase])
    bpms_list_y = intersect([model_propagation, model_cor, model_back_propagation, model_back_cor, measured_ver_phase])

    _write_phase_for_plane(file_phase_x, element_name, "X", bpms_list_x, measured_hor_phase, measured_hor_beta, model_propagation, model_cor, model_back_propagation, model_back_cor)

    _write_phase_for_plane(file_phase_y, element_name, "Y", bpms_list_y, measured_ver_phase, measured_ver_beta, model_propagation, model_cor, model_back_propagation, model_back_cor)


def _write_phase_for_plane(file_phase, element_name, plane, bpms_list, measured_phase, measured_beta, model_propagation, model_cor, model_back_propagation, model_back_cor):
    first_bpm = bpms_list[0][1]
    last_bpm = bpms_list[-1][1]
    (beta_start, err_beta_start, alfa_start, err_alfa_start,
     beta_end, err_beta_end, alfa_end, err_alfa_end) = sbs_beta_writer._get_start_end_betas(bpms_list, measured_beta, plane)

    # Fix for phase jump at the start of the ring
    tune, fix_start_s = {}, None
    first_bpm_on_ring = measured_phase.NAME[np.argmin(measured_phase.S)]
    if first_bpm_on_ring in model_propagation.NAME:
        fix_start_s = model_propagation.S[model_propagation.indx[first_bpm_on_ring]]
        tune["X"], tune["Y"] = measured_phase.Q1, measured_phase.Q2

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        model_s = measured_phase.S[measured_phase.indx[bpm_name]]

        meas_phase = (getattr(measured_phase, "PHASE" + plane)[measured_phase.indx[bpm_name]] -
                      getattr(measured_phase, "PHASE" + plane)[measured_phase.indx[first_bpm]]) % 1

        std_err_phase = getattr(measured_phase, "STDPH" + plane)[measured_phase.indx[bpm_name]]

        model_prop_phase = (getattr(model_propagation, "MU" + plane)[model_propagation.indx[bpm_name]] -
                            getattr(model_propagation, "MU" + plane)[model_propagation.indx[first_bpm]]) % 1

        model_cor_phase = (getattr(model_cor, "MU" + plane)[model_cor.indx[bpm_name]] -
                           getattr(model_cor, "MU" + plane)[model_cor.indx[first_bpm]]) % 1

        meas_phase_back = (getattr(measured_phase, "PHASE" + plane)[measured_phase.indx[bpm_name]] -
                           getattr(measured_phase, "PHASE" + plane)[measured_phase.indx[last_bpm]]) % 1

        model_back_propagation_phase = (getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[last_bpm]] -
                                        getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[bpm_name]]) % 1

        model_back_cor_phase = (getattr(model_back_cor, "MU" + plane)[model_back_cor.indx[last_bpm]] -
                                getattr(model_back_cor, "MU" + plane)[model_back_cor.indx[bpm_name]]) % 1

        prop_phase_difference = (meas_phase - model_prop_phase) % 1
        if prop_phase_difference > 0.5:
            prop_phase_difference = prop_phase_difference - 1

        back_prop_phase_difference = (meas_phase_back - model_back_propagation_phase) % 1
        if back_prop_phase_difference > 0.5:
            back_prop_phase_difference = back_prop_phase_difference - 1

        if fix_start_s is not None:
            if bpm_s >= fix_start_s:
                prop_phase_difference += tune[plane]
            elif bpm_s <= fix_start_s:
                back_prop_phase_difference -= tune[plane]

        prop_cor_phase = -getattr(model_propagation, "MU" + plane)[model_propagation.indx[bpm_name]] + getattr(model_cor, "MU" + plane)[model_cor.indx[bpm_name]]
        back_cor_phase = getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[bpm_name]] - getattr(model_back_cor, "MU" + plane)[model_back_cor.indx[bpm_name]]

        prop_phase_error = _propagate_error_phase(err_beta_start, err_alfa_start, model_prop_phase, beta_start, alfa_start)
        cor_phase_error = _propagate_error_phase(err_beta_start, err_alfa_start, model_cor_phase, beta_start, alfa_start)

        back_phase_error = _propagate_error_phase(err_beta_end, err_alfa_end, model_back_propagation_phase, beta_end, alfa_end)
        back_cor_phase_error = _propagate_error_phase(err_beta_end, err_alfa_end, model_back_cor_phase, beta_end, alfa_end)

        # Error for the phase difference -> sqrt(e1**2 + e2**2 - 2cov(e1, e2)) and assuming no covariance
        prop_meas_diff_error = np.sqrt(prop_phase_error ** 2 + std_err_phase ** 2)
        prop_cor_diff_error = np.sqrt(prop_meas_diff_error ** 2 + cor_phase_error ** 2)
        back_meas_diff_error = np.sqrt(back_phase_error ** 2 + std_err_phase ** 2)
        back_cor_diff_error = np.sqrt(back_meas_diff_error ** 2 + back_cor_phase_error ** 2)

        file_phase.add_table_row([bpm_name, bpm_s, meas_phase, std_err_phase, prop_phase_difference, prop_meas_diff_error, prop_cor_phase, prop_cor_diff_error, back_prop_phase_difference, back_meas_diff_error, back_cor_phase, back_cor_diff_error, model_s])

    file_phase.write_to_file()


def _get_phase_tfs_files(element_name, save_path):
    file_phase_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsphasext_" + element_name + ".out"))
    file_phase_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsphaseyt_" + element_name + ".out"))
    file_phase_x.add_string_descriptor("TYPE", "USER")  # Needed for sbs match
    file_phase_y.add_string_descriptor("TYPE", "USER")  # Needed for sbs match

    file_phase_x.add_column_names(["NAME", "S", "MEASPHASEX", "STDERRPHASEX", "PROPPHASEX", "ERRPROPPHASEX", "CORPHASEX", "ERRCORPHASEX", "BACKPHASEX", "ERRBACKPHASEX", "BACKCORPHASEX", "ERRBACKCORPHASEX", "MODEL_S"])
    file_phase_x.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    file_phase_y.add_column_names(["NAME", "S", "MEASPHASEY", "STDERRPHASEY", "PROPPHASEY", "ERRPROPPHASEY", "CORPHASEY", "ERRCORPHASEY", "BACKPHASEY", "ERRBACKPHASEY", "BACKCORPHASEY", "ERRBACKCORPHASEY", "MODEL_S"])
    file_phase_y.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    return file_phase_x, file_phase_y


def _propagate_error_phase(errb0, erra0, dphi, bet0, alf0):
    return math.sqrt((((1/2.*np.cos(4*np.pi*dphi)*alf0/bet0)-(1/2.*np.sin(4*np.pi*dphi)/bet0)-(1/2.*alf0/bet0))*errb0)**2+((-(1/2.*np.cos(4*np.pi*dphi))+(1/2.))*erra0)**2)/(2*np.pi)  # @IgnorePep8
