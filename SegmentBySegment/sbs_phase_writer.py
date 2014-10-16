import __init__  # @UnusedImport
import os
import SegmentBySegment
import sbs_beta_writer
import numpy as np
import math

from Utilities import tfs_file_writer


def write_phase(element_name, measured_hor_phase, measured_ver_phase, measured_hor_beta, measured_ver_beta,
                   input_model, model_propagation, model_cor, save_path):

    file_phase_x, file_phase_y = _get_phase_tfs_files(element_name, save_path)

    bpms_list = SegmentBySegment.intersect([measured_hor_phase, input_model, model_cor, model_propagation])

    _write_phase_for_plane(file_phase_x, element_name, "X", bpms_list, measured_hor_phase, measured_hor_beta, model_propagation, model_cor)

    _write_phase_for_plane(file_phase_y, element_name, "Y", bpms_list, measured_ver_phase, measured_ver_beta, model_propagation, model_cor)


def _write_phase_for_plane(file_phase, element_name, plane, bpms_list, measured_phase, measured_beta, model_propagation, model_cor):
    first_bpm = bpms_list[0][1]
    (beta_start, err_beta_start, alfa_start, err_alfa_start,
     _, _, _, _) = sbs_beta_writer._get_start_end_betas(bpms_list, measured_beta, plane)

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        model_index = model_propagation.indx[bpm_name]
        meas_index = measured_phase.indx[bpm_name]

        model_phase = (getattr(model_propagation, "MU" + plane)[model_index] -
                       getattr(model_propagation, "MU" + plane)[model_propagation.indx[first_bpm]]) % 1

        meas_phase = (getattr(measured_phase, "PHASE" + plane)[meas_index] -
                      getattr(measured_phase, "PHASE" + plane)[measured_phase.indx[first_bpm]]) % 1

        phase_difference = (meas_phase - model_phase) % 1
        if phase_difference > 0.5:
            phase_difference = phase_difference - 1

        cor_phase = -getattr(model_propagation, "MU" + plane)[model_index] + getattr(model_cor, "MU" + plane)[model_index]

        prop_error = _propagate_error_phase(err_beta_start, err_alfa_start, meas_phase, beta_start, alfa_start)

        std_error = getattr(measured_phase, "STDPH" + plane)[meas_index]

        file_phase.add_table_row([bpm_name, measured_phase.S[meas_index], meas_phase, phase_difference, std_error, prop_error, cor_phase, bpm_s])

    file_phase.write_to_file()


def _get_phase_tfs_files(element_name, save_path):
    file_phase_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsphasext_" + element_name + ".out"))
    file_phase_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsphaseyt_" + element_name + ".out"))

    file_phase_x.add_column_names(["NAME", "S", "PHASEX", "PHASEXT", "ERRORX", "PHSTDX", "PHASE_PLAY", "MODEL_S"])
    file_phase_x.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    file_phase_y.add_column_names(["NAME", "S", "PHASEY", "PHASEYT", "ERRORY", "PHSTDY", "PHASE_PLAY", "MODEL_S"])
    file_phase_y.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

    return file_phase_x, file_phase_y


def _propagate_error_phase(errb0, erra0, dphi, bet0, alf0):
    return math.sqrt((((1/2.*np.cos(4*np.pi*dphi)*alf0/bet0)-(1/2.*np.sin(4*np.pi*dphi)/bet0)-(1/2.*alf0/bet0))*errb0)**2+((-(1/2.*np.cos(4*np.pi*dphi))+(1/2.))*erra0)**2)/(2*np.pi)  # @IgnorePep8
