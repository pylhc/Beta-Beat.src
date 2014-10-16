import __init__  # @UnusedImport
import os
import SegmentBySegment
import numpy as np
import math

from math import sqrt
from Utilities import tfs_file_writer


def write_beta(element_name, is_element, measured_hor_beta, measured_ver_beta, input_model, propagated_models, save_path, filesum_b):
    file_alfa_x, file_beta_x, file_alfa_y, file_beta_y = _get_beta_tfs_files(element_name, save_path, is_element)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected

    if not is_element:
        bpms_list = SegmentBySegment.intersect([measured_hor_beta, input_model, model_cor, model_propagation, model_back_propagation])
    else:
        bpms_list = SegmentBySegment.intersect([input_model, model_cor, model_propagation, model_back_propagation])

    _write_beta_for_plane(file_alfa_x, file_beta_x, "X",
                          element_name, bpms_list, measured_hor_beta,
                          input_model, model_cor, model_propagation, model_back_propagation,
                          save_path, is_element)

    _write_beta_for_plane(file_alfa_y, file_beta_y, "Y",
                          element_name, bpms_list, measured_ver_beta,
                          input_model, model_cor, model_propagation, model_back_propagation,
                          save_path, is_element)


def _get_beta_tfs_files(element_name, save_path, is_element):
    file_beta_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetax_" + element_name + ".out"))
    file_alfa_x = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsalfax_" + element_name + ".out"))
    file_beta_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetay_" + element_name + ".out"))
    file_alfa_y = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsalfay_" + element_name + ".out"))

    if not is_element:
        file_beta_x.add_column_names(["NAME", "S", "BETPROPX", "ERRBETPROPX", "BETCORX", "ERRBETCORX", "BETXMDL", "MODEL_S"])
        file_beta_x.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_x.add_column_names(["NAME", "S", "ALFPROPX", "ERRALFPROPX", "ALFCORX", "ERRALFCORX", "ALFXMDL", "MODEL_S"])
        file_alfa_x.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])

        file_beta_y.add_column_names(["NAME", "S", "BETPROPY", "ERRBETPROPY", "BETCORY", "ERRBETCORY", "BETYMDL", "MODEL_S"])
        file_beta_y.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_y.add_column_names(["NAME", "S", "ALFPROPY", "ERRALFPROPY", "ALFCORY", "ERRALFCORY", "ALFYMDL", "MODEL_S"])
        file_alfa_y.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    else:
        file_beta_x.add_column_names(["NAME", "S", "BETPROPX", "ERRBETPROPX", "BETXMDL", "MODEL_S"])
        file_beta_x.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_x.add_column_names(["NAME", "S", "ALFPROPX", "ERRALFPROPX", "ALFXMDL", "MODEL_S"])
        file_alfa_x.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le"])

        file_beta_y.add_column_names(["NAME", "S", "BETPROPY", "ERRBETPROPY", "BETYMDL", "MODEL_S"])
        file_beta_y.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le"])
        file_alfa_y.add_column_names(["NAME", "S", "ALFPROPY", "ERRALFPROPY", "ALFYMDL", "MODEL_S"])
        file_alfa_y.add_column_datatypes(["%bpm_s", "%le", "%le", "%le", "%le", "%le"])

    return file_alfa_x, file_beta_x, file_alfa_y, file_beta_y


def _write_beta_for_plane(file_alfa, file_beta, plane, element_name, bpms_list, measured_beta, model, model_cor, model_propagation, model_back_propagation, output_path, is_element):
    first_bpm = bpms_list[0][1]

    beta_start = getattr(measured_beta, "BET" + plane)[measured_beta.indx[first_bpm]]
    alfa_start = getattr(measured_beta, "ALF" + plane)[measured_beta.indx[first_bpm]]

    err_beta_start = sqrt(getattr(measured_beta, "ERRBET" + plane)[measured_beta.indx[first_bpm]] ** 2 + getattr(measured_beta, "STDBET" + plane)[measured_beta.indx[first_bpm]] ** 2)
    err_alfa_start = sqrt(getattr(measured_beta, "ERRALF" + plane)[measured_beta.indx[first_bpm]] ** 2 + getattr(measured_beta, "STDALF" + plane)[measured_beta.indx[first_bpm]] ** 2)

    last_bpm = bpms_list[-1][1]

    beta_end = getattr(measured_beta, "BET" + plane)[measured_beta.indx[last_bpm]]
    alfa_end = -getattr(measured_beta, "ALF" + plane)[measured_beta.indx[last_bpm]]

    err_beta_end = sqrt(getattr(measured_beta, "ERRBET" + plane)[measured_beta.indx[last_bpm]] ** 2 + getattr(measured_beta, "STDBET" + plane)[measured_beta.indx[last_bpm]] ** 2)
    err_alfa_end = sqrt(getattr(measured_beta, "ERRALF" + plane)[measured_beta.indx[last_bpm]] ** 2 + getattr(measured_beta, "STDALF" + plane)[measured_beta.indx[last_bpm]] ** 2)

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        model_s = model.S[model.indx[bpm_name]]
        beta_model = getattr(model, "BET" + plane)[model.indx[bpm_name]]
        alfa_model = getattr(model, "ALF" + plane)[model.indx[bpm_name]]

        beta_propagation = getattr(model_propagation, "BET" + plane)[model_propagation.indx[bpm_name]]
        beta_back_propagation = getattr(model_back_propagation, "BET" + plane)[model_back_propagation.indx[bpm_name]]

        alfa_propagation = getattr(model_propagation, "ALF" + plane)[model_propagation.indx[bpm_name]]
        alfa_back_propagation = getattr(model_back_propagation, "ALF" + plane)[model_back_propagation.indx[bpm_name]]

        delta_phase = (getattr(model_propagation, "MU" + plane)[model_propagation.indx[bpm_name]]) % 1

        if not is_element:
            meas_beta = getattr(measured_beta, "BET" + plane)[measured_beta.indx[bpm_name]]
            err_beta_prop = _propagate_error_beta(err_beta_start, err_alfa_start, delta_phase, meas_beta, beta_start, alfa_start)

            meas_alfa = getattr(measured_beta, "ALF" + plane)[measured_beta.indx[bpm_name]]
            err_alfa_prop = _propagate_error_alfa(err_beta_start, err_alfa_start, delta_phase, meas_alfa, beta_start, alfa_start)

            beta_cor = getattr(model_cor, "BET" + plane)[model_cor.indx[bpm_name]]
            err_beta_cor = 0  # TODO: Propagate?

            alfa_cor = getattr(model_cor, "ALF" + plane)[model_cor.indx[bpm_name]]
            err_alfa_cor = 0  # TODO: Propagate?

            file_beta.add_table_row([bpm_name, bpm_s, meas_beta, err_beta_prop, beta_cor, err_beta_cor, beta_model, model_s])
            file_alfa.add_table_row([bpm_name, bpm_s, meas_alfa, err_alfa_prop, alfa_cor, err_alfa_cor, alfa_model, model_s])
        else:
            err_beta_prop = _propagate_error_beta(err_beta_start, err_alfa_start, delta_phase, beta_propagation, beta_start, alfa_start)
            err_alfa_prop = _propagate_error_alfa(err_beta_start, err_alfa_start, delta_phase, alfa_propagation, beta_start, alfa_start)

            delta_phase_back = (getattr(model_back_propagation, "MU" + plane)[model_back_propagation.indx[bpm_name]]) % 1
            err_beta_back = _propagate_error_beta(err_beta_end, err_alfa_end, delta_phase_back, beta_back_propagation, beta_end, alfa_end)
            err_alfa_back = _propagate_error_beta(err_beta_end, err_alfa_end, delta_phase_back, beta_back_propagation, beta_end, alfa_end)

            beta_f = (1 / err_beta_prop ** 2 * beta_propagation + 1 / err_beta_back ** 2 * beta_back_propagation) / (1 / err_beta_prop ** 2 + 1 / err_beta_back ** 2)
            err_beta_f = sqrt(1 / (1 / err_beta_prop ** 2 + 1 / err_beta_back ** 2))

            alfa_f = (1 / err_alfa_prop ** 2 * alfa_propagation + 1 / err_alfa_back ** 2 * alfa_back_propagation) / (1 / err_alfa_prop ** 2 + 1 / err_alfa_back ** 2)
            err_alfa_f = sqrt(1 / (1 / err_alfa_prop ** 2 + 1 / err_alfa_back ** 2))

            std_wght = sqrt(2 * (1 / err_beta_prop ** 2 * (beta_propagation - beta_f) ** 2 + 1 / err_beta_back ** 2 * (beta_back_propagation - beta_f) ** 2) / (1 / err_beta_prop ** 2 + 1 / err_beta_back ** 2))
            err_beta_f = sqrt(err_beta_f ** 2 + std_wght ** 2)

            std_wght = sqrt(2 * (1 / err_alfa_prop ** 2 * (alfa_propagation - alfa_f) ** 2 + 1 / err_alfa_back ** 2 * (alfa_back_propagation - alfa_f) ** 2) / (1 / err_alfa_prop ** 2 + 1 / err_alfa_back ** 2))
            err_alfa_f = sqrt(err_alfa_f ** 2 + std_wght ** 2)

            file_alfa.add_table_row([bpm_name, bpm_s, alfa_f, err_alfa_f, alfa_model, model_s])
            file_beta.add_table_row([bpm_name, bpm_s, beta_f, err_beta_f, beta_model, model_s])

    file_beta.write_to_file()
    file_alfa.write_to_file()


def _propagate_error_beta(errb0, erra0, dphi, bets, bet0, alf0):
    return math.sqrt((bets*np.sin(4*np.pi*dphi)*alf0/bet0 + bets*np.cos(4*np.pi*dphi)/bet0)**2*errb0**2 + (bets*np.sin(4*np.pi*dphi))**2*erra0**2)  # @IgnorePep8


def _propagate_error_alfa(errb0, erra0, dphi, alfs, bet0, alf0):
    return math.sqrt(((alfs*((np.sin(4*np.pi*dphi)*alf0/bet0) + (np.cos(4*np.pi*dphi)/bet0))) - (np.cos(4*np.pi*dphi)*alf0/bet0) + (np.sin(4*np.pi*dphi)/bet0))**2*errb0**2 + ((np.cos(4*np.pi*dphi)) - (alfs*np.sin(4*np.pi*dphi)))**2*erra0**2)  # @IgnorePep8
