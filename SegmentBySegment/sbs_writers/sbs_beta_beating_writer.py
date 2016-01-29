import __init__  # @UnusedImport
import os
import sys

import numpy as np
from math import sqrt
from Utilities import tfs_file_writer


def write_beta_beat(element_name,
                    measured_hor_beta_phase, measured_ver_beta_phase,
                    measured_hor_beta_amp, measured_ver_beta_amp,
                    measured_kmod_hor, measured_kmod_ver,
                    propagated_models, save_path):
    file_beta_phase_x, file_beta_phase_y = _get_beta_beat_from_phase_tfs_files(element_name, save_path)
    file_beta_amp_x, file_beta_amp_y = _get_beta_beat_from_amp_tfs_files(element_name, save_path)
    file_avg_beta_beat_x, file_avg_beta_beat_y = _get_q_avg_beta_beat_tfs_files(element_name, save_path)

    model_propagation = propagated_models.propagation
    model_back_propagation = propagated_models.back_propagation
    model_cor = propagated_models.corrected
    model_back_cor = propagated_models.corrected_back_propagation

    phase_bpms_list = intersect([model_cor, model_propagation,
                                 model_back_propagation, model_back_cor,
                                 measured_hor_beta_phase, measured_ver_beta_phase])

    _write_beta_beat_from_phase_for_plane(
        file_beta_phase_x, "X", phase_bpms_list, measured_hor_beta_phase,
        model_propagation, model_cor, model_back_propagation, model_back_cor,
        save_path)
    _write_beta_beat_from_phase_for_plane(
        file_beta_phase_y, "Y", phase_bpms_list, measured_ver_beta_phase,
        model_propagation, model_cor, model_back_propagation, model_back_cor,
        save_path)

    if measured_hor_beta_amp is not None and measured_ver_beta_amp is not None:
        amp_bpms_list = intersect([model_cor, model_propagation,
                                   model_back_propagation, model_back_cor,
                                   measured_hor_beta_amp, measured_ver_beta_amp])
        _write_beta_beat_from_amp_for_plane(
            file_beta_amp_x, "X", amp_bpms_list, measured_hor_beta_amp,
            model_propagation, model_cor, model_back_propagation, model_back_cor,
            save_path)
        _write_beta_beat_from_amp_for_plane(
            file_beta_amp_y, "Y", amp_bpms_list, measured_ver_beta_amp,
            model_propagation, model_cor, model_back_propagation, model_back_cor,
            save_path)

    if measured_kmod_hor is not None and measured_kmod_ver is not None:
        kmod_bpms_list = intersect([model_cor, model_propagation,
                                    model_back_propagation, model_back_cor,
                                    measured_kmod_hor, measured_kmod_ver])
        _write_q_avg_beta_beat_for_plane(
            file_avg_beta_beat_x, "X", kmod_bpms_list, measured_kmod_hor,
            model_propagation, model_cor, model_back_propagation, model_back_cor,
            save_path
        )
        _write_q_avg_beta_beat_for_plane(
            file_avg_beta_beat_y, "Y", kmod_bpms_list, measured_kmod_ver,
            model_propagation, model_cor, model_back_propagation, model_back_cor,
            save_path
        )


def _get_beta_beat_from_phase_tfs_files(element_name, save_path):
    file_beta_beat = {}
    file_beta_beat["X"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetabeatx_" + element_name + ".out"))
    file_beta_beat["Y"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsbetabeaty_" + element_name + ".out"))

    for plane in ["X", "Y"]:
        column_names = [
            "NAME", "S",
            "BETABEATPHASE" + plane, "ERRBETABEATPHASE" + plane,
            "BETABEATCOR" + plane,
            "BETABEATPHASEBACK" + plane, "ERRBETABEATPHASEBACK" + plane,
            "BETABEATCORBACK" + plane,
            "BET" + plane + "MDL", "MODEL_S"
        ]
        file_beta_beat[plane].add_column_names(column_names)
        file_beta_beat[plane].add_column_datatypes(["%bpm_s"] + ["%le"] * (len(column_names) - 1))

    return file_beta_beat["X"], file_beta_beat["Y"]


def _get_beta_beat_from_amp_tfs_files(element_name, save_path):
    file_beta_beat = {}
    file_beta_beat["X"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsampbetabeatx_" + element_name + ".out"))
    file_beta_beat["Y"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsampbetabeaty_" + element_name + ".out"))

    for plane in ["X", "Y"]:
        column_names = [
            "NAME", "S",
            "BETABEATAMP" + plane, "ERRBETABEATAMP" + plane,
            "BETABEATAMPBACK" + plane, "ERRBETABEATAMPBACK" + plane,
            "BET" + plane + "MDL", "MODEL_S"
        ]
        file_beta_beat[plane].add_column_names(column_names)
        file_beta_beat[plane].add_column_datatypes(["%bpm_s"] + ["%le"] * (len(column_names) - 1))

    return file_beta_beat["X"], file_beta_beat["Y"]


def _get_q_avg_beta_beat_tfs_files(element_name, save_path):
    file_avg_beta_beat = {}
    file_avg_beta_beat["X"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsqavgbetabeatx_" + element_name + ".out"))
    file_avg_beta_beat["Y"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbsqavgbetabeaty_" + element_name + ".out"))

    for plane in ["X", "Y"]:
        column_names = [
            "NAME", "S",
            "AVGBETABEAT" + plane, "ERRAVGBETABEAT" + plane,
            "AVGBETABEATCOR" + plane,
            "AVGBETABEATBACK" + plane, "ERRAVGBETABEATBACK" + plane,
            "AVGBETABEATBACKCOR" + plane
        ]
        file_avg_beta_beat[plane].add_column_names(column_names)
        file_avg_beta_beat[plane].add_column_datatypes(["%bpm_s"] + ["%le"] * (len(column_names) - 1))

    return file_avg_beta_beat["X"], file_avg_beta_beat["Y"]


def _write_beta_beat_from_phase_for_plane(file_beta_beat, plane, bpms_list,
                                          measured_beta_phase, model_propagation, model_cor,
                                          model_back_propagation, model_back_cor,
                                          output_path):

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        beta_propagation = _get_from_twiss(model_propagation, bpm_name, "BET", plane)
        beta_back_propagation = _get_from_twiss(model_back_propagation, bpm_name, "BET", plane)

        # Beta from phase beating (front)
        beta_phase = _get_from_twiss(measured_beta_phase, bpm_name, "BET", plane)
        beta_beat_phase = (beta_phase - beta_propagation) / beta_propagation
        syst_err_beta_phase = _get_from_twiss(measured_beta_phase, bpm_name, "ERRBET", plane)
        rand_err_beta_phase = _get_from_twiss(measured_beta_phase, bpm_name, "STDBET", plane)
        err_beta_beat_phase = sqrt(syst_err_beta_phase ** 2 + rand_err_beta_phase ** 2) / beta_propagation

        # Beta from corrected model beating (front)
        beta_cor = _get_from_twiss(model_cor, bpm_name, "BET", plane)
        beta_beat_cor = (beta_cor - beta_propagation) / beta_propagation

        # Beta from phase beating (back)
        beta_beat_phase_back = (beta_phase - beta_back_propagation) / beta_back_propagation
        syst_err_beta_phase_back = _get_from_twiss(measured_beta_phase, bpm_name, "ERRBET", plane)
        rand_err_beta_phase_back = _get_from_twiss(measured_beta_phase, bpm_name, "STDBET", plane)
        err_beta_beat_phase_back = sqrt(syst_err_beta_phase_back ** 2 + rand_err_beta_phase_back ** 2) / beta_back_propagation

        # Beta from corrected model beating (back)
        beta_back_cor = _get_from_twiss(model_back_cor, bpm_name, "BET", plane)
        beta_beat_back_cor = (beta_back_cor - beta_back_propagation) / beta_back_propagation

        model_s = measured_beta_phase.S[measured_beta_phase.indx[bpm_name]]
        beta_model = _get_from_twiss(measured_beta_phase, bpm_name, "BET", plane, "MDL")
        file_beta_beat.add_table_row([
            bpm_name, bpm_s,
            beta_beat_phase, err_beta_beat_phase,
            beta_beat_cor,
            beta_beat_phase_back, err_beta_beat_phase_back,
            beta_beat_back_cor,
            beta_model, model_s
        ])
    file_beta_beat.write_to_file()


def _write_beta_beat_from_amp_for_plane(file_beta_beat, plane, bpms_list,
                                        measured_beta_amp, model_propagation, model_cor,
                                        model_back_propagation, model_back_cor,
                                        output_path):

    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        beta_propagation = _get_from_twiss(model_propagation, bpm_name, "BET", plane)
        beta_back_propagation = _get_from_twiss(model_back_propagation, bpm_name, "BET", plane)

        # Beta from amplitude beating (front)
        beta_amp = _get_from_twiss(measured_beta_amp, bpm_name, "BET", plane)
        beta_beat_amp = (beta_amp - beta_propagation) / beta_propagation
        err_beta_beat_amp = _get_from_twiss(measured_beta_amp, bpm_name, "BET", plane, "STD") / beta_propagation

        # Beta from amplitude beating (back)
        beta_amp_back = _get_from_twiss(measured_beta_amp, bpm_name, "BET", plane)
        beta_beat_amp_back = (beta_amp_back - beta_back_propagation) / beta_back_propagation
        err_beta_beat_amp_back = _get_from_twiss(measured_beta_amp, bpm_name, "BET", plane, "STD") / beta_back_propagation

        model_s = measured_beta_amp.S[measured_beta_amp.indx[bpm_name]]
        beta_model = _get_from_twiss(measured_beta_amp, bpm_name, "BET", plane, "MDL")
        file_beta_beat.add_table_row([
            bpm_name, bpm_s,
            beta_beat_amp, err_beta_beat_amp,
            beta_beat_amp_back, err_beta_beat_amp_back,
            beta_model, model_s
        ])
    file_beta_beat.write_to_file()


def _write_q_avg_beta_beat_for_plane(file_avg_beta_beat, plane, element_list,
                                      measured_kmod, model_propagation, model_cor,
                                      model_back_propagation, model_back_cor,
                                      output_path):
    for element in element_list:
        element_name = element[1]
        if "MQ" in element_name.upper():  # Only compute average beta for quads
            element_s = element[0]
            element_length = _get_from_twiss(model_propagation, element, "L", "")

            measured_avg_beta = _get_from_twiss(measured_kmod, element, "AVGBET", plane)
            err_measured_avg_beta = _get_from_twiss(measured_kmod, element, "ERRAVGBET", plane)

            avg_beta_prop = _get_average_beta_for_model_element(model_propagation, element, plane, element_length)
            avg_beta_beating_prop = (measured_avg_beta - avg_beta_prop) / avg_beta_prop
            err_avg_beta_beating_prop = err_measured_avg_beta / avg_beta_prop

            avg_beta_cor = _get_average_beta_for_model_element(model_cor, element, plane, element_length)
            avg_beta_beating_cor = (avg_beta_cor - avg_beta_prop) / avg_beta_prop

            avg_beta_back = _get_average_beta_for_model_element(model_back_propagation, element, plane, element_length)
            avg_beta_beating_back = (measured_avg_beta - avg_beta_back) / avg_beta_back
            err_avg_beta_beating_back = err_measured_avg_beta / avg_beta_back

            avg_beta_back_cor = _get_average_beta_for_model_element(model_back_cor, element, plane, element_length)
            avg_beta_beating_back_cor = (avg_beta_back_cor - avg_beta_back) / avg_beta_back

            file_avg_beta_beat.add_table_row([
                element_name, element_s,
                avg_beta_beating_prop, err_avg_beta_beating_prop,
                avg_beta_beating_cor,
                avg_beta_beating_back, err_avg_beta_beating_back,
                avg_beta_beating_back_cor
            ])
    file_avg_beta_beat.write_to_file()


def _get_average_beta_for_model_element(model, element, plane, length):
    beta = _get_from_twiss(model, element, "BET", plane)
    alpha = _get_from_twiss(model, element, "ALF", plane)
    k1l = _get_from_twiss(model, element, "K1L", "")
    avg_beta, _ = _average_beta_in_quad(beta, alpha, length, k1l, 0., 0.)
    return avg_beta


def _get_from_twiss(twiss_data, element, column, plane, suffix=""):
    return getattr(twiss_data, column + plane + suffix)[twiss_data.indx[element]]


def intersect(list_of_files):
    '''Pure intersection of all bpm names in all files '''
    if len(list_of_files) == 0:
        print "Nothing to intersect!!!!"
        sys.exit()
    z = list_of_files[0].NAME
    for b in list_of_files:
        z = filter(lambda x: x in z, b.NAME)
    #SORT by S
    result = []
    x0 = list_of_files[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))
    result.sort()
    return result


def _average_beta_in_quad(bet1, alf1, L, K1L, beterr, alferr):
    K = np.sqrt(K1L / L)
    KL = K * L
    KL2 = 2 * KL
    if K1L < 0:
        sinhc = np.sinh(KL2) / KL2
        avg_beta = 0.5 * bet1 * (1 + sinhc) + alf1 * np.sinh(KL) ** 2 / KL / K + (sinhc - 1) * (1 + alf1 ** 2) / (2. * bet1 * K ** 2)
        avg_beta_err = np.sqrt((0.5 * (1 + sinhc) - (sinhc - 1) * (1 + alf1 ** 2) / (2. * bet1 ** 2 * K ** 2)) ** 2 * beterr ** 2 + (np.sinh(KL) ** 2 / KL / K + (sinhc - 1) * (2 * alf1) / (2. * bet1 * K ** 2)) ** 2 * alferr ** 2)
    else:
        sinc = np.sin(KL2) / KL2
        avg_beta = 0.5 * bet1 * (1 + sinc) + alf1 * np.sin(KL) ** 2 / KL / K + (-sinc + 1) * (1 + alf1 ** 2) / (2. * bet1 * K ** 2)
        avg_beta_err = np.sqrt((0.5 * (1 + sinc) - (sinc - 1) * (1 + alf1 ** 2) / (2. * bet1 ** 2 * K ** 2)) ** 2 * beterr ** 2 + (np.sin(KL) ** 2 / KL / K + (sinc - 1) * (2 * alf1) / (2. * bet1 * K ** 2)) ** 2 * alferr ** 2)
    return avg_beta, avg_beta_err
