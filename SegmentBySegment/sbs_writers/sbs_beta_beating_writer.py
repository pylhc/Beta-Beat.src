import __init__  # @UnusedImport
import os
import sys

from math import sqrt
from Utilities import tfs_file_writer


def write_beta_beat(element_name,
                    measured_hor_beta_phase, measured_ver_beta_phase,
                    measured_hor_beta_amp, measured_ver_beta_amp,
                    measured_kmod_hor, measured_kmod_ver,
                    propagated_models, save_path):
    file_beta_phase_x, file_beta_phase_y = _get_beta_beat_from_phase_tfs_files(element_name, save_path)
    file_beta_amp_x, file_beta_amp_y = _get_beta_beat_from_amp_tfs_files(element_name, save_path)
    file_kmod_beta_beat_x, file_kmod_beta_beat_y = _get_kmod_beta_beat_tfs_files(element_name, save_path)

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
        _write_kmod_beta_beat_for_plane(
            file_kmod_beta_beat_x, "X", kmod_bpms_list, measured_kmod_hor,
            model_propagation, model_cor, model_back_propagation, model_back_cor,
            save_path
        )
        _write_kmod_beta_beat_for_plane(
            file_kmod_beta_beat_y, "Y", kmod_bpms_list, measured_kmod_ver,
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


def _get_kmod_beta_beat_tfs_files(element_name, save_path):
    file_kmod_beta_beat = {}
    file_kmod_beta_beat["X"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbskmodbetabeatx_" + element_name + ".out"))
    file_kmod_beta_beat["Y"] = tfs_file_writer.TfsFileWriter.open(os.path.join(save_path, "sbskmodbetabeaty_" + element_name + ".out"))

    for plane in ["X", "Y"]:
        column_names = [
            "NAME", "S",
            "BETABEAT" + plane, "ERRBETABEAT" + plane,
            "BETABEATCOR" + plane,
            "BETABEATBACK" + plane, "ERRBETABEATBACK" + plane,
            "BETABEATBACKCOR" + plane
        ]
        file_kmod_beta_beat[plane].add_column_names(column_names)
        file_kmod_beta_beat[plane].add_column_datatypes(["%bpm_s"] + ["%le"] * (len(column_names) - 1))

    return file_kmod_beta_beat["X"], file_kmod_beta_beat["Y"]


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


def _write_kmod_beta_beat_for_plane(file_kmod_beta_beat, plane, bpms_list,
                                    measured_kmod, model_propagation, model_cor,
                                    model_back_propagation, model_back_cor,
                                    output_path):
    for bpm in bpms_list:
        bpm_s = bpm[0]
        bpm_name = bpm[1]

        beta_propagation = _get_from_twiss(model_propagation, bpm_name, "BET", plane)
        beta_back_propagation = _get_from_twiss(model_back_propagation, bpm_name, "BET", plane)

        # Beta from kmod beating (front)
        beta_kmod = _get_from_twiss(measured_kmod, bpm_name, "BET", plane)
        beta_beat_kmod = (beta_kmod - beta_propagation) / beta_propagation
        err_beta_beat_kmod = _get_from_twiss(measured_kmod, bpm_name, "STDBET", plane) / beta_propagation

        # Beta from kmod beating (back)
        beta_kmod_back = _get_from_twiss(measured_kmod, bpm_name, "BET", plane)
        beta_beat_kmod_back = (beta_kmod_back - beta_back_propagation) / beta_back_propagation
        err_beta_beat_kmod_back = _get_from_twiss(measured_kmod, bpm_name, "STDBET", plane) / beta_back_propagation

        model_s = measured_kmod.S[measured_kmod.indx[bpm_name]]
        beta_model = _get_from_twiss(measured_kmod, bpm_name, "BET", plane, "MDL")
        file_kmod_beta_beat.add_table_row([
            bpm_name, bpm_s,
            beta_beat_kmod, err_beta_beat_kmod,
            beta_beat_kmod_back, err_beta_beat_kmod_back,
            beta_model, model_s
        ])
    file_kmod_beta_beat.write_to_file()


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
    # SORT by S
    result = []
    x0 = list_of_files[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))
    result.sort()
    return result
