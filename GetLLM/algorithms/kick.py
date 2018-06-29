"""
.. module: kick

Created on 29/06/18

:author: Lukas Malina

It computes kick actions.
"""
from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas
from os.path import join
from utils import logging_tools
import logging
import pandas as pd
import numpy as np
import compensate_excitation


def _calculate_kick(model, mad_ac, getllm_d, files_x, files_y, beta_d, phase_d, beta_from_phase, model_cut, error_cut, output, header_dict, files_dict):
    """
    Fills the following TfsFiles:
     - getkick.out getkickphase.out getkickac.out
    """
    phase_beta = {}
    model_beta = {}

    for plane in ["X", "Y"]:
        model_beta[plane] = pd.DataFrame(model).loc[:, ['S', 'BET' + plane]]
        phase_beta[plane] = pd.merge(model_beta, beta_from_phase[plane].loc[:, ['BET' + plane, 'ERRBET' + plane]],
                                     how='inner', left_index=True, right_index=True, suffixes=('_mdl', '_phase'))
        # TODO filter on modelcut and error cut of beta-beating method,plane-wise
        phase_beta[plane].rename(columns={'BET' + plane + '_phase': 'BET' + plane}, inplace=True)
    try:
        tunes_actions, numbers_bpms = getkick(files_x, files_y, model_beta, phase_beta)
    except IndexError:  # occurs if either no x or no y files exist
        return pd.DataFrame, pd.DataFrame
    column_names = ["DPP", "QX", "QXRMS", "QY", "QYRMS", "NATQX", "NATQXRMS", "NATQY",
                         "NATQYRMS", "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX",
                         "2JXSTD", "2JY", "2JYSTD"]
    kick_model = pd.DataFrame(data=tunes_actions[:, :18], columns=column_names)
    kick_phase = pd.DataFrame(data=np.concatenate((tunes_actions[:, :9], tunes_actions[:, 18:])), columns=column_names)
    model_header = _get_header(header_dict, model_cut, error_cut, from_model=True)
    phase_header = _get_header(header_dict, model_cut, error_cut)
    tfs_pandas.write_tfs(join(output, model_header.FILENAME), kick_model, model_header)
    tfs_pandas.write_tfs(join(output, phase_header.FILENAME), kick_phase, phase_header)
    actions_x, actions_y = tunes_actions[:, 18:20], tunes_actions[:, 20:22] # sqrt2jx, sqrt2Jy

    if getllm_d.accelerator.excitation != AccExcitationMode.FREE:
        tfs_file = files_dict['getkickac.out']
        tfs_file.add_float_descriptor("RescalingFactor_for_X", beta_d.x_ratio_f)
        tfs_file.add_float_descriptor("RescalingFactor_for_Y", beta_d.y_ratio_f)
        tfs_file.add_column_names(column_names + ["sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES", "2JXSTDRES", "2JYRES", "2JYSTDRES"])
        tfs_file.add_column_datatypes(column_types_list + ["%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        [inv_jx, inv_jy, tunes, dpp] = compensate_excitation.getkickac(
            mad_ac, [files_x, files_y], phase_d.ac2bpmac_x, phase_d.ac2bpmac_y, getllm_d.accelerator.get_beam_direction(), getllm_d.lhc_phase)
        for i in range(0, len(dpp)):
            #TODO: in table will be the ratio without f(beta_d.x_ratio) used but rescaling factor is f version(beta_d.x_ratio_f). Check it (vimaier)
            list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i], tunes[4][i], tunes[5][i], tunes[6][i], tunes[7][i], inv_jx[i][0], inv_jx[i][1], inv_jy[i][0], inv_jy[i][1], (inv_jx[i][0] ** 2), (2 * inv_jx[i][0] * inv_jx[i][1]), (inv_jy[i][0] ** 2), (2 * inv_jy[i][0] * inv_jy[i][1]), (inv_jx[i][0] / math.sqrt(beta_d.x_ratio)), (inv_jx[i][1] / math.sqrt(beta_d.x_ratio)), (inv_jy[i][0] / math.sqrt(beta_d.y_ratio)), (inv_jy[i][1] / math.sqrt(beta_d.y_ratio)), (inv_jx[i][0] ** 2 / beta_d.x_ratio), (2 * inv_jx[i][0] * inv_jx[i][1] / beta_d.x_ratio), (inv_jy[i][0] ** 2 / beta_d.y_ratio), (2 * inv_jy[i][0] * inv_jy[i][1] / beta_d.y_ratio)]
            tfs_file.add_table_row(list_row_entries)
            actions_x, actions_y = inv_jx, inv_jx

    return actions_x, actions_y


def _get_header(header_dict, model_cut, error_cut, from_model=False):
    header = header_dict.copy()
    if from_model:
        header['COMMENT'] = "Calculates the kick from the model beta function"
        header['FILENAME'] = "getkick.out"
    else:
        header['COMMENT'] = "Calculates the kick from the beta function from phase"
        header['FILENAME'] = "getkickphase.out"
        header["Threshold_for_abs(beta_d-beta_m)/beta_m"] = model_cut
        header["Threshold_for_uncert(beta_d)/beta_d"] = error_cut
    return header


def getkick(files_x, files_y, model_beta, phase_beta):
    out = np.zeros([len(files_x), 25])
    for i in range(len(files_x)):
        files_x[i] = 0
        files_y[i] =0
        action_x_model, nbpms_x_model = gen_kick_calc(files_x[i], model_beta["X"], "X")
        action_y_model, nbpms_y_model = gen_kick_calc(files_y[i], model_beta["Y"], "Y")
        action_x_phase, nbpms_x_phase = gen_kick_calc(files_x[i], phase_beta["X"], "X")
        action_y_phase, nbpms_y_phase = gen_kick_calc(files_y[i], phase_beta["Y"], "Y")
        # what if the following is not there - except KeyError?
        out[i, 0] = files_x[i].DPP
        out[i, 1] = files_x[i].Q1
        out[i, 2] = files_x[i].Q1RMS
        out[i, 3] = files_y[i].Q2
        out[i, 4] = files_y[i].Q2RMS
        out[i, 5] = files_x[i].NATQ1
        out[i, 6] = files_x[i].NATQ1RMS
        out[i, 7] = files_y[i].NATQ2
        out[i, 8] = files_y[i].NATQ2RMS
        out[i, 9:17] = np.ravel(np.concatenate((action_x_model, action_y_model), axis=1))
        out[i, 17:] = np.ravel(np.concatenate((action_x_phase, action_y_phase), axis=1))

    return out, np.array[nbpms_x_model, nbpms_y_model, nbpms_x_phase, nbpms_y_phase] # numbers should from all files, now only last number


def gen_kick_calc(lin, beta, plane):
    frame = pd.DataFrame(beta).loc[:, ['S', 'BET' + plane]]
    frame = pd.merge(frame, lin.loc[:, ['PK2PK']], how='inner', left_index=True, right_index=True)
    meansqrt2j = (frame.loc[:, 'PK2PK'].values / 2) / np.sqrt(frame.loc[:, 'BET' + plane].values)
    mean2j = np.square(frame.loc[:, 'PK2PK'].values / 2) / frame.loc[:, 'BET' + plane].values
    return np.array[[np.mean(meansqrt2j), np.std(meansqrt2j)], [np.mean(mean2j), np.std(mean2j)]], len(mean2j)
