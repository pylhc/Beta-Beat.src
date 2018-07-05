"""
.. module: kick

Created on 29/06/18

:author: Lukas Malina

It computes kick actions.
# TODO use only arc BPMs
# TODO
"""
from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas
from os.path import join
import pandas as pd
import numpy as np
import compensate_excitation


def calculate_kick(model, mad_ac, getllm_d, files_x, files_y, beta_d, phase_d, output, header_dict,
                   files_dict):
    """
    Fills the following TfsFiles:
     - getkick.out getkickac.out
    """
    model_betas = pd.DataFrame(model).loc[:, ['S', 'BETX', 'BETY']]
    try:
        tunes_actions = _getkick(files_x, files_y, model_betas)
    except IndexError:  # occurs if either no x or no y files exist
        return files_dict, pd.DataFrame, pd.DataFrame
    column_names = ["DPP", "QX", "QXRMS", "QY", "QYRMS", "NATQX", "NATQXRMS", "NATQY", "NATQYRMS",
                    "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX", "2JXSTD", "2JY",
                    "2JYSTD"]
    kick_frame = pd.DataFrame(data=tunes_actions, columns=column_names)
    header = _get_header(header_dict)
    tfs_pandas.write_tfs(join(output, header['FILENAME']), kick_frame, header)
    actions_x, actions_y = tunes_actions[:, 9:11], tunes_actions[:, 11:13]  # sqrt2jx, sqrt2Jy

    if getllm_d.accelerator.excitation != AccExcitationMode.FREE:
        tfs_file = files_dict['getkickac.out']
        tfs_file.add_float_descriptor("RescalingFactor_for_X", beta_d.x_ratio_f)
        tfs_file.add_float_descriptor("RescalingFactor_for_Y", beta_d.y_ratio_f)
        tfs_file.add_column_names(
            column_names + ["sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES",
                            "2JXSTDRES", "2JYRES", "2JYSTDRES"])
        tfs_file.add_column_datatypes(25 * ["%le"])
        [inv_jx, inv_jy, tunes, dpp] = compensate_excitation.getkickac(
            mad_ac, [files_x, files_y], phase_d.ac2bpmac_x, phase_d.ac2bpmac_y,
            getllm_d.accelerator.get_beam_direction(), getllm_d.lhc_phase)
        for i in range(0, len(dpp)):
            # TODO: in table will be the ratio without f(beta_d.x_ratio) used but rescaling factor
            #  is f version(beta_d.x_ratio_f). Check it (vimaier)
            list_row_entries = [dpp[i], tunes[0][i], tunes[1][i], tunes[2][i], tunes[3][i],
                                tunes[4][i], tunes[5][i], tunes[6][i], tunes[7][i], inv_jx[i][0],
                                inv_jx[i][1], inv_jy[i][0], inv_jy[i][1], (inv_jx[i][0] ** 2),
                                (2 * inv_jx[i][0] * inv_jx[i][1]), (inv_jy[i][0] ** 2),
                                (2 * inv_jy[i][0] * inv_jy[i][1]),
                                (inv_jx[i][0] / np.sqrt(beta_d.x_ratio)),
                                (inv_jx[i][1] / np.sqrt(beta_d.x_ratio)),
                                (inv_jy[i][0] / np.sqrt(beta_d.y_ratio)),
                                (inv_jy[i][1] / np.sqrt(beta_d.y_ratio)),
                                (inv_jx[i][0] ** 2 / beta_d.x_ratio),
                                (2 * inv_jx[i][0] * inv_jx[i][1] / beta_d.x_ratio),
                                (inv_jy[i][0] ** 2 / beta_d.y_ratio),
                                (2 * inv_jy[i][0] * inv_jy[i][1] / beta_d.y_ratio)]

            tfs_file.add_table_row(list_row_entries)
            actions_x, actions_y = inv_jx, inv_jx

    return files_dict, actions_x, actions_y


def _get_header(header_dict):
    header = header_dict.copy()
    header['COMMENT'] = "Calculates the kick from the model beta function"
    header['FILENAME'] = "getkick.out"
    return header


def _getkick(files_x, files_y, model_beta):
    out = np.zeros([len(files_x), 17])
    for i in range(len(files_x)):
        action_x_model = _gen_kick_calc(files_x[i], model_beta, "X")
        action_y_model = _gen_kick_calc(files_y[i], model_beta, "Y")
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
        out[i, 9:] = np.ravel(np.concatenate((action_x_model, action_y_model), axis=1))
    return out


def _gen_kick_calc(lin, beta, plane):
    frame = pd.merge(beta, lin.loc[:, ['PK2PK']], how='inner', left_index=True, right_index=True)
    meansqrt2j = (frame.loc[:, 'PK2PK'].values / 2) / np.sqrt(frame.loc[:, 'BET' + plane].values)
    mean2j = np.square(frame.loc[:, 'PK2PK'].values / 2) / frame.loc[:, 'BET' + plane].values
    return np.array([[np.mean(meansqrt2j), np.std(meansqrt2j)], [np.mean(mean2j), np.std(mean2j)]])
