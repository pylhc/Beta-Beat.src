"""
.. module: kick

Created on 29/06/18

:author: Lukas Malina

It computes kick actions.
# TODO use only arc BPMs
"""
from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas, bpm
from os.path import join
import pandas as pd
import numpy as np
import compensate_excitation


def calculate_kick(measure_input, input_files, model, mad_ac, beta_d, header_dict):
    """
    Fills the following TfsFiles:
     - getkick.out getkickac.out
    """
    model_betas = pd.DataFrame(model).loc[:, ['S', 'BETX', 'BETY']]
    try:
        tunes_actions = _getkick(input_files["X"], input_files["Y"], model_betas)
    except IndexError:  # occurs if either no x or no y files exist
        return pd.DataFrame, pd.DataFrame
    column_names = ["DPP", "QX", "QXRMS", "QY", "QYRMS", "NATQX", "NATQXRMS", "NATQY", "NATQYRMS",
                    "sqrt2JX", "sqrt2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JX", "2JXSTD", "2JY",
                    "2JYSTD"]
    kick_frame = pd.DataFrame(data=tunes_actions, columns=column_names)
    header = _get_header(header_dict, beta_d)
    tfs_pandas.write_tfs(join(measure_input.outputdir, header['FILENAME']), kick_frame, header)
    actions_x, actions_y = tunes_actions[:, 9:11], tunes_actions[:, 11:13]  # sqrt2jx, sqrt2Jy

    if measure_input.accelerator.excitation != AccExcitationMode.FREE:
        column_names_ac = column_names + ["sqrt2JXRES", "sqrt2JXSTDRES", "sqrt2JYRES", "sqrt2JYSTDRES", "2JXRES",
                            "2JXSTDRES", "2JYRES", "2JYSTDRES"]
        [inv_jx, inv_jy, tunes, dpp] = getkickac([input_files["X"], input_files["Y"]], mad_ac)
        datas=[]
        for i in range(0, len(dpp)):
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
            datas.append(list_row_entries)
        kick_frame_ac = pd.DataFrame(data=np.array(datas), columns=column_names_ac)
        header_ac = _get_header(header_dict, beta_d, ac=True)
        tfs_pandas.write_tfs(join(measure_input.outputdir, header['FILENAME']), kick_frame_ac, header_ac)
        actions_x, actions_y = inv_jx, inv_jx
    return actions_x, actions_y


def _get_header(header_dict, beta_d, ac=False):
    header = header_dict.copy()
    header['COMMENT'] = "Calculates the kick from the model beta function"
    header['FILENAME'] = 'getkick' + ac * 'ac' + '.out'
    if ac:
        header["RescalingFactor_for_X"] = beta_d.x_ratio_f
        header["RescalingFactor_for_Y"] = beta_d.y_ratio_f

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


def getkickac(MADTwiss_ac, files, accel):
    invarianceJx = []
    invarianceJy = []

    tunex = []
    tuney = []
    tunexRMS = []
    tuneyRMS = []

    nat_tunex = []
    nat_tuney = []
    nat_tunexRMS = []
    nat_tuneyRMS = []

    dpp = []

    all_bpms_x = bpm.model_intersect(bpm.intersect(files[0]), MADTwiss_ac)
    all_bpms_y = bpm.model_intersect(bpm.intersect(files[1]), MADTwiss_ac)

    good_bpms_for_kick_x = all_bpms_x[accel.get_element_types_mask(all_bpms_x, "arc_bpm")]
    good_bpms_for_kick_y = all_bpms_y[accel.get_element_types_mask(all_bpms_y, "arc_bpm")]
    Jx2sq, Jx2sq_std = compensate_excitation.get_kick_from_bpm_list_w_acdipole(MADTwiss_ac, good_bpms_for_kick_x,
                                                                               files[0], 'X')
    Jy2sq, Jy2sq_std = compensate_excitation.get_kick_from_bpm_list_w_acdipole(MADTwiss_ac, good_bpms_for_kick_y,
                                                                               files[1], 'Y')

    for j in range(len(files[0])):
        tw_x = files[0][j]
        tw_y = files[1][j]
        invarianceJx.append([Jx2sq[j], Jx2sq_std[j]])
        invarianceJy.append([Jy2sq[j], Jy2sq_std[j]])

        dpp.append(getattr(tw_x, "DPP", 0.0))

        tunex.append(getattr(tw_x, "Q1", 0.0))
        tuney.append(getattr(tw_y, "Q2", 0.0))
        tunexRMS.append(getattr(tw_x, "Q1RMS", 0.0))
        tuneyRMS.append(getattr(tw_y, "Q2RMS", 0.0))

        nat_tunex.append(getattr(tw_x, "NATQ1", 0.0))
        nat_tuney.append(getattr(tw_y, "NATQ2", 0.0))
        nat_tunexRMS.append(getattr(tw_x, "NATQ1RMS", 0.0))
        nat_tuneyRMS.append(getattr(tw_y, "NATQ2RMS", 0.0))

    tune_values_list = [tunex, tunexRMS, tuney, tuneyRMS, nat_tunex, nat_tunexRMS, nat_tuney,
                        nat_tuneyRMS]
    return [invarianceJx, invarianceJy, tune_values_list, dpp]