"""
.. module: beta_from_amplitude

Created on 05/07/18

:author: Lukas Malina

It computes beta from amplitude.
"""

from os.path import join
import numpy as np
import pandas as pd
import compensate_excitation
from utils import tfs_pandas
from model.accelerators.accelerator import AccExcitationMode
# TODO all action scaling should be done with arc BPMs


def calculate_beta_from_amplitude(getllm_d, twisses, tune_d, phase_d, accelerator, header_dict, output):  # twises - dictionary by plane
    """
    Calculates beta and fills the following TfsFiles:
        getampbetax.out        getampbetax_free.out        getampbetax_free2.out
        getampbetay.out        getampbetay_free.out        getampbetay_free2.out

    Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            accel and beam_direction are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances
        'phase_d': _PhaseData (In-param, values will only be read)
            Holds results from get_phases

    Returns:
    """
    mad_twiss = accelerator.get_model_tfs()
    if accelerator.excitation != AccExcitationMode.FREE:
        mad_ac = accelerator.get_driven_tfs()
    else:
        mad_ac = mad_twiss

    for plane in ['X', 'Y']:
        # column_names = ["NAME", "S", "COUNT", "BET" + plane, "BET" + plane + "STD", "BET" + plane + "MDL", "MU" + plane + "MDL", "BET" + plane + "RES", "BET" + plane + "STDRES"]
        if not twiss_d.has_zero_dpp_x():
            beta_amp = beta_from_amplitude(mad_ac, twisses[plane], 'X')
            beta_amp['DPP'] = 0
            x_ratio = np.mean(beta_phase / beta_amp)  # over good arc bpms : both betas positive abs(beta_phase / beta_amp)<100
            beta_amp['BET' + plane + 'RES'] = beta_amp.loc[:, 'BET' + plane] * x_ratio
            beta_amp['BET' + plane + 'STDRES'] = beta_amp.loc[:, 'BET' + plane + 'STD'] * x_ratio
            header_d = _get_header(header_dict, tune_d, np.std(beta_amp.loc[:, 'DELTABET' + plane].values), x_ratio, 'getampbeta' + plane.lower() + '.out', free=False)
            tfs_pandas.write_tfs(join(output, header_d['FILENAME']), beta_amp, header_d, save_index='NAME')
            # -- ac to free amp beta
            if getllm_d.accelerator.excitation is not AccExcitationMode.FREE:
                beta_amp_f = compensate_excitation.get_free_beta_from_amp_eq(mad_ac, twiss_d.zero_dpp_x, tune_d.q1, tune_d.q1f, phase_d.ac2bpmac_x, 'H', getllm_d, arcbpms)
                x_ratio_f = np.mean(beta_phase / beta_amp)  # over good arc bpms : both betas positive, 0.1 < abs(beta_phase / beta_amp)<10
                header_f = _get_header(header_dict, tune_d, np.std(beta_amp_f.loc[:, 'DELTABET' + plane].values), x_ratio_f, 'getampbeta' + plane.lower() + '_free.out', free=True)

                beta_amp_f['BET' + plane + 'RES'] = beta_amp_f.loc[:, 'BET' + plane] * x_ratio_f
                beta_amp_f['BET' + plane + 'STDRES'] = beta_amp_f.loc[:, 'BET' + plane + 'STD'] * x_ratio_f
                tfs_pandas.write_tfs(join(output, header_f['FILENAME']), beta_amp_f, header_f, save_index='NAME')

                # FREE2 calculation
                beta_amp_f2 = pd.DataFrame(beta_amp)
                beta_amp_f2['BET' + plane] = _get_free_amp_beta(beta_amp_f2.loc[:, 'BET' + plane], mad_ac, mad_twiss, plane)
                beta_amp_f2['BET' + plane + 'RES'] = _get_free_amp_beta(beta_amp_f2.loc[:, 'BET' + plane + 'RES'], mad_ac, mad_twiss, plane)
                header_f2 = _get_header(header_dict, tune_d, np.std(beta_amp_f2.loc[:, 'DELTABET' + plane].values), x_ratio,
                                           'getampbeta' + plane.lower() + '_free2.out', free=True)
                tfs_pandas.write_tfs(join(output, header_f2['FILENAME']), beta_amp_f2, header_f2, save_index='NAME')


def _get_free_amp_beta(df_meas,  mad_ac, mad_twiss, plane):
    df = pd.merge(pd.DataFrame(df_meas), mad_ac.loc[:, 'BET' + plane], how='inner',
                       left_index=True, right_index=True, suffixes=('', 'ac'))
    df = pd.merge(df, mad_twiss.loc[:, 'BET' + plane], how='inner', left_index=True,
                       right_index=True, suffixes=('', 'f'))
    return df.loc[:, "BET" + plane] * df.loc[:, "BET" + plane + 'f'] / df.loc[:, "BET" + plane + 'ac']


def beta_from_amplitude(model, list_of_df, plane):
    df_amp_beta = pd.DataFrame(model).loc[:, ['S', 'MU' + plane, 'BET' + plane]]
    df_amp_beta.rename(columns={'MU' + plane: 'MU' + plane + 'MDL',
                                'BET' + plane: 'BET' + plane + 'MDL'}, inplace=True)
    df_amp_beta = df_amp_beta.assign(AMPX=0.0, AMPY=0.0, COUNT=len(list_of_df))
    measured_amps_columns = []
    for i, df in enumerate(list_of_df):
        df_amp_beta = pd.merge(df_amp_beta, df.loc[:, 'AMP' + plane], how='inner', left_index=True,
                               right_index=True, suffixes=('', str(i + 1)))
        measured_amps_columns.append('AMP' + plane + str(i + 1))
    df_amp_beta['AMP' + plane] = np.mean(df_amp_beta.loc[:, measured_amps_columns].values, axis=1)
    # amplitudes are first averaged over files then squared and averaged over BPMs
    kick = np.mean(np.square(df_amp_beta.loc[:,'AMP' + plane].values) / df_amp_beta.loc[:, 'BET' + plane + 'MDL'].values)
    # amplitudes are first squared then averaged
    kick2 = np.mean(np.square(df_amp_beta.loc[:, measured_amps_columns].values) / df_amp_beta.loc[:, 'BET' + plane + 'MDL'].values, axis=0)

    df_amp_beta['BET' + plane] = np.square(df_amp_beta.loc[:, 'AMP' + plane].values) / kick
    df_amp_beta['BET' + plane + 'STD'] = np.std((np.square(df_amp_beta.loc[:, measured_amps_columns].values).T / kick2).T, axis=1)
    df_amp_beta['DELTABET' + plane] = (df_amp_beta.loc[:, 'BET' + plane] - df_amp_beta.loc[:, 'BET' + plane + 'MDL']) / df_amp_beta.loc[:, 'BET' + plane + 'MDL']
    return df_amp_beta.loc[:, ['S', 'COUNT', 'BET' + plane, 'BET' + plane + 'STD',
                    'BET' + plane + 'MDL', 'MU' + plane + 'MDL', 'DELTABET' + plane]]


def _get_header(header_dict, tune_d, rmsbbeat, scaling_factor, file_name, free=False):
    header = header_dict.copy()
    if free:
        header['Q1'] = tune_d.q1f
        header['Q2'] = tune_d.q2f
    else:
        header['Q1'] = tune_d.q1
        header['Q2'] = tune_d.q2
    header['RMSbetabeat'] = rmsbbeat
    header['RescalingFactor'] = scaling_factor
    header['FILENAME'] = file_name
    return header
