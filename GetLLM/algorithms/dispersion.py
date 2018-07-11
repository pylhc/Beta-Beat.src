"""
.. module: dispersion

Created on 28/06/18

:author: Lukas Malina

It computes orbit, dispersion and normalised dispersion.
"""
from os.path import join
import pandas as pd
import numpy as np
from utils import tfs_pandas
from utils import stats
SCALES = {'um': 1.0e-6, 'mm': 1.0e-3, 'cm': 1.0e-2, 'm': 1.0}


def calculate_orbit_and_dispersion(input_files, tune_d, model, header_dict, unit, cut, beta_from_phase, output):
    """
    Calculates orbit and dispersion, fills the following TfsFiles:
       getCOx.out        getCOy.out
       getNDx.out        getDx.out        getDy.out
    Args:
        twiss_d: Holds tfs pandas of the linx/y files. Later to be replaced by two lists of tfs_pandas
        tune_d: Holds tunes and phase advances
        model:  Model tfs panda
        header_dict: OrderedDict containing information about the analysis
        unit: Unit of orbit measurements, e.g. m, mm
        cut: orbit cut
        beta_from_phase:
        output: output folder

    Returns:

    """
    for plane in PLANES:
        orbit_header = _get_header(header_dict, tune_d, 'getCO' + plane.lower() + '.out', orbit=True)
        dispersion_header = _get_header(header_dict, tune_d, 'getD' + plane.lower() + '.out')
        _calculate_orbit(model, input_files, plane, orbit_header, output)
        _calculate_dispersion(model, input_files, plane, dispersion_header, unit, cut, output)
        if plane == 'X':
            ndx_header = _get_header(header_dict, tune_d, 'getNDx.out', orbit=False)
            _calculate_normalised_dispersion(model, input_files, plane, beta_from_phase, ndx_header,
                                             unit, cut, output)


def _get_header(header_dict, tune_d, filename, orbit=False):
    header = header_dict.copy()
    if orbit:
        header['TABLE'] = 'ORBIT'
        header['TYPE'] = 'ORBIT'
    # TODO: ['SEQUENCE'] getllm_d.accel
    header['Q1'] = tune_d.q1
    header['Q2'] = tune_d.q2
    header['FILENAME'] = filename
    return header


def _calculate_orbit(model, input_files, plane, header, output):
    df_orbit = pd.DataFrame(model).loc[:, ['S', 'MU' + plane, plane]]
    df_orbit.rename(columns={'MU' + plane: 'MU' + plane + 'MDL', plane: plane + 'MDL'}, inplace=True)
    df_orbit = pd.merge(df_orbit, input_files.get_joined_frame(plane, ['CO', 'CORMS']), how='inner',
                        left_index=True, right_index=True)
    df_orbit['COUNT'] = len(input_files.get_columns(df_orbit, ['CO']))
    df_orbit[plane] = stats.weighted_mean(input_files.get_data(df_orbit, 'CO'), axis=1)
    df_orbit['STD' + plane] = stats.weighted_error(input_files.get_data(df_orbit, 'CO'), axis=1)
    output_df = df_orbit.loc[:, ['S', 'COUNT', plane, 'STD' + plane, plane + 'MDL', 'MU' + plane + 'MDL']]
    tfs_pandas.write_tfs(join(output, header['FILENAME']), output_df, header, save_index='NAME')
    return output_df


def _calculate_dispersion(model, input_files, plane, header, unit, cut, output, order=1):
    df_orbit = pd.DataFrame(model).loc[:, ['S', 'MU' + plane, 'DP' + plane, 'D' + plane, plane]]
    df_orbit.rename(columns={'MU' + plane: 'MU' + plane + 'MDL', 'DP' + plane: 'DP' + plane + 'MDL',
                             'D' + plane: 'D' + plane + 'MDL', plane: plane + 'MDL'}, inplace=True)
    df_orbit = pd.merge(df_orbit, input_files.get_joined_frame(plane, ['CO', 'CORMS']), how='inner',
                        left_index=True, right_index=True)
    df_orbit['COUNT'] = len(input_files.get_columns(df_orbit, ['CO']))
    dpps = input_files.get_dpps(plane)
    if np.max(dpps) - np.min(dpps) == 0.0:
        return  # temporary solution
        # raise ValueError('Cannot calculate dispersion, only a single momentum data')
    fit = np.polyfit(dpps, SCALES[unit] * input_files.get_data(df_orbit, 'CO').T, order, cov=True)
    # in the fit results the coefficients are sorted by power in decreasing order
    df_orbit['D' + plane] = fit[0][-2, :].T
    df_orbit['STDD' + plane] = np.sqrt(fit[1][-2, -2, :].T)
    df_orbit[plane] = fit[0][-1, :].T
    df_orbit['STD' + plane] = np.sqrt(fit[1][-1, -1, :].T)
    # since we get variances from the fit, maybe we can include the variances of fitted points
    df_orbit = df_orbit.loc[np.abs(df_orbit.loc[:, plane]) < cut*SCALES[unit], :]
    df_orbit['DP' + plane] = _calculate_dp(model,
                                           df_orbit.loc[:, ['D' + plane, 'STDD' + plane]], plane)
    output_df = df_orbit.loc[:,
                             ['S', 'COUNT', 'D' + plane, 'STDD' + plane, plane, 'STD' + plane, 'DP' + plane,
                              'D' + plane + 'MDL', 'DP' + plane + 'MDL', 'MU' + plane + 'MDL']]
    tfs_pandas.write_tfs(join(output, header['FILENAME']), output_df, header, save_index='NAME')
    return output_df


def _calculate_normalised_dispersion(model, list_of_df, beta, header, unit, cut, output):
    #TODO there are no errors from orbit
    df_orbit = pd.DataFrame(model).loc[:, ['S', 'MUX', 'DPX', 'DX', 'X', 'BETX']]
    df_orbit['NDXMDL'] = df_orbit.loc[:, 'DX'] / np.sqrt(df_orbit.loc[:, 'BETX'])
    df_orbit.rename(columns={'MUX': 'MUXMDL', 'DPX': 'DPXMDL', 'DX': 'DXMDL', 'X': 'XMDL'}, inplace=True)
    df_orbit = df_orbit.assign(CO=0.0, CORMS=0.0, AMPX=0.0, COUNT=len(list_of_df))
    orbit_columns = []
    dpps = np.empty(len(list_of_df))
    amps = []
    for i, df in enumerate(list_of_df):
        df_orbit = pd.merge(df_orbit, df.loc[:, ['CO', 'CORMS', 'AMPX']], how='inner', left_index=True,
                            right_index=True, suffixes=('', str(i + 1)))
        orbit_columns.append('CO' + str(i+1))
        dpps[i] = df.DPP
        if df.DPP == 0.0:
            amps.append('AMPX' + str(i+1))
    df_orbit = pd.merge(df_orbit, beta.loc[:, ['BETX', 'ERRBETX']], how='inner', left_index=True,
                        right_index=True, suffixes=('', '_phase'))
    if np.max(dpps) - np.min(dpps) == 0.0:
        return  # temporary solution
        # raise ValueError('Cannot calculate dispersion, only a single dpoverp')
    fit = np.polyfit(np.array(dpps), SCALES[unit] * df_orbit.loc[:, orbit_columns].values.T, 1, cov=True)
    amps = np.array(amps)
    df_orbit['NDX_unscaled'] = fit[0][-2, :].T / np.mean(df_orbit.loc[:, amps].values, axis=1) # TODO there is no error from AMPX
    df_orbit['STDNDX_unscaled'] = np.sqrt(fit[1][-2,-2, :].T) / np.mean(df_orbit.loc[:, amps].values, axis=1)
    df_orbit = df_orbit.loc[np.abs(fit[0][-1, :].T) < cut * SCALES[unit], :]
    global_factor = np.sum(df_orbit.loc[:, 'NDXMDL'].values) / np.sum(df_orbit.loc[:, 'NDX_unscaled'].values) # TODO should be arc BPMS
    df_orbit['NDX'] = global_factor * df_orbit.loc[:, 'NDX_unscaled']
    df_orbit['STDNDX'] = global_factor * df_orbit.loc[:, 'STDNDX_unscaled']
    df_orbit['DX'] = df_orbit.loc[:, 'NDX'] * np.sqrt(df_orbit.loc[:, 'BETX_phase'])
    df_orbit['STDDX'] = df_orbit.loc[:, 'STDNDX'] * np.sqrt(df_orbit.loc[:, 'BETX_phase'])
    df_orbit['DPX'] = _calculate_dp(model, df_orbit.loc[:, ['DX', 'STDDX']], "X")
    output_df = df_orbit.loc[:, ['S', 'COUNT', 'NDX', 'STDNDX', 'DX', 'DPX',
                                 'NDXMDL', 'DXMDL', 'DPXMDL', 'MUXMDL']]
    tfs_pandas.write_tfs(join(output, header['FILENAME']), output_df, header, save_index='NAME')
    return output_df


def _calculate_dp(model, disp, plane):
    df = pd.DataFrame(model).loc[:, ['S', 'MU' + plane, 'DP' + plane, 'D' + plane,
                                         'BET' + plane, 'ALF' + plane]]
    df = pd.merge(df, disp.loc[:, ['D' + plane, 'STDD' + plane]], how='inner', left_index=True,
                  right_index=True, suffixes=('', 'meas'))
    shifted = np.roll(df.index.values, 1)
    print(shifted[0])
    print(shifted[-1])
    p_mdl_12 = df.loc[shifted, 'MU' + plane].values - df.loc[:, 'MU' + plane].values
    p_mdl_12[-1] = p_mdl_12[-1] + model['Q' + str(1+(plane == "Y"))]
    phi_12 = p_mdl_12 * 2 * np.pi
    m11 = np.sqrt(df.loc[shifted, 'BET' + plane] / df.loc[:, 'BET' + plane]) * (np.cos(phi_12)
                  + df.loc[:, 'ALF' + plane] * np.sin(phi_12))
    m12 = np.sqrt(df.loc[shifted, 'BET' + plane] * df.loc[:, 'BET' + plane]) * np.sin(phi_12)
    m13 = df.loc[shifted, 'D' + plane] - m11 * df.loc[:, 'D' + plane] - m12 * df.loc[:, 'DP' + plane]
    return (-m13 + df.loc[shifted, 'D' + plane + 'meas'] - m11 * df.loc[:, 'D' + plane + 'meas']) / m12
