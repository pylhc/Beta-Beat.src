"""
:module: correction.getdiff

Created on 24/02/18

:author: Lukas Malina

Calculates the difference between GetLLM output and correction plugged in the model.
Provide as first argument the path to the output files of GetLLM.

model inputs:
    twiss_cor.dat and twiss_no.dat

outputs in measurement directory:
    phasex.out and phasey.out
    bbx.out and bby.out
    dx.out, dy.out and ndx.out
    couple.out and chromatic_coupling.out

TODOs and Notes:
    GetLlmMeasurement: possibly extend and use with measurement filters from global correction
                        to be used in sbs, new corrections, getdiff, plot_export

    Expected values after correction to be put in, little tricky with phase column names
    No coupling in twiss_no.dat? not used
"""
from __future__ import print_function

import sys

import numpy as np
import pandas as pd
from os.path import join, isdir, exists, split

import __init__
from segment_by_segment.segment_by_segment import GetLlmMeasurement
from twiss_optics.optics_class import TwissOptics
from utils.tfs_pandas import read_tfs, write_tfs


def getdiff(meas_path=None):
    if meas_path is None:
        meas_path = sys.argv[1]

    if not isdir(meas_path):
        raise IOError("No valid measurement directory:" + meas_path)
    corrected_model_path = join(meas_path, 'twiss_cor.dat')
    uncorrected_model_path = join(meas_path, 'twiss_no.dat')

    meas = GetLlmMeasurement(meas_path)
    twiss_cor = read_tfs(corrected_model_path).set_index('NAME', drop=False)
    twiss_no = read_tfs(uncorrected_model_path).set_index('NAME', drop=False)
    coup_cor = TwissOptics(twiss_cor, quick_init=True).get_coupling(method='cmatrix')
    coup_no = TwissOptics(twiss_no, quick_init=True).get_coupling(method='cmatrix')
    model = pd.merge(twiss_cor, twiss_no, how='inner', on='NAME', suffixes=('_c', '_n'))
    coupling_model = pd.merge(coup_cor, coup_no, how='inner', left_index=True, right_index=True,
                              suffixes=('_c', '_n'))
    coupling_model['NAME'] = coupling_model.index.values

    for plane in ['x', 'y']:
        _write_beta_diff_file(meas_path, meas, model, plane)
        _write_phase_diff_file(meas_path, meas, model, plane)
        _write_disp_diff_file(meas_path, meas, model, plane)
    _write_coupling_diff_file(meas_path, meas, coupling_model)
    _write_norm_disp_diff_file(meas_path, meas, model)
    _write_chromatic_coupling_files(meas_path, corrected_model_path)


def _write_beta_diff_file(meas_path, meas, model, plane):
    up = plane.upper()
    tw = pd.merge(meas.beta[plane], model, how='inner', on='NAME')
    tw['MEA'] = ((tw.loc[:, 'BET' + up] - tw.loc[:, 'BET' + up + 'MDL'])
                 / tw.loc[:, 'BET' + up + 'MDL'])
    tw['ERROR'] = tw.loc[:, 'ERRBET' + up] / tw.loc[:, 'BET' + up + 'MDL']
    tw['MODEL'] = ((tw.loc[:, 'BET' + up + '_c'] - tw.loc[:, 'BET' + up + '_n'])
                   / tw.loc[:, 'BET' + up + '_n'])
    write_tfs(join(meas_path, 'bb' + plane + '.out'),
              tw.loc[:, ['NAME', 'S', 'MEA', 'ERROR', 'MODEL']])


def _write_phase_diff_file(meas_path, meas, model, plane):
    up = plane.upper()
    tw = pd.merge(meas.phase[plane], model, how='inner', on='NAME')
    tw['MEA'] = tw.loc[:, 'PHASE' + up]
    tw['ERROR'] = tw.loc[:, 'STDPH' + up]
    tw['MODEL'] = np.concatenate((np.diff(tw.loc[:, 'MU' + up + '_c']), np.array([0.0])))
    tw['DIFF'] = tw.loc[:, 'PHASE' + up] - tw.loc[:, 'PH' + up + 'MDL']
    tw['DIFF_MDL'] = tw.loc[:, 'MODEL'] - tw.loc[:, 'PH' + up + 'MDL']
    write_tfs(join(meas_path, 'phase' + plane + '.out'),
              tw.loc[tw.index[:-1], ['NAME', 'S', 'MEA', 'ERROR', 'MODEL', 'DIFF', 'DIFF_MDL']])


def _write_disp_diff_file(meas_path, meas, model, plane):
    try:
        up = plane.upper()
        tw = pd.merge(meas.disp[plane], model, how='inner', on='NAME')
        tw['MEA'] = tw.loc[:, 'D' + up] - tw.loc[:, 'D' + up + 'MDL']
        tw['ERROR'] = tw.loc[:, 'STDD' + up]
        tw['MODEL'] = tw.loc[:, 'D' + up + '_c'] - tw.loc[:, 'D' + up + '_n']
        write_tfs(join(meas_path, 'd' + plane + '.out'),
                  tw.loc[:, ['NAME', 'S', 'MEA', 'ERROR', 'MODEL']])
    except IOError:
        pass


def _write_norm_disp_diff_file(meas_path, meas, model):
    try:
        tw = pd.merge(meas.norm_disp, model, how='inner', on='NAME')
        tw['MEA'] = tw.loc[:, 'NDX'] - tw.loc[:, 'NDXMDL']
        tw['ERROR'] = tw.loc[:, 'STDNDX']
        tw['MODEL'] = (tw.loc[:, 'DX_c'] / np.sqrt(tw.loc[:, 'BETX_c'])
                       - tw.loc[:, 'DX_n'] / np.sqrt(tw.loc[:, 'BETX_n']))
        write_tfs(join(meas_path, 'ndx.out'), tw.loc[:, ['NAME', 'S', 'MEA', 'ERROR', 'MODEL']])
    except IOError:
        pass


def _write_coupling_diff_file(meas_path, meas, model):
    tw = pd.merge(meas.coupling, model, how='inner', on='NAME')
    tw['F1001re'] = tw.loc[:, 'F1001R']
    tw['F1001im'] = tw.loc[:, 'F1001I']
    tw['F1001e'] = tw.loc[:, 'FWSTD1']
    tw['F1001re_m'] = np.real(tw.loc[:, 'F1001_c'])
    tw['F1001im_m'] = np.imag(tw.loc[:, 'F1001_c'])
    tw['F1001W_prediction'] = np.sqrt(np.square(tw.loc[:, 'F1001re'] - tw.loc[:, 'F1001re_m'])
                                      + np.square(tw.loc[:, 'F1001im'] - tw.loc[:, 'F1001im_m']))
    tw['in_use'] = 1
    write_tfs(join(meas_path, 'couple.out'),
              tw.loc[:, ['NAME', 'S', 'F1001re', 'F1001im', 'F1001e', 'F1001re_m', 'F1001im_m',
                         'F1001W', 'F1001W_prediction', 'in_use']])


def _write_chromatic_coupling_files(meas_path, cor_path):
    try:
        twiss_plus = read_tfs(join(split(cor_path)[0], "twiss_cor_dpp.dat"), index='NAME')
        twiss_min = read_tfs(join(split(cor_path)[0], "twiss_cor_dpm.dat"), index='NAME')
        deltap = np.abs(twiss_plus.DELTAP - twiss_min.DELTAP)
        plus = TwissOptics(twiss_plus, quick_init=True).get_coupling(method='cmatrix')
        minus = TwissOptics(twiss_min, quick_init=True).get_coupling(method='cmatrix')
        model = pd.merge(plus, minus, how='inner', left_index=True, right_index=True,
                         suffixes=('_p', '_m'))
        model['NAME'] = model.index.values
        if exists(join(meas_path, "chromcoupling_free.out")):
            meas = read_tfs(join(meas_path, "chromcoupling_free.out"))
        else:
            meas = read_tfs(join(meas_path, "chromcoupling.out"))
        tw = pd.merge(meas, model, how='inner', on='NAME')
        cf1001 = (tw.loc[:, 'F1001_p'] - tw.loc[:, 'F1001_m']) / deltap
        tw['Cf1001r_model'] = np.real(cf1001)
        tw['Cf1001i_model'] = np.imag(cf1001)
        write_tfs(join(meas_path, 'chromatic_coupling.out'),
                  tw.loc[:, ['NAME', 'S', 'Cf1001r', 'Cf1001rERR', 'Cf1001i', 'Cf1001iERR',
                             'Cf1001r_model', 'Cf1001i_model']])
    except IOError:
        pass


if __name__ == "__main__":
    getdiff()
