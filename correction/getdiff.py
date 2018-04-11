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
from utils import logging_tools

LOG = logging_tools.get_logger(__name__)


# Main invocation ############################################################


def getdiff(meas_path=None):
    """ Calculates the differences between measurement, corrected and uncorrected model.

    After running madx and creating the model with (twiss_cor) and without (twiss_no)
    corrections, run this functions to get tfs-files with the differences between measurements
    and models.

    Args:
        meas_path (str): Path to the measurement folder.
        Needs to contain twiss_cor.dat and twiss_no.dat.
    """
    if meas_path is None:
        meas_path = sys.argv[1]

    LOG.debug("Started 'getdiff' for measurment dir '{:s}'".format(meas_path))

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
    LOG.debug("Finished 'getdiff'.")


# Writing Functions ##########################################################


def _write_beta_diff_file(meas_path, meas, model, plane):
    LOG.debug("Calculating beta diff.")
    up = plane.upper()
    tw = pd.merge(meas.beta[plane], model, how='inner', on='NAME')
    tw['MEA'] = ((tw.loc[:, 'BET' + up] - tw.loc[:, 'BET' + up + 'MDL'])
                 / tw.loc[:, 'BET' + up + 'MDL'])
    tw['ERROR'] = tw.loc[:, 'ERRBET' + up] / tw.loc[:, 'BET' + up + 'MDL']
    tw['MODEL'] = ((tw.loc[:, 'BET' + up + '_c'] - tw.loc[:, 'BET' + up + '_n'])
                   / tw.loc[:, 'BET' + up + '_n'])
    tw['EXPECT'] = tw['MEA'] - tw['MODEL']
    write_tfs(join(meas_path, 'bb' + plane + '.out'),
              tw.loc[:, ['NAME', 'S', 'MEA', 'ERROR', 'MODEL', 'EXPECT']])


def _write_phase_diff_file(meas_path, meas, model, plane):
    LOG.debug("Calculating phase diff.")
    up = plane.upper()
    tw = pd.merge(meas.phase[plane], model, how='inner', on='NAME')
    tw['MEA'] = tw.loc[:, 'PHASE' + up]
    tw['ERROR'] = tw.loc[:, 'STDPH' + up]
    tw['MODEL'] = np.concatenate((np.diff(tw.loc[:, 'MU' + up + '_c']), np.array([0.0])))
    tw['DIFF'] = tw.loc[:, 'PHASE' + up] - tw.loc[:, 'PH' + up + 'MDL']
    tw['DIFF_MDL'] = tw.loc[:, 'MODEL'] - tw.loc[:, 'PH' + up + 'MDL']
    tw['EXPECT'] = tw['DIFF'] - tw['DIFF_MDL']
    write_tfs(join(meas_path, 'phase' + plane + '.out'),
              tw.loc[tw.index[:-1],
                     ['NAME', 'S', 'MEA', 'ERROR', 'MODEL', 'DIFF', 'DIFF_MDL', 'EXPECT']])


def _write_disp_diff_file(meas_path, meas, model, plane):
    LOG.debug("Calculating dispersion diff.")
    try:
        up = plane.upper()
        tw = pd.merge(meas.disp[plane], model, how='inner', on='NAME')
        tw['MEA'] = tw.loc[:, 'D' + up] - tw.loc[:, 'D' + up + 'MDL']
        tw['ERROR'] = tw.loc[:, 'STDD' + up]
        tw['MODEL'] = tw.loc[:, 'D' + up + '_c'] - tw.loc[:, 'D' + up + '_n']
        tw['EXPECT'] = tw['MEA'] - tw['MODEL']
        write_tfs(join(meas_path, 'd' + plane + '.out'),
                  tw.loc[:, ['NAME', 'S', 'MEA', 'ERROR', 'MODEL', 'EXPECT']])
    except IOError:
        LOG.debug("Dispersion measurements not found. Skipped.")


def _write_norm_disp_diff_file(meas_path, meas, model):
    LOG.debug("Calculating normalized dispersion diff.")
    try:
        tw = pd.merge(meas.norm_disp, model, how='inner', on='NAME')
        tw['MEA'] = tw.loc[:, 'NDX'] - tw.loc[:, 'NDXMDL']
        tw['ERROR'] = tw.loc[:, 'STDNDX']
        tw['MODEL'] = (tw.loc[:, 'DX_c'] / np.sqrt(tw.loc[:, 'BETX_c'])
                       - tw.loc[:, 'DX_n'] / np.sqrt(tw.loc[:, 'BETX_n']))
        tw['EXPECT'] = tw['MEA'] - tw['MODEL']
        write_tfs(join(meas_path, 'ndx.out'),
                  tw.loc[:, ['NAME', 'S', 'MEA', 'ERROR', 'MODEL', 'EXPECT']])
    except IOError:
        LOG.debug("Normalized dispersion measurements not found. Skipped.")


def _write_coupling_diff_file(meas_path, meas, model):
    LOG.debug("Calculating coupling diff.")
    tw = pd.merge(meas.coupling, model, how='inner', on='NAME')
    out_columns = ['NAME', 'S']
    for rdt in ['F1001', 'F1010']:
        tw[rdt+'re'] = tw.loc[:, rdt+'R']
        tw[rdt+'im'] = tw.loc[:, rdt+'I']
        tw[rdt+'e'] = tw.loc[:, 'FWSTD1']
        tw[rdt+'re_m'] = np.real(tw.loc[:, rdt+'_c'])
        tw[rdt+'im_m'] = np.imag(tw.loc[:, rdt+'_c'])
        tw[rdt+'re_prediction'] = tw.loc[:, rdt+'re'] - tw.loc[:, rdt+'re_m']
        tw[rdt+'im_prediction'] = tw.loc[:, rdt+'im'] - tw.loc[:, rdt+'im_m']
        tw[rdt+'W_prediction'] = np.sqrt(np.square(tw[rdt+'re_prediction'])
                                         + np.square(tw[rdt+'im_prediction']))

        out_columns += [rdt+'re', rdt+'im', rdt+'e',
                        rdt+'re_m', rdt+'im_m',
                        rdt+'W', rdt+'W_prediction',
                        rdt+'re_prediction', rdt+'im_prediction']
        
    tw['in_use'] = 1
    out_columns += ['in_use']
    write_tfs(join(meas_path, 'couple.out'), tw.loc[:, out_columns])


def _write_chromatic_coupling_files(meas_path, cor_path):
    LOG.debug("Calculating chromatic coupling diff.")
    # TODO: Add Cf1010
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
        tw['Cf1001r_prediction'] = tw.loc[:, 'Cf1001r'] - tw.loc[:, 'Cf1001r_model']
        tw['Cf1001i_prediction'] = tw.loc[:, 'Cf1001i'] - tw.loc[:, 'Cf1001i_model']
        write_tfs(join(meas_path, 'chromatic_coupling.out'),
                  tw.loc[:, ['NAME', 'S',
                             'Cf1001r', 'Cf1001rERR',
                             'Cf1001i', 'Cf1001iERR',
                             'Cf1001r_model', 'Cf1001i_model',
                             'Cf1001r_prediction', 'Cf1001i_prediction']])
    except IOError:
        LOG.debug("Chromatic coupling measurements not found. Skipped.")


# Script Mode ################################################################


if __name__ == "__main__":
    getdiff()
