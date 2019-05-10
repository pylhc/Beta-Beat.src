import numpy as np
import pandas as pd
from tfs_files import tfs_pandas
from utils import logging_tools, iotools
from utils.entrypoint import entrypoint, EntryPointParameters, ArgumentError
import os
import time
import calendar

LOG = logging_tools.get_logger(__name__)

PLANES = ['X', 'Y']

SUFFIX_TFS = '_corrCoup'


def _get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--kickfiles",
        help="Directory where the SDDS file are stored",
        name="kickfiles",
        type=str,
        required=True
    )
    params.add_parameter(
        flags="--analysisdir",
        help="Directory of the analysis",
        name="analysisdir",
        type=str,
        required=True
    )
    params.add_parameter(
        flags="--suffix",
        help="Suffix of file, _free or _free2",
        name="suffix",
        type=str,
        default=''
    )
    params.add_parameter(
        flags="--ac",
        help="Use of AC-dipole",
        name="ac",
        action="store_true"
    )

    return params


def _calc_cosh2P(coupling_df):

    f1010 = coupling_df['F1010R'] + 1j * coupling_df['F1010I']
    f1001 = coupling_df['F1001R'] + 1j * coupling_df['F1001I']

    twoP = np.sqrt(np.abs(
        (np.abs(2*f1010))**2 -
        (np.abs(2*f1001))**2))

    cosh2P = np.cosh(twoP)

    coupling_df['COSH2P'] = cosh2P

    return coupling_df


def _load_lin_files(lindir, plane):

    filename = os.path.split(lindir)[-1]
    lin_df = tfs_pandas.read_tfs(os.path.join(
        lindir, '{:}.sdds.lin{:}'.format(filename, plane.lower())), index='NAME')

    return lin_df


def _load_beta_files(opt, plane):

    beta_df = tfs_pandas.read_tfs(os.path.join(
        opt.analysisdir, 'getbeta{:}{:}.out'.format(plane.lower(), opt.suffix)), index='NAME')

    return beta_df


def _get_time(lindir):
    
    filename = os.path.split(lindir)[-1]
    kickdate, kicktime = filename.split('@')[-2:]
    kick_timestamp = time.strptime(
        kickdate+' '+kicktime, '%Y_%m_%d %H_%M_%S_%f')
    utc_timestamp = calendar.timegm(kick_timestamp)
    return float(utc_timestamp)


def _calc_corrected_action(opt, lindir, coupling_df, corr_kick_df):
    
    timestamp = _get_time(lindir)
    kicks = np.zeros((8))

    for i, plane in enumerate(PLANES):

        lin_df = _load_lin_files(lindir, plane)
        beta_df = _load_beta_files(opt, plane)

        if opt.ac:
            frame = pd.merge(beta_df, lin_df.loc[:, [
                             'AMP'+plane.lower()]], how='inner', left_index=True, right_index=True)
            frame = pd.merge(frame, coupling_df.loc[:, [
                             'COSH2P']], how='inner', left_index=True, right_index=True)
            amps = 2.0 * \
                frame.loc[:, 'AMP'+plane.lower()].values / \
                (frame.loc[:, 'COSH2P'].values)
        else:
            frame = pd.merge(beta_df, lin_df.loc[:, [
                             'PK2PK']], how='inner', left_index=True, right_index=True)
            frame = pd.merge(frame, coupling_df.loc[:, [
                             'COSH2P']], how='inner', left_index=True, right_index=True)
            amps = frame.loc[:, 'PK2PK'].values / \
                (2.0 * frame.loc[:, 'COSH2P'].values)

        meansqrt2j = amps / np.sqrt(frame.loc[:, 'BET' + plane + 'MDL'].values)
        mean2j = np.square(amps) / frame.loc[:, 'BET' + plane + 'MDL'].values
        kicks[i*4:(i+1)*4] = np.mean(meansqrt2j), np.std(
            meansqrt2j), np.mean(mean2j), np.std(mean2j)

    corr_kick_df.loc[timestamp] = kicks

    return corr_kick_df


def _get_kickfiles(opt):
    kickfiles = opt.kickfiles.split(',')
    return kickfiles


@entrypoint(_get_params(), strict=True)
def correct_action_for_coupling(opt):
    LOG.info("Starting the dark magic")

    coupling_df = tfs_pandas.read_tfs(os.path.join(
        opt.analysisdir, 'getcouple{:}.out'.format(opt.suffix)), index='NAME')
    coupling_df = _calc_cosh2P(coupling_df)

    kick_df = tfs_pandas.read_tfs(os.path.join(
        opt.analysisdir, 'getkick{:}.out'.format('ac' if opt.ac else '')), index='TIME')

    column_names = [s + SUFFIX_TFS
                    for s in ["sqrt2JX", "sqrt2JXSTD", "2JX", "2JXSTD", "sqrt2JY", "sqrt2JYSTD", "2JY", "2JYSTD"]]

    corr_kicks_df = pd.DataFrame(columns=column_names)

    kickfiles = _get_kickfiles(opt)

    for kick in kickfiles:
        corr_kick_df = _calc_corrected_action(
            opt, kick, coupling_df, corr_kicks_df)

    kick_df = kick_df.join( corr_kick_df )
    
    LOG.info("Finishing the dark magic")

    return coupling_df, kick_df


# Script Mode ##################################################################


if __name__ == '__main__':
    correct_action_for_coupling()
