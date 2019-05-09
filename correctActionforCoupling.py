import numpy as np 
import pandas as pd 
from tfs_files import tfs_pandas
from utils import logging_tools, iotools
from utils.entrypoint import entrypoint, EntryPointParameters, ArgumentError

LOG = logging_tools.get_logger(__name__)

PLANES = ['X', 'Y']


def _get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--bpmfile",
        help="Directory where the SDDS file are stored",
        name="bpmfile",
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

    return params




def _calc_cosh2P( coupling_df ):

    f1010 = coupling_df['F1010R'] + 1j * coupling_df['F1010I']
    f1001 = coupling_df['F1001R'] + 1j * coupling_df['F1001I']

    twoP = np.sqrt( 
        (np.abs(2*f1010))**2 -
        (np.abs(2*f1001))**2 )

    cosh2P = np.cosh( twoP )

    return cosh2P


def _load_lin_files(filename, plane):

    lin_df = tfs_pandas.read_tfs( '{:}.lin{:}'.format( filename, plane.lower() ) )

    return lin_df


def _load_beta_files(analysisdir, plane):

    beta_df = tfs_pandas.read_tfs( 'getbeta{:}'.format( plane.lower() ) )

    return beta_df

@entrypoint( _get_params(), strict=True )
def correct_action_for_coupling(opt):
    LOG.info( "Starting the dark magic" )
    _calc_cosh2P(  )

    for plane in PLANES:

        # lin_df = _load_lin_files( , plane )

        # beta_df =  _load_beta_files(plane)

        frame = beta_df.join( [lin_df, cosh2P], how='inner' )
        # frame = pd.merge(beta_df, lin_df.loc[:, ['PK2PK']], how='inner', left_index=True, right_index=True)
        amps = frame.loc[:, 'PK2PK'].values / (2.0 * frame.loc[:, 'cosh2P'].values )
        meansqrt2j = amps / np.sqrt(frame.loc[:, 'BET' + plane].values)
        mean2j = np.square(amps) / frame.loc[:, 'BET' + plane].values

        return np.array([[np.mean(meansqrt2j), np.std(meansqrt2j)], [np.mean(mean2j), np.std(mean2j)]])

    LOG.info( "Starting the dark magic" )

    return coupling_df


# Script Mode ##################################################################


if __name__ == '__main__':
    correct_action_for_coupling()