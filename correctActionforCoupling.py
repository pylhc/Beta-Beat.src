import numpy as np 
import pandas as pd 
from tfs_files import tfs_pandas

PLANES = ['X', 'Y']


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


def _load_beta_files( plane):

    beta_df = tfs_pandas.read_tfs( 'getbeta{:}'.format( plane.lower() ) )

    return beta_df

def correct_action_for_coupling():

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

    return 

coupl = tfs_pandas.read_tfs( '/afs/cern.ch/work/m/mihofer/public/Beta-Beat.src_mod_CouplingAction/TestData/2019-04-26/LHCB1/Results/nominal_injection/getcouple.out' )

print( _calc_cosh2P( coupling_df=coupl ) )
