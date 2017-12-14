import os
import re
import json
import numpy as np
import pandas as pd
import itertools
from math import factorial
import matplotlib.pyplot as plt
from Utilities.plotting import plot_style as pstyle
from Utilities import logging_tools as logtool
from Utilities import tfs_pandas as tfs
from Utilities.contexts import timeit
from twiss_optics import sequence_parser

LOG = logtool.get_logger(__name__)


def assertion(condition, exception):
    """ Raise Exception if condition is not fulfilled """
    if not condition:
        raise exception


def regex_in(regex, lst):
    """ Return boolean array of length lst, determining if that element starts with regex """
    return np.array([re.search(regex, element) is not None for element in lst])


def get_all_rdts(n):
    """ Returns list of all valid RDTs of order 2 to n """
    assertion(n > 1, ValueError("'n' must be greater 1 for resonance driving terms."))
    permut = [x for x in itertools.product(range(n + 1), repeat=4)
           if 1 < sum(x) <= n and not (x[0] == x[1] and x[2] == x[3])]
    return ['f{:d}{:d}{:d}{:d}'.format(j, k, l, m) for j, k, l, m in sorted(permut, key=sum)]


################################
#   Phase Advance Functions
################################

def get_phase_advances(twiss_df):
    """
    Calculate phase advances between all elements
    :return: Matrices similar to DPhi(i,j) = Phi(j) - Phi(i)
    """
    LOG.debug("Calculating Phase Advances:")
    phase_advance_dict = dict.fromkeys(['X', 'Y'])
    with timeit(lambda t:
                LOG.debug("  Phase Advances calculated in {:f}s".format(t))):
        for plane in ['X', 'Y']:
            colmn_phase = "MU" + plane

            phases_mdl = twiss_df.loc[twiss_df.index, colmn_phase]
            phase_advances = pd.DataFrame((phases_mdl[:, None] - phases_mdl[None, :]),
                                          index=twiss_df.index,
                                          columns=twiss_df.index)
            # Do not calculate dphi and tau here.
            # only slices of phase_advances as otherwise super slow
            phase_advance_dict[plane] = phase_advances
    return phase_advance_dict


def dphi(data, q):
    """ Return dphi from phase advances in data """
    return data + np.where(data <= 0, q, 0)  # '<=' seems to be what MAD-X does


def tau(data, q):
    """ Return tau from phase advances in data """
    return data + np.where(data <= 0, q / 2, -q / 2)  # '<=' seems to be what MAD-X does