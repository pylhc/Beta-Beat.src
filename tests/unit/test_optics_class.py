import os
import sys

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis.extra.pandas import range_indexes, columns, data_frames
from hypothesis.strategies import integers, floats, tuples

from twiss_optics.optics_class import TwissOptics
from twiss_optics.twiss_functions import get_all_rdts
from utils.tfs_pandas import TfsDataFrame

sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
))

MAX_NRES = 40
CURRENT_DIR = os.path.dirname(__file__)

# Input Data #################################################################

cmatrix_dataframes = data_frames(
    columns=columns(["S", "ALFX", "ALFY", "BETX", "BETY", "R11", "R12", "R21", "R22"],
                    elements=floats(allow_infinity=False, allow_nan=False),
                    dtype=np.float64),
    index=range_indexes(min_size=1, max_size=MAX_NRES)
)

full_dataframes = data_frames(
    columns=columns(["S", "BETX", "BETY", "MUX", "MUY", "DX", "DY",
                     "K0L", "K0SL", "K1L", "K1SL", "K2L", "K2SL", "K3L", "K3SL"],
                    elements=floats(allow_infinity=False, allow_nan=False),
                    dtype=np.float64),
    index=range_indexes(min_size=1, max_size=MAX_NRES)
)

tunes = tuples(floats(min_value=0, max_value=100, allow_infinity=False, allow_nan=False),
               floats(min_value=0, max_value=100, allow_infinity=False, allow_nan=False))


# Tests ######################################################################

def test_empty_dataframe():
    with pytest.raises(IndexError):
        TwissOptics(pd.DataFrame())


def test_no_index_given():
    df = pd.DataFrame(np.zeros([3,3]),
                      columns=["S", "BETX", "BETY"])
    with pytest.raises(IndexError):
        TwissOptics(df)


def test_import_from_file():
    twiss_file = os.path.join(CURRENT_DIR, "..", "inputs", "models", "flat_beam1", "twiss.dat")
    TwissOptics(twiss_file)


def test_missing_columns():
    df = pd.DataFrame(index=["BPM" + str(i) for i in range(10)],
                      columns=["S", "BETX", "BETY"])
    to = TwissOptics(df)
    with pytest.raises(KeyError):
        to.get_coupling()


@given(df=cmatrix_dataframes)
def test_calculate_cmatrix(df):
    df = _pd_to_tfs(df, (0, 0))
    to = TwissOptics(df.copy())
    cmat = to.get_cmatrix()
    couple = to.get_coupling(method="cmatrix")
    gamma = to.get_gamma()

    assert all([c in cmat for c in ['S', 'C11', 'C12', 'C21', 'C22']])
    assert all([c in couple for c in ['S', 'F1001', 'F1010']])
    assert all([c in gamma for c in ['S', 'GAMMA']])
    assert len(cmat.index.values) == len(df.index.values)


@given(df=full_dataframes,
       q=tunes,
       order=integers(min_value=2, max_value=4))
def test_rdt_calculation(df, q, order):
    df = _pd_to_tfs(df, q)
    rdt_names = get_all_rdts(order)
    to = TwissOptics(df, quick_init=False)
    to.get_rdts(rdt_names)


@given(df=full_dataframes, q=tunes)
def test_linear_dispersion(df, q):
    df = _pd_to_tfs(df, q)
    to = TwissOptics(df, quick_init=False)
    to.get_linear_dispersion()


@given(df=full_dataframes, q=tunes)
def test_linear_chromaticity(df, q):
    df = _pd_to_tfs(df, q)
    to = TwissOptics(df, quick_init=False)
    to.get_linear_chromaticity()


@given(df=full_dataframes, q=tunes)
def test_chromatic_beating(df, q):
    df = _pd_to_tfs(df, q)
    to = TwissOptics(df, quick_init=False)
    to.get_chromatic_beating()

# Utilities ##################################################################


def _pd_to_tfs(df, q):
    df = TfsDataFrame(df)
    df.index = ["BPM{}".format(i) for i in df.index.values]
    df.headers["Q1"] = q[0]
    df.headers["Q2"] = q[1]
    return df
