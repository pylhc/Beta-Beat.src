import os
import sys

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis.extra.pandas import range_indexes, data_frames, column
from hypothesis.strategies import integers, floats, tuples

from twiss_optics.optics_class import TwissOptics
from twiss_optics.twiss_functions import get_all_rdts
from tfs_files.tfs_pandas import TfsDataFrame
from utils.contexts import suppress_warnings

sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
))

MAX_NRES = 40
CURRENT_DIR = os.path.dirname(__file__)

# Input Data #################################################################


def generic_column(name):
    return column(name,
                  elements=floats(allow_infinity=False, allow_nan=False),
                  dtype=np.float64)


def generic_column_pos(name):
    return column(name,
                  elements=floats(min_value=0,
                                  allow_infinity=False, allow_nan=False),
                  dtype=np.float64)


def s_column():
    return column("S",
                  elements=floats(min_value=0, max_value=27000,
                                  allow_infinity=False, allow_nan=False),
                  dtype=np.float64)


def alf_column(plane):
    return column("ALF" + plane,
                  elements=floats(allow_infinity=False, allow_nan=False),
                  dtype=np.float64)


def bet_column(plane):
    return column("BET" + plane,
                  elements=floats(min_value=1e-7,
                                  allow_infinity=False, allow_nan=False),
                  dtype=np.float64)


def mu_column(plane):
    return generic_column_pos("MU" + plane)


def d_column(plane):
    return generic_column("D" + plane)


def cmatrix_dataframes():
    df = data_frames(
        columns=[s_column(),
                 alf_column("X"), alf_column("Y"),
                 bet_column("X"), bet_column("Y"),
                 generic_column("R")],
        index=range_indexes(min_size=2, max_size=MAX_NRES)
    )
    return df


def full_dataframes():
    df = data_frames(
        columns=[s_column(),
                 bet_column("X"), bet_column("Y"),
                 mu_column("X"), mu_column("Y"),
                 d_column("X"), d_column("Y"),
                 generic_column("K0L"), generic_column("K0SL"),
                 generic_column("K1L"), generic_column("K1SL"),
                 generic_column("K2L"), generic_column("K2SL"),
                 generic_column("K3L"), generic_column("K3SL")],
        index=range_indexes(min_size=2, max_size=MAX_NRES)
    )
    return df


def tunes():
    return tuples(floats(min_value=0, max_value=100, allow_infinity=False, allow_nan=False),
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


@given(df=cmatrix_dataframes())
@settings(deadline=1000)
def test_calculate_cmatrix(df):
    df.loc[:, "R11"] = np.sin(df["R"])
    df.loc[:, "R22"] = df["R11"]
    df.loc[:, "R21"] = np.cos(df["R"])
    df.loc[:, "R12"] = -df["R21"]
    df = _pd_to_tfs(df, (0, 0))
    to = TwissOptics(df.copy())
    with suppress_warnings(RuntimeWarning):
        cmat = to.get_cmatrix()
        couple = to.get_coupling(method="cmatrix")
        gamma = to.get_gamma()

    assert all([c in cmat for c in ['S', 'C11', 'C12', 'C21', 'C22']])
    assert all([c in couple for c in ['S', 'F1001', 'F1010']])
    assert all([c in gamma for c in ['S', 'GAMMA']])
    assert len(cmat.index.values) == len(df.index.values)


@given(df=full_dataframes(),
       q=tunes(),
       order=integers(min_value=2, max_value=3))
@settings(deadline=1000)
def test_rdt_calculation(df, q, order):
    df = _pd_to_tfs(df, q)
    rdt_names = get_all_rdts(order)
    to = TwissOptics(df, quick_init=False)
    with suppress_warnings(RuntimeWarning):
        rdts = to.get_rdts(rdt_names)
    assert all([c in rdts for c in ["S"] + rdt_names])


@given(df=full_dataframes(), q=tunes())
@settings(deadline=1000)
def test_linear_dispersion(df, q):
    df = _pd_to_tfs(df, q)
    to = TwissOptics(df, quick_init=False)
    with suppress_warnings(RuntimeWarning):
        disp = to.get_linear_dispersion()
    assert all([c in disp for c in ["S", "DX", "DY"]])


@given(df=full_dataframes(), q=tunes())
def test_linear_chromaticity(df, q):
    df = _pd_to_tfs(df, q)
    to = TwissOptics(df, quick_init=False)
    with suppress_warnings(RuntimeWarning):
        lchrom = to.get_linear_chromaticity()
    assert len(lchrom) == 2


@given(df=full_dataframes(), q=tunes())
def test_chromatic_beating(df, q):
    df = _pd_to_tfs(df, q)
    to = TwissOptics(df, quick_init=False)
    with suppress_warnings(RuntimeWarning):
        chrom_beat = to.get_chromatic_beating()
    assert all([c in chrom_beat for c in ["S", "DBEATX", "DBEATY"]])


# Utilities ##################################################################


def _pd_to_tfs(df, q):
    df = TfsDataFrame(df)
    df.index = ["BPM{}".format(i) for i in df.index.values]
    df.headers["Q1"] = q[0]
    df.headers["Q2"] = q[1]
    return df
