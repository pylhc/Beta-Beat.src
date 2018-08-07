from __future__ import print_function
import sys
import os
import pandas as pd
import numpy as np
sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
))
from harmonic_analysis import clean


def test_detect_known_bad_bpms():
    df = _get_fake_df(3, 3)

    cleaned = clean.detect_known_bad_bpms(df, [])
    assert not len(cleaned.values)

    cleaned = clean.detect_known_bad_bpms(df, ["BPM1"])
    assert len(cleaned.values) == 1
    assert cleaned[0] == "BPM1"

    cleaned = clean.detect_known_bad_bpms(df, ["BPM0", "BPM2"])
    assert len(cleaned.values) == 2
    assert cleaned[0] == "BPM0"
    assert cleaned[1] == "BPM2"


def test_detect_bpms_with_exact_zeros():
    df = _get_fake_df(3, 3)

    df.loc[:, :] = 1.
    assert clean.detect_bpms_with_exact_zeros(df).size == 0

    df.loc["BPM1", 1] = 0.
    assert clean.detect_bpms_with_exact_zeros(df).size == 1
    assert clean.detect_bpms_with_exact_zeros(df)[0] == "BPM1"


def test_detect_bpms_with_spikes():
    df = _get_fake_df(3, 3)
    df.loc["BPM1", 1] = 10.
    df.loc["BPM2", 0:1] = -10.

    cleaned = clean.detect_bpms_with_spikes(df, 5)
    assert len(cleaned.values) == 2
    assert cleaned[0] == "BPM1"
    assert cleaned[1] == "BPM2"


def test_detect_flat_bpms():
    df = _get_fake_df(3, 3)
    df.loc["BPM1", 1] = 10.
    df.loc["BPM2", 0:1] = -10.

    cleaned = clean.detect_flat_bpms(df, 5)
    assert len(cleaned.values) == 1
    assert cleaned[0] == "BPM0"


def test_index_union():
    indx0 = pd.Index([])
    indx1 = pd.Index("BPM1 BPM3".split())
    indx2 = pd.Index("BPM2".split())
    indx3 = pd.Index("BPM3".split())
    assert clean._index_union(indx0, indx0).size == 0

    union = clean._index_union(indx1, indx1)
    assert union.size == indx1.size
    assert union[0] == indx1[0]

    union = clean._index_union(indx1, indx2)
    assert union.size == indx1.size + indx2.size
    assert union.values[0] == "BPM1"
    assert union.values[1] == "BPM2"
    assert union.values[2] == "BPM3"

    union = clean._index_union(indx1, indx3)
    assert union.size == 2


def _get_fake_df(n_bpms, n_samples):
    index = ["BPM{}".format(i) for i in range(n_bpms)]
    df = pd.DataFrame(index=index, data=np.zeros((n_bpms, n_samples)))
    return df
