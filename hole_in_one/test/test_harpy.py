import sys
import os
import pandas as pd
import numpy as np
from hypothesis import given
from hypothesis.strategies import integers
sys.path.append(os.path.abspath(".."))
import harpy


@given(integers(min_value=1, max_value=2000),
       integers(min_value=1, max_value=2000),
       integers(min_value=1, max_value=2000),
       integers(min_value=1, max_value=2000))
def test_highest_coefs_output_shape(nbpm, nres, target, toll):
    freqs = _get_fake_df(nbpm, nres)
    coefs = _get_fake_df(nbpm, nres)
    max_coefs, max_freqs = harpy._search_highest_coefs(
        target, toll, freqs, coefs
    )
    assert max_coefs.shape == (coefs.shape[0], )
    assert max_freqs.shape == (freqs.shape[0], )


@given(integers(min_value=1, max_value=2000),
       integers(min_value=1, max_value=2000),
       integers(min_value=1, max_value=2000),
       integers(min_value=1, max_value=2000))
def test_bpm_names_are_passed_correctly_to_output(nbpm, nres, target, toll):
    freqs = _get_fake_df(nbpm, nres)
    coefs = _get_fake_df(nbpm, nres)
    max_coefs, max_freqs = harpy._search_highest_coefs(
        target, toll, freqs, coefs
    )
    assert (max_coefs.index == coefs.index).all()
    assert (max_freqs.index == freqs.index).all()


@given(integers(min_value=1, max_value=2000))
def test_a_line_is_found_within_the_tollerance_window(nres):
    freqs = _get_fake_df(1, nres)
    coefs = _get_fake_df(1, nres)
    freqs.loc["BPM0", :] = np.arange(nres)
    coefs.loc[:, :] = 1
    max_coefs, max_freqs = harpy._search_highest_coefs(nres/2, 3, freqs, coefs)
    assert nres/2 - 3 <= max_freqs["BPM0"]
    assert max_freqs["BPM0"] <= nres/2 + 3


def test_search_highest_coefs():
    freqs = _get_fake_df(3, 10)
    freqs.loc["BPM0", :] = 0.
    freqs.loc["BPM1", :] = np.arange(10)
    freqs.loc["BPM2", :] = 10.
    coefs = _get_fake_df(3, 10)
    coefs.loc[:, :] = np.arange(10)
    coefs.loc["BPM1", 6] = -6  # Max negative value
    max_coefs, max_freqs = harpy._search_highest_coefs(4, 2, freqs, coefs)

    assert max_coefs.shape == (coefs.shape[0], )
    assert max_freqs.shape == (freqs.shape[0], )

    assert not max_coefs.loc["BPM0"]
    assert not max_coefs.loc["BPM2"]
    assert not max_freqs.loc["BPM0"]
    assert not max_freqs.loc["BPM2"]

    assert max_coefs.loc["BPM1"] == -6.
    assert max_freqs.loc["BPM1"] == 6.


def _get_fake_df(n_bpms, n_samples):
    index = ["BPM{}".format(i) for i in range(n_bpms)]
    df = pd.DataFrame(index=index, data=np.zeros((n_bpms, n_samples)))
    return df
