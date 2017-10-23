"""
This module is the actual implementation of the resonance search for
turn-by-turn data.
It uses a combination of the laskar method with SVD decomposition to speed up
the search of resonances.

TODOs:
- Mising nattunez - not needed yet
- Improve search of higher order resonances smaller tolerance
- Other option of noise estimate is to take average of to take average
   amplitude of fft of few last taken or first not taken modes
"""
from __future__ import print_function
import multiprocessing
import logging
from functools import partial
import numpy as np
import pandas as pd
from Utilities import outliers

try:
    from scipy.fftpack import fft as _fft
except ImportError:
    from numpy.fft import fft as _fft

NUMBA_AVAILABLE = True
try:
    from numba import jit
except ImportError:
    NUMBA_AVAILABLE = False


LOGGER = logging.getLogger(__name__)

PI2I = 2 * np.pi * complex(0, 1)

RESONANCE_LISTS = {
    "X": ((0, 1, 0), (-2, 0, 0), (0, 2, 0), (-3, 0, 0), (-1, -1, 0),
          (2, -2, 0), (0, -2, 0), (1, -2, 0), (-1, 3, 0), (1, 2, 0),
          (-2, 1, 0), (1, 1, 0), (2, 0, 0), (-1, -2, 0), (3, 0, 0)),
    "Y": ((1, 0, 0), (-1, 1, 0), (-2, 0, 0), (1, -1, 0), (0, -2, 0),
          (0, -3, 0), (2, 1, 0), (-1, 3, 0), (1, 1, 0), (-1, 2, 0)),
    "Z": ((1, 0, 1), (0, 1, 1))
}

MAIN_LINES = {"X": (1, 0, 0), "Y": (0, 1, 0), "Z": (0, 0, 1)}
Z_TOLERANCE = 0.0001
NUM_HARMS = 300
NUM_HARMS_SVD = 100

PROCESSES = multiprocessing.cpu_count()


def harpy(bpm_matrix, usv, plane, harpy_input, panda, model_tfs):
    """Searches for the strongest resonances in each row of the bpm_matrix.

    TODO: More elaborate description here.

    Args:
        bpm_matrix: Pandas DataFrame that containing the samples with shape
            (n_bpms, n_turns) and the BPM names as row index.
        usv: SVD decomposition of bpm_matrix, shoud be an iterable such that
            u, s, v = usv is possible.
        plane: A string containing either X or Y.
        harpy_input: HarpyInput instance containing user provided info.
        panda: A pandas DataFrame where the resonance results will be put.
        model_tfs: A pandas DataFrame to get the model dispersion.

    Returns:
        A tuple containing the resulting DataFrame object filled with the
            resonance information (Representing the _linx or _liny files),
            a dictionary of DataFrames with form: bpm_name -> bpm_spectrum and
            a list of string descibing which BPMs are marked as bad.
    """
    frequencies, bpm_coefficients = harmonic_analysis(
        bpm_matrix,
        usv=usv,
        mode=harpy_input.harpy_mode,
        sequential=harpy_input.sequential,
    )

    all_bpms_coefs = pd.DataFrame(data=bpm_coefficients,
                                  index=bpm_matrix.index)
    all_bpms_freqs = pd.DataFrame(data=frequencies,
                                  index=bpm_matrix.index)
    all_bpms_spectr = {"COEFS": all_bpms_coefs, "FREQS": all_bpms_freqs}

    panda, not_tune_bpms = _get_main_resonances(
        harpy_input, frequencies, bpm_coefficients,
        plane, harpy_input.tolerance, panda
    )
    if harpy_input.tunez > 0.0:
        panda, _, _ = _get_main_resonances(
            harpy_input, frequencies, bpm_coefficients,
            "Z", Z_TOLERANCE, panda
        )

    cleaned_by_tune_bpms = clean_by_tune(panda.loc[:, 'TUNE' + plane],
                                         harpy_input.tune_clean_limit)
    panda = panda.loc[panda.index.difference(cleaned_by_tune_bpms)]
    bpm_matrix = bpm_matrix.loc[panda.index]
    bpm_coefficients = bpm_coefficients.loc[panda.index]
    frequencies = frequencies.loc[panda.index]

    bad_bpms_summary = _get_bad_bpms_summary(not_tune_bpms, cleaned_by_tune_bpms)
    panda = _amp_and_mu_from_avg(bpm_matrix, plane, panda)
    # TODO: I removed the old noise!
    panda = _get_natural_tunes(frequencies, bpm_coefficients, harpy_input, plane, panda)

    resonances_freqs = _compute_resonances_freqs(plane, harpy_input)
    if harpy_input.tunez > 0.0:
        resonances_freqs.update(_compute_resonances_freqs("Z", harpy_input))
 
    panda = _resonance_search(frequencies, bpm_coefficients,
                              harpy_input.tolerance, resonances_freqs, panda)

    panda['BINDEX'] = 0
    panda['SLABEL'] = 0
    return (panda,
            all_bpms_spectr,
            bad_bpms_summary)


def harmonic_analysis(bpm_matrix=None, usv=None, mode="bpm", sequential=False):
    """
    Performs the laskar method on every of the BPMs signals contained in each
    row of the bpm_matrix pandas DataFame.

    Args:
        bpm_matrix: Pandas DataFrame containing the signals of each bpm in
            each row. In bpm mode this parameter is obligatory.
        usv: Truncated SVD decomposition of the bpm matrix. It must contain a
            tuple (U, S, V), where U must be a DataFrame with the bpm names
            as index.
        mode: one of 'bpm', 'svd', or 'fast'. Check 'harmonic_analysis_bpm'
            documentation for 'bpm' mode and 'harmonic_analysis_svd" for 'svd'
            and 'fast'.
        sequential: If true, it will run all the computations in a single
            core.
    Returns:
        frequencies: A numpy array with the frequencies found per BPM.
        bpm_coefficients: A numpy array containing the complex coefficients
            found per BPM.
    """
    if mode == "bpm":
        if bpm_matrix is None:
            raise ValueError("bpm_matrix has to be provided "
                             "for the bpm mode")
        frequencies, bpm_coefficients = harmonic_analysis_bpm(
            bpm_matrix,
            sequential=sequential,
            num_harms=NUM_HARMS,
        )
    elif mode in ("svd", "fast"):
        if usv is None:
            raise ValueError("SVD decomposition has to be provided "
                             "for fast or svd modes")
        num_harms = NUM_HARMS if mode == "svd" else NUM_HARMS_SVD
        frequencies, bpm_coefficients = harmonic_analysis_svd(
            usv,
            fast=mode == "fast",
            sequential=sequential,
            num_harms=num_harms,
        )
    else:
        raise ValueError("Invalid harpy mode: {}".format(mode))
    return frequencies, bpm_coefficients


def harmonic_analysis_bpm(bpm_matrix,
                          sequential=False, num_harms=NUM_HARMS):
    """
    Performs the laskar method on every of the BPMs signals contained in each
    row of the bpm_matrix pandas DataFame. This method will run the full
    laskar analysis in each BPM in parallel (slow).

    Args:
        bpm_matrix: Pandas DataFrame containing the signals of each bpm in
            each row. In bpm mode this parameter is obligatory.
        sequential: If true, it will run all the computations in a single
            core.
        num_harms: Number of harmonics to compute per BPM.
    Returns:
        frequencies: A numpy array with the frequencies found per BPM.
        bpm_coefficients: A numpy array containing the complex coefficients
            found per BPM.
    """
    frequencies, bpm_coefficients = _parallel_laskar(
        bpm_matrix, sequential, num_harms,
    )
    return frequencies, bpm_coefficients


def harmonic_analysis_svd(usv, fast=False,
                          sequential=False, num_harms=NUM_HARMS):
    """
    Performs the laskar method on every of the BPMs signals contained in each
    row of the bpm_matrix pandas DataFame. It takes advantage of the
    truncation of the V matrix in the SVD decomposition to speed up the
    calculations. Mode 'svd' is slower but more precise, 'bpm' mode is the
    fastest but makes more assumptions.

    Args:
        usv: Truncated SVD decomposition of the bpm matrix. It must contain a
            tuple (U, S, V), where U must be a DataFrame with the bpm names
            as index.
        fast: If true, it will use the 'fast' calculation mode ('svd'
            otherwise)
        sequential: If true, it will run all the computations in a single
            core.
        num_harms: Number of harmonics to compute per BPM.
    Returns:
        frequencies: A numpy array with the frequencies found per BPM.
        bpm_coefficients: A numpy array containing the complex coefficients
            found per BPM.
    """
    sv = np.dot(np.diag(usv[1]), usv[2])
    fake_bpms = ["MODE{}".format(i) for i in range(sv.shape[0])]
    sv = pd.DataFrame(index=fake_bpms, data=sv)
    if fast:
        frequencies, _ = _laskar_per_mode(np.mean(sv, axis=0), num_harms)
    else:
        frequencies, _ = _parallel_laskar(
            sv, sequential, num_harms,
        )
        frequencies = np.ravel(frequencies)
    svd_coefficients = _compute_coefs_for_freqs(sv, frequencies)
    bpm_coefficients = usv[0].dot(svd_coefficients.values)
    frequencies = pd.DataFrame(
        index=bpm_coefficients.index,
        data=np.outer(np.ones(bpm_coefficients.shape[0]), frequencies)
    )
    return frequencies, bpm_coefficients


def clean_by_tune(tunes, tune_clean_limit):
    """
    This function looks for outliers in the tunes pandas Series and returns
    their indices.

    Args:
        tunes: Pandas series with the tunes per BPM and the BPM names as
            index.
        tune_clean_limit: No BPM will find as oulier if its distance to the
            average is lower than this limit.
    """
    bad_bpms_mask = outliers.get_filter_mask(
        tunes,
        limit=tune_clean_limit,
    )
    bad_bpms_names = tunes[~bad_bpms_mask].index
    return bad_bpms_names


def _get_bad_bpms_summary(not_tune_bpms, cleaned_by_tune_bpms):
    bad_bpms_summary = []
    for bpm_name in not_tune_bpms:
        bad_bpms_summary.append(bpm_name +
                                " The main resonance has not been found.")
    for bpm_name in cleaned_by_tune_bpms:
        bad_bpms_summary.append(
            "{} tune is too far from average"
            .format(bpm_name)
        )
    return bad_bpms_summary


def _parallel_laskar(samples, sequential, num_harms):
    freqs = pd.DataFrame(
        index=samples.index,
        data=np.zeros((samples.shape[0], num_harms), dtype=np.float)
    )
    coefs = pd.DataFrame(
        index=samples.index,
        data=np.zeros((samples.shape[0], num_harms), dtype=np.complex128)
    )

    def _collect_results(bpm_name, freq_and_coef):
        freq, coef = freq_and_coef
        freqs.loc[bpm_name] = freq
        coefs.loc[bpm_name] = coef

    pool = multiprocessing.Pool(np.min([PROCESSES, samples.shape[0]]))
    for bpm_name in samples.index:
        args = (samples.loc[bpm_name, :], num_harms)
        callback = partial(_collect_results, bpm_name)
        if sequential:
            callback(_laskar_per_mode(*args))
        else:
            pool.apply_async(_laskar_per_mode, args,
                             callback=callback)
    pool.close()
    pool.join()
    return freqs, coefs


def _laskar_per_mode(samples, number_of_harmonics):
    freqs, coefs = laskar_method(samples.values, number_of_harmonics)
    return np.array(freqs), np.array(coefs)


def _compute_coefs_for_freqs(samples, freqs):  # Samples is a matrix
    n = samples.shape[1]
    coefficients = samples.dot(
        np.exp(-PI2I * np.outer(np.arange(n), freqs))
    ) / n
    return coefficients


def _get_main_resonances(harpy_input, frequencies, coefficients,
                         plane, tolerance, panda):
    h, v, l = MAIN_LINES[plane]
    freq = ((h * harpy_input.tunex) +
            (v * harpy_input.tuney) +
            (l * harpy_input.tunez))
    max_coefs, max_freqs = _search_highest_coefs(
        freq, tolerance, frequencies, coefficients
    )
    if not np.any(max_coefs):
        raise ValueError(
            "No main " + plane + " resonances found, "
            "try to increase the tolerance or adjust the tunes"
        )
    bad_bpms_by_tune = coefficients.loc[max_coefs == 0.].index
    panda = panda.loc[panda.index.difference(bad_bpms_by_tune)]
    panda['TUNE'+plane] = max_freqs
    panda['AMP'+plane] = np.abs(max_coefs)
    panda['MU'+plane] = np.angle(max_coefs) / (2 * np.pi)
    return panda, bad_bpms_by_tune


def _search_highest_coefs(freq, tolerance, frequencies, coefficients):
    min = freq - tolerance
    max = freq + tolerance
    on_window_mask = (frequencies.values >= min) & (frequencies.values <= max)
    filtered_amps = coefficients.abs().where(on_window_mask, other=0.)
    max_indices = filtered_amps.idxmax(axis=1)
    max_coefs = coefficients.lookup(max_indices.index, max_indices.values)
    # TODO: Better solution?:
    max_coefs = pd.Series(index=coefficients.index, data=max_coefs)
    max_freqs = frequencies.lookup(max_indices.index, max_indices.values)
    max_freqs = (max_coefs != 0) * max_freqs  # Zero bad freqs
    return max_coefs, max_freqs


def _amp_and_mu_from_avg(bpm_matrix, plane, panda):
    avg_coefs = _compute_coefs_for_freqs(
        bpm_matrix,
        np.mean(panda.loc[:, 'TUNE' + plane])
    )
    panda['AVG_AMP'+plane] = np.abs(avg_coefs)
    panda['AVG_MU'+plane] = np.angle(avg_coefs) / (2 * np.pi)
    return panda


def _get_natural_tunes(frequencies, coefficients, harpy_input, plane, panda):
    if harpy_input.nattunex is None or harpy_input.nattuney is None:
        return panda
    if plane == "X":
        freq = harpy_input.nattunex
    if plane == "Y":
        freq = harpy_input.nattunex
    max_coefs, max_freqs = _search_highest_coefs(
        freq, harpy_input.tolerance, frequencies, coefficients
    )
    panda['NATTUNE'+plane] = max_freqs
    panda['NATAMP'+plane] = np.abs(max_coefs)
    return panda


# Vectorized - coefficiens are matrix:
def _resonance_search(frequencies, coefficients, tolerance,
                      resonances_freqs, panda):
    results = panda
    for resonance, resonance_freq in resonances_freqs.iteritems():
        max_coefs, max_freqs = _search_highest_coefs(
            resonance_freq, tolerance, frequencies, coefficients
        )
        resstr = _get_resonance_suffix(resonance)
        results['AMP' + resstr] = np.abs(max_coefs)
        results['PHASE' + resstr] = np.angle(max_coefs) / (2 * np.pi)
        results['FREQ' + resstr] = max_freqs
    return results


def _get_resonance_suffix(resonance):
    x, y, z = resonance
    if z == 0:
        resstr = (str(x) + str(y)).replace("-", "_")
    else:
        resstr = (str(x) + str(y) + str(z)).replace("-", "_")
    return resstr


def _compute_resonances_freqs(plane, harpy_input):
    """
    Computes the frequencies for all the resonances listed in the
    constante RESONANCE_LISTS, together with the natural tunes
    frequencies if given.
    """
    resonances_freqs = {}
    freqs = [(resonance_h * harpy_input.tunex) +
             (resonance_v * harpy_input.tuney) +
             (resonance_l * harpy_input.tunez)
             for (resonance_h, resonance_v, resonance_l) in RESONANCE_LISTS[plane]]
    # Move to [0, 1] domain.
    freqs = [freq + 1. if freq < 0. else freq for freq in freqs]
    resonances_freqs = dict(zip(RESONANCE_LISTS[plane], freqs))
    return resonances_freqs


def laskar_method(tbt, num_harmonics):
    samples = tbt[:]  # Copy the samples array.
    n = len(samples)
    int_range = np.arange(n)
    coefficients = []
    frequencies = []
    for _ in range(num_harmonics):
        # Compute this harmonic frequency and coefficient.
        dft_data = _fft(samples)
        frequency = _jacobsen(dft_data,n)
        coefficient = _compute_coef(samples, frequency * n) / n

        # Store frequency and amplitude
        coefficients.append(coefficient)
        frequencies.append(frequency)

        # Subtract the found pure tune from the signal
        new_signal = coefficient * np.exp(PI2I * frequency * int_range)
        samples = samples - new_signal

    coefficients, frequencies = zip(*sorted(zip(coefficients, frequencies),
                                            key=lambda tuple: np.abs(tuple[0]),
                                            reverse=True))
    return frequencies, coefficients


def fft_method(tbt, num_harmonics):
    samples = tbt[:]  # Copy the samples array.
    n = float(len(samples))
    dft_data = _fft(samples)
    indices = np.argsort(np.abs(dft_data))[::-1][:num_harmonics]
    frequencies = indices / n
    coefficients = dft_data[indices] / n
    return frequencies, coefficients


def _jacobsen(dft_values,n):
    """
    This method interpolates the real frequency of the
    signal using the three highest peaks in the FFT.
    """
    k = np.argmax(np.abs(dft_values))
    r = dft_values
    delta = np.tan(np.pi / n) / (np.pi / n)
    kp = (k + 1) % n
    km = (k - 1) % n
    delta = delta * np.real((r[km] - r[kp]) / (2 * r[k] - r[km] - r[kp]))
    return (k + delta) / n
 

def _compute_coef_simple(samples, kprime):
    """
    Computes the coefficient of the Discrete Time Fourier
    Transform corresponding to the given frequency (kprime).
    """
    n = len(samples)
    freq = kprime / n
    exponents = np.exp(-PI2I * freq * np.arange(n))
    coef = np.sum(exponents * samples)
    return coef


def _compute_coef_goertzel(samples, kprime):
    """
    Computes the coefficient of the Discrete Time Fourier
    Transform corresponding to the given frequency (kprime).
    This function is faster than the previous one if compiled
    with Numba.
    """
    n = len(samples)
    a = 2 * np.pi * (kprime / n)
    b = 2 * np.cos(a)
    c = np.exp(-complex(0, 1) * a)
    d = np.exp(-complex(0, 1) * ((2 * np.pi * kprime) / n) * (n - 1))
    s0 = 0.
    s1 = 0.
    s2 = 0.
    for i in range(n - 1):
        s0 = samples[i] + b * s1 - s2
        s2 = s1
        s1 = s0
    s0 = samples[n - 1] + b * s1 - s2
    y = s0 - s1 * c
    return y * d


def _which_compute_coef():
    if not NUMBA_AVAILABLE:
        LOGGER.warn("Cannot import numba, using numpy functions.")
        return _compute_coef_simple
    LOGGER.debug("Using compiled Numba functions.")
    return jit(_compute_coef_goertzel, nopython=True, nogil=True)

_compute_coef = _which_compute_coef()
