"""
This module is the actual implementation of the resonance search for
turn-by-turn data.
It uses a combination of the laskar method with SVD decomposition to speed up
the search of resonances.
The only public function and entry point is svd_harpy(...).

TODOs:
- Mising nattunez - not needed yet
- Improve search of higher order resonances smaller tolerance
- If non-zero tunez try chroma_calculation
- The amplitudes are not rescaled by ampx/y
- Add header info, tunes, nattunes, dpp
- Separate resonance for non-zero tunez computed if tunez > 0
- Other option of noise estimate is to take average of to take average
   amplitude of fft of few last taken or first not taken modes
"""
from __future__ import print_function
import multiprocessing
import logging
import numpy as np
import pandas as pd
from harpy.harmonic_analysis import HarmonicAnalysis

from harpy import _python_path_manager
_python_path_manager.append_betabeat()
from Utilities import outliers
from model import manager

LOGGER = logging.getLogger(__name__)

PI2I = 2 * np.pi * complex(0, 1)

SPECTR_COLUMN_NAMES = ["FREQ", "AMP"]

RESONANCE_LISTS = {
    "X": ((0, 1, 0), (-2, 0, 0), (0, 2, 0), (-3, 0, 0), (-1, -1, 0),
          (2, -2, 0), (0, -2, 0), (1, -2, 0), (-1, 3, 0), (1, 2, 0),
          (-2, 1, 0), (1, 1, 0), (2, 0, 0), (-1, -2, 0), (3, 0, 0)),
    "Y": ((1, 0, 0), (-1, 1, 0), (-2, 0, 0), (1, -1, 0), (0, -2, 0),
          (0, -3, 0), (2, 1, 0), (-1, 3, 0), (1, 1, 0), (-1, 2, 0)),
    "Z": ((1, 0, 1), (0, 1, 1))
}

MAIN_LINES = {"X": (1, 0, 0),
              "Y": (0, 1, 0),
              "Z": (0, 0, 1)}
Z_TOLERANCE = 0.0001

DEFAULT_DFT_METHOD = "laskar"

NUM_HARMS = 300
NUM_HARMS_SVD = 100

PROCESSES = multiprocessing.cpu_count()


def svd_harpy(bpm_names, bpm_matrix, usv, plane,
              harpy_input, panda, model_tfs):
    """Searches for the strongest resonances in each row of the bpm_matrix.

    TODO: More elaborate description here.

    Args:
        bpm_names: List of bpm names with the same order as bpm_matrix.
        bpm_matrix: Matrix containing the samples with shape (n_bpms, n_turns).
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
    frequencies, bpm_coefficients = _harmonic_analysis(
        harpy_input, bpm_matrix, usv
    )
    if frequencies.ndim < 2:
        frequencies = np.outer(np.ones(bpm_coefficients.shape[0]), frequencies)
    all_bpms_coefs = pd.DataFrame(
        data=bpm_coefficients,
        index=bpm_names
    )
    all_bpms_freqs = pd.DataFrame(
        data=frequencies,
        index=bpm_names
    )
    all_bpms_spectr = {"COEFS": all_bpms_coefs, "FREQS": all_bpms_freqs}

    panda, bad_bpms_mask = _get_main_resonances(
        harpy_input, frequencies, bpm_coefficients,
        plane, harpy_input.tolerance, panda
    )
    if harpy_input.tunez > 0.0:
        panda, bad_bpms_mask = _get_main_resonances(
            harpy_input, frequencies, bpm_coefficients,
            "Z", Z_TOLERANCE, panda
        )
    panda, bad_bpms, bad_bpms_mask = _clean_by_tune(harpy_input, plane, panda)
    panda = _amp_and_mu_from_avg(bpm_matrix, bad_bpms_mask, plane, panda)
    panda = _get_noise(bpm_matrix, panda)
    panda = _get_natural_tunes(frequencies, bpm_coefficients,
                               harpy_input, plane, panda)

    resonances_freqs = _compute_resonances_freqs(plane, harpy_input)
    if harpy_input.tunez > 0.0:
        resonances_freqs.update(_compute_resonances_freqs("Z", harpy_input))
        if plane == "X":
            dpoverp_amp, pos, neg = get_dpoverp_amp(
                model_tfs, panda, bpm_names, bpm_matrix,
            )
    panda = _resonance_search(frequencies, bpm_coefficients,
                              harpy_input.tolerance, resonances_freqs, panda)

    panda['BINDEX'] = 0
    panda['SLABEL'] = 0
    return (panda.loc[bad_bpms_mask, :],
            all_bpms_spectr,
            bad_bpms)


def _harmonic_analysis(harpy_input, bpm_matrix, usv):
    sv = np.dot(np.diag(usv[1]), usv[2])
    if harpy_input.harpy_mode == "bpm":
        frequencies, bpm_coefficients = _parallel_laskar(
            bpm_matrix, harpy_input.sequential, NUM_HARMS,
        )
    elif harpy_input.harpy_mode == "fast":
        frequencies, _ = _laskar_per_mode(np.mean(sv, axis=0), NUM_HARMS)
        svd_coefficients = _compute_coefs_for_freqs(sv, frequencies)
        bpm_coefficients = np.dot(usv[0], svd_coefficients)
    elif harpy_input.harpy_mode == "svd":
        frequencies, _ = _parallel_laskar(
            sv, harpy_input.sequential, NUM_HARMS_SVD
        )
        frequencies = np.ravel(frequencies)
        svd_coefficients = _compute_coefs_for_freqs(sv, frequencies)
        bpm_coefficients = np.dot(usv[0], svd_coefficients)
    else:
        raise ValueError("Wrong harpy mode: " + harpy_input.harpy_mode)
    return frequencies, bpm_coefficients


def _parallel_laskar(samples, sequential, num_harms):
    freqs = []
    coefs = []

    def _collect_results(freq, coef):
        freqs.append(freq)
        coefs.append(coef)

    pool = multiprocessing.Pool(np.min([PROCESSES, samples.shape[0]]))
    for i in range(samples.shape[0]):
        args = (samples[i, :], num_harms)
        if sequential:
            _collect_results(*_laskar_per_mode(*args))
        else:
            pool.apply_async(_laskar_per_mode, args,
                             callback=lambda res: _collect_results(*res))
    pool.close()
    pool.join()
    return np.array(freqs), np.array(coefs)


def _laskar_per_mode(samples, number_of_harmonics):
    har_analysis = HarmonicAnalysis(samples)
    freqs, coefs = har_analysis.laskar_method(number_of_harmonics)
    return np.array(freqs), np.array(coefs)


def _compute_coefs_for_freqs(samples, freqs):  # Samples is a matrix
    n = samples.shape[1]
    coefficients = np.dot(
        samples,
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
    bad_bpms_mask = max_coefs != 0.
    panda['TUNE'+plane] = max_freqs
    panda['AMP'+plane] = np.abs(max_coefs)
    panda['MU'+plane] = np.angle(max_coefs) / (2 * np.pi)
    return panda, bad_bpms_mask


def _search_highest_coefs(freq, tolerance, frequencies, coefficients):
    min = freq - tolerance
    max = freq + tolerance
    on_window_mask = (frequencies >= min) & (frequencies <= max)
    filtered_amps = np.where(on_window_mask,
                              np.abs(coefficients),
                              np.zeros_like(coefficients))
    filtered_coefs = np.where(on_window_mask,
                              coefficients,
                              np.zeros_like(coefficients))
    max_indices = np.argmax(filtered_amps, axis=1)
    max_coefs = filtered_coefs[np.arange(coefficients.shape[0]), max_indices]
    max_freqs = frequencies[np.arange(frequencies.shape[0]), max_indices]
    max_freqs = (max_coefs != 0) * max_freqs  # Zero bad freqs
    return max_coefs, max_freqs


def _clean_by_tune(harpy_input, plane, panda):
    bad_bpms = []
    bad_bpms_mask = outliers.get_filter_mask(
        panda.loc[:, 'TUNE' + plane],
        limit=harpy_input.tune_clean_limit,
    )
    for i in np.arange(len(bad_bpms_mask)):
        if not bad_bpms_mask[i]:
            bad_bpms.append(panda.at[i, 'NAME'] +
                            " tune is too far from average")
    return panda, bad_bpms, bad_bpms_mask


def _amp_and_mu_from_avg(bpm_matrix, bad_bpm_mask, plane, panda):
    avg_coefs = _compute_coefs_for_freqs(
        bpm_matrix,
        np.mean(panda.loc[bad_bpm_mask, 'TUNE' + plane])
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


def _get_noise(bpm_matrix, panda):
    result = panda
    leng = 2 ** (int(np.log2(bpm_matrix.shape[1])) - 2)
    leng = 512 if leng > 512 else leng
    pos = int(leng / 64)
    noise = np.empty([bpm_matrix.shape[0],
                      int(np.floor(bpm_matrix.shape[1] / leng))])
    for i in range(bpm_matrix.shape[1] / leng):
        noise[:, i] = np.partition(
            np.abs(np.fft.fft(bpm_matrix[:, :leng])), pos, axis=1
        )[:, pos] / leng
    result['NOISE'] = np.mean(noise, axis=1)
    return result


# TODO: The amplitudes are not rescaled by ampx/y
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


# # full period is one, phase is of the arc bpms in horizontal plane
def get_dpoverp_amp(model_tfs, panda, bpm_names, bpm_samples):
    # Interval is around integer phases for positive dpoverp
    # and around halfes for the negative dpoverp
    arc_bpm_names, arc_bpm_samples = _get_lhc_arc_bpms(bpm_names, bpm_samples)
    tunez = np.mean(panda.set_index("NAME").loc[arc_bpm_names, "TUNEZ"])
    tunez_phase = _phase_mean(panda.set_index("NAME").loc[arc_bpm_names, "MUZ"])
    turns = bpm_samples.shape[1]
    pos = _get_positive_dpoverp_intervals(tunez, tunez_phase, turns)
    neg = _get_negative_dpoverp_intervals(tunez, tunez_phase, turns)
    # We assume 3D kicks only happen in LHC for now.
    co = np.zeros_like(arc_bpm_names)
    model_dx = model_tfs.set_index("NAME").loc[arc_bpm_names, "DX"]
    length = 0
    for bin in pos:
        mini, maxi = bin
        co = co + np.sum(arc_bpm_samples[:, mini:maxi], axis=1)
        length = length + maxi - mini
    for bin in neg:
        mini, maxi = bin
        co = co - np.sum(arc_bpm_samples[:, mini:maxi], axis=1)
        length = length + maxi - mini
    co = co / 1e3  # Going from mm to m.
    codx = co * model_dx / length
    dx2 = model_dx ** 2
    return np.sum(codx) / np.sum(dx2), pos, neg


def _get_lhc_arc_bpms(bpm_names, bpm_samples):
    accel_cls = manager.get_accel_class("lhc")
    arc_bpms_mask = accel_cls.get_arc_bpms_mask(bpm_names)
    arc_bpm_samples = bpm_samples[arc_bpms_mask]
    arc_bpm_names = bpm_names[arc_bpms_mask]
    return arc_bpm_names, arc_bpm_samples


def _phase_mean(phases):
    return np.angle(np.sum(np.exp(PI2I * phases))) / (2 * np.pi)


def get_chroma(harpy_input, plane, bpm_names, bpm_samples, panda):
    arc_bpm_names, arc_bpm_samples = _get_lhc_arc_bpms(bpm_names, bpm_samples)
    tunez = np.mean(panda.set_index("NAME").loc[arc_bpm_names, "TUNEZ"])
    tunez_phase = _phase_mean(panda.set_index("NAME").loc[arc_bpm_names, "MUZ"])
    dpoverp_amp, pos, neg = get_dpoverp_amp(tunez, tunez_phase, bpm_samples)
    for bin in pos:
        mini, maxi = bin
        coefficients = panda[:, "AMP" + plane] * np.exp(PI2I * panda[:, "MU" + plane])
        frequencies = panda[:, "TUNE" + plane]
        new_signal = coefficients * np.exp(PI2I * np.outer(frequencies, np.arange(bpm_samples.shape[1])))
        # Remove synchro-betatron line?
        bpm_int = (bpm_samples - new_signal)[:, mini:maxi]
        # or subtract the synchrotron line? should be about the same
        bpm_int = (bpm_int.T-np.mean(bpm_int, axis=1)).T
        length = length + maxi - mini
    for bin in neg:
        mini, maxi = bin
        co = co - np.sum(bpm_samples[:, mini:maxi], axis=1)
        length = length + maxi - mini


def _get_positive_dpoverp_intervals(tunez, tunez_phase, turns):
    halfwidth = 0.1  # at 90 % of the maximum
    intervals = []
    start = (-halfwidth - tunez_phase) / tunez
    end = (halfwidth - tunez_phase) / tunez
    for i in range(int(turns * tunez) + 1):
        start_i = start + float(i) / tunez
        end_i = end + float(i) / tunez
        if start_i > 0 and end_i < turns:
            intervals.append((int(start_i), int(end_i)))
    return intervals


def _get_negative_dpoverp_intervals(tunez, tunez_phase, turns):
    return _get_positive_dpoverp_intervals(tunez, tunez_phase + 0.5, turns)


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
             for (resonance_h,
                  resonance_v,
                  resonance_l) in RESONANCE_LISTS[plane]]
    # Move to [0, 1] domain.
    freqs = [freq + 1. if freq < 0. else freq for freq in freqs]
    resonances_freqs = dict(zip(RESONANCE_LISTS[plane], freqs))
    return resonances_freqs
