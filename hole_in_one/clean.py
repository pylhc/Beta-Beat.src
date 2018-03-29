from __future__ import print_function
import logging
from datetime import datetime
import numpy as np
import pandas as pd

SPARSE_AVAILABLE = True
try:
    from scipy.sparse.linalg.eigen.arpack.arpack import svds
except ImportError:
    SPARSE_AVAILABLE = False

PI2I = 2 * np.pi * complex(0, 1)
LOGGER = logging.getLogger(__name__)

# piotr: for data 20160408 01:45+ (40cm) B1
LIST_OF_KNOWN_BAD_BPMS = ["BPM.17L8.B1", "BPM.16L8.B1", "BPM.8R8.B1", "BPM.9R8.B1", # H B1 (big ones)
                          "BPM.26L8.B1","BPM.24L8.B1", "BPM.10R6.B1","BPM.8R6.B1","BPM.32R1.B1","BPM.33R1.B1", # V B1 (big ones)
                          "BPM.12R2.B1","BPM.13R2.B1","BPM.15R6.B1","BPM.16R6.B1","BPM.19L7.B1","BPM.18L7.B1", # H B1 (small ones)
                          "BPM.21R7.B1","BPM.22R7.B1","BPM.20R8.B1","BPM.21R8.B1","BPM.19L2.B1","BPM.18L2.B1",  # H B1 (small ones)
                          "BPMR.7L5.B1","BPM.6L5.B1","BPM.8L1.B1","BPM.6L1.B1",
                          "BPM.22R3.B2", "BPM.23R3.B2", "BPMR.6L7.B2", "BPMWC.6L7.B2"] # New from ATS MD 27-07-2016
LIST_OF_WRONG_POLARITY_BPMS = []

# Noise to signal limit
NTS_LIMIT = 8.


def clean(bpm_data, clean_input, file_date):
    
    LOGGER.debug ("clean: number of BPMs in the input %s ", bpm_data.index.size)
    known_bad_bpm_names = clean_input.bad_bpms + LIST_OF_KNOWN_BAD_BPMS
    known_bad_bpms = detect_known_bad_bpms(bpm_data, known_bad_bpm_names)
    bpm_flatness = detect_flat_bpms(bpm_data, clean_input.peak_to_peak)
    bpm_spikes = detect_bpms_with_spikes(bpm_data, clean_input.max_peak)
    if  not clean_input.no_exact_zeros :
      exact_zeros = detect_bpms_with_exact_zeros(bpm_data)
    else :
      exact_zeros = pd.DataFrame()
      LOGGER.debug ("clean: Skipped exact zero check %d",exact_zeros.size)
      

    original_bpms = bpm_data.index
    
    all_bad_bpms = _index_union(known_bad_bpms, bpm_flatness,
                                bpm_spikes, exact_zeros)
    bpm_data = bpm_data.loc[bpm_data.index.difference(all_bad_bpms)]

    bad_bpms_with_reasons = _get_bad_bpms_summary(
        clean_input, known_bad_bpms, bpm_flatness, bpm_spikes, exact_zeros
    )
    _report_clean_stats(original_bpms.size, bpm_data.index.size)

    wrong_plrty_names = (LIST_OF_WRONG_POLARITY_BPMS +
                         clean_input.wrong_polarity_bpms)
    bpm_data = fix_polarity(wrong_plrty_names, bpm_data)

    if not clean_input.noresync:
        bpm_data = resync_bpms(bpm_data, file_date)

    return bpm_data, bad_bpms_with_reasons


def fix_polarity(wrong_polarity_names, bpm_data):
    """
    Fixes wrong polarity
    """
    bpm_data.loc[wrong_polarity_names] *= -1
    return bpm_data


def resync_bpms(bpm_data, file_date):
    """
    Resynchronizes BPMs between the injection point and start of the lattice if
    the acquisition date is after 2016-04-01.
    """
    if file_date > datetime(2016, 4, 1):
        LOGGER.debug("Will resynchronize BPMs")
        bpm_data.loc[bpm_data.index.str.endswith("L2.B1")] =\
            np.roll(bpm_data.loc[bpm_data.index.str.endswith("L2.B1")], -1, axis=1)
        bpm_data.loc[bpm_data.index.str.endswith("R8.B2")] =\
            np.roll(bpm_data.loc[bpm_data.index.str.endswith("R8.B2")], -1, axis=1)
        bpm_data = bpm_data.iloc[:, :-1]
    return bpm_data


def detect_known_bad_bpms(bpm_data, list_of_bad_bpms):
    """
    Searches for known bad BPMs
    """
    known_bad_bpms = bpm_data.index.intersection(list_of_bad_bpms)
    return known_bad_bpms


def detect_flat_bpms(bpm_data, min_peak_to_peak):
    """
    Detects BPMs with the same values for all turns
    """
    cond = ((bpm_data.max(axis=1) - bpm_data.min(axis=1)).abs() <
            min_peak_to_peak)
    bpm_flatness = bpm_data[cond].index
    if bpm_flatness.size:
        LOGGER.debug(
            "Flat BPMS detected (diff min/max <= %s. BPMs removed: %s",
            min_peak_to_peak,
            bpm_flatness.size,
        )
    return bpm_flatness


def detect_bpms_with_spikes(bpm_data, max_peak_cut):
    """
    Detects BPMs with spikes > max_peak_cut
    """
    too_high = bpm_data[bpm_data.max(axis=1) > max_peak_cut].index
    too_low = bpm_data[bpm_data.min(axis=1) < -max_peak_cut].index
    bpm_spikes = too_high.union(too_low)
    if bpm_spikes.size:
        LOGGER.debug("Spikes > %s detected. BPMs removed: %s",
                     max_peak_cut,
                     bpm_spikes.size)
    return bpm_spikes


def detect_bpms_with_exact_zeros(bpm_data):
    """
    Detects BPMs with exact zeros due to OP workaround
    """
    exact_zeros = bpm_data[~np.all(bpm_data, axis=1)].index
    if exact_zeros.size:
        LOGGER.debug("Exact zeros detected. BPMs removed: %s",
                     exact_zeros.size)
    return exact_zeros


def svd_decomposition(clean_input, bpm_data):
    # Parameters for matrix normalisation
    sqrt_number_of_turns = np.sqrt(bpm_data.shape[1])
    bpm_data_mean = bpm_data.values.mean()
    normalized_data = (bpm_data - bpm_data_mean) / sqrt_number_of_turns
    svd_functs = {
        "NUM": _get_singular_value_decomposition,
        "SPA": _get_singular_value_decomposition_sparse,
        "RAN": _get_singular_value_decomposition_random,
    }
    U, S, V = svd_functs[clean_input.svd_mode.upper()[:3]](
        normalized_data,
        clean_input.sing_val
    )
    num = np.sum(S > 0.)
    U = pd.DataFrame(index=bpm_data.index, data=U)
    USV = U.loc[:, :num], S[:num], V[:num, :]
    if num < USV[1].shape[0]:
        LOGGER.warn("Zero singular values detected.")
    return USV, bpm_data_mean, sqrt_number_of_turns


def svd_clean(bpm_data, clean_input):
    USV, bpm_data_mean, sqrt_number_of_turns = svd_decomposition(clean_input, bpm_data)
    clean_U, dominance_summary = _clean_dominant_bpms(
        USV[0],
        clean_input.single_svd_bpm_threshold
    )

    # Reconstruct the SVD-cleaned data
    USV = clean_U, USV[1], USV[2]
    A = USV[0].dot(np.dot(np.diag(USV[1]), USV[2]))
    A = A * sqrt_number_of_turns + bpm_data_mean

    # BPM resolution computation
    bpm_res = (A - bpm_data.loc[clean_U.index]).std(axis=1)
    LOGGER.debug("Average BPM resolution: %s", str(np.mean(bpm_res)))
    LOGGER.debug("np.mean(np.std(A, axis=1): %s", np.mean(np.std(A, axis=1)))
    if np.mean(bpm_res) > NTS_LIMIT * np.mean(np.std(A, axis=1)):
        raise ValueError(
            "The data is too noisy. The most probable explanation"
            " is that there was no excitation or it was very low.")

    # Denormalize the data
    V = (USV[2].T - np.mean(USV[2], axis=1)).T
    USV = (USV[0], USV[1] * sqrt_number_of_turns, V)

    good_bpm_data = A
    return good_bpm_data, bpm_res, dominance_summary, USV


# HELPER FUNCTIONS #########################


def _get_bad_bpms_summary(clean_input, known_bad_bpms,
                          bpm_flatness, bpm_spikes, exact_zeros):
    bad_bpms_summary = []
    for bpm_name in known_bad_bpms:
        bad_bpms_summary.append("{} Known bad BPM".format(bpm_name))
    for bpm_name in bpm_flatness:
        bad_bpms_summary.append(
            "{} Flat BPM, the difference "
            "between min/max is smaller than {}"
            .format(bpm_name, clean_input.peak_to_peak)
        )
    for bpm_name in bpm_spikes:
        bad_bpms_summary.append(
            "{} Spiky BPM, found spike higher than {}"
            .format(bpm_name, clean_input.max_peak)
        )
    for bpm_name in exact_zeros:
        bad_bpms_summary.append(
            "{} Found an exact zero"
            .format(bpm_name)
        )
    return bad_bpms_summary


def _report_clean_stats(n_total_bpms, n_good_bpms):
     
    LOGGER.debug("Filtering done:")

    if n_total_bpms == 0:
        raise ValueError("Total Number of BPMs after filtering is zero " )
    
    n_bad_bpms = n_total_bpms - n_good_bpms
    LOGGER.debug(
        "(Statistics for file reading) Total BPMs: {0}, "
        "Good BPMs: {1} ({2:2.2f}%), bad BPMs: {3} ({4:2.2f}%)"
        .format(n_total_bpms,
                n_good_bpms, 100.0 * (n_good_bpms / float(n_total_bpms)),
                n_bad_bpms, 100.0 * (n_bad_bpms / float(n_total_bpms))))

    if (n_good_bpms / float(n_total_bpms)) < 0.5:
        raise ValueError(
            "More than half of BPMs are bad. "
            "This could be cause a bunch not present in the machine has been "
            "selected or because a problem with the phasing of the BPMs."
        )


def _index_union(*indices):
    new_index = pd.Index([])
    for index in indices:
        new_index = new_index.union(index)
    return new_index


# Removes noise floor: having MxN matrix
# and requiring K singular values we get SVD ((MxK) x diag(K) x (K,N))
def _get_singular_value_decomposition(matrix, num):
    LOGGER.debug("Using NumPy SVD")
    if num > min(matrix.shape):
        LOGGER.warn("Requested more singular values than available(={0})"
                    .format(min(matrix.shape)))
        return np.linalg.svd(matrix, full_matrices=False)
    LOGGER.debug("Amount of singular values to keep: {0}".format(num))
    U, S, V = np.linalg.svd(matrix, full_matrices=False)
    return (U[:, :num], S[:num], V[:num, :])


def _get_singular_value_decomposition_sparse(matrix, num):
    if not SPARSE_AVAILABLE:
        LOGGER.warn("Cannot import scipy sparse SVD, falling back to Numpy.")
        return _get_singular_value_decomposition(matrix, num)
    LOGGER.debug("Using Sparse SVD")
    return svds(matrix, k=num)


# SVD of a normalized matrix using random matrix to
# get lower rank approximation
def _get_singular_value_decomposition_random(matrix, num):
    LOGGER.debug("Using Randomized SVD")
    Q = np.linalg.qr(np.dot(matrix,
                            np.random.randn(matrix.shape[1], num + 6)))[0]
    U, S, V = np.linalg.svd(np.dot(np.transpose(Q), matrix),
                            full_matrices=False)
    return (np.dot(Q, U)[:, :num], S[:num], V[:num, :])


def _clean_dominant_bpms(U, single_svd_bpm_threshold):
    if single_svd_bpm_threshold < 1 / np.sqrt(2):
        LOGGER.warn("Careful, the single_svd_bpm_threshold looks too low %s.",
                    single_svd_bpm_threshold)
    dominant_bpms = U[np.max(U.abs(), axis=1) > single_svd_bpm_threshold].index
    dominance_summary = []
    if dominant_bpms.size > 0:
        for bpm_name in dominant_bpms:
            dominance_summary.append(
                "{} Detected from SVD, single peak value is greater then {}"
                .format(bpm_name, str(single_svd_bpm_threshold))
            )
        LOGGER.debug("Bad BPMs from SVD detected. Number of BPMs removed: %s",
                     dominant_bpms.size)
    clean_U = U.loc[U.index.difference(dominant_bpms)]
    return clean_U, dominance_summary
