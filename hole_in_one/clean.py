import logging
from datetime import datetime
import numpy as np
from scipy.fftpack import fft as scipy_fft

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
LIST_OF_WRONG_POLARITY_BPMS_BOTH_PLANES = []

########################
# Somewhere we should dump the parameters used for the analysis to a file - main, or some inputdata class?

def clean(bpm_names, bpm_data, clean_input, file_date):
    # Loops through BPM names
    known_bad_bpms = np.zeros([len(bpm_names)], dtype=bool)
    for i in range(len(bpm_names)):
        # Searches for known bad BPMs
        if bpm_names[i] in (LIST_OF_KNOWN_BAD_BPMS + clean_input.bad_bpms):
            known_bad_bpms[i] = True
        # Fixes wrong polarity
        if bpm_names[i] in (LIST_OF_WRONG_POLARITY_BPMS_BOTH_PLANES + clean_input.wrong_polarity_bpms):
            bpm_data[i, :] = -1. * bpm_data[i, :]
        # Resynchronizes BPMs between the injection point and start of the lattice
        if (not clean_input.noresync) and _resync_from_date(file_date):
            LOGGER.debug("Will resynchronize BPMs")
            if (bpm_names[i][-5:].upper() == "L2.B1" or
                    bpm_names[i][-5:].upper() == "R8.B2"):
                bpm_data[i, :] = np.roll(bpm_data[i, :], -1)
    maxima = _get_max_bpm_readings(bpm_data)
    minima = _get_min_bpm_readings(bpm_data)
    bpm_flatness = _detect_flat_bpms(maxima, minima, clean_input.peak_to_peak)
    bpm_spikes = _detect_bpms_with_spikes(maxima, minima, clean_input.max_peak)
    exact_zeros = _detect_bpms_with_exact_zeros(bpm_data)
    bad_bpms = np.zeros_like(exact_zeros)
    # Merges all bad BPMs
    bad_bpms = ((exact_zeros) | (known_bad_bpms) | (bpm_flatness) | (bpm_spikes))
    bad_bpm_indices = np.nonzero(bad_bpms)[0]
    bad_bpms_with_reasons = []
    for i in bad_bpm_indices:
        reason_for_bad_bpm = ("Found an exact zero," * int(exact_zeros[i]) +
                              "Known bad BPM," * int(known_bad_bpms[i]) +
                              ("Flat BPM, the difference between min/max is smaller than " +
                               str(clean_input.peak_to_peak)) * int(bpm_flatness[i]) +
                              ("Spiky BPM, found spike higher than " +
                               str(clean_input.max_peak)) * int(bpm_spikes[i]))
        bad_bpms_with_reasons.append(bpm_names[i] + " " + reason_for_bad_bpm)
    good_bpms = np.logical_not(bad_bpms)

    LOGGER.debug(
        "(Statistics for file reading) Total BPMs: {0}, Good BPMs: {1} ({2:2.2f}%), bad BPMs: {3} ({4:2.2f}%)"
        .format(len(good_bpms), np.sum(good_bpms), 100.0 * np.mean(good_bpms),
                len(bad_bpm_indices), 100.0 * (1 - np.mean(good_bpms))))
    LOGGER.info("Filtering done")
    
    if np.mean(good_bpms)<0.5:
        raise ValueError("Message to user: More than half of BPMs are bad. -> This could be cause a bunch not present in the machine has been selected or because a problem with the phasing of the BPMs.")

    if not clean_input.noresync:
        new_bpm_data = bpm_data[good_bpms, :-1]
    else:
        new_bpm_data = bpm_data[good_bpms, :]
    # TODO check if this is the best way to cut turns, some checks?
    if clean_input.startturn > 0 or clean_input.endturn < new_bpm_data.shape[1]:
        start=max(0,clean_input.startturn)
        end=min(clean_input.endturn, new_bpm_data.shape[1])
        new_bpm_data = new_bpm_data[:,start:end]
        
    return bpm_names[good_bpms], new_bpm_data, bad_bpms_with_reasons


def svd_decomposition(clean_input, bpm_data):
    # Parameters for matrix normalisation
    sqrt_number_of_turns = np.sqrt(bpm_data.shape[1])
    bpm_data_mean = np.mean(bpm_data)
    normalized_data = (bpm_data - bpm_data_mean) / sqrt_number_of_turns
    if clean_input is not None:
        sv_to_keep = clean_input.sing_val
    else:
        # If not clean keep all sing. values.
        sv_to_keep = bpm_data.shape[0]
    svd_functs = {
        "NUM": _get_singular_value_decomposition,
        "SPA": _get_singular_value_decomposition_sparse,
        "RAN": _get_singular_value_decomposition_random,
    }
    USV = svd_functs[clean_input.svd_mode.upper()[:3]](
        normalized_data,
        sv_to_keep
    )
    num = np.sum(USV[1] > 0.)
    USV = USV[0][:, :num], USV[1][:num], USV[2][:num, :]
    if num < USV[1].shape[0]:
        LOGGER.warn("Zero singular values detected.")
    return USV, bpm_data_mean, sqrt_number_of_turns


def svd_clean(bpm_names, bpm_data, clean_input):
    USV, bpm_data_mean, sqrt_number_of_turns = svd_decomposition(clean_input, bpm_data)

    if clean_input is not None:
        bpm_dominance, bad_bpms_with_reasons = _detect_dominant_bpms_with_reasons(
            bpm_names,
            USV[0],
            clean_input.single_svd_bpm_threshold
        )
        good_bpms = np.logical_not(bpm_dominance)
        LOGGER.debug(">> Values in GOOD BPMs:  %s", np.sum(good_bpms))

        # Reconstruct the SVD-cleaned data
        A = np.dot(USV[0][good_bpms],np.dot(np.diag(USV[1]), USV[2])) * sqrt_number_of_turns + bpm_data_mean
        bpm_res = np.std(A - bpm_data[good_bpms, :], axis=1)
        LOGGER.debug("Average BPM resolution: %s", str(np.mean(bpm_res)))
        good_bpm_names = bpm_names[good_bpms]
        good_bpm_data = A
        usv2 = (USV[2].T - np.mean(USV[2], axis=1)).T
        usv = (USV[0][good_bpms], USV[1] * sqrt_number_of_turns, usv2)
        print("np.mean(bpm_res)", np.mean(bpm_res), "np.mean(np.std(A, axis=1)", np.mean(np.std(A, axis=1)))
        if np.mean(bpm_res) > 8*np.mean(np.std(A, axis=1)): #This is factor 8 here is somewhat arbitrary.
           raise ValueError("Message to user: The data is too noisy. The most explanation is that the excitation was No excitation very low or that there ")
    else:
        usv2 = (USV[2].T - np.mean(USV[2], axis=1)).T
        usv = (USV[0], USV[1] * sqrt_number_of_turns, usv2)
        bad_bpms_with_reasons = []
        bpm_res = np.zeros_like(bpm_names)
        good_bpm_names = bpm_names
        good_bpm_data = bpm_data
    return good_bpm_names, good_bpm_data, bpm_res, bad_bpms_with_reasons, usv



# HELPER FUNCTIONS #########################

# Detects BPMs with exact zeros due to OP workaround
def _detect_bpms_with_exact_zeros(bpm_data):
    exact_zeros = np.logical_not(np.all(bpm_data, axis=1))
    if np.any(exact_zeros):
        LOGGER.debug("Exact zeros detected. BPMs removed: " +
                     str(np.sum(exact_zeros)))
    return exact_zeros


def _get_max_bpm_readings(bpm_data):
    return np.max(bpm_data, axis=1)


def _get_min_bpm_readings(bpm_data):
    return np.min(bpm_data, axis=1)


# Detects BPMs with the same values for all turns
def _detect_flat_bpms(maxima, minima, min_peak_to_peak):
    bpm_flatness = np.abs(maxima - minima) <= min_peak_to_peak
    if np.any(bpm_flatness):
        LOGGER.debug("Flat BPMS detected (diff min/max <= %s. BPMs removed: %s",
                     min_peak_to_peak,
                     np.sum(bpm_flatness))
    return bpm_flatness


# Detects BPMs with spikes > max_peak_cut
def _detect_bpms_with_spikes(maxima, minima, max_peak_cut):
    m3 = np.empty([len(maxima), 2])
    m3[:, 0] = np.abs(maxima)
    m3[:, 1] = np.abs(minima)
    bpm_spikes = np.max(m3, axis=1) > max_peak_cut
    if np.any(bpm_spikes):
        LOGGER.debug("Spikes > %s detected. BPMs removed: %s",
                     max_peak_cut,
                     np.sum(bpm_spikes))
    return bpm_spikes


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


# In every column of U (left-singular vectors) identifies the largest value,
# if it exceeds the limit, it finds the BPMs index (row)
def _detect_dominant_bpms_with_reasons(bpm_names, U, single_svd_bpm_threshold):
    dominant = np.max(np.abs(U), axis=0) > single_svd_bpm_threshold
    if np.any(dominant):
        dominance = np.zeros(U.shape[0], dtype=bool)
        bad_bpm_indices = set(np.argmax(np.abs(U[:, dominant]), axis=0))
        LOGGER.debug("Bad BPM indices: " + str(bad_bpm_indices))
        LOGGER.debug("Bad BPMs from SVD detected. Number of BPMs removed: {0}"
                     .format(len(bad_bpm_indices)))
        bad_bpms_with_reasons = []
        for i in bad_bpm_indices:
            bad_bpms_with_reasons.append(
                bpm_names[i] +
                " Detected from SVD, single peak value is greater then " +
                str(single_svd_bpm_threshold)
            )
            dominance[i] = True
        return dominance, bad_bpms_with_reasons
    else:
        return np.zeros(U.shape[0], dtype=bool), []


def _resync_from_date(acqdate):
    return acqdate > datetime(2016, 4, 1)
