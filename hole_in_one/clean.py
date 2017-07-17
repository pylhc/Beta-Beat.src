import time
import numpy as np
import logging
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
LIST_OF_OUT_OF_SYNC_BPMS_B1 = ["BPM.33L2.B1", "BPM.32L2.B1", "BPM.31L2.B1", "BPM.30L2.B1",
                               "BPM.29L2.B1", "BPM.28L2.B1", "BPM.27L2.B1", "BPM.26L2.B1", "BPM.25L2.B1", "BPM.24L2.B1", "BPM.23L2.B1", "BPM.22L2.B1", "BPM.21L2.B1", "BPM.20L2.B1",
                               "BPM.19L2.B1", "BPM.18L2.B1", "BPM.17L2.B1", "BPM.16L2.B1", "BPM.15L2.B1", "BPM.14L2.B1", "BPM.13L2.B1", "BPM.12L2.B1", "BPM.11L2.B1", "BPM.10L2.B1",
                               "BPM.9L2.B1", "BPM.8L2.B1", "BPM.7L2.B1", "BPMR.6L2.B1", "BPMYB.5L2.B1", "BPMYB.4L2.B1", "BPMWI.4L2.B1", "BPMS.2L2.B1", "BPMSW.1L2.B1"]
LIST_OF_OUT_OF_SYNC_BPMS_B2 = ["BPM.34R8.B2", "BPM.33R8.B2", "BPM.32R8.B2", "BPM.31R8.B2", "BPM.30R8.B2",
                               "BPM.29R8.B2", "BPM.28R8.B2", "BPM.27R8.B2", "BPM.26R8.B2", "BPM.25R8.B2", "BPM.24R8.B2", "BPM.23R8.B2", "BPM.22R8.B2", "BPM.21R8.B2", "BPM.20R8.B2",
                               "BPM.19R8.B2", "BPM.18R8.B2", "BPM.17R8.B2", "BPM.16R8.B2", "BPM.15R8.B2", "BPM.14R8.B2", "BPM.13R8.B2", "BPM.12R8.B2", "BPM.11R8.B2", "BPM.10R8.B2",
                               "BPM.9R8.B2", "BPM.8R8.B2", "BPM_A.7R8.B2", "BPMR.6R8.B2", "BPMYB.5R8.B2", "BPMYB.4R8.B2", "BPMWI.4R8.B2", "BPMS.2R8.B2", "BPMSW.1R8.B2"]
LIST_OF_WRONG_POLARITY_BPMS_BOTH_PLANES = []

########################
# Somewhere we should dump the parameters used for the analysis to a file - main, or some inputdata class?
# FFT of SVD decomposed data is much faster, but imprecise, we would need to do some super-clever integration of recomposed spectral lines
# TO DO in main: cut the turns if requested start:end
# Other method for BPM swapping and yet another one for dpoverp
#===================================================================================================


def clean(bpm_names, bpm_data, clean_input):
    time_start = time.time()
    # Loops through BPM names
    known_bad_bpms = np.zeros([len(bpm_names)], dtype=bool)
    for i in range(len(bpm_names)):
        # Searches for known bad BPMs
        if bpm_names[i] in LIST_OF_KNOWN_BAD_BPMS:
            known_bad_bpms[i] = True
        # Fixes wrong polarity
        if bpm_names[i] in LIST_OF_WRONG_POLARITY_BPMS_BOTH_PLANES:
            bpm_data[i, :] = -1. * bpm_data[i, :]
        # Resynchronizes BPMs between the injection point and start of the lattice
        if not clean_input.noresync:
            if (bpm_names[i] in LIST_OF_OUT_OF_SYNC_BPMS_B1 or
                    bpm_names[i] in LIST_OF_OUT_OF_SYNC_BPMS_B2):
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

    if np.any(exact_zeros):
        LOGGER.debug("Exact zeros detected. BPMs removed: " +
                     str(np.sum(exact_zeros)))
    if np.any(known_bad_bpms):
        LOGGER.debug("Known bad BPMs removed: " +
                     str(np.sum(known_bad_bpms)))
    if np.any(bpm_flatness):
        LOGGER.debug("Flat BPMS detected (diff min/max <= {0}. BPMs removed: {1}"
                     .format(clean_input.peak_to_peak, np.sum(bpm_flatness)))
    if np.any(bpm_spikes):
        LOGGER.debug("Spikes > {0}mm detected. BPMs removed: {1}"
                     .format(clean_input.max_peak, np.sum(bpm_spikes)))
    LOGGER.debug(
        "(Statistics for file reading) Total BPMs: {0}, Good BPMs: {1} ({2:2.2f}%), bad BPMs: {3} ({4:2.2f}%)"
        .format(len(good_bpms), np.sum(good_bpms), 100.0 * np.mean(good_bpms),
                len(bad_bpm_indices), 100.0 * (1 - np.mean(good_bpms))))
    LOGGER.info("Filtering done")
    LOGGER.debug(">>Time for filtering: {0}s".format(time.time() - time_start))

    if not clean_input.noresync:
        new_bpm_data = bpm_data[good_bpms, :-1]
    else:
        new_bpm_data = bpm_data[good_bpms, :]
    return bpm_names[good_bpms], new_bpm_data, bad_bpms_with_reasons


def svd_clean(bpm_names, bpm_data, clean_input):
    # Parameters for matrix normalisation
    sqrt_number_of_turns = np.sqrt(bpm_data.shape[1])
    bpm_data_mean = np.mean(bpm_data)
    normalized_data = (bpm_data - bpm_data_mean) / sqrt_number_of_turns
    sv_to_keep = clean_input.sing_val
    svd_functs = {
        "NUM": _get_singular_value_decomposition,
        "SPA": _get_singular_value_decomposition_sparse,
        "RAN": _get_singular_value_decomposition_random,
    }

    time_start = time.time()
    USV = svd_functs[clean_input.svd_mode.upper()[:3]](
        normalized_data,
        sv_to_keep
    )
    LOGGER.debug(">> Time for SVD: {0}s".format(time.time() - time_start))

    bpm_dominance, bad_bpms_with_reasons = _detect_dominant_bpms_with_reasons(
        bpm_names,
        USV[0],
        clean_input.single_svd_bpm_threshold
    )
    good_bpms = np.logical_not(bpm_dominance)
    LOGGER.debug(">> Values in GOOD BPMs:  {0}".format(np.sum(good_bpms)))

    # Reconstruct the SVD-cleaned data
    A = (np.dot(USV[0][good_bpms],
         np.dot(np.diag(USV[1]), USV[2])) * sqrt_number_of_turns) + bpm_data_mean
       
    bpm_res = np.std(A - bpm_data[good_bpms, :], axis=1)
    LOGGER.debug("Average BPM resolution: " + str(np.mean(bpm_res)))
    LOGGER.debug(">> Time for svd_clean: {0}s"
                 .format(time.time() - time_start))
    return bpm_names[good_bpms], A, bpm_res, bad_bpms_with_reasons, (USV[0][good_bpms], USV[1], USV[2])


# HELPER FUNCTIONS #########################

# Detects BPMs with exact zeros due to OP workaround
def _detect_bpms_with_exact_zeros(bpm_data):
    return np.logical_not(np.all(bpm_data, axis=1))


def _get_max_bpm_readings(bpm_data):
    return np.max(bpm_data, axis=1)


def _get_min_bpm_readings(bpm_data):
    return np.min(bpm_data, axis=1)


# Detects BPMs with the same values for all turns
def _detect_flat_bpms(maxima, minima, min_peak_to_peak):
    return np.abs(maxima - minima) <= min_peak_to_peak


# Detects BPMs with spikes > max_peak_cut
def _detect_bpms_with_spikes(maxima, minima, max_peak_cut):
    m3 = np.empty([len(maxima), 2])
    m3[:, 0] = np.abs(maxima)
    m3[:, 1] = np.abs(minima)
    return np.max(m3, axis=1) > max_peak_cut


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


# For one bunch, one plane, all BPMs after previous filtering(optional)
# returns the decomposition to USV matrices and parameters needed for later
# recomposition, only noise cleaning is done here
def svd_for_fft(bpm_names, bpm_data, singular_values_amount_to_keep=12,
                single_svd_bpm_threshold=0.925):
    time_start = time.time()
    # SVD of normalised matrix
    sqrt_number_of_turns = np.sqrt(bpm_data.shape[1])
    bpm_data_mean = np.mean(bpm_data)
    USV = get_singular_value_decomposition_section(
        (bpm_data - bpm_data_mean) /
        sqrt_number_of_turns, singular_values_amount_to_keep
    )
    LOGGER.debug(">> Time for SVD: {0}s".format(time.time() - time_start))
    return USV, sqrt_number_of_turns, bpm_data_mean


###############################################################

def laskar_method(num_harmonics, sample):
    samples = sample  # Copy the samples array.
    n = len(samples)
    ints = np.arange(n)
    coefficients = []
    frequencies = []
    for _ in range(num_harmonics):
        # Compute this harmonic frequency and coefficient.
        dft_data = scipy_fft(samples)
        frequency = jacobsen(dft_data)
        coefficient = compute_coef_simple(samples, frequency * n) / n

        # Store frequency and amplitude
        coefficients.append(coefficient)
        frequencies.append(frequency)

        # Subtract the found pure tune from the signal
        new_signal = coefficient * np.exp(PI2I * frequency * ints)
        samples = samples - new_signal

    coefficients, frequencies = zip(*sorted(zip(coefficients, frequencies),
                                                key=lambda tuple: np.abs(tuple[0]),
                                                reverse=True))
    return frequencies, coefficients


def laskar_method_modes_freqs(samples, freqs):
    n = samples.shape[1]
    ints = np.arange(n)
    frequencies = np.empty([len(freqs)])
    coefficients = np.empty([samples.shape[0],len(freqs)], dtype=complex)
    i = 0
    for i, frequency in enumerate(freqs):
        frequencies[i] = frequency
        coefficients[:,i] = np.sum(np.exp(-PI2I * frequency * ints) * samples, axis=1)/ n
        
    #coefficients, frequencies = zip(*sorted(zip(coefficients, frequencies), key=lambda tuple: np.abs(tuple[0]), reverse=True))
    return frequencies, coefficients



def jacobsen(dft_values):
    """
    This method interpolates the real frequency of the
    signal using the three highest peaks in the FFT.
    """
    k = np.argmax(np.abs(dft_values))
    n = len(dft_values)
    r = dft_values
    delta = np.tan(np.pi / n) / (np.pi / n)
    kp = (k + 1) % n
    km = (k - 1) % n
    delta = delta * np.real((r[km] - r[kp]) / (2 * r[k] - r[km] - r[kp]))
    return (k + delta) / n


def compute_coef_simple(samples, kprime):
    """
    Computes the coefficient of the Discrete Time Fourier
    Transform corresponding to the given frequency (kprime).
    """
    n = len(samples)
    freq = kprime / n
    exponents = np.exp(-PI2I * freq * np.arange(n))
    coef = np.sum(exponents * samples)
    return coef

def _get_allowed_length(rang=[300, 10000], p2max=14, p3max=9, p5max=6):
    ind = np.indices((p2max, p3max, p5max))
    nums = (np.power(2, ind[0]) *
            np.power(3, ind[1]) *
            np.power(5, ind[2])).reshape(p2max * p3max * p5max)
    nums = nums[(nums > rang[0]) & (nums <= rang[1])]
    return np.sort(nums)
 

###########################################################
