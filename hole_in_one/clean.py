import sys
import time
import numpy as np

# internal options
PRINT_TIMES = False  # If True, (more) execution times will be printed
PRINT_DEBUG = False  # If True, internal debug information will be printed

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
# FFT of SVD decomposed data is much faster, but imprecise, we would need to do some super-clever integration of recomposed spectral lines
# TO DO in main: cut the turns if requested start:end
# Other method for BPM swapping and yet another one for dpoverp
#===================================================================================================


def clean(bpm_names, bpm_data, min_peak_to_peak=0.000001, max_peak_cut=10.0,
          do_resync=False):   # For one bunch, one plane, all BPMs
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
        if do_resync:
            if (bpm_names[i][-5:].upper() == "L2.B1" or
                    bpm_names[i][-5:].upper() == "R8.B2"):
                bpm_data[i, :] = np.roll(bpm_data[i, :], -1)

    maxima = get_max_bpm_readings(bpm_data)
    minima = get_min_bpm_readings(bpm_data)
    bpm_flatness = detect_flat_bpms(maxima, minima, min_peak_to_peak)
    bpm_spikes = detect_bpms_with_spikes(maxima, minima, max_peak_cut)
    exact_zeros = detect_bpms_with_exact_zeros(bpm_data)
    bad_bpms = np.zeros_like(exact_zeros)
    # Merges all bad BPMs
    bad_bpms = ((exact_zeros) | (known_bad_bpms) | (bpm_flatness) | (bpm_spikes))
    bad_bpm_indices = np.nonzero(bad_bpms)[0]
    bad_bpms_with_reasons = []
    for i in bad_bpm_indices:
        reason_for_bad_bpm = ("Found an exact zero," * int(exact_zeros[i]) +
                              "Known bad BPM," * int(known_bad_bpms[i]) +
                              ("Flat BPM, the difference between min/max is smaller than " + str(min_peak_to_peak)) * int(bpm_flatness[i]) +
                              ("Spiky BPM, found spike higher than " + str(max_peak_cut)) * int(bpm_spikes[i]))
        bad_bpms_with_reasons.append(bpm_names[i] + " " + reason_for_bad_bpm)
    good_bpms = np.logical_not(bad_bpms)
    
    if np.any(exact_zeros):
        print "Exact zeros detected. BPMs removed: " + str(np.sum(exact_zeros))
    if np.any(known_bad_bpms):
        print "Known bad BPMs removed: " + str(np.sum(known_bad_bpms))
    if np.any(bpm_flatness):
        print "Flat BPMS detected (diff min/max <= {0}. BPMs removed: {1}".format(min_peak_to_peak, np.sum(bpm_flatness))
    if np.any(bpm_spikes):
        print "Spikes > {0}mm detected. BPMs removed: {1}".format(max_peak_cut, np.sum(bpm_spikes))
    print "(Statistics for file reading) Total BPMs: {0}, Good BPMs: {1} ({2:2.2f}%), bad BPMs: {3} ({4:2.2f}%)"\
            .format(len(good_bpms), np.sum(good_bpms), 100.0 * np.mean(good_bpms), len(bad_bpm_indices), 100.0 * (1-np.mean(good_bpms)))
    print "Filtering done"
    if PRINT_TIMES:
        print ">>Time for filtering: {0}s".format(time.time() - time_start)

    if do_resync:
        return bpm_names[good_bpms], bpm_data[good_bpms,:-1], bad_bpms_with_reasons
    else:
        return bpm_names[good_bpms], bpm_data[good_bpms,:], bad_bpms_with_reasons  


def svd_clean(bpm_names, bpm_data, singular_values_amount_to_keep=12, single_svd_bpm_threshold=0.925, svd_mode='NUM'):
    time_start = time.time()
    # Parameters for matrix normalisation
    sqrt_number_of_turns = np.sqrt(bpm_data.shape[1])
    bpm_data_mean = np.mean(bpm_data)
    if svd_mode.upper()[:3]=="RAN":
        print "Using Randomized SVD"
        USV = get_singular_value_decomposition_random((bpm_data - bpm_data_mean) / sqrt_number_of_turns, singular_values_amount_to_keep)
    elif svd_mode.upper()[:3]=="SPA":
        try:
            from scipy.sparse.linalg.eigen.arpack.arpack import svds as get_singular_value_decomposition_sparse
            print "Using Sparse SVD"
            USV = get_singular_value_decomposition_sparse((bpm_data - bpm_data_mean) / sqrt_number_of_turns, k = singular_values_amount_to_keep)
        except ImportError:
            print "Failed to import SciPy, using NumPy SVD"
            USV = get_singular_value_decomposition((bpm_data - bpm_data_mean) / sqrt_number_of_turns, singular_values_amount_to_keep)
    else:
        print "Using NumPy SVD"
        USV = get_singular_value_decomposition((bpm_data - bpm_data_mean) / sqrt_number_of_turns,singular_values_amount_to_keep)
    if PRINT_TIMES:
        print ">> Time for SVD: {0}s".format(time.time() - time_start)
    
    bpm_dominance , bad_bpms_with_reasons = detect_dominant_bpms_with_reasons(bpm_names, USV[0], single_svd_bpm_threshold)
    good_bpms = np.logical_not(bpm_dominance)
    print ">> Values in GOOD BPMs:  {0}".format(np.sum(good_bpms)) 
   
    #reconstruct the SVD-cleaned data 
    A = (np.dot(USV[0][good_bpms], np.dot(np.diag(USV[1]), USV[2])) * sqrt_number_of_turns) + bpm_data_mean
    
    bpm_res = np.std(A - bpm_data[good_bpms,:], axis=1)
    print "Average BPM resolution: " + str(np.mean(bpm_res))
    if PRINT_TIMES:
        print ">> Time for svd_clean: {0}s".format(time.time() - time_start)

    return bpm_names[good_bpms], A, bpm_res, bad_bpms_with_reasons


#################################   HELPER FUNCTIONS #########################

# Detects BPMs with exact zeros due to OP workaround
def detect_bpms_with_exact_zeros(bpm_data):
    return np.logical_not(np.all(bpm_data, axis=1))


def get_max_bpm_readings(bpm_data):
    return np.max(bpm_data, axis=1)


def get_min_bpm_readings(bpm_data):
    return np.min(bpm_data, axis=1)


# Detects BPMs with the same values for all turns
def detect_flat_bpms(maxima, minima, min_peak_to_peak):
    return np.abs(maxima - minima) <= min_peak_to_peak


# Detects BPMs with spikes > max_peak_cut
def detect_bpms_with_spikes(maxima, minima, max_peak_cut):
    m3=np.empty([len(maxima),2])    
    m3[:,0]=np.abs(maxima)
    m3[:,1]=np.abs(minima)
    return np.max(m3, axis=1) > max_peak_cut


#Removes noise floor: having MxN matrix 
#and requiring K singular values we get SVD ((MxK) x diag(K) x (K,N))
def get_singular_value_decomposition(matrix,num):
    if num > min(matrix.shape):
        print "Requested more singular values than available(={0})".format(min(matrix.shape))
        return np.linalg.svd(matrix, full_matrices=False)
    print "Amount of singular values to keep: {0}".format(num)     
    U, S, V= np.linalg.svd(matrix, full_matrices=False)
    return (U[:,:num],S[:num],V[:num,:])


# In every column of U (left-singular vectors) identifies the largest value, 
# if it exceeds the limit, it finds the BPMs index (row)
def detect_dominant_bpms_with_reasons(bpm_names, U, single_svd_bpm_threshold):
    dominant = np.max(np.abs(U), axis=0) > single_svd_bpm_threshold
    if np.any(dominant):
        dominance=np.zeros(U.shape[0],dtype=bool) 
        bad_bpm_indices = set(np.argmax(np.abs(U[:,dominant]), axis=0))
        if PRINT_DEBUG:
            print "Bad BPM indices: " + str(bad_bpm_indices)
        print "Bad BPMs from SVD detected. Number of BPMs removed: {0}".format(len(bad_bpm_indices))
        bad_bpms_with_reasons = []
        for i in bad_bpm_indices:
            bad_bpms_with_reasons.append(bpm_names[i] + " Detected from SVD, single peak value is greater then " + str(single_svd_bpm_threshold))
            dominance[i]=True    
        return dominance, bad_bpms_with_reasons
    else:
        return np.zeros(U.shape[0],dtype=bool), []

#SVD of a normalized matrix using random matrix to get lower rank approximation
def get_singular_value_decomposition_random(matrix,num):
    Q = np.linalg.qr(np.dot(matrix, np.random.randn(matrix.shape[1], num+6)))[0]
    U, S, V = np.linalg.svd(np.dot(np.transpose(Q), matrix), full_matrices=False)
    return (np.dot(Q,U)[:,:num], S[:num], V[:num,:])


# For one bunch, one plane, all BPMs after previous filtering(optional) returns the decomposition to USV matrices
# and parameters needed for later recomposition, only noise cleaning is done here
def svd_for_fft(bpm_names, bpm_data, singular_values_amount_to_keep=12, single_svd_bpm_threshold=0.925):
    time_start = time.time()
    # SVD of normalised matrix
    sqrt_number_of_turns = np.sqrt(bpm_data.shape[1])
    bpm_data_mean = np.mean(bpm_data)
    USV = get_singular_value_decomposition_section((bpm_data - bpm_data_mean) / sqrt_number_of_turns, singular_values_amount_to_keep)
    if PRINT_TIMES:
        print ">> Time for SVD: {0}s".format(time.time() - time_start)
    return USV, sqrt_number_of_turns, bpm_data_mean
