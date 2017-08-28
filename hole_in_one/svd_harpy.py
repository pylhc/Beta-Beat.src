from __future__ import print_function
import os
import multiprocessing
import logging
import numpy as np
import pandas as pd
from harpy.harmonic_analysis import HarmonicAnalysis

from harpy import _python_path_manager
_python_path_manager.append_betabeat()
from Utilities import outliers  # noqa
LOGGER = logging.getLogger(__name__)

PI2I = 2 * np.pi * complex(0, 1)

SPECTR_COLUMN_NAMES = ["FREQ", "AMP"]

RESONANCE_LISTS = {"X": ((0, 1, 0), (-2, 0, 0), (0, 2, 0), (-3, 0, 0), (-1, -1, 0),
                         (2, -2, 0), (0, -2, 0), (1, -2, 0), (-1, 3, 0), (1, 2, 0), 
                         (-2, 1, 0), (1, 1, 0), (2, 0, 0), (-1, -2, 0), (3, 0, 0)),
                   "Y": ((1, 0, 0), (-1, 1, 0), (-2, 0, 0), (1, -1, 0), (0, -2, 0),
                         (0, -3, 0), (2, 1, 0), (-1, 3, 0), (1, 1, 0), (-1, 2, 0)),
                   "Z": ((1, 0, 1), (0, 1 ,1))
                    }

MAIN_LINES = {"X": (1, 0, 0),
              "Y": (0, 1, 0),
              "Z": (0, 0, 1)}

DEFAULT_DFT_METHOD = "laskar"

NUM_HARMS = 300

NUM_HARMS_SVD = 100

PROCESSES = multiprocessing.cpu_count()

#TODO Mising nattunez - not needed yet
#TODO Improve search of higher order resonances smaller tolerance
# TODO if non-zero tunez try chroma_calculation 
# TODO the amplitudes are not rescaled by ampx/y
# TODO add header info, tunes, nattunes, dpp
#TODO separate resonance for non-zero tunez computed if tunez >0
#other option of noise estimate is to take average of to take average amplitude of fft of few last taken or first not taken modes  
def svd_harpy(bpm_names, bpm_matrix, usv, plane, harpy_input, panda):
    SV = np.dot(np.diag(usv[1]), usv[2])
    if harpy_input.harpy_mode == 'fast':
        freqs = laskar_per_mode(np.mean(SV,axis=0), NUM_HARMS)
    else:
        pool = multiprocessing.Pool(np.min([PROCESSES, SV.shape[0]]))
        freqs = []
        for i in range(SV.shape[0]):
            args = (SV[i, :], NUM_HARMS_SVD)
            if harpy_input.sequential:
                freqs.extend(laskar_per_mode(*args))
            else:
                pool.apply_async(laskar_per_mode, args, callback=freqs.extend)
        pool.close()
        pool.join()
    frequencies = np.array(freqs)
    svd_coefficients = compute_coefs_for_freqs(SV, frequencies)
    bpm_coefficients = np.dot(usv[0], svd_coefficients)
    
    spectrum_phases=pd.DataFrame(data=np.angle(bpm_coefficients) / (2 * np.pi), columns=frequencies, index=bpm_names, dtype=float) 
    spectrum_amplitudes=pd.DataFrame(data=np.abs(bpm_coefficients), columns=frequencies, index=bpm_names, dtype=float)
    panda, mask, bad_bpms = get_main_resonances(frequencies, bpm_coefficients, harpy_input, plane, panda, bpm_matrix)
    panda = get_noise(bpm_matrix, panda)
    panda = get_natural_tunes(frequencies, bpm_coefficients, harpy_input, plane, panda)
    resonances_freqs = compute_resonances_freqs(plane, harpy_input)
    if harpy_input.tunez > 0.0:
        resonances_freqs.update(compute_resonances_freqs("Z", harpy_input))
    panda = resonance_search(frequencies, bpm_coefficients, harpy_input.tolerance, resonances_freqs, panda)
    panda['BINDEX']=0
    panda['SLABEL']=0
    return panda.loc[mask,:], {'PHASE': spectrum_phases, 'AMP': spectrum_amplitudes}, bad_bpms


def laskar_per_mode(sv, number_of_harmonics):
    har_analysis = HarmonicAnalysis(sv)
    freqs, _ = har_analysis.laskar_method(number_of_harmonics)
    return freqs


def compute_coefs_for_freqs(samples, freqs): # samples is a matrix
    n = samples.shape[1]
    coefficients = np.dot(samples, np.exp(-PI2I * np.outer(np.arange(n), freqs))) / n
    return coefficients


def get_main_resonances(frequencies, coefficients, harpy_input, plane, panda, bpm_matrix):
    results = panda
    h, v, l = MAIN_LINES[plane]
    freq = (h * harpy_input.tunex) + (v * harpy_input.tuney) + (l * harpy_input.tunez)
    min = freq - harpy_input.tolerance
    max = freq + harpy_input.tolerance
    indices = np.where((frequencies >= min) & (frequencies <= max))[0]
    bad_bpms = []
    if len(indices) == 0:
        raise ValueError("No main resonances found, try increase the tolerance or adjust the tunes")
    elif len(indices) == 1:
        LOGGER.info("Only a single frequency found within tolerance for main resonance")
        results['TUNE'+plane] = frequencies[indices[0]]
        results['AMP'+plane] = np.abs(coefficients[:,indices[0]])
        results['MU'+plane] = np.angle(coefficients[:,indices[0]]) / (2 * np.pi)
        results['AVG_AMP'+plane] = results.loc[:,'AMP'+plane]
        results['AVG_MU'+plane] = results.loc[:,'MU'+plane]
        mask = np.ones(bpm_matrix.shape[0], dtype=bool)
    else:
        max_indices = indices[np.argmax(np.abs(coefficients[:,indices]),axis=1)]
        results['TUNE'+plane] = frequencies[max_indices]
        results['AMP'+plane] = np.abs(coefficients[np.arange(coefficients.shape[0]), max_indices])
        results['MU'+plane] = np.angle(coefficients[np.arange(coefficients.shape[0]), max_indices])/ (2 * np.pi)
        mask = outliers.get_filter_mask(results.loc[:,'TUNE'+plane], limit=1e-5)
        for i in np.arange(len(mask)):
            if not mask[i]:
                bad_bpms.append((results.at[i,'NAME'] + " Tune is too far from average"))
        avg_coefs = compute_coefs_for_freqs(bpm_matrix, np.mean(results.loc[mask,'TUNE'+plane]))
        results['AVG_AMP'+plane] = np.abs(avg_coefs)
        results['AVG_MU'+plane] = np.angle(avg_coefs) / (2 * np.pi)
    if harpy_input.tunez > 0.0:
        z_tolerance = 0.0001
        h, v, l = MAIN_LINES["Z"]
        freq = (h * harpy_input.tunex) + (v * harpy_input.tuney) + (l * harpy_input.tunez)
        min = freq - z_tolerance
        max = freq + z_tolerance
        indices = np.where((frequencies >= min) & (frequencies <= max))[0]
        if len(indices) == 0:
            raise ValueError("No longitudinal resonances found, try adjust the tune")
        elif len(indices) == 1:
            LOGGER.info("Only a single frequency found within tolerance for longitudinal resonance")
            results['TUNEZ'] = frequencies[indices[0]]
            results['AMPZ'] = np.abs(coefficients[:,indices[0]])
            results['MUZ'] = np.angle(coefficients[:,indices[0]]) / (2 * np.pi)
        else:
            max_indices = indices[np.argmax(np.abs(coefficients[:,indices]),axis=1)]
            results['TUNEZ'] = frequencies[max_indices]
            results['AMPZ'] = np.abs(coefficients[np.arange(coefficients.shape[0]), max_indices])
            results['MUZ'] = np.angle(coefficients[np.arange(coefficients.shape[0]), max_indices])/ (2 * np.pi)

    return results, mask, bad_bpms


def get_natural_tunes(frequencies, coefficients, harpy_input, plane, panda):
    results = panda
    if harpy_input.nattunex is not None and harpy_input.nattuney is not None:
        if plane == "X":
            freq = harpy_input.nattunex
        if plane == "Y":
            freq = harpy_input.nattunex
        min = freq - harpy_input.tolerance
        max = freq + harpy_input.tolerance
        indices = np.where((frequencies >= min) & (frequencies <= max))[0]
        if len(indices) == 0:
            results['NATTUNE'+plane]=0.0
            results['NATAMP'+plane]=0.0
        elif len(indices) == 1:
            results['NATTUNE'+plane]=frequencies[indices[0]]
            results['NATAMP'+plane]=np.abs(coefficients[:,indices[0]])
        else:
            max_indices = indices[np.argmax(np.abs(coefficients[:,indices]),axis=1)]
            results['NATTUNE'+plane] = frequencies[max_indices]
            results['NATAMP'+plane] = np.abs(coefficients[np.arange(coefficients.shape[0]), max_indices])
    else:
        results['NATTUNE'+plane]=0.0
        results['NATAMP'+plane]=0.0
    return results


def get_noise(bpm_matrix, panda):
    result= panda
    leng=512
    pos = 8
    ns=np.empty([bpm_matrix.shape[0],int(np.floor(bpm_matrix.shape[1]/leng))])
    for i in range(bpm_matrix.shape[1]/leng):
        ns[:,i] = np.partition(np.abs(np.fft.fft(bpm_matrix[:,:leng])), pos, axis=1)[:,pos]/leng
    result['NOISE'] = np.mean(ns, axis=1)
    return result


    # the amplitudes are not rescaled by ampx/y
def resonance_search(frequencies, coefficients, tolerance, resonances_freqs, panda): # vectorized - coefficiens are matrix
    results = panda
    bins = [((resonance_freq - tolerance, resonance_freq + tolerance), resonance)
            for resonance, resonance_freq in resonances_freqs.iteritems()]
    for bin, resonance in bins:
        min, max = bin
        resstr = get_resonance_suffix(resonance)
        indices = np.where((frequencies >= min) & (frequencies <= max))[0]
        if len(indices) == 0:
            results['AMP'+resstr] = 0.0
            results['PHASE'+resstr] = 0.0
            results['FREQ'+resstr] = 0.0
            continue
        elif len(indices) == 1:
            results['AMP'+resstr] = np.abs(coefficients[:,indices[0]])
            results['PHASE'+resstr] = np.angle(coefficients[:,indices[0]]) / (2 * np.pi)
            results['FREQ'+resstr] = frequencies[indices[0]]
        else:
            max_indices = indices[np.argmax(np.abs(coefficients[:,indices]),axis=1)]
            results['AMP'+resstr] = np.abs(coefficients[np.arange(coefficients.shape[0]), max_indices])
            results['PHASE'+resstr] = np.angle(coefficients[np.arange(coefficients.shape[0]), max_indices])/ (2 * np.pi)
            results['FREQ'+resstr] = frequencies[max_indices]
            # TODO: Is it right to remove already used lines? I dont think so...:
    return results

def get_dpoverp_amp(tunez, tunez_phase, bpm_samples):# full period is one, phase is of the arc bpms in horizontal plane
    # assumes non-zero tunez
    # interval is around integer phases for positive dpoverp
    # and around halfes for the negative dpoverp
    turns=bpm_samples.shape[1]
    pos=get_positive_dpoverp_intervals(tunez, tunez_phase, turns)
    neg=get_negative_dpoverp_intervals(tunez, tunez_phase, turns)
    co=np.zeros_like(arc_bpm_names)
    length=0
    for bin in pos:
        mini, maxi = bin
        co = co + np.sum(bpm_samples[:,mini:maxi], axis=1)
        length = length + maxi - mini
    for bin in neg:
        mini, maxi = bin
        co = co - np.sum(bpm_samples[:,mini:maxi], axis=1)
        length = length + maxi - mini
    codx = co * dx / length #dispersion at the same BPMs TODO check units
    dx2 = dx * dx
    return np.sum(codx) / np.sum(dx2), pos, neg


def get_chroma(tunez, tunez_phase, harpy.input, plane, bpm_samlpes, panda):
    dpoverp_amp, pos, neg = get_dpoverp_amp(tunez, tunez_phase, bpm_samples)
    for bin in pos:
        mini, maxi = bin
        bpm_int=bpm_samples[:,mini:maxi]
        bpm_int = (bpm_int.T-np.mean(bpm_int,axis=1)).T # or subtract the synchrotron line? should be about the same
        main_coefs = compute_coefs_for_freqs(bpm_int, panda.loc[:,'TUNE'+plane])
        bpm_int = bpm_int - (main_coefs * np.exp(PI2I * np.outer(np.arange(n), panda.loc[:,'TUNE'+plane])))
        #we subtracted the main line - now we have signals to look for natural tunes
        length = length + maxi - mini
    for bin in neg:
        mini, maxi = bin
        co = co - np.sum(bpm_samples[:,mini:maxi], axis=1)
        length = length + maxi - mini    

def get_positive_dpoverp_intervals(tunez, tunez_phase, turns):
    halfwidth=0.1 # at 90 % of the maximum
    intervals=[]
    start = (-halfwidth - tunez_phase) / tunez
    end = (halfwidth - tunez_phase) / tunez
    for i in range(int(turns*tunez)+1):
        start_i = start + float(i) / tunez   
        end_i = end + float(i) / tunez     
        if start_i > 0 and end_i < turns:
            intervals.append((int(start_i),int(end_i)))
    return intervals


def get_negative_dpoverp_intervals(tunez, tunez_phase, turns):
    return get_positive_dpoverp_intervals(tunez, tunez_phase + 0.5, turns)


def get_resonance_suffix(resonance):           
    x, y, z = resonance
    if z == 0:
        resstr = (str(x) + str(y)).replace("-", "_")
    else:
        resstr = (str(x) + str(y) + str(z)).replace("-", "_")
    return resstr


def compute_resonances_freqs(plane, harpy_input):
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


#def create_lin_header(plane, tune, rms_tune):
#    plane_number = "1" if plane == "X" else "2"
#    lin_header={}
#    lin_header["DPP"] = 0.0
#    lin_header["Q" + plane_number] = tune
#    lin_header["Q" + plane_number + "RMS"] = rms_tune
#    return lin_header
