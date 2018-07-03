import os
from utils import outliers
import numpy as np

PI2 = 2 * np.pi
PI2I = PI2 * 1j


def get_start_and_end_turns_kicker(matrix, num):
    matrix = matrix.sub(np.mean(matrix, axis=1), axis=1)
    u, s, v = np.linalg.svd(matrix.sub(np.mean(matrix, axis=1), axis=1), full_matrices=False)
    main_mode = v[0, :]
    start = np.argmax(np.abs(main_mode))
    envelope = np.maximum.accumulate(np.log10(np.abs(main_mode[::-1])))[::-1]
    end = len(main_mode) - 1
    indx = np.argmax(np.linspace(envelope[start], envelope[-1], num=(len(main_mode)-start)) - envelope[start:])
    if envelope[indx] > envelope[start]-1:
        end = start + indx
    return start, end


def get_phase_advances_and_unscaled_betas(u, s, v, tune, tol): # u is indexed by BPM names, right?
    n = v.shape[1]
    coefs_in_tol = np.dot(u, np.fft.fft2(np.dot(s, v), axis=1)[:, (tune - tol) * n:(tune + tol) * n])
    tune_ind = np.argmax(np.abs(coefs_in_tol), axis = 1)
    coefs = coefs_in_tol[np.arange(coefs_in_tol.shape[0]), tune_ind] / n
    tune_mask = outliers.get_filter_mask(tune_ind, limit=2.5)
    phase_advances = np.remainder(np.diff(np.angle(coefs[tune_mask]) / PI2), 1.0) # on (0,1)
    betas_unscaled = np.square(np.abs(coefs[tune_mask]))
    return phase_advances, betas_unscaled, tune_mask

from os import listdir


# onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
def _find_measurement_files(meas_path):
    onlyfiles = [f for f in listdir(meas_path) if os.path.isfile(os.path.join(meas_path, f))]
