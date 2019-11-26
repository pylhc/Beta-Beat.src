import sys
import math
import copy
import numpy as np
from numpy import sin, cos, tan
from collections import OrderedDict
import pandas as pd
import utils.bpm
import phase
from constants import PI, TWOPI, kEPSILON
from tfs_files import tfs_pandas
import numpy as np
import matplotlib.pyplot as plt
DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


def get_beta_amplitude(amplitude,action,action_error,calibration,calibration_error):
    action_matrix = np.matrix(action)
    action_matrix_error = np.matrix(action_error)
    action_matrix_transpose = np.transpose(action_matrix)
    action_matrix_error_transpose = np.transpose(action_matrix_error)
    amplitude_normalized = np.array(amplitude)**2*np.array(1/action_matrix_transpose)
    beta = np.array(amplitude_normalized)
    beta_error = np.array(beta)*np.array(action_matrix_error_transpose)*(np.array(1/action_matrix_transpose))
    beta_error_calibration = np.array(beta) * np.array(calibration_error) * np.array(1/calibration)
    beta_average = beta.mean(0)
    beta_std = beta.std(0)
    beta_error_average = (beta_error).mean(0)
    beta_error_calibration_average = beta_error_calibration.mean(0)
    beta_error_total = (beta_error_average **2 + beta_error_calibration_average **2 + (beta_std/len(beta)**0.5)**2)**0.5
    results = [beta_average,beta_error_total]
    return results


def get_beta_beating_rms(beta_phase,beta_phase_error,beta_model):
    beta_beating = np.array((beta_phase-beta_model)/beta_model)
    beta_beating_error = np.array(beta_phase_error/beta_model)
    beta_beating_average = np.mean(beta_beating)
    beta_beating_rms = np.sqrt(np.average(np.array(beta_beating)**2))
    beta_beating_std = np.std(beta_beating)
    beta_beating_average_error = np.sqrt(np.sum(np.array(beta_phase_error/beta_model)**2))/len(beta_beating)
    return beta_beating_average,beta_beating_average_error,beta_beating_rms,beta_beating



def get_beta_from_amp_eq(Psid,amp,action_phase,action_phase_error,action_model,action_model_error,calibration,calibration_error,accel):
    beta_result_phase_average,beta_result_model_average = get_beta_values(Psid,amp,action_phase,action_phase_error,action_model,action_model_error,calibration,calibration_error)
    beta_beating_average,beta_beating_average_error,beta_beating_rms,beta_beating = get_beta_beating_rms(beta_result_phase_average[0],beta_result_phase_average[1],beta_result_model_average[0])
    return beta_result_phase_average,beta_result_model_average,beta_beating_average


def get_beta_values(Psid,amp,action_phase,action_phase_error,action_model,action_model_error,calibration,calibration_error):
    beta_result_phase_average =  get_beta_amplitude(amp,action_phase,action_phase_error,calibration,calibration_error)
    beta_result_model_average =  get_beta_amplitude(amp,action_model,action_model_error,calibration,calibration_error)
    return beta_result_phase_average,beta_result_model_average




