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

def get_beta_model_values(all_bpms,MADTwiss_ac,plane):
    if plane == 'H':
        beta_model = np.array(
            [MADTwiss_ac.BETX[MADTwiss_ac.indx[b]] for b in all_bpms]
        )

    if plane == 'V':
        beta_model = np.array(
            [MADTwiss_ac.BETY[MADTwiss_ac.indx[b]] for b in all_bpms]
        )
    return beta_model

def get_free_beta_values(psid,psid_ac2bpmac,s_lastbpm,bpmac,k_bpmac,amp,action_phase,action_phase_error,action_model,action_model_error,calibration,calibration_error,tune_values_list,all_bpms):
    driven_tune = np.average(tune_values_list[0])
    natural_tune = np.average((tune_values_list[2]))
    natural_tune_error = np.average(tune_values_list[1])
    Psid = fix_phase_jump(psid,s_lastbpm,all_bpms,psid_ac2bpmac,bpmac,k_bpmac,driven_tune)
    r = sin(np.pi * (np.array(driven_tune) - np.array(natural_tune))) / sin(np.pi * (np.array(driven_tune) + np.array(natural_tune)))
    beta_result_phase =  get_beta_amplitude_compensated(Psid,amp,action_phase,action_phase_error,driven_tune,natural_tune,natural_tune_error,r,calibration,calibration_error)
    beta_result_model =  get_beta_amplitude_compensated(Psid,amp,action_model,action_model_error,driven_tune,natural_tune,natural_tune_error,r,calibration,calibration_error)
    return beta_result_phase,beta_result_model


def get_beta_amplitude_compensated(Psid,amplitude,action,action_error,driven_tune,natural_tune,natural_tune_error,r,calibration,calibration_error):
    action_matrix = np.matrix(action)
    action_matrix_error = np.matrix(action_error)
    action_matrix_transpose = np.transpose(action_matrix)
    action_matrix_error_transpose = np.transpose(action_matrix_error)
    compensation_matrix =  (1 + r ** 2 + 2 * r * np.cos(2 * np.transpose(Psid))) / (1 - r ** 2)
    amplitude_normalized = np.array(amplitude)**2*np.array(1/action_matrix_transpose)
    beta = np.array(amplitude_normalized)*np.array(compensation_matrix)
    beta_error = np.array(beta)*np.array(action_matrix_error_transpose)*(np.array(1/action_matrix_transpose))
    beta_error_calibration = np.array(beta) * np.array(calibration_error) * np.array(1/calibration)
    beta_average = beta.mean(0)
    beta_std = beta.std(0)
    beta_error_average = (beta_error).mean(0)
    beta_error_calibration_average = beta_error_calibration.mean(0)
    lambda_value = np.sin(np.pi*(np.array(driven_tune)-np.array(natural_tune)))/np.sin(np.pi*(np.array(driven_tune)+np.array(natural_tune)))
    beta_error_tune = (2*np.pi**2*(2*np.cos(np.pi*np.array(driven_tune))*np.sin(np.pi*np.array(driven_tune))/(np.sin(np.pi*(np.array(driven_tune)+np.array(natural_tune)))**2)**2)*np.array(natural_tune_error))**2
    beta_error_tune_complete = np.average(np.transpose(np.array(beta_error_tune))*beta*2/(1-lambda_value)**2)
    beta_error_total = (beta_error_average **2 + beta_error_calibration_average **2 + (beta_std/len(beta)**0.5)**2 + beta_error_tune_complete**2)**0.5
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

def fix_phase_jump(psid_all,s_lastbpm,all_bpms,psid_ac2bpmac,bpmac,k_bpmac,Qd):
    Psid_tot = []
    i = 0
    for psid in psid_all.transpose():
        for k in range(len(all_bpms)):
               try:
                   if all_bpms[k][0] > s_lastbpm:
                        psid[k] += 2 * np.pi * np.average(Qd)
               except:
                   pass
        psid = psid - (psid[k_bpmac] - psid_ac2bpmac[bpmac])
        Psid = psid + np.pi * np.average(Qd)
        Psid [k_bpmac:] = (Psid[k_bpmac:] - 2 * np.pi * np.average(Qd))
        Psid_tot.append(Psid)
        i = i+1
    return np.matrix(np.array(Psid_tot)).transpose()


def get_beta_from_amp_eq(MADTwiss_ac,Psid,amp,action_phase,action_phase_error,action_model,action_model_error,calibration,calibration_error,psid_ac2bpmac,tune_values_list,bd,op,good_bpms,plane):
    if bd ==  1 and op == 1:
        s_lastbpm = MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L2.B1']]
    elif bd == -1 and op == 1 :
        s_lastbpm = MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L8.B2']]
    else:
        s_lastbpm = MADTwiss_ac.S[-1]
    #-- Determine the BPM closest to the AC dipole and its position
    bpm_ac1 = psid_ac2bpmac.keys()[0]
    bpm_ac2 = psid_ac2bpmac.keys()[1]
    if bpm_ac1 in good_bpms:
        k_bpmac = good_bpms.index(bpm_ac1)
        bpmac = bpm_ac1
    elif bpm_ac2 in good_bpms:
        k_bpmac = good_bpms.index(bpm_ac2)
        bpmac = bpm_ac2
    else:
        print >> sys.stderr, 'WARN: BPMs next to AC dipoles missing. Was looking for: bpm_ac1="{0}" and bpm_ac2="{1}"'.format(bpm_ac1, bpm_ac2)
        return {},{},0.0, []
    beta_model = get_beta_model_values(good_bpms,MADTwiss_ac,plane)
    beta_result_phase,beta_result_model = get_free_beta_values(Psid,psid_ac2bpmac,s_lastbpm,bpmac,k_bpmac,amp,action_phase,action_phase_error,action_model,action_model_error,calibration,calibration_error,tune_values_list,good_bpms)
    beta_beating_average,beta_beating_average_error,beta_beating_rms,beta_beating = get_beta_beating_rms(beta_result_phase[0],beta_result_phase[1],beta_model)
    return beta_result_phase,beta_result_model,beta_beating_average

