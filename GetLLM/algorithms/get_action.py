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


def get_bpm_names_kick(MADTwiss_ac,Files,plane,accel):
    all_bpms = utils.bpm.model_intersect(
        utils.bpm.intersect(Files),
        MADTwiss_ac,
    )
    all_bpms = [(b[0], str.upper(b[1])) for b in all_bpms]
    all_bpms_names = [str.upper(b[1]) for b in all_bpms]
    # select BPMs needed to compute action
    if accel == "LHCB1" or accel == "LHCB2":
        good_bpms_for_kick_all,ir_bpms = intersect_bpm_list_with_arc_bpms(all_bpms)
        good_bpms_for_kick = [str.upper(b[1]) for b in good_bpms_for_kick_all]
    elif accel == "PSB":
        good_bpms_for_kick_all = intersect_bpm_list_psb(all_bpms)
        good_bpms_for_kick = [str.upper(b[1]) for b in good_bpms_for_kick_all]
    else:
        good_bpms_for_kick = all_bpms_names
    return good_bpms_for_kick


def intersect_bpm_list_with_arc_bpms(bpms_list):
    bpm_arcs = []
    bpm_irs = []
    # Selecting ARC BPMs
    for b in bpms_list:
        if (((b[1][4]) == '1' and (b[1][5]) >= '4') or
                (b[1][4]) == '2' or
                (b[1][4]) == '3'):
            bpm_arcs.append(b)
        else:
            bpm_irs.append(b)
    return bpm_arcs,bpm_irs


def intersect_bpm_list_psb(bpm_list, accel):
    all_bpms_filter = []
    for i in range(len(bpm_list)):
             bpm_end_list = ["4L3","6L3","11L3","12L3","14L3"]
             if bpm_list[i][1][-3:] not in bpm_end_list:
                  all_bpms_filter.append(bpm_list[i])
    ''' For non LHC it is called instead of intersect_bpm_list_with_arc_bpms '''
    return all_bpms_filter


BAD_BPM_LIST = ['BPM.15R8.B1', 'BPM.16R3.B1', 'BPM.31L5.B1', 'BPM.23L6.B1',
                'BPM.22R8.B1', 'BPM.11R6.B1', 'BPM.18R6.B1', 'BPM.34R5.B1']


def intersect_bpms_list_with_bad_known_bpms(bpms_list):
    bpm_arcs_clean = []
    for b in bpms_list:
        if b[1] not in BAD_BPM_LIST:
             bpm_arcs_clean.append(b)
    import sys
    return bpm_arcs_clean

def get_beta_phase_values(good_bpms_for_kick,beta_d):
    #-- Beta phase
    # Selecting beta phase in the arcs (not _free)
    beta_phase = np.array(
        [beta_d[bpm][0] for bpm in good_bpms_for_kick]
    )
    beta_phase_error = np.array(
        [beta_d[bpm][3] for bpm in good_bpms_for_kick]
    )
    return beta_phase,beta_phase_error
    # Selecting beta model for the arcs and for all the BPMS
    #-- Model beta and phase advance


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


def get_amplitude_calibration_and_beta(MADTwiss_ac,beta_d, Files, plane, bd,accel,all_bpms):
    [beta_phase,beta_phase_errors] = get_beta_phase_values(all_bpms,beta_d)
    beta_model = get_beta_model_values(all_bpms,MADTwiss_ac,plane)
    amp = get_amplitude(Files,all_bpms,plane,bd)
    #amp,calibration,calibration_error,psid = get_amplitude_calibration_and_phase(Files,all_bpms,plane,bd)
    #return amp,calibration,calibration_error,psid,beta_phase,beta_phase_errors,beta_model
    return amp,beta_phase,beta_phase_errors,beta_model




def get_amplitude(Files,all_bpms,plane,bd):
    amp = []
    calibration = []
    calibration_error = []
    psid = []
    for i in range(len(Files)):
        if plane == 'H':
            amp.append(np.array(
                    [2 * Files[i].AMPX[Files[i].indx[b]] for b in all_bpms]
            ))
        if plane == 'V':
            amp.append(np.array(
                    [2 * Files[i].AMPY[Files[i].indx[b]] for b in all_bpms]
            ))
    amp = np.array(amp)
    amp_transpose = amp.transpose()
    return amp
    

def get_kick_from_bpm_list_w_ACdipole_phase(amp, beta_phase,beta_phase_error):
    '''
    @author: F Carlier, A Garcia-Tabares
    Function calculates kick from measurements with AC dipole using the amplitude of the main line. The main line
    amplitude is obtained from Drive/SUSSIX and is normalized with the model beta-function.

    '''
    beta_phase_filter = np.array(beta_phase)[(np.array(beta_phase) > -1)]
    actions = (((amp).transpose() * (1/(beta_phase.transpose()**0.5))[:, np.newaxis]).mean(0))**2
    actions_error_spread = ((amp).transpose() * (1/(beta_phase.transpose()**0.5))[:, np.newaxis]).std(0)
    actions_error_std = (actions_error_spread/len(amp[0])**0.5)
    actions_error_phase = ((((amp**2).transpose() * 1/(beta_phase.transpose()**2)[:, np.newaxis] * beta_phase_error.transpose()[:, np.newaxis])**2).sum(0))**0.5/len(amp[0])
    return actions,actions_error_spread,actions_error_std,actions_error_phase


def get_kick_from_bpm_list_w_ACdipole_model(amp, beta_model):
    actions =  (((amp).transpose() * (1/(beta_model.transpose()**0.5))[:, np.newaxis]).mean(0))**2
    actions_error_spread = (((amp).transpose() * (1/(beta_model.transpose()**0.5))[:, np.newaxis]).std(0))**2
    actions_error_std = (actions_error_spread/len(amp[0])**0.5)
    return actions,actions_error_spread,actions_error_std


def get_action_values(
       amp_arcs,beta_phase_arcs,beta_phase_error_arcs,beta_model_arcs):
    action_phase,action_phase_std,action_phase_error,actions_phase_error_beta = get_kick_from_bpm_list_w_ACdipole_phase(
       amp_arcs,beta_phase_arcs,beta_phase_error_arcs)
    action_model,action_model_std,action_model_error = get_kick_from_bpm_list_w_ACdipole_model(
       amp_arcs,beta_model_arcs)
    return action_phase,action_phase_error,action_model,action_model_error


def get_action(MADTwiss_ac,beta_d, Files, plane, bd,accel):
    good_bpms_for_kick = get_bpm_names_kick(MADTwiss_ac,Files,plane,accel)
    amp_kick,beta_phase,beta_phase_errors,beta_model = get_amplitude_calibration_and_beta(MADTwiss_ac,beta_d, Files, plane, bd,accel,good_bpms_for_kick)
    action_phase,action_phase_error,action_model,action_model_error = get_action_values(
       amp_kick,beta_phase,beta_phase_errors,beta_model)
    action_summary_phase = [action_phase,action_phase_error]
    action_summary_model = [action_model,action_model_error]
    return action_summary_phase,action_summary_model


