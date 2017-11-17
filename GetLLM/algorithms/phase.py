'''
.. module: phase
Created on 27 May 2013

@author: ?, vimaier

@version: 0.0.1

GetLLM.algorithms.phase.py stores helper functions for phase calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import traceback
import math
import numpy as np
import compensate_ac_effect
from Utilities import tfs_file_writer
from model.accelerators.accelerator import AccExcitationMode
from constants import PI, TWOPI, kEPSILON
from Utilities import logging_tools

import pandas as pd
from time import time

DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)
UNION = True
LOGGER = logging_tools.get_logger(__name__)

#===================================================================================================
# main part
#===================================================================================================

class PhaseData(object):
    ''' File for storing results from get_phases.
        Storing results from getphases_tot are not needed since values not used again.
    '''

    def __init__(self):
        self.ac2bpmac_x = None
        self.ac2bpmac_y = None

        self.phase_advances_x = None # horizontal phase
        self.phase_advances_free_x = None
        self.phase_advances_free2_x = None
        self.phase_advances_y = None # horizontal phase
        self.phase_advances_free_y = None
        self.phase_advances_free2_y = None


def calculate_phase(getllm_d, twiss_d, tune_d, model, model_driven, elements, files_dict):
    '''
    Calculates phase and fills the following TfsFiles:
        ``getphasex.out        getphasex_free.out        getphasex_free2.out``
        ``getphasey.out        getphasey_free.out        getphasey_free2.out``

    :Parameters:
        'getllm_d': GetllmData (In-param, values will only be read)
            lhc_phase, accel and beam_direction are used.
        'twiss_d': TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': TuneData (In/Out-param, values will be read and set)
            Holds tunes and phase advances

    :Return: PhaseData, _TuneData
        an instance of PhaseData with the result of this function
        the same instance as param tune_d to indicate changes in the instance.
        
        The phase data will be a pandas.Panel with 3 dataframes ``MEAS``, ``MODEL``, ``ERRMEAS``
        
        ``phase_d.phase_advances_free_x[MEAS]``:
            
        +------++--------+--------+--------+--------+
        |      ||  BPM1  |  BPM2  |  BPM3  |  BPM4  | 
        +======++========+========+========+========+
        | BPM1 ||   0    | phi_21 | phi_31 | phi_41 | 
        +------++--------+--------+--------+--------+
        | BPM2 || phi_12 |    0   | phi_32 | phi_42 | 
        +------++--------+--------+--------+--------+
        | BPM3 || phi_13 | phi_23 |   0    | phi_43 | 
        +------++--------+--------+--------+--------+
        
        The phase advance between BPM_i and BPM_j can be obtained via::
            
            phi_ij = phase_advances.loc["MEAS", "BPM_i", "BPM_j"]
    '''
    # get common bpms
    phase_d = PhaseData()
    if UNION:
        bpmsx = twiss_d.zero_dpp_unionbpms_x
        bpmsy = twiss_d.zero_dpp_unionbpms_y
    else:
        bpmsx = twiss_d.zero_dpp_commonbpms_x
        bpmsy = twiss_d.zero_dpp_commonbpms_y



    print 'Calculating phase'
    
    # ============= Calculate tunes ===================================================================================

    if twiss_d.has_zero_dpp_x():
        #-- Calculate tune_x from files, weighted average based on rms
        q1_files = np.zeros(len(twiss_d.zero_dpp_x))
        q1_inv_Var = np.zeros(len(twiss_d.zero_dpp_x))
        for i, twiss_file in enumerate(twiss_d.zero_dpp_x):
            q1_files[i] = twiss_file.headers["Q1"]
            q1rms = twiss_file.headers["Q1RMS"]
            if q1rms == 0:
                q1rms = 1000
            q1_inv_Var[i] = 1.0 / float(q1rms) ** 2
        q1 = np.sum(q1_files * q1_inv_Var) / np.sum(q1_inv_Var)
        tune_d.q1 = q1
        tune_d.q1f = q1
        
        phase_d.phase_advances_free_x, tune_d.mux = get_phases(getllm_d, model_driven, twiss_d.zero_dpp_x, bpmsx, q1, 'H')
        if not twiss_d.has_zero_dpp_y():
            print 'liny missing and output x only ...'

    
    if twiss_d.has_zero_dpp_y():
        #-- Calculate tune_x from files
        q2_files = np.zeros(len(twiss_d.zero_dpp_y))
        q2_inv_Var = np.zeros(len(twiss_d.zero_dpp_y))
        for i, twiss_file in enumerate(twiss_d.zero_dpp_y):
            q2_files[i] = twiss_file.headers["Q2"]
            q2rms = twiss_file.headers["Q2RMS"]
            if q2rms == 0:
                q2rms = 1000
            q2_inv_Var[i] = 1.0 / float(q2rms) ** 2
        q2 = np.sum(q2_files * q2_inv_Var) / np.sum(q2_inv_Var)
        tune_d.q2 = q2
        tune_d.q2f = q2

        phase_d.phase_advances_free_y, tune_d.muy = get_phases(getllm_d, model_driven, twiss_d.zero_dpp_y, bpmsy, q2, 'V')
        if not twiss_d.has_zero_dpp_x():
            print 'linx missing and output y only ...'

    # ============= Calculate the phases ==============================================================================

    #---- ac to free phase from eq and the model
    if getllm_d.accelerator.excitation != AccExcitationMode.FREE:
        if twiss_d.has_zero_dpp_x():
            tune_d.q1f = tune_d.q1 - getllm_d.accelerator.drv_tune_x + getllm_d.accelerator.nat_tune_x#-- Free H-tune
            phase_d.phase_advances_x = phase_d.phase_advances_free_x
            phase_d.ac2bpmac_x = compensate_ac_effect.GetACPhase_AC2BPMAC(model_driven, bpmsx, tune_d.q1, tune_d.q1f, 'H', getllm_d)
            [phase_d.phase_advances_free_x, tune_d.muxf] = compensate_ac_effect.get_free_phase_eq(model, twiss_d.zero_dpp_x, twiss_d.zero_dpp_commonbpms_x,
                                                                                                    tune_d.q1, tune_d.q1f, phase_d.ac2bpmac_x, 'H', model.Q1 % 1.0, getllm_d)
#            [phase_d.phase_advances_free2_x, tune_d.muxf2] = _get_free_phase(phase_d.phase_advances_free_x, tune_d.q1, tune_d.q1f, bpmsx, model_driven, model, "H")
        if twiss_d.has_zero_dpp_y():
            phase_d.phase_advances_y = phase_d.phase_advances_free_y
            tune_d.q2f =  tune_d.q2 - getllm_d.accelerator.drv_tune_y + getllm_d.accelerator.nat_tune_y #-- Free V-tune
            phase_d.ac2bpmac_y = compensate_ac_effect.GetACPhase_AC2BPMAC(elements, bpmsy, tune_d.q2, tune_d.q2f, 'V', getllm_d)
            [phase_d.phase_advances_free_y, tune_d.muyf] = compensate_ac_effect.get_free_phase_eq(model, twiss_d.zero_dpp_y, twiss_d.zero_dpp_commonbpms_y,
                                                                                                    tune_d.q2, tune_d.q2f, phase_d.ac2bpmac_y, 'V',
                                                                                                    model.Q2%1, getllm_d)
#            [phase_d.phase_advances_free2_y, tune_d.muyf2] = _get_free_phase(phase_d.phase_advances_free_y, tune_d.q2, tune_d.q2f, bpmsy, model_driven, model, "V")


    # ============= Write the phases to file ==========================================================================

    #---- H plane result
    print "output files"
    if twiss_d.has_zero_dpp_x():
        print "x ouptut"
        files_dict["getphasex_free.out"] = write_phase_file(files_dict["getphasex_free.out"], "H", phase_d.phase_advances_free_x, model, elements, tune_d.q1f, tune_d.q2f, getllm_d.accelerator)
        files_dict["getphasetotx_free.out"] = write_phasetot_file(files_dict["getphasetotx_free.out"], "H", phase_d.phase_advances_free_x, model, elements, tune_d.q1f, tune_d.q2f, getllm_d.accelerator)
        #-- ac to free phase
        if getllm_d.accelerator.excitation != AccExcitationMode.FREE:
            #-- from eq
            files_dict["getphasex.out"] = write_phase_file(files_dict["getphasex.out"], "H", phase_d.phase_advances_x, model, elements, tune_d.q1, tune_d.q2, getllm_d.accelerator)
            files_dict["getphasetotx.out"] = write_phasetot_file(files_dict["getphasetotx.out"], "H", phase_d.phase_advances_x, model, elements, tune_d.q1, tune_d.q2, getllm_d.accelerator)

    #---- V plane result
    if twiss_d.has_zero_dpp_y():
        files_dict["getphasey_free.out"] = write_phase_file(files_dict["getphasey_free.out"], "V", phase_d.phase_advances_free_y, model, elements, tune_d.q1f, tune_d.q2f, getllm_d.accelerator)
        files_dict["getphasetoty_free.out"] = write_phasetot_file(files_dict["getphasetoty_free.out"], "V", phase_d.phase_advances_free_y, model, elements, tune_d.q1f, tune_d.q2f, getllm_d.accelerator)
        #-- ac to free phase
        if getllm_d.accelerator.excitation != AccExcitationMode.FREE:
            #-- from eq
            files_dict["getphasey.out"] = write_phase_file(files_dict["getphasey.out"], "V", phase_d.phase_advances_y, model, elements, tune_d.q1, tune_d.q2, getllm_d.accelerator)
            files_dict["getphasetoty.out"] = write_phasetot_file(files_dict["getphasetoty.out"], "V", phase_d.phase_advances_y, model, elements, tune_d.q1, tune_d.q2, getllm_d.accelerator)

    return phase_d, tune_d
# END calculate_phase ------------------------------------------------------------------------------

#===================================================================================================
# helper-functions
#===================================================================================================
#TODO: awful name! what does this function??? (vimaier)
def _phi_last_and_last_but_one(phi, ftune):
    if ftune <= 0:
        ftune += 1
    phi += ftune
    if phi > 1:
        phi -= 1
    return phi

def t_value_correction(_num):
    ''' Calculations are based on Hill, G. W. (1970)
    Algorithm 396: Student's t-quantiles. Communications of the ACM, 
    13(10), 619-620.

    http://en.wikipedia.org/wiki/Quantile_function#The_Student.27s_t-distribution

    It is not implemented directly here because a library for the erfinv() function, the inverse error function
    cannot be accessed from our servers in their current python installation (Jan-2015).
    (http://en.wikipedia.org/wiki/Error_function#Inverse_function)
    '''
    num = int(_num)
    correction_dict = {2:1.8394733927562799, 3:1.3224035682262103, 4:1.1978046912864673, 
                       5:1.1424650980932523, 6:1.1112993008590089, 7:1.0913332519214189, 
                       8:1.0774580800762166, 9:1.0672589736833817, 10:1.0594474783177483,
                       11:1.053273802733051, 12:1.0482721313740653, 13:1.0441378866779087,
                       14:1.0406635564353071, 15:1.0377028976401199, 16:1.0351498875115406,
                       17:1.0329257912610941, 18:1.0309709166064416, 19:1.029239186837585, 
                       20:1.0276944692596461}
    if num > 1 and num <=20:
        t_factor = correction_dict[num]
    else:
        t_factor = 1
    return t_factor

vec_t_value_correction = np.vectorize(t_value_correction, otypes=[int])

def calc_phase_mean(phase0, norm):
    ''' phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi '''
    phase0 = np.array(phase0)%norm
    phase1 = (phase0 + .5*norm) % norm - .5*norm
    phase0ave = np.mean(phase0)
    phase1ave = np.mean(phase1)
    # Since phase0std and phase1std are only used for comparing, I modified the expressions to avoid
    # math.sqrt(), np.mean() and **2.
    # Old expressions:
    #     phase0std = math.sqrt(np.mean((phase0-phase0ave)**2))
    #     phase1std = math.sqrt(np.mean((phase1-phase1ave)**2))
    # -- vimaier
    mod_phase0std = sum(abs(phase0-phase0ave))
    mod_phase1std = sum(abs(phase1-phase1ave))
    if mod_phase0std < mod_phase1std:
        return phase0ave
    else:
        return phase1ave % norm

def calc_phase_std(phase0, norm):
    ''' phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi '''
    phase0 = np.array(phase0)%norm
    phase1 = (phase0 + .5*norm) % norm - .5*norm
    phase0ave = np.mean(phase0)
    phase1ave = np.mean(phase1)

    # Omitted unnecessary computations. Old expressions:
    #     phase0std=sqrt(mean((phase0-phase0ave)**2))
    #     phase1std=sqrt(mean((phase1-phase1ave)**2))
    #     return min(phase0std,phase1std)
    # -- vimaier
    phase0std_sq = np.sum((phase0-phase0ave)**2)
    phase1std_sq = np.sum((phase1-phase1ave)**2)

    min_phase_std = min(phase0std_sq, phase1std_sq)
    if len(phase0) > 1:
        phase_std = math.sqrt(min_phase_std/(len(phase0)-1))
        phase_std = phase_std * t_value_correction(len(phase0)-1)
    else:
        phase_std = 0
    return phase_std

def _get_phases_total(mad_twiss, src_files, commonbpms, tune, plane, beam_direction, accel, lhc_phase):
    #-- Last BPM on the same turn to fix the phase shift by tune for exp data of LHC
    s_lastbpm = None
    if lhc_phase == "1":
        print "phase jump correction"
        if accel == "JPARC":
            print "-- no total phase jump correction with JPARC"
            #s_lastbpm = mad_twiss.S[mad_twiss.indx['MOH.3']]
        elif accel == "LHCB1":
            s_lastbpm = mad_twiss.S[mad_twiss.indx['BPMSW.1L2.B1']]
            print "-> for LHC"
        elif accel == "LHCB2":
            s_lastbpm = mad_twiss.S[mad_twiss.indx['BPMSW.1L8.B2']]
            print "-> for LHC"
        elif accel == "PETRA":
            s_lastbpm = mad_twiss.S[mad_twiss.indx['BPM_SOR_46']]
            print "-> for PETRA {0:f}".format(tune)
    
    bn1 = str.upper(commonbpms[0][1])
    phase_t = {}
    if DEBUG:
        print "Reference BPM:", bn1, "Plane:", plane
    for i in range(0, len(commonbpms)):
        bn2 = str.upper(commonbpms[i][1])
        if plane == 'H':
            phmdl12 = (mad_twiss.MUX[mad_twiss.indx[bn2]]-mad_twiss.MUX[mad_twiss.indx[bn1]]) % 1
        if plane == 'V':
            phmdl12 = (mad_twiss.MUY[mad_twiss.indx[bn2]]-mad_twiss.MUY[mad_twiss.indx[bn1]]) % 1

        phi12 = []
        for twiss_file in src_files:
            # Phase is in units of 2pi
            if plane == 'H':
                phm12 = (twiss_file.MUX[twiss_file.indx[bn2]]-twiss_file.MUX[twiss_file.indx[bn1]]) % 1
            if plane == 'V':
                phm12 = (twiss_file.MUY[twiss_file.indx[bn2]]-twiss_file.MUY[twiss_file.indx[bn1]]) % 1
            #-- To fix the phase shift by tune in LHC
            try:
                if s_lastbpm is not None and commonbpms[i][0] > s_lastbpm:
                    phm12 += beam_direction*tune
            except:
                traceback.print_exc()
            phi12.append(phm12)
        phi12 = np.array(phi12)
        # for the beam circulating reversely to the model
        if beam_direction == -1:
            phi12 = 1 - phi12

        #phstd12=math.sqrt(np.average(phi12*phi12)-(np.average(phi12))**2.0+2.2e-15)
        #phi12=np.average(phi12)
        phstd12 = calc_phase_std(phi12, 1.)
        phi12   = calc_phase_mean(phi12, 1.)
        phase_t[bn2] = [phi12, phstd12, phmdl12, bn1]

    return phase_t


#IMPORTANT_PAIRS = {"BPMYA.5R6.B2": ["BPMWB.4R5.B2", "BPMWB.4R1.B2"]}


def get_phases(getllm_d, mad_twiss, Files, bpm, tune_q, plane):
    """
    Calculates phase.
    tune_q will be used to fix the phase shift in LHC.
    
    ``phase_advances["MEAS"]`` contains the measured phase advances.
    
    ``phase_advances["MODEL"]`` contains the model phase advances.
    
    ``phase_advances["ERRMEAS"]`` contains the error of the measured phase advances as deterined by the standard
    deviation scaled by sqrt(number_of_files).::
    
        phase_advances.loc["MEAS", bpm_namei, bpm_namej]
    
    yields the phase advance ``phi_ij`` between BPMi and BPMj 
    """
    acc = getllm_d.accelerator
    plane_mu = "MUX" if plane == "H" else "MUY"
    bd = acc.get_beam_direction()
    number_commonbpms = bpm.shape[0]
    muave= 0  # TODO: find out what this is and who needs it
    
    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if getllm_d.lhc_phase == "1":
        print "correcting phase jump"
        k_lastbpm = acc.get_k_first_BPM(bpm.index)
    else:
        print "phase jump will not be corrected"
        k_lastbpm = len(bpm.index)

    if UNION:
        phase_advances = _get_phases_union(bpm, number_commonbpms, bd, plane_mu, mad_twiss, Files, k_lastbpm)
    else:
        phase_advances = _get_phases_intersection(bpm, number_commonbpms, bd, plane_mu, mad_twiss, Files, k_lastbpm)
    
    return phase_advances, muave

def _get_phases_intersection(bpm, number_commonbpms, bd, plane_mu, mad_twiss, Files, k_lastbpm): 

    # pandas panel that stores the model phase advances, measurement phase advances and measurement errors
    phase_advances = pd.Panel(items=["MODEL", "MEAS", "ERRMEAS"], major_axis=bpm.index, minor_axis=bpm.index)

    phases_mdl = np.array(mad_twiss.loc[bpm.index, plane_mu])
    phase_advances["MODEL"] = (phases_mdl[np.newaxis,:] - phases_mdl[:,np.newaxis]) % 1.0

    # loop over the measurement files
    phase_matr_meas = np.empty((len(Files), number_commonbpms, number_commonbpms))  # temporary 3D matrix that stores the phase advances
    for i in range(len(Files)):
        file_tfs = Files[i]
        phases_meas = bd * np.array(file_tfs.loc[bpm.index, plane_mu]) #-- bd flips B2 phase to B1 direction
        #phases_meas[k_lastbpm:] += tune_q  * bd
        meas_matr = (phases_meas[np.newaxis,:] - phases_meas[:,np.newaxis]) 
        phase_matr_meas[i] = np.where(meas_matr > 0, meas_matr, meas_matr + 1.0)
        
    phase_advances["MEAS"] = np.mean(phase_matr_meas, axis=0) % 1.0
    phase_advances["ERRMEAS"] = np.std(phase_matr_meas, axis=0) * t_value_correction(len(Files)) / np.sqrt(len(Files))

    return phase_advances

def _get_phases_union(bpm, number_commonbpms, bd, plane_mu, mad_twiss, Files, k_lastbpm):

    LOGGER.debug("calculating phases with union of measurement files")
    LOGGER.debug("maximum {:d} measurements per BPM".format(len(Files)))
    # pandas panel that stores the model phase advances, measurement phase advances and measurement errors
    phase_advances = pd.Panel(items=["MODEL", "MEAS", "ERRMEAS", "NFILES"], major_axis=bpm.index, minor_axis=bpm.index)

    phases_mdl = np.array(mad_twiss.loc[bpm.index, plane_mu])
    phase_advances["MODEL"] = (phases_mdl[np.newaxis,:] - phases_mdl[:,np.newaxis]) % 1.0
    
    # loop over the measurement files
    phase_matr_meas = pd.Panel(items=range(len(Files)), major_axis=bpm.index, minor_axis=bpm.index)  # temporary 3D matrix that stores the phase advances
    #phase_matr_count = pd.Panel(items=range(len(Files)), major_axis=bpm.index, minor_axis=bpm.index)
    for i in range(len(Files)):
        file_tfs = Files[i]
        phases_meas = bd * np.array(file_tfs.loc[:, plane_mu]) #-- bd flips B2 phase to B1 direction
        #phases_meas[k_lastbpm:] += tune_q  * bd
        
        meas_matr = (phases_meas[np.newaxis,:] - phases_meas[:,np.newaxis]) 
        phase_matr_meas.loc[i] = pd.DataFrame(data=np.where(meas_matr > 0, meas_matr, meas_matr + 1.0),
                                              index=file_tfs.index, columns=file_tfs.index)
        
    phase_matr_meas = phase_matr_meas.values
    mean = np.nanmean(phase_matr_meas, axis=0) % 1.0
    phase_advances["MEAS"] = mean
    nfiles = np.sum(~np.isnan(phase_matr_meas), axis=0)
    
    phase_advances["NFILES"] = nfiles
    phase_advances["ERRMEAS"] = np.nanstd(phase_matr_meas, axis=0) / np.sqrt(nfiles) * vec_t_value_correction(nfiles)

    LOGGER.debug("t_value correction.................[YES]")
    LOGGER.debug("optimistic errorbars...............[YES]")

    return phase_advances

#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def _get_free_phase(phase, tune_ac, tune, bpms, model_driven, model, plane):
    '''
    :Parameters:
        'phase': dict
            (bpm_name:string) --> (phase_list:[phi12,phstd12,phi13,phstd13,phmdl12,phmdl13,bn2])
            phi13, phstd13, phmdl12 and phmdl13 are note used.
    '''
    if DEBUG:
        print "Calculating free phase using model"
    raise NotImplementedError("new phase table algortihm not yet implemented")
    phasef = {}
    phi = []

    for bpm in bpms["NAME"]:
        bn1 = bpm.upper()

        phase_list = phase[bn1]
        phi12 = phase_list[0]
        phstd12 = phase_list[1]
        bn2 = phase_list[6]
        bn2s = model.S[model.indx[bn2]]
        #model ac
        if plane == "H":
            ph_ac_m = model_driven.MUX[model_driven.indx[bn2]]-model_driven.MUX[model_driven.indx[bn1]]
            ph_m = model.MUX[model.indx[bn2]]-model.MUX[model.indx[bn1]]
        else:
            ph_ac_m = model_driven.MUY[model_driven.indx[bn2]]-model_driven.MUY[mad_ac.indx[bn1]]
            ph_m = model.MUY[model.indx[bn2]]-model.MUY[model.indx[bn1]]

        # take care the last BPM
        if bn1 == bpms[-1][1].upper():
            ph_ac_m += tune_ac
            ph_ac_m = ph_ac_m % 1
            ph_m += tune
            ph_m = ph_m % 1

        phi12f = phi12-(ph_ac_m-ph_m)
        phi.append(phi12f)
        phstd12f = phstd12
        phmdl12f = ph_m

        phasef[bn1] = phi12f, phstd12f, phmdl12f, bn2, bn2s

    mu = sum(phi)

    return phasef, mu

def get_free_phase_total(phase, bpms, plane, mad_twiss, mad_ac):
    '''
    :Parameters:
        'phase': dict
            (bpm_name:string) --> (phase_list:[phi12,phstd12,phmdl12,bn1])
            phmdl12 and bn1 are note used.
    '''
    if DEBUG:
        print "Calculating free total phase using model"

    first = bpms[0][1]

    phasef = {}

    for bpm in bpms:
        bn2 = bpm[1].upper()

        if plane == "H":
            ph_ac_m = (mad_ac.MUX[mad_ac.indx[bn2]] - mad_ac.MUX[mad_ac.indx[first]]) % 1
            ph_m = (mad_twiss.MUX[mad_twiss.indx[bn2]] - mad_twiss.MUX[mad_twiss.indx[first]]) % 1
        else:
            ph_ac_m = (mad_ac.MUY[mad_ac.indx[bn2]] - mad_ac.MUY[mad_ac.indx[first]]) % 1
            ph_m = (mad_twiss.MUY[mad_twiss.indx[bn2]] - mad_twiss.MUY[mad_twiss.indx[first]]) % 1

        phase_list = phase[bn2]
        phi12 = phase_list[0]
        phstd12 = phase_list[1]

        phi12 = phi12-(ph_ac_m-ph_m)
        phstd12 = phstd12

        phasef[bn2] = phi12, phstd12, ph_m

    return phasef

#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def write_phase_file(tfs_file, plane, phase_advances, model, elements, tune_x, tune_y, accel):
    
    plane_char = "X" if plane == "H" else "Y"
    plane_mu = "MU" + plane_char
    plane_tune = tune_x if plane == "H" else tune_y
    
    tfs_file.add_float_descriptor("Q1", tune_x)
    tfs_file.add_float_descriptor("Q2", tune_y)
    tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "PHASE" + plane_char, "STDPH" + plane_char,
                               "PH{}MDL".format(plane_char), "MU{}MDL".format(plane_char), "NFILES"])
    if UNION:
        tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    else:
        tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le", "%s"])


    meas = phase_advances["MEAS"]
    mod = phase_advances["MODEL"]
    err = phase_advances["ERRMEAS"]

    if UNION:
        nfiles = phase_advances["NFILES"]
    bd = accel.get_beam_direction()
    
    intersected_model = model.loc[meas.index]

    for elem1, elem2 in accel.get_important_phase_advances():

        mus1 = elements.loc[elem1, plane_mu] - intersected_model.loc[:, plane_mu]
        minmu1 = abs(mus1).idxmin()
        
        mus2 = intersected_model.loc[:, plane_mu] - elements.loc[elem2, plane_mu]
        minmu2 = abs(mus2).idxmin()
        
        try:
            bpm_phase_advance = meas.loc[minmu1, minmu2]
            model_value = elements.loc[elem2, plane_mu] - elements.loc[elem1, plane_mu]

            if (elements.loc[elem2, "S"] - elements.loc[elem1, "S"]) * bd < 0.0:
                bpm_phase_advance += plane_tune
                model_value += plane_tune
            bpm_err = err.loc[minmu1, minmu2]
            phase_to_first = -mus1.loc[minmu1]
            phase_to_second = -mus2.loc[minmu2]
            
            ph_result = ((bpm_phase_advance + phase_to_first + phase_to_second) * bd)% 1.0
            
            model_value = (model_value * bd) % 1.0
            
            tfs_file.add_string_descriptor(elem1 + "__to__" + elem2 + "___MODL", 
                                          "{:8.4f}     {:6s} = {:6.2f} deg".format(model_value, "", (model_value) * 360))
            tfs_file.add_string_descriptor(elem1 + "__to__" + elem2 + "___MEAS", 
                                          "{:8.4f}  +- {:6.4f} = {:6.2f} +- {:3.2f} deg ({:8.4f} + {:8.4f} [{}, {}])".format(
                                                  ph_result, bpm_err, ph_result * 360, bpm_err * 360,
                                                  bpm_phase_advance,
                                                  phase_to_first + phase_to_second,
                                                  minmu1, minmu2) )
        except KeyError as e:
            LOGGER.error("Couldn't calculate the phase advance because " + e)
            
    for i in range(len(meas.index)-1):
        
        if UNION:
            nf = nfiles[meas.index[i+1]][meas.index[i]]
        else:
            nf = "ALL"
        
        tfs_file.add_table_row([
                 meas.index[i],
                 meas.index[i+1],
                 model.loc[meas.index[i], "S"],
                 model.loc[meas.index[i+1], "S"],
                 meas[meas.index[i+1]][meas.index[i]],
                 err[meas.index[i+1]][meas.index[i]],
                 mod[meas.index[i+1]][meas.index[i]],
                 model.loc[meas.index[i], plane_mu],
                 nf
                ])
    return tfs_file

def write_phasetot_file(tfs_file, plane, phase_advances, model, elements, tune_x, tune_y, accel):
    
    plane_char = "X" if plane == "H" else "Y"
    plane_mu = "MU" + plane_char
    plane_tune = tune_x if plane == "H" else tune_y
    meas = phase_advances["MEAS"]
    mod = phase_advances["MODEL"]
    err = phase_advances["ERRMEAS"]
    
    tfs_file.add_float_descriptor("Q1", tune_x)
    tfs_file.add_float_descriptor("Q2", tune_y)
    tfs_file.add_column_names(["NAME", "NAME2", "S", "S1", "PHASE" + plane_char, "STDPH" + plane_char, "PH{}MDL".format(plane_char), "MU{}MDL".format(plane_char)])
    tfs_file.add_column_datatypes(["%s", "%s", "%le", "%le", "%le", "%le", "%le", "%le"])

    
    for i in range(len(meas.index)-1):
        tfs_file.add_table_row([
                 meas.index[i+1],
                 meas.index[0],
                 model.loc[meas.index[i+1], "S"],
                 model.loc[meas.index[0], "S"],
                 meas.loc[meas.index[i+1]][meas.index[0]],
                 err.loc[meas.index[i+1]][meas.index[0]],
                 mod.loc[meas.index[i+1]][meas.index[0]],
                 model.loc[meas.index[i], plane_mu]
                ])
    return tfs_file
