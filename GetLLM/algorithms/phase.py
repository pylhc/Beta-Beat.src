'''
Created on 27 May 2013

@author: vimaier

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

import Utilities.bpm


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


#===================================================================================================
# helper-functions
#===================================================================================================
#------------ Get phases
#TODO: awful name! what does this function??? (vimaier)
def phiLastAndLastButOne(phi, ftune):
    if ftune > 0.0:
        phit = phi+ftune
        if phit > 1.0:
            phit = phit-1.0
    elif ftune <= 0.0:
        phit = phi+(1.0+ftune)
        if phit > 1.0:
            phit = phit-1.0
    return phit

def phase_mean(phase0, norm):  #-- phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi
    phase0 = np.array(phase0)%norm
    phase1 = (phase0+0.5*norm)%norm - 0.5*norm
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

def calc_phase_std(phase0, norm):  #-- phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi
    phase0 = np.array(phase0)%norm
    phase1 = (phase0+0.5*norm)%norm-0.5*norm
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
    phase_std = math.sqrt(min_phase_std/len(phase0))
    
    return phase_std


def get_phases_total(mad_twiss, src_files, tune, plane, beam_direction, accel, lhc_phase):
    commonbpms = Utilities.bpm.intersect(src_files)
    commonbpms = Utilities.bpm.model_intersect(commonbpms, mad_twiss)
    #-- Last BPM on the same turn to fix the phase shift by tune for exp data of LHC
    s_lastbpm = None
    if lhc_phase == "1" and accel == "LHCB1": 
        s_lastbpm = mad_twiss.S[mad_twiss.indx['BPMSW.1L2.B1']]
    if lhc_phase == "1" and accel == "LHCB2": 
        s_lastbpm = mad_twiss.S[mad_twiss.indx['BPMSW.1L8.B2']]

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
            phi12 = 1.0-phi12

        #phstd12=math.sqrt(np.average(phi12*phi12)-(np.average(phi12))**2.0+2.2e-15)
        #phi12=np.average(phi12)
        phstd12 = calc_phase_std(phi12, 1.0)
        phi12   = phase_mean(phi12, 1.0)
        phase_t[bn2] = [phi12, phstd12, phmdl12, bn1]

    return [phase_t, commonbpms]


def get_phases(getllm_d, mad_twiss, list_of_files, tune_q, plane):
    commonbpms = Utilities.bpm.intersect(list_of_files)
    commonbpms = Utilities.bpm.model_intersect(commonbpms, mad_twiss)
    length_commonbpms = len(commonbpms)

    #-- Last BPM on the same turn to fix the phase shift by tune_q for exp data of LHC
    if getllm_d.lhc_phase == "1" and getllm_d.accel == "LHCB1": 
        s_lastbpm = mad_twiss.S[mad_twiss.indx['BPMSW.1L2.B1']]
    if getllm_d.lhc_phase == "1" and getllm_d.accel == "LHCB2": 
        s_lastbpm = mad_twiss.S[mad_twiss.indx['BPMSW.1L8.B2']]

    mu = 0.0
    tunem = []
    phase = {} # Dictionary for the output containing [average phase, rms error]
    
    p_mdl_12 = 0.0 # phase model between BPM1 and BPM2
    p_mdl_13 = 0.0 # phase model between BPM1 and BPM3
    
    for i in range(0, length_commonbpms): # To find the integer part of tune as well, the loop is up to the last monitor
        bn1 = str.upper(commonbpms[i%length_commonbpms][1])
        bn2 = str.upper(commonbpms[(i+1)%length_commonbpms][1])
        bn3 = str.upper(commonbpms[(i+2)%length_commonbpms][1])
        
        if bn1 == bn2 :
            print >> sys.stderr, "There seem two lines with the same BPM name "+bn1+" in linx/y file."
            print >> sys.stderr, "Please check your input data....leaving GetLLM."
            sys.exit(1)
            
        if plane == 'H':
            p_mdl_12 = mad_twiss.MUX[mad_twiss.indx[bn2]] - mad_twiss.MUX[mad_twiss.indx[bn1]]
            p_mdl_13 = mad_twiss.MUX[mad_twiss.indx[bn3]] - mad_twiss.MUX[mad_twiss.indx[bn1]]
        elif plane == 'V':
            p_mdl_12 = mad_twiss.MUY[mad_twiss.indx[bn2]] - mad_twiss.MUY[mad_twiss.indx[bn1]]
            p_mdl_13 = mad_twiss.MUY[mad_twiss.indx[bn3]] - mad_twiss.MUY[mad_twiss.indx[bn1]]
            
        if i == length_commonbpms-2:
            if plane == 'H':
                madtune = mad_twiss.Q1 % 1.0
            elif plane == 'V':
                madtune = mad_twiss.Q2 % 1.0
            
            if madtune > 0.5:
                madtune -= 1.0
            
            p_mdl_13 = p_mdl_13 % 1.0
            p_mdl_13 = phiLastAndLastButOne(p_mdl_13, madtune)
        elif i == length_commonbpms-1:
            if plane == 'H':
                madtune = mad_twiss.Q1 % 1.0
            elif plane == 'V':
                madtune = mad_twiss.Q2 % 1.0
            
            if madtune > 0.5:
                madtune -= 1.0
            
            p_mdl_12 = p_mdl_12 % 1.0
            p_mdl_13 = p_mdl_13 % 1.0
            p_mdl_12 = phiLastAndLastButOne(p_mdl_12, madtune)
            p_mdl_13 = phiLastAndLastButOne(p_mdl_13, madtune)




        p_i_12 = [] # list with phases between BPM1 and BPM2 from src files
        p_i_13 = [] # list with phases between BPM1 and BPM3 from src files
        tunemi = []
        for src_twiss in list_of_files:
            # Phase is in units of 2pi
            if plane == 'H':
                p_m_12 = (src_twiss.MUX[src_twiss.indx[bn2]]-src_twiss.MUX[src_twiss.indx[bn1]]) # the phase advance between BPM1 and BPM2
                p_m_13 = (src_twiss.MUX[src_twiss.indx[bn3]]-src_twiss.MUX[src_twiss.indx[bn1]]) # the phase advance between BPM1 and BPM3
                tunemi.append(src_twiss.TUNEX[src_twiss.indx[bn1]])
            elif plane == 'V':
                p_m_12 = (src_twiss.MUY[src_twiss.indx[bn2]]-src_twiss.MUY[src_twiss.indx[bn1]]) # the phase advance between BPM1 and BPM2
                p_m_13 = (src_twiss.MUY[src_twiss.indx[bn3]]-src_twiss.MUY[src_twiss.indx[bn1]]) # the phase advance between BPM1 and BPM3
                tunemi.append(src_twiss.TUNEY[src_twiss.indx[bn1]])
            
            #-- To fix the phase shift by tune_q in LHC
            try:
                if mad_twiss.S[mad_twiss.indx[bn1]] <= s_lastbpm:
                    if mad_twiss.S[mad_twiss.indx[bn2]] > s_lastbpm: 
                        p_m_12 += getllm_d.beam_direction*tune_q
                    if mad_twiss.S[mad_twiss.indx[bn3]] > s_lastbpm: 
                        p_m_13 += getllm_d.beam_direction*tune_q
                if mad_twiss.S[mad_twiss.indx[bn1]] > s_lastbpm:
                    if mad_twiss.S[mad_twiss.indx[bn2]] <= s_lastbpm: 
                        p_m_12 += -getllm_d.beam_direction*tune_q
                    if mad_twiss.S[mad_twiss.indx[bn3]] <= s_lastbpm: 
                        p_m_13 += -getllm_d.beam_direction*tune_q
            except UnboundLocalError: 
                pass # s_lastbpm is not always defined
            
            if p_m_12 < 0:
                p_m_12 += 1
            if p_m_13 < 0:
                p_m_13 += 1
            p_i_12.append(p_m_12)
            p_i_13.append(p_m_13)

        p_i_12 = np.array(p_i_12)
        p_i_13 = np.array(p_i_13)
        if getllm_d.beam_direction == -1: # for the beam circulating reversely to the model
            p_i_12 = 1.0-p_i_12
            p_i_13 = 1.0-p_i_13

        #if any(p_i_12)>0.9 and i !=len(commonbpms): # Very small phase advance could result in larger than 0.9 due to measurement error
        #       print 'Warning: there seems too large phase advance! '+bn1+' to '+bn2+' = '+str(p_i_12)+'in plane '+plane+', recommended to check.'
        p_std_12 = calc_phase_std(p_i_12, 1.0)
        p_std_13 = calc_phase_std(p_i_13, 1.0)
        p_i_12 = phase_mean(p_i_12, 1.0)
        p_i_13 = phase_mean(p_i_13, 1.0)
        tunemi = np.array(tunemi)
        if i < length_commonbpms-1 :
            tunem.append(np.average(tunemi))

        # Note that the phase advance between the last monitor and the first monitor should be find by taking into account the fractional part of tune.
        if i == length_commonbpms-2:
            tunem = np.array(tunem)
            tune = np.average(tunem)
            p_i_13 = phiLastAndLastButOne(p_i_13, tune)
        elif i == length_commonbpms-1:
            p_i_12 = phiLastAndLastButOne(p_i_12, tune)
            p_i_13 = phiLastAndLastButOne(p_i_13, tune)
        mu = mu+p_i_12

        small = 0.0000001
        if (abs(p_mdl_12) < small):
            p_mdl_12 = small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn2+" in MAD model is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        if (abs(p_mdl_13) < small):
            p_mdl_13 = small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn3+" in MAD model is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        if (abs(p_i_12) < small ):
            p_i_12 = small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn2+" in measurement is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        if (abs(p_i_13) < small):
            p_i_13 = small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn3+" in measurement is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        phase[bn1] = [p_i_12, p_std_12, p_i_13, p_std_13, p_mdl_12, p_mdl_13, bn2]

    return [phase, tune, mu, commonbpms]


#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def getfreephase(phase,Qac,Q,bpms,MADTwiss_ac,MADTwiss,plane):
    '''
    :Parameters:
        'phase': dict
            (bpm_name:string) --> (phase_list:[phi12,phstd12,phi13,phstd13,phmdl12,phmdl13,bn2])
            phi13, phstd13, phmdl12 and phmdl13 are note used.
        
    '''

    if DEBUG:
        print "Calculating free phase using model"

    phasef={}
    phi=[]

    for bpm in bpms:
        bn1=bpm[1].upper()

        phase_list = phase[bn1]
        phi12 = phase_list[0]
        phstd12 = phase_list[1]
        bn2 =phase_list[6]
        bn2s=MADTwiss.S[MADTwiss.indx[bn2]]
        #model ac
        if plane=="H":
            ph_ac_m=MADTwiss_ac.MUX[MADTwiss_ac.indx[bn2]]-MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]]
            ph_m=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
        else:
            ph_ac_m=MADTwiss_ac.MUY[MADTwiss_ac.indx[bn2]]-MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]]
            ph_m=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]

        # take care the last BPM
        if bn1==bpms[-1][1].upper():
            ph_ac_m+=Qac; ph_ac_m=ph_ac_m%1
            ph_m   +=Q  ; ph_m   =ph_m%1

        phi12f=phi12-(ph_ac_m-ph_m)
        phi.append(phi12f)
        phstd12f=phstd12
        phmdl12f=ph_m

        phasef[bn1]=phi12f,phstd12f,phmdl12f,bn2,bn2s

    mu=sum(phi)


    return phasef,mu,bpms

def get_free_phase_total(phase,bpms,plane,MADTwiss,MADTwiss_ac):
    '''
    :Parameters:
        'phase': dict
            (bpm_name:string) --> (phase_list:[phi12,phstd12,phmdl12,bn1])
            phmdl12 and bn1 are note used.
        
    '''
    if DEBUG:
        print "Calculating free total phase using model"

    first=bpms[0][1]

    phasef={}

    for bpm in bpms:
        bn2=bpm[1].upper()

        if plane=="H":

            ph_ac_m=(MADTwiss_ac.MUX[MADTwiss_ac.indx[bn2]]-MADTwiss_ac.MUX[MADTwiss_ac.indx[first]])%1
            ph_m=(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[first]])%1

        else:
            ph_ac_m=(MADTwiss_ac.MUY[MADTwiss_ac.indx[bn2]]-MADTwiss_ac.MUY[MADTwiss_ac.indx[first]])%1
            ph_m=(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[first]])%1

        phase_list = phase[bn2]
        phi12 = phase_list[0]
        phstd12 = phase_list[1]


        phi12=phi12-(ph_ac_m-ph_m)
        phstd12=phstd12

        phasef[bn2]=phi12,phstd12,ph_m

    return phasef,bpms