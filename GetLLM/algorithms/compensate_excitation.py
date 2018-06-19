'''
Created on 27 May 2013

@author: ?

@version: 0.0.1

GetLLM.algorithms.compensate_ac_effect.py stores helper functions to compensate the AC dipole effect
based on analytic formulae (by R. Miyamoto).
This module is not intended to be executed. It stores only functions for GetLLM.


Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import math
import copy
import numpy as np
from numpy import sin, cos, tan
import re

import utils.bpm
import phase
from SegmentBySegment.SegmentBySegment import get_good_bpms
from __builtin__ import raw_input
from constants import PI, TWOPI, kEPSILON
from SegmentBySegment.sbs_writers.sbs_phase_writer import FIRST_BPM_B1
import pandas as pd
from utils import logging_tools

LOGGER = logging_tools.get_logger(__name__)

DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# helper-functions
#===================================================================================================
#---------  The following is functions

def GetACPhase_AC2BPMAC(commonbpms, Qd, Q, plane, acc):
    """Returns the necessary values for the exciter compensation.

    Args:
        model (pandas.DataFrame): model twiss as DataFrame.
        commonbpms (pandas.DataFrame): commonbpms (see GetLLM._get_commonbpms)
        Qd, Q (float): Driven and natural fractional tunes.
        plane (char): H, V
        accelerator: accelerator class instance.

    Returns tupel(a,b,c,d):
        a (string): name of the nearest BPM.
        b (float): compensated phase advance between the exciter and the nearest BPM.
        c (int): k of the nearest BPM.
        d (string): name of the exciter element.
    """
    bd = acc.get_beam_direction()
    r = sin(PI * (Qd - Q)) / sin(PI * (Qd + Q))
    plane_mu = "MUX" if plane == "H" else "MUY"
    [k, bpmac1], exciter = acc.get_exciter_bpm(plane, commonbpms)
    model_driven = acc.get_driven_tfs()
    ##return bpmac1, np.arctan((1 + r) / (1 - r) * tan(TWOPI * model.loc[bpmac1, plane_mu] + PI * Q)) % PI - PI * Qd, k
    try:
        psi = (
            np.arctan((1 + r) / (1 - r) *
                      tan(TWOPI * (model_driven.loc[bpmac1, plane_mu] - model_driven.loc[exciter, plane_mu]) )
                     ) / TWOPI
        ) % .5 - .5
    except:
        psi = acc.get_elements_tfs().loc[bpmac1, plane_mu] - acc.get_elements_tfs().loc[exciter, plane_mu]
    return bpmac1, psi, k, exciter


def get_free_phase_eq(model, Files, bpm, Qd, Q, ac2bpmac, plane, Qmdl, getllm_d):
    """
       Calculates phase free from excitation effect.

       Args:
           model (pandas.DataFrame): the model twiss file. Must have the names as index.
           Files (list of pandas.DataFrame): the input files.
           bpm (pandas.DataFrame): commonbpms (see GetLLM._get_commonbpms)
           Qd (float): driven fractional tune.
           Q (float): natural fractional tune.
           ac2bpmac: result of GetACPhase_AC2BPMAC.
           plane (char): H or V.
           Qmdl: natural model tune.
           getllm_d (GetLLM_Data): GetLLM_data created by GetLLM.

        Returns:
            pandas.Panel: A pandas.Panel with 3 dataframes:
            MODEL: model phase advances between each element
            MEAS: measured phase advances
            ERRMEAS: the standard deviation of the measured values.

    """
    LOGGER.info(
        "Compensating excitation for plane {2:s}. Q = {0:f}, Qd = {1:f}".format(Q, Qd, plane))

    acc = getllm_d.accelerator
    psid_ac2bpmac = ac2bpmac[1]
    k_bpmac = ac2bpmac[2]
    plane_mu = "MUX" if plane == "H" else "MUY"
    number_commonbpms = bpm.shape[0]
    bd = getllm_d.accelerator.get_beam_direction()

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if getllm_d.lhc_phase == "1":
        LOGGER.info("correcting phase jump")
        k_lastbpm = acc.get_k_first_BPM(bpm.index)
    else:
        k_lastbpm = len(bpm)
        LOGGER.info("phase jump will not be corrected")

    # pandas panel that stores the model phase advances, measurement phase advances and measurement errors
    phase_advances = pd.Panel(
        items=["MODEL", "MEAS", "ERRMEAS"],
        major_axis=bpm.index, minor_axis=bpm.index)

    phases_mdl = np.array(model.loc[bpm.index, plane_mu])
    phase_advances["MODEL"] = ((phases_mdl[np.newaxis, :] - phases_mdl[:, np.newaxis])) % 1.0

    #-- Global parameters of the driven motion
    r = sin(PI * (Qd - Q)) / sin(PI * (Qd + Q))
    LOGGER.debug(plane + " compensation lambda = {}".format(r))
    LOGGER.debug(plane + " k_bpmac = {}".format(k_bpmac))
    LOGGER.debug(plane + " psid_ac2bpmac = {}".format(psid_ac2bpmac))
    LOGGER.debug(plane + " bpmac = {}".format(ac2bpmac[0]))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    # temporary 3D matrix that stores the phase advances
    phase_matr_meas = np.empty((len(Files), number_commonbpms, number_commonbpms))
    sin_phase_matr_meas = np.zeros((number_commonbpms, number_commonbpms))
    cos_phase_matr_meas = np.zeros((number_commonbpms, number_commonbpms))
    for i, file_tfs in enumerate(Files):
        phases_meas = bd * np.array(file_tfs.loc[bpm.index, plane_mu]) #-- bd flips B2 phase to B1 direction
        phases_meas[k_lastbpm+1:] += Qd# * bd
        psid = phases_meas - (phases_meas[k_bpmac] - psid_ac2bpmac) # OK, untill here, it is Psi(s, s_ac)
        psid += .5 * Qd# * bd
        psid[k_bpmac:] -= Qd

        Psi = (np.arctan((1 - r) / (1 + r) * np.tan(TWOPI * psid)) / TWOPI) % 0.5
        Psi = np.where(psid % 1.0 > 0.5, Psi + .5, Psi)
        Psi -= Psi[0]
        Psi[k_bpmac:] += Q
        LOGGER.debug(psid[np.isnan(Psi)])

        meas_matr = (Psi[np.newaxis,:] - Psi[:,np.newaxis])
        phase_matr_meas[i] = meas_matr
        sin_phase_matr_meas += np.sin(meas_matr * TWOPI)
        cos_phase_matr_meas += np.cos(meas_matr * TWOPI)

    phase_advances["NFILES"] = len(Files)

    phase_advances["MEAS"] = (np.arctan2(sin_phase_matr_meas / len(Files), cos_phase_matr_meas /
                                         len(Files)) / TWOPI) % 1

    if phase.OPTIMISTIC:
        phase_advances["ERRMEAS"] = np.std(phase_matr_meas, axis=0) * t_value_correction(len(Files)) / np.sqrt(len(Files))
    else:
        R = np.sqrt(
            (sin_phase_matr_meas * sin_phase_matr_meas + cos_phase_matr_meas * cos_phase_matr_meas)
            ) / len(Files)
        phase_advances["ERRMEAS"] = np.sqrt(-2.0 * np.log(R)) / np.sqrt(len(Files))

    return phase_advances, 0

def get_free_beta_from_amp_eq(MADTwiss_ac, Files, Qd, Q, ac2bpmac, plane, getllm_d, commonbpms):
    # TODO: check if this function excludes additional BPMs, right now commonbpms is return unchanged but maybe this has
    # to be different

    #-- Select common BPMs
    bd = getllm_d.accelerator.get_beam_direction()
    

    good_bpms_for_kick = intersect_bpm_list_with_arc_bpms(
        intersect_bpms_list_with_bad_known_bpms(commonbpms)
    )

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if getllm_d.lhc_phase == "1":
        LOGGER.info("correcting phase jump")
        s_lastbpm = getllm_d.accelerator.get_s_first_BPM()

    #-- Determine the BPM closest to the AC dipole and its position
    k_bpmac = ac2bpmac[2]
    bpmac = ac2bpmac[0]
    psid_ac2bpmac = ac2bpmac[1]

    #-- Model beta and phase advance
    if plane == 'H':
        betmdl = MADTwiss_ac.loc[commonbpms.index, "BETX"]
    if plane == 'V':
        betmdl = MADTwiss_ac.loc[commonbpms.index, "BETY"]

    #-- Global parameters of the driven motion
    r = sin(np.pi * (Qd - Q)) / sin(np.pi * (Qd + Q))

    # TODO: Use std to compute errorbars.
    sqrt2j, sqrt2j_std = get_kick_from_bpm_list_w_ACdipole(
        MADTwiss_ac, good_bpms_for_kick, Files, plane
    )

    #-- Loop for files
    betall = np.zeros((len(commonbpms.index), len(Files)))
    for i in range(len(Files)):
        if plane == 'H':
            amp = 2 * Files[i].loc[commonbpms.index, "AMPX"].values
            psid = bd * 2 * np.pi * Files[i].loc[commonbpms.index, "MUX"]
        if plane == 'V':
            amp = 2 * Files[i].loc[commonbpms.index, "AMPY"].values
            psid = bd * 2 * np.pi * Files[i].loc[commonbpms.index, "MUY"]

        # This loop is just to fix the phase jump at the beginning of the ring.
        for k in range(len(commonbpms.index)):
            try:
                if all_bpms[k][0] > s_lastbpm:
                    psid[k] += 2 * np.pi * Qd
            except:
                pass

        psid = psid - (psid[k_bpmac] - psid_ac2bpmac)
        Psid = psid + np.pi * Qd
        Psid[k_bpmac:] = Psid[k_bpmac:] - 2 * np.pi * Qd
        bet = ((amp / sqrt2j[i]) ** 2 *
               (1 + r ** 2 + 2 * r * np.cos(2 * Psid)) / (1 - r ** 2))
        for bpm_index in range(len(commonbpms.index)):
            betall[bpm_index][i] = bet[bpm_index]

    #-- Output
    #result = {}
    #bb = []
    #for k in range(len(commonbpms.index)):
    #   betave = np.mean(betall[k])
    #   betstd = np.std(betall[k])
    #   bb.append((betave - betmdl[k]) / betmdl[k])
    #   result[all_bpms[k][1]] = [betave, betstd, all_bpms[k][0]]
    #bb = math.sqrt(np.mean(np.array(bb) ** 2))

    betave = np.mean(betall, 1)
    betstd = np.std(betall, 1)
    bb = math.sqrt(np.mean((betave / betmdl - 1.0) ** 2))
    result = dict(zip(commonbpms.index, zip(betave, betstd, commonbpms.loc[:, "S"].values)))
    return result, bb, commonbpms


def intersect_bpm_list_with_arc_bpms(bpms_list):
    """Intersects bpms_list with LHC arc BPMs.

    TODO: move this to accelerator class!!!!!
    """

    BPM_distance = bpms_list.index.str.extract("BPM[_,A-Z]*\\.([0-9]+)[R,L].*").values.astype(int)
    return bpms_list.loc[BPM_distance > 13]


BAD_BPM_LIST = ['BPM.15R8.B1', 'BPM.16R3.B1', 'BPM.31L5.B1', 'BPM.23L6.B1',
                'BPM.22R8.B1', 'BPM.11R6.B1', 'BPM.18R6.B1', 'BPM.34R5.B1']


def intersect_bpms_list_with_bad_known_bpms(bpms_list):
    return bpms_list.loc[bpms_list.index.difference(BAD_BPM_LIST)]


def get_kick_from_bpm_list_w_ACdipole(MADTwiss_ac, bpm_list, measurements, plane):
    '''
    @author: F Carlier
    Function calculates kick from measurements with AC dipole using the amplitude of the main line. The main line
    amplitude is obtained from Drive/SUSSIX and is normalized with the model beta-function.

    Input:
        bpm_list:     Can be any list of bpms. Preferably only arc bpms for kick calculations, but other bad bpms may be
                      included as well.
        measurements: List of measurements when analyzing multiple measurements at once
        plane:        Either H or V
    Output:
        actions_sqrt:       is a list containing the actions of each measurement. Notice this is the square root of the
                            action, so sqrt(2JX) or sqrt(2JY) depending on the plane
        actions_sqrt_err:   is the list containing the errors for sqrt(2Jx/y) for each measurement.
    '''
    if plane == 'H':
        betmdl = MADTwiss_ac.loc[bpm_list.index, "BETX"].values
    if plane == 'V':
        betmdl = MADTwiss_ac.loc[bpm_list.index, "BETY"].values

    actions_sqrt = []
    actions_sqrt_err = []

    for meas in measurements:
        if plane == 'H':
            amp = 2 * meas.loc[bpm_list.index, "AMPX"].values
        if plane == 'V':
            amp = 2 * meas.loc[bpm_list.index, "AMPY"].values
        actions_sqrt.append(np.average(amp / np.sqrt(betmdl)))
        actions_sqrt_err.append(np.std(amp / np.sqrt(betmdl)))

    return actions_sqrt, actions_sqrt_err


def GetFreeCoupling_Eq(MADTwiss, FilesX, FilesY, bpms, Qh, Qv, Qx, Qy, accelerator):
    """Calculates coupling using Ryoichi's formula for AC dipole compensation.
       Details of this algorithms is in http://www.agsrhichome.bnl.gov/AP/ap_notes/ap_note_410.pdf

    Args:
        MADTwiss: model twiss file
        FilesX: horizontal measurement files.
        FilesY: vertical measurement files.
        bpms: list of commonbpms
        Qh, Qv: natural tunes
        Qx, Qy: driven tunes
        psih_ac2bpmac, psiv_ac2bpmac: result of GetACPhase_AC2BPMAC()
        accelerator: accelerator class instance.

    """

    #-- Check linx/liny files, may be redundant
    if len(FilesX)!=len(FilesY): return [{},[]]

    ac2bpmac_h = GetACPhase_AC2BPMAC(bpms, Qx, Qh, "H", accelerator)
    ac2bpmac_v = GetACPhase_AC2BPMAC(bpms, Qy, Qv, "V", accelerator)

    horBPMsCopensation =[]
    verBPMsCopensation = []
    psid_ac2bpmac_h = ac2bpmac_h[1]
    k_bpmac_h = ac2bpmac_h[2]
    psid_ac2bpmac_v = ac2bpmac_v[1]
    k_bpmac_v = ac2bpmac_v[2]

    bd = accelerator.get_beam_direction()
    fqwList = []
    if True:

       #-- Global parameters of the driven motion
        dh =Qh-Qx
        dv =Qv-Qy
        rh =sin(np.pi*(Qh-Qx))/sin(np.pi*(Qh+Qx))
        rv =sin(np.pi*(Qv-Qy))/sin(np.pi*(Qv+Qy))
        rch=sin(np.pi*(Qh-Qy))/sin(np.pi*(Qh+Qy))
        rcv=sin(np.pi*(Qx-Qv))/sin(np.pi*(Qx+Qv))

        #-- Loop for files
        f1001Abs =np.zeros((len(bpms),len(FilesX)))
        f1010Abs =np.zeros((len(bpms),len(FilesX)))
        f1001xArg=np.zeros((len(bpms),len(FilesX)))
        f1001yArg=np.zeros((len(bpms),len(FilesX)))
        f1010xArg=np.zeros((len(bpms),len(FilesX)))
        f1010yArg=np.zeros((len(bpms),len(FilesX)))
        for i in range(len(FilesX)):

            #-- Read amplitudes and phases
            amph  =     np.array([FilesX[i].loc[b, "AMPX"]    for b in bpms.index])
            ampv  =     np.array([FilesY[i].loc[b, "AMPY"]    for b in bpms.index])
            amph01=     np.array([FilesX[i].loc[b, "AMP01"]   for b in bpms.index])
            ampv10=     np.array([FilesY[i].loc[b, "AMP10"]   for b in bpms.index])
            psih  =2*np.pi*np.array([FilesX[i].loc[b, "MUX"]     for b in bpms.index])
            psiv  =2*np.pi*np.array([FilesY[i].loc[b, "MUY"]     for b in bpms.index])
            psih01=2*np.pi*np.array([FilesX[i].loc[b, "PHASE01"] for b in bpms.index])
            psiv10=2*np.pi*np.array([FilesY[i].loc[b, "PHASE10"] for b in bpms.index])
            #-- I'm not sure this is correct for the coupling so I comment out this part for now (by RM 9/30/11).
            #for k in range(len(bpm)):
            #       try:
            #               if bpm[k][0]>s_lastbpm:
            #                       psih[k]  +=bd*2*np.pi*Qh  #-- To fix the phase shift by Qh
            #                       psiv[k]  +=bd*2*np.pi*Qv  #-- To fix the phase shift by Qv
            #                       psih01[k]+=bd*2*np.pi*Qv  #-- To fix the phase shift by Qv
            #                       psiv10[k]+=bd*2*np.pi*Qh  #-- To fix the phase shift by Qh
            #       except: pass

            #-- Construct Fourier components
            #   * be careful for that the note is based on x+i(alf*x*bet*x')).
            #   * Calculating Eqs (87)-(92) by using Eqs (47) & (48) (but in the Fourier space) in the note.
            #   * Note that amph(v)01 is normalized by amph(v) and it is un-normalized in the following.
            dpsih  =np.append(psih[1:]  ,2*np.pi*Qh+psih[0]  )-psih
            dpsiv  =np.append(psiv[1:]  ,2*np.pi*Qv+psiv[0]  )-psiv
            dpsih01=np.append(psih01[1:],2*np.pi*Qv+psih01[0])-psih01
            dpsiv10=np.append(psiv10[1:],2*np.pi*Qh+psiv10[0])-psiv10

            X_m10=2*amph*np.exp(-1j*psih)
            Y_0m1=2*ampv*np.exp(-1j*psiv)
            X_0m1=amph*np.exp(-1j*psih01)/(1j*sin(dpsih))*(amph01*np.exp(1j*dpsih)-np.append(amph01[1:],amph01[0])*np.exp(-1j*dpsih01))
            X_0p1=amph*np.exp( 1j*psih01)/(1j*sin(dpsih))*(amph01*np.exp(1j*dpsih)-np.append(amph01[1:],amph01[0])*np.exp( 1j*dpsih01))
            Y_m10=ampv*np.exp(-1j*psiv10)/(1j*sin(dpsiv))*(ampv10*np.exp(1j*dpsiv)-np.append(ampv10[1:],ampv10[0])*np.exp(-1j*dpsiv10))
            Y_p10=ampv*np.exp( 1j*psiv10)/(1j*sin(dpsiv))*(ampv10*np.exp(1j*dpsiv)-np.append(ampv10[1:],ampv10[0])*np.exp( 1j*dpsiv10))

            #-- Construct f1001hv, f1001vh, f1010hv (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
            f1001hv=-np.conjugate(1/(2j)*Y_m10/X_m10)  #-- - sign from the different def
            f1001vh=-1/(2j)*X_0m1/Y_0m1             #-- - sign from the different def
            f1010hv=-1/(2j)*Y_p10/np.conjugate(X_m10)  #-- - sign from the different def
            f1010vh=-1/(2j)*X_0p1/np.conjugate(Y_0m1)  #-- - sign from the different def
    ##              f1001hv=conjugate(1/(2j)*Y_m10/X_m10)
    ##              f1001vh=1/(2j)*X_0m1/Y_0m1
    ##              f1010hv=1/(2j)*Y_p10/conjugate(X_m10)
    ##              f1010vh=1/(2j)*X_0p1/conjugate(Y_0m1)

            #-- Construct phases psih, psiv, Psih, Psiv w.r.t. the AC dipole
            psih = psih - (psih[k_bpmac_h] - psid_ac2bpmac_h) # OK, untill here, it is Psi(s, s_ac)
            psiv = psiv - (psiv[k_bpmac_v] - psid_ac2bpmac_v) # OK, untill here, it is Psi(s, s_ac)

            Psih=psih-np.pi*Qh
            Psih[:k_bpmac_h]=Psih[:k_bpmac_h]+2*np.pi*Qh
            Psiv=psiv-np.pi*Qv
            Psiv[:k_bpmac_v]=Psiv[:k_bpmac_v]+2*np.pi*Qv

            Psix=np.arctan((1-rh)/(1+rh)*np.tan(Psih))%np.pi
            Psiy=np.arctan((1-rv)/(1+rv)*np.tan(Psiv))%np.pi
            for k in range(len(bpms)):
                if Psih[k]%(2*np.pi)>np.pi: Psix[k]=Psix[k]+np.pi
                if Psiv[k]%(2*np.pi)>np.pi: Psiy[k]=Psiy[k]+np.pi

            psix=Psix-np.pi*Qx
            psix[k_bpmac_h:]=psix[k_bpmac_h:]+2*np.pi*Qx
            psiy=Psiy-np.pi*Qy
            psiy[k_bpmac_v:]=psiy[k_bpmac_v:]+2*np.pi*Qy

            #-- Construct f1001h, f1001v, f1010h, f1010v (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
            f1001h=1/math.sqrt(1-rv**2)*(np.exp(-1j*(Psiv-Psiy))*f1001hv+rv*np.exp( 1j*(Psiv+Psiy))*f1010hv)
            f1010h=1/math.sqrt(1-rv**2)*(np.exp( 1j*(Psiv-Psiy))*f1010hv+rv*np.exp(-1j*(Psiv+Psiy))*f1001hv)
            f1001v=1/math.sqrt(1-rh**2)*(np.exp( 1j*(Psih-Psix))*f1001vh+rh*np.exp(-1j*(Psih+Psix))*np.conjugate(f1010vh))
            f1010v=1/math.sqrt(1-rh**2)*(np.exp( 1j*(Psih-Psix))*f1010vh+rh*np.exp(-1j*(Psih+Psix))*np.conjugate(f1001vh))

            #-- Construct f1001 and f1010 from h and v BPMs (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
            g1001h          =np.exp(-1j*((psih-psih[k_bpmac_h])-(psiy-psiy[k_bpmac_v])))*(ampv/amph*amph[k_bpmac_h]/ampv[k_bpmac_v])*f1001h[k_bpmac_h]
            g1001h[:k_bpmac_h]=1/(np.exp(2*np.pi*1j*(Qh-Qy))-1)*(f1001h-g1001h)[:k_bpmac_h]
            g1001h[k_bpmac_h:]=1/(1-np.exp(-2*np.pi*1j*(Qh-Qy)))*(f1001h-g1001h)[k_bpmac_h:]

            g1010h          =np.exp(-1j*((psih-psih[k_bpmac_h])+(psiy-psiy[k_bpmac_v])))*(ampv/amph*amph[k_bpmac_h]/ampv[k_bpmac_v])*f1010h[k_bpmac_h]
            g1010h[:k_bpmac_h]=1/(np.exp(2*np.pi*1j*(Qh+Qy))-1)*(f1010h-g1010h)[:k_bpmac_h]
            g1010h[k_bpmac_h:]=1/(1-np.exp(-2*np.pi*1j*(Qh+Qy)))*(f1010h-g1010h)[k_bpmac_h:]

            g1001v          =np.exp(-1j*((psix-psix[k_bpmac_h])-(psiv-psiv[k_bpmac_v])))*(amph/ampv*ampv[k_bpmac_v]/amph[k_bpmac_h])*f1001v[k_bpmac_v]
            g1001v[:k_bpmac_v]=1/(np.exp(2*np.pi*1j*(Qx-Qv))-1)*(f1001v-g1001v)[:k_bpmac_v]
            g1001v[k_bpmac_v:]=1/(1-np.exp(-2*np.pi*1j*(Qx-Qv)))*(f1001v-g1001v)[k_bpmac_v:]

            g1010v          =np.exp(-1j*((psix-psix[k_bpmac_h])+(psiv-psiv[k_bpmac_v])))*(amph/ampv*ampv[k_bpmac_v]/amph[k_bpmac_h])*f1010v[k_bpmac_v]
            g1010v[:k_bpmac_v]=1/(np.exp(2*np.pi*1j*(Qx+Qv))-1)*(f1010v-g1010v)[:k_bpmac_v]
            g1010v[k_bpmac_v:]=1/(1-np.exp(-2*np.pi*1j*(Qx+Qv)))*(f1010v-g1010v)[k_bpmac_v:]
    
            f1001x=np.exp(1j*(psih-psix))*f1001h
            f1001x=f1001x-rh*np.exp(-1j*(psih+psix))/rch*np.conjugate(f1010h)
            f1001x=f1001x-2j*sin(np.pi*dh)*np.exp(1j*(Psih-Psix))*g1001h
            f1001x=f1001x-2j*sin(np.pi*dh)*np.exp(-1j*(Psih+Psix))/rch*np.conjugate(g1010h)
            f1001x=1/math.sqrt(1-rh**2)*sin(np.pi*(Qh-Qy))/sin(np.pi*(Qx-Qy))*f1001x
    
            f1010x=np.exp(1j*(psih-psix))*f1010h
            f1010x=f1010x-rh*np.exp(-1j*(psih+psix))*rch*np.conjugate(f1001h)
            f1010x=f1010x-2j*sin(np.pi*dh)*np.exp(1j*(Psih-Psix))*g1010h
            f1010x=f1010x-2j*sin(np.pi*dh)*np.exp(-1j*(Psih+Psix))*rch*np.conjugate(g1001h)
            f1010x=1/math.sqrt(1-rh**2)*sin(np.pi*(Qh+Qy))/sin(np.pi*(Qx+Qy))*f1010x
    
            f1001y=np.exp(-1j*(psiv-psiy))*f1001v
            f1001y=f1001y+rv*np.exp(1j*(psiv+psiy))/rcv*f1010v
            f1001y=f1001y+2j*sin(np.pi*dv)*np.exp(-1j*(Psiv-Psiy))*g1001v
            f1001y=f1001y-2j*sin(np.pi*dv)*np.exp(1j*(Psiv+Psiy))/rcv*g1010v
            f1001y=1/math.sqrt(1-rv**2)*sin(np.pi*(Qx-Qv))/sin(np.pi*(Qx-Qy))*f1001y
    
            f1010y=np.exp(1j*(psiv-psiy))*f1010v
            f1010y=f1010y+rv*np.exp(-1j*(psiv+psiy))*rcv*f1001v
            f1010y=f1010y-2j*sin(np.pi*dv)*np.exp(1j*(Psiv-Psiy))*g1010v
            f1010y=f1010y+2j*sin(np.pi*dv)*np.exp(-1j*(Psiv+Psiy))*rcv*g1001v
            f1010y=1/math.sqrt(1-rv**2)*sin(np.pi*(Qx+Qv))/sin(np.pi*(Qx+Qy))*f1010y
    
            #-- For B2, must be double checked
            if bd == -1:
                f1001x=-np.conjugate(f1001x)
                f1001y=-np.conjugate(f1001y)
                f1010x=-np.conjugate(f1010x)
                f1010y=-np.conjugate(f1010y)
    
            #-- Separate to amplitudes and phases, amplitudes averaged to cancel math.sqrt(betv/beth) and math.sqrt(beth/betv)
            for k in range(len(bpms)):
                f1001Abs[k][i] =math.sqrt(abs(f1001x[k]*f1001y[k]))
                f1010Abs[k][i] =math.sqrt(abs(f1010x[k]*f1010y[k]))
                f1001xArg[k][i]=np.angle(f1001x[k])%(2*np.pi)
                f1001yArg[k][i]=np.angle(f1001y[k])%(2*np.pi)
                f1010xArg[k][i]=np.angle(f1010x[k])%(2*np.pi)
                f1010yArg[k][i]=np.angle(f1010y[k])%(2*np.pi)
    
        #-- Output
        fwqw={}
        goodbpm=[]
        for k in range(len(bpms)):
            bname = bpms.index[k]
            #-- Bad BPM flag based on phase
            badbpm=0
            f1001xArgAve = phase.calc_phase_mean(f1001xArg[k],2*np.pi)
            f1001yArgAve = phase.calc_phase_mean(f1001yArg[k],2*np.pi)
            f1010xArgAve = phase.calc_phase_mean(f1010xArg[k],2*np.pi)
            f1010yArgAve = phase.calc_phase_mean(f1010yArg[k],2*np.pi)
            #This seems to be to conservative or somethings...
            if min(abs(f1001xArgAve-f1001yArgAve),2*np.pi-abs(f1001xArgAve-f1001yArgAve))>np.pi/2: badbpm=1
            if min(abs(f1010xArgAve-f1010yArgAve),2*np.pi-abs(f1010xArgAve-f1010yArgAve))>np.pi/2: badbpm=1
            #-- Output
	    badbpm=0
            if badbpm==0:
                f1001AbsAve = np.mean(f1001Abs[k])
                f1010AbsAve = np.mean(f1010Abs[k])
                f1001ArgAve = phase.calc_phase_mean(np.append(f1001xArg[k],f1001yArg[k]),2*np.pi)
                f1010ArgAve = phase.calc_phase_mean(np.append(f1010xArg[k],f1010yArg[k]),2*np.pi)
                f1001Ave = f1001AbsAve*np.exp(1j*f1001ArgAve)
                f1010Ave = f1010AbsAve*np.exp(1j*f1010ArgAve)
                f1001AbsStd = math.sqrt(np.mean((f1001Abs[k]-f1001AbsAve)**2))
                f1010AbsStd = math.sqrt(np.mean((f1010Abs[k]-f1010AbsAve)**2))
                f1001ArgStd = phase.calc_phase_std(np.append(f1001xArg[k],f1001yArg[k]),2*np.pi)
                f1010ArgStd = phase.calc_phase_std(np.append(f1010xArg[k],f1010yArg[k]),2*np.pi)
                fwqw[bname] = [[f1001Ave          ,f1001AbsStd       ,f1010Ave          ,f1010AbsStd       ],
                                 [f1001ArgAve/(2*np.pi),f1001ArgStd/(2*np.pi),f1010ArgAve/(2*np.pi),f1010ArgStd/(2*np.pi)]]  #-- Phases renormalized to [0,1)
                goodbpm.append(bname)
                
    #-- Global parameters not implemented yet
        
        fqwList.append(fwqw)
        
    
    fwqw = copy.deepcopy(fqwList[0])
    for key in fwqw:
        for a in range(1, len(fqwList)):
            tmp = fqwList[a]
            fwqw[key][0][0] = fwqw[key][0][0] + tmp[key][0][0]
            fwqw[key][0][1] = fwqw[key][0][1] + tmp[key][0][1]
            fwqw[key][0][2] = fwqw[key][0][2] + tmp[key][0][2]
            fwqw[key][0][3] = fwqw[key][0][3] + tmp[key][0][3]
            if(key is 'BPMWB.4R5.B1'):
                print fwqw[key][1]
        fwqw[key][0][0]=fwqw[key][0][0]/len(fqwList)
        fwqw[key][0][1]=fwqw[key][0][1]/len(fqwList)
        fwqw[key][0][2]=fwqw[key][0][2]/len(fqwList)
        fwqw[key][0][3]=fwqw[key][0][3]/len(fqwList)
    fwqw['Global']=['"null"','"null"']        
    return [fwqw,bpms.loc[goodbpm]]

def GetFreeIP2_Eq(MADTwiss, Files, Qd, Q, ac2bpmac, plane, accelerator, bpms, op):

    if ac2bpmac is None:
        return [{}, []]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    s_lastbpm = accelerator.get_s_first_BPM()
    bd = accelerator.get_beam_direction()

    #-- Determine the BPM closest to the AC dipole and its position
    bpmac = ac2bpmac[0]
    psid_bpmac = ac2bpmac[1]
    k_bpmac = ac2bpmac[2]

    #-- Global parameters of the driven motion
    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))

    #-- Determine Psid (w.r.t the AC dipole) for each file
    Psidall=[]
    for i in range(len(Files)):
        if plane=='H': psid = bd * 2 * np.pi * Files[i].loc[:, "MUX"]  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid = bd * 2 * np.pi * Files[i].loc[:, "MUY"]   #-- bd flips B2 phase to B1 direction
        for k in range(len(bpms)):
            try:
                if bpm.iloc[k]["S"] > s_lastbpm: psid[k] += 2 * np.pi * Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac] - psid_bpmac)
        Psid=psid+np.pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*np.pi*Qd
        Psidall.append(Psid)

    #-- Loop for IPs
    result={}
    for ip in ('1','2','5','8'):

        bpml = 'BPMSW.1L' + ip + '.' + accelerator.__class__.__name__[-2:]
        bpmr = 'BPMSW.1R' + ip + '.' + accelerator.__class__.__name__[-2:]
        if bpml in bpms.index and bpmr in bpms.index:

            #-- Model values
            L = 0.5 * (MADTwiss.loc[bpmr, "S"] - MADTwiss.loc[bpml, "S"])
            if L < 0: L += 0.5 * MADTwiss.LENGTH
            if plane=='H':
                betlmdl = MADTwiss.loc[bpml, "BETX"]
                alflmdl = MADTwiss.loc[bpml, "ALFX"]
            if plane=='V':
                betlmdl = MADTwiss.loc[bpml, "BETY"]
                alflmdl = MADTwiss.loc[bpml, "ALFY"]
            betsmdl=betlmdl/(1+alflmdl**2)
            betmdl =betlmdl-2*alflmdl*L+L**2/betsmdl
            alfmdl =alflmdl-L/betsmdl
            dsmdl  =alfmdl*betsmdl

            #-- Measurement for each file
            betall=[]
            alfall=[]
            betsall=[]
            dsall=[]
            rt2Jall=[]
            for i in range(len(Files)):
                try:    #-- Maybe not needed, to avoid like math.sqrt(-...)
                    if plane=='H':
                        al = Files[i].loc[bpml, "AMPX"]
                        ar = Files[i].loc[bpmr, "AMPX"]
                    if plane=='V':
                        al = Files[i].loc[bpml, "AMPY"]
                        ar = Files[i].loc[bpmr, "AMPY"]

                    # OMG, what is this?
                    #Psidl=Psidall[i][list(zip(*bpm)[1]).index(bpml)]
                    #Psidr=Psidall[i][list(zip(*bpm)[1]).index(bpmr)]
                    # try with:
                    Psidl = Psidall[i][bpms.index.get_loc(bpml)]
                    Psidr = Psidall[i][bpms.index.get_loc(bpmr)]

                    dpsid=Psidr-Psidl

                    #-- betd, alfd, and math.sqrt(2Jd) at BPM_left from amp and phase advance
                    betdl=2*L*al/(ar*sin(dpsid))
                    alfdl=(al-ar*cos(dpsid))/(ar*sin(dpsid))
                    rt2J =math.sqrt(al*ar*sin(dpsid)/(2*L))
                    #-- Convert to free bet and alf
                    betl=(1+r**2+2*r*np.cos(2*Psidl))/(1-r**2)*betdl
                    alfl=((1+r**2+2*r*np.cos(2*Psidl))*alfdl+2*r*sin(2*Psidl))/(1-r**2)
                    #-- Calculate IP parameters
                    bets=betl/(1+alfl**2)
                    bet =betl-2*alfl*L+L**2/bets
                    alf =alfl-L/bets
                    ds  =alf*bets
                    betall.append(bet)
                    alfall.append(alf)
                    betsall.append(bets)
                    dsall.append(ds)
                    rt2Jall.append(rt2J)
                except:
                    pass

            #-- Ave and Std
            betall =np.array(betall) ; betave =np.mean(betall) ; betstd =math.sqrt(np.mean((betall-betave)**2))
            alfall =np.array(alfall) ; alfave =np.mean(alfall) ; alfstd =math.sqrt(np.mean((alfall-alfave)**2))
            betsall=np.array(betsall); betsave=np.mean(betsall); betsstd=math.sqrt(np.mean((betsall-betsave)**2))
            dsall  =np.array(dsall)  ; dsave  =np.mean(dsall)  ; dsstd  =math.sqrt(np.mean((dsall-dsave)**2))
            rt2Jall=np.array(rt2Jall); rt2Jave=np.mean(rt2Jall); rt2Jstd=math.sqrt(np.mean((rt2Jall-rt2Jave)**2))
            result['IP'+ip]=[betave,betstd,betmdl,alfave,alfstd,alfmdl,betsave,betsstd,betsmdl,dsave,dsstd,dsmdl,rt2Jave,rt2Jstd]

    return result

def getkickac(MADTwiss_ac,files,psih_ac2bpmac,psiv_ac2bpmac,bd,op):

    invarianceJx = []
    invarianceJy = []

    tunex = []
    tuney = []
    tunexRMS = []
    tuneyRMS = []

    nat_tunex = []
    nat_tuney = []
    nat_tunexRMS = []
    nat_tuneyRMS = []

    dpp = []

    all_bpms_x = utils.bpm.model_intersect(utils.bpm.intersect(files[0]), MADTwiss_ac )
    all_bpms_y = utils.bpm.model_intersect(utils.bpm.intersect(files[1]), MADTwiss_ac )

    good_bpms_for_kick_x = intersect_bpm_list_with_arc_bpms( intersect_bpms_list_with_bad_known_bpms(all_bpms_x) )
    good_bpms_for_kick_y = intersect_bpm_list_with_arc_bpms( intersect_bpms_list_with_bad_known_bpms(all_bpms_y) )
    Jx2sq, Jx2sq_std = get_kick_from_bpm_list_w_ACdipole(MADTwiss_ac, good_bpms_for_kick_x, files[0], 'H')
    Jy2sq, Jy2sq_std = get_kick_from_bpm_list_w_ACdipole(MADTwiss_ac, good_bpms_for_kick_y, files[1], 'V')

    for j in range(len(files[0])):

        tw_x = files[0][j]
        tw_y = files[1][j]
        invarianceJx.append([Jx2sq[j], Jx2sq_std[j]])
        invarianceJy.append([Jy2sq[j], Jy2sq_std[j]])
        
        dpp.append(getattr(tw_x, "DPP", 0.0))

        tunex.append(getattr(tw_x, "Q1", 0.0))
        tuney.append(getattr(tw_y, "Q2", 0.0))
        tunexRMS.append(getattr(tw_x, "Q1RMS", 0.0))
        tuneyRMS.append(getattr(tw_y, "Q2RMS", 0.0))

        nat_tunex.append(getattr(tw_x, "NATQ1", 0.0))
        nat_tuney.append(getattr(tw_y, "NATQ2", 0.0))
        nat_tunexRMS.append(getattr(tw_x, "NATQ1RMS", 0.0))
        nat_tuneyRMS.append(getattr(tw_y, "NATQ2RMS", 0.0))

    tune_values_list = [tunex, tunexRMS, tuney, tuneyRMS, nat_tunex, nat_tunexRMS, nat_tuney, nat_tuneyRMS]
    return [invarianceJx, invarianceJy, tune_values_list, dpp]

######### end ac-dipole stuff

