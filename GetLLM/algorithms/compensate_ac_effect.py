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

import numpy as np
from numpy import sin, cos, tan

import Utilities.bpm
import phase
from SegmentBySegment.SegmentBySegment import get_good_bpms
from __builtin__ import raw_input
from constants import PI, TWOPI, kEPSILON

DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# constants
#===================================================================================================




#===================================================================================================
# helper-functions
#===================================================================================================
#---------  The following is functions
def GetACPhase_AC2BPMAC(MADTwiss,Qd,Q,plane,oa, acdipole):

    if acdipole == "ACD":
        dipole_nameH = 'MKQA.6L4.'+oa[3:]
        dipole_nameV= 'MKQA.6L4.' + oa[3:]
        if   oa=='LHCB1':
            bpmac1='BPMYA.5L4.B1'
            bpmac2='BPMYB.6L4.B1'
        elif oa=='LHCB2':
            bpmac1='BPMYB.5L4.B2'
            bpmac2='BPMYA.6L4.B2'
        else:
            return {}
    elif acdipole == "ADT":
        if oa == 'LHCB1':
            dipole_nameH = "ADTKH.C5L4.B1"
            dipole_nameV = "ADTKV.B5R4.B1"
            if plane=="H":
                bpmac1='BPMWA.B5L4.B1'
                bpmac2='BPMWA.A5L4.B1'
            else:
                bpmac1 = 'BPMWA.B5R4.B1'
                bpmac2 = 'BPMWA.A5R4.B1'
                
        elif oa == 'LHCB2':
            dipole_nameH = "ADTKV.C5L4.B2"
            dipole_nameV = "ADTKH.B5R4.B2"
            if plane == "H":
                bpmac1 = 'BPMWA.B5R4.B2'
                bpmac2 = 'BPMWA.A5R4.B2'
            else:
                bpmac1 = 'BPMWA.B5L4.B2'
                bpmac2 = 'BPMWA.A5L4.B2'
        else:
            return {}

    if plane=='H':
        psi_ac2bpmac1=MADTwiss.MUX[MADTwiss.indx[bpmac1]]-MADTwiss.MUX[MADTwiss.indx[dipole_nameH]]  #-- B1 direction for B2
        psi_ac2bpmac2=MADTwiss.MUX[MADTwiss.indx[bpmac2]]-MADTwiss.MUX[MADTwiss.indx[dipole_nameH]]  #-- B1 direction for B2
    if plane=='V':
        psi_ac2bpmac1=MADTwiss.MUY[MADTwiss.indx[bpmac1]]-MADTwiss.MUY[MADTwiss.indx[dipole_nameV]]  #-- B1 direction for B2
        psi_ac2bpmac2=MADTwiss.MUY[MADTwiss.indx[bpmac2]]-MADTwiss.MUY[MADTwiss.indx[dipole_nameV]]  #-- B1 direction for B2

    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))
    psid_ac2bpmac1=np.arctan((1+r)/(1-r)*tan(2*np.pi*psi_ac2bpmac1-np.pi*Q))%np.pi-np.pi+np.pi*Qd
    psid_ac2bpmac2=np.arctan((1+r)/(1-r)*tan(2*np.pi*psi_ac2bpmac2+np.pi*Q))%np.pi-np.pi*Qd

    return {bpmac1:psid_ac2bpmac1,bpmac2:psid_ac2bpmac2}


def get_free_phase_total_eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,op):

    #-- Select common BPMs
    bpm=Utilities.bpm.model_intersect(Utilities.bpm.intersect(Files),MADTwiss)
    bpm=[(b[0],str.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op == "1":
        if bd== 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
        if bd==-1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]
    else:
        print "LHC phase will not be corrected [total phase]"
    #-- Determine the BPM closest to the AC dipole and its position
    
    # WHY does this code exist?
#     for b in psid_ac2bpmac.keys():
#         if '5L4' in b: bpmac1=b
#         if '6L4' in b: bpmac2=b
    bpmac1 = psid_ac2bpmac.keys()[0]
    bpmac2 = psid_ac2bpmac.keys()[1]
    try:
        k_bpmac=list(zip(*bpm)[1]).index(bpmac1)
        bpmac=bpmac1
    except:
        try:
            k_bpmac=list(zip(*bpm)[1]).index(bpmac2)
            bpmac=bpmac2
        except:
            return [{},[]]

    # -- Model phase advances
    if plane == 'H':
        psimdl = np.array([(MADTwiss.MUX[MADTwiss.indx[b[1]]]-MADTwiss.MUX[MADTwiss.indx[bpm[0][1]]])%1 for b in bpm])
    if plane == 'V':
        psimdl=np.array([(MADTwiss.MUY[MADTwiss.indx[b[1]]]-MADTwiss.MUY[MADTwiss.indx[bpm[0][1]]])%1 for b in bpm])

    # -- Global parameters of the driven motion
    r = sin(np.pi * (Qd - Q)) / sin(np.pi * (Qd + Q))

    # -- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    psiall=np.zeros((len(bpm), len(Files)))
    for i in range(len(Files)):
        if plane == 'H':
            psid = bd * 2 * np.pi * np.array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid=bd*2*np.pi*np.array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*np.pi*Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+np.pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*np.pi*Qd
        Psi=np.arctan((1-r)/(1+r)*np.tan(Psid))%np.pi
        for k in range(len(bpm)):
            if Psid[k]%(2*np.pi)>np.pi: Psi[k]=Psi[k]+np.pi
        psi=Psi-Psi[0]
        psi[k_bpmac:]=psi[k_bpmac:]+2*np.pi*Q
        for k in range(len(bpm)): psiall[k][i]=psi[k]/(2*np.pi)  #-- phase range back to [0,1)

    #-- Output
    result={}
    for k in range(len(bpm)):
        psiave = phase.calc_phase_mean(psiall[k],1)
        psistd = phase.calc_phase_std(psiall[k],1)
        result[bpm[k][1]]=[psiave,psistd,psimdl[k],bpm[0][1]]

    return [result,bpm]


def get_free_phase_eq(MADTwiss, Files, Qd, Q, psid_ac2bpmac, plane, bd, op, Qmdl, acdipole):


# 
#     print "\33[38;2;255;200;0m"
#     print "         /\\                     "
#     print "        //\\\\                   "
#     print "       //  \\\\                  "
#     print "      // || \\\\                 "
#     print "     //  ||  \\\\                "
#     print "    //   ||   \\\\               "
#     print "   //    ||    \\\\              "
#     print "  //            \\\\           "
#     print " //      ()      \\\\           "
#     print "//________________\\\\            "
#     print "--------------------               "
#     print "                               "
#     print "INFO: using changed code for ac compensation\n DO NOT TRUST THIS CODE\33[0m"

    print "Compensating {3:s} effect for plane {2:s}. Q = {0:f}, Qd = {1:f}".format(Q, Qd, plane, acdipole)

    #-- Select common BPMs
    bpm = Utilities.bpm.model_intersect(Utilities.bpm.intersect(Files), MADTwiss)
    bpm = [(b[0], str.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    
    if op == "1":
        print "correcting phase jump"
        if 'MOH_3' in MADTwiss.NAME:
            print "--> for JPARC"
            s_lastbpm = MADTwiss.S[MADTwiss.indx['MOH_3']]
        else:
            print "--> for LHC"
            if bd == 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
            if bd == -1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]
    else:
        print "phase jump will not be corrected"

    #-- Determine the position of the AC dipole BPM
    bpmac1 = psid_ac2bpmac.keys()[0]
    bpmac2 = psid_ac2bpmac.keys()[1]

    try:
        k_bpmac = list(zip(*bpm)[1]).index(bpmac1)
        bpmac = bpmac1
    except:
        try:
            k_bpmac = list(zip(*bpm)[1]).index(bpmac2)
            bpmac = bpmac2
        except:
            print >> sys.stderr,'WARN: BPMs next to AC dipoles missing. AC dipole effects not calculated for '+plane+' with eqs !'
            return [{}, 0.0, []]

    #-- Model phase advances
    if plane=='H': psimdl=np.array([MADTwiss.MUX[MADTwiss.indx[b[1]]] for b in bpm])
    if plane=='V': psimdl=np.array([MADTwiss.MUY[MADTwiss.indx[b[1]]] for b in bpm])
    
    # <<<<<<<<<< ICH
    psiijmdl = [None] * 10
    psiijall = []


    for which_psi in range(1,11):
        psiijmdl[which_psi-1] = (np.append(psimdl[which_psi:], psimdl[:which_psi] + Qmdl) - psimdl) % 1
        psiijall.append(np.zeros((len(bpm), len(Files))))

    #-- Global parameters of the driven motion
    r=sin(PI * (Qd - Q)) / sin(PI * (Qd + Q))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    
    for i in range(len(Files)):
        psid = []
        if plane == 'H':
            psid = bd * TWOPI * np.array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane == 'V':
            psid = bd * TWOPI * np.array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0] > s_lastbpm: psid[k] += TWOPI * Qd  #-- To fix the phase shift by Q
            except: pass
        psid = psid - (psid[k_bpmac] - psid_ac2bpmac[bpmac])  # OK, untill here, it is Psi(s, s_ac)
        Psid=psid + PI * Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-TWOPI * Qd
        
        gamma = Psid*2
        alpha = psid + PI * Qd
        alpha[k_bpmac:] = alpha[k_bpmac:] + TWOPI * (Q - Qd) + kEPSILON
        
        Psi = np.arctan((1 - r) / (1 + r) * np.tan(Psid)) % PI  # Ryoichi
        Psi_a = np.arctan((1 + r * np.sin(alpha - gamma) / np.sin(alpha)) / (1 + r * np.cos(alpha - gamma)/np.cos(alpha)))
        for k in range(len(bpm)):
            if Psid[k] % TWOPI > PI: Psi[k] = Psi[k] + PI
        psi = Psi - Psi[0]
        psi[k_bpmac:] = psi[k_bpmac:] + TWOPI * Q
        
        # <<<<<<<<<< ICH
        psiij = [None] * 10
        for j in range(1, 11):
            psiij[j-1] = (np.append(psi[j:], psi[:j] + TWOPI * Q) - psi)/ TWOPI
        for k in range(len(bpm)):
            # <<<<<<<<<< ICH
            for j in range(0, 10):
                psiijall[j][k][i] = psiij[j][k]
            
    #-- Output
   
    result={}
    result_={}
    muave=0.0  #-- mu is the same as psi but w/o mod
    muave_ = 0.0
    for k in range(len(bpm)):
        # <<<<<<<<<< ICH
        psiijave = [None] * 10
        psiijstd = [None] * 10
        bnj = [None] * 11
        
        for j in range(0,11):
            bnj[j] = str.upper(bpm[(k + j) % len(bpm)][1])
        
        for j in range(0,10):
            psiijave[j] = phase.calc_phase_mean(psiijall[j][k],1)
            psiijstd[j] = phase.calc_phase_std(psiijall[j][k],1)
            result_["".join([plane, bnj[0], bnj[j + 1]])] = [psiijave[j],psiijstd[j],psiijmdl[j][k]]
            
        muave_ += psiijave[0]
        try:    result_[bpm[k][1]]=[psiijave[0],psiijstd[0],psiijave[1],psiijstd[1],psiijmdl[0][k],psiijmdl[1][k],bpm[k+1][1]]
        except: result_[bpm[k][1]]=[psiijave[0],psiijstd[0],psiijave[1],psiijstd[1],psiijmdl[0][k],psiijmdl[1][k],bpm[0][1]]    #-- The last BPM
        

    return [result_,muave_,bpm]

def get_free_phase_eq_intermediat(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,op,Qmdl):

    #-- Select common BPMs
    bpm=Utilities.bpm.model_intersect(Utilities.bpm.intersect(Files),MADTwiss)
    bpm=[(b[0],str.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and bd== 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    if op=="1" and bd==-1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    #-- Determine the position of the AC dipole BPM
    for b in psid_ac2bpmac.keys():
        if '5L4' in b: bpmac1=b
        if '6L4' in b: bpmac2=b
    try:
        k_bpmac=list(zip(*bpm)[1]).index(bpmac1)
        bpmac=bpmac1
    except:
        try:
            k_bpmac=list(zip(*bpm)[1]).index(bpmac2)
            bpmac=bpmac2
        except:
            print >> sys.stderr,'WARN: BPMs next to AC dipoles missing. AC dipole effects not calculated for '+plane+' with eqs !'
            return [{}, 0.0, []]

    #-- Model phase advances
    if plane=='H': psimdl=np.array([MADTwiss.MUX[MADTwiss.indx[b[1]]] for b in bpm])
    if plane=='V': psimdl=np.array([MADTwiss.MUY[MADTwiss.indx[b[1]]] for b in bpm])
    
    # <<<<<<<<<< ICH
    psiijmdl = [None] * 10
    psiijall = []
    # ========== OLD
    # >>>>>>>>>>
    
    
    # <<<<<<<<<< ICH
    for which_psi in range(1,11):
        psiijmdl[which_psi-1] = (np.append(psimdl[which_psi:], psimdl[:which_psi] + Qmdl) - psimdl) % 1
        psiijall.append(np.zeros((len(bpm),len(Files))))
     # ========== OLD   
    psi12mdl= (np.append(psimdl[1:], psimdl[0] + Qmdl)-psimdl)%1
    psi13mdl=(np.append(psimdl[2:],psimdl[:2]+Qmdl)-psimdl)%1
    psi14mdl=(np.append(psimdl[3:],psimdl[:3] +Qmdl)-psimdl)%1
    psi15mdl=(np.append(psimdl[4:],psimdl[:4] +Qmdl)-psimdl)%1
    psi16mdl=(np.append(psimdl[5:],psimdl[:5] +Qmdl)-psimdl)%1
    psi17mdl=(np.append(psimdl[6:],psimdl[:6] +Qmdl)-psimdl)%1
    psi18mdl=(np.append(psimdl[7:],psimdl[:7] +Qmdl)-psimdl)%1
    psi19mdl=(np.append(psimdl[8:],psimdl[:8] +Qmdl)-psimdl)%1
    psi110mdl=(np.append(psimdl[9:],psimdl[:9] +Qmdl)-psimdl)%1
    psi111mdl=(np.append(psimdl[10:],psimdl[:10] +Qmdl)-psimdl)%1
    # >>>>>>>>>>

    #-- Global parameters of the driven motion
    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    
    # <<<<<<<<<< ICH
    # ========== OLD    
    psi12all=np.zeros((len(bpm),len(Files)))
    psi13all=np.zeros((len(bpm),len(Files)))
    psi14all=np.zeros((len(bpm),len(Files)))
    psi15all=np.zeros((len(bpm),len(Files)))
    psi16all=np.zeros((len(bpm),len(Files)))
    psi17all=np.zeros((len(bpm),len(Files)))
    psi18all=np.zeros((len(bpm),len(Files)))
    psi19all=np.zeros((len(bpm),len(Files)))
    psi110all=np.zeros((len(bpm),len(Files)))
    psi111all=np.zeros((len(bpm),len(Files)))
    # >>>>>>>>>>
    
    for i in range(len(Files)):
        if plane=='H': psid=bd*2*np.pi*np.array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid=bd*2*np.pi*np.array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*np.pi*Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+np.pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*np.pi*Qd
        Psi=np.arctan((1-r)/(1+r)*np.tan(Psid))%np.pi
        for k in range(len(bpm)):
            if Psid[k]%(2*np.pi)>np.pi: Psi[k]=Psi[k]+np.pi
        psi=Psi-Psi[0]
        psi[k_bpmac:]=psi[k_bpmac:]+2*np.pi*Q
        
        # <<<<<<<<<< ICH
        psiij = [None] * 10
        for j in range(1,11):
            psiij[j-1] = (np.append(psi[j:],psi[:j] +2*np.pi*Q)-psi)/(2*np.pi)
        # ========== OLD 
        psi12=(np.append(psi[1:],psi[0] +2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        #psi12=(np.append(psi[1:],psi[0])-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi13=(np.append(psi[2:],psi[:2]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi14=(np.append(psi[3:],psi[:3]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi15=(np.append(psi[4:],psi[:4]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi16=(np.append(psi[5:],psi[:5]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi17=(np.append(psi[6:],psi[:6]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi18=(np.append(psi[7:],psi[:7]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi19=(np.append(psi[8:],psi[:8]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi110=(np.append(psi[9:],psi[:9]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi111=(np.append(psi[10:],psi[:10]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        # >>>>>>>>>>
        
       
        for k in range(len(bpm)):
            # <<<<<<<<<< ICH
            for j in range(0,10):
                psiijall[j][k][i] = psiij[j][k]
            
            # ========== OLD     
            psi12all[k][i]=psi12[k]
            psi13all[k][i]=psi13[k]
            psi14all[k][i]=psi14[k]
            psi15all[k][i]=psi15[k]
            psi16all[k][i]=psi16[k]
            psi17all[k][i]=psi17[k]
            psi18all[k][i]=psi18[k]
            psi19all[k][i]=psi19[k]
            psi110all[k][i]=psi110[k]
            psi111all[k][i]=psi111[k]
            
            
           
         
           # >>>>>>>>>>
    #-- Output
   
    result={}
    result_={}
    muave=0.0  #-- mu is the same as psi but w/o mod
    muave_ = 0.0
    for k in range(len(bpm)):
        # <<<<<<<<<< ICH
        psiijave = [None] * 10
        psiijstd = [None] * 10
        bnj = [None] * 11
        
        for j in range(0,11):
            bnj[j] = str.upper(bpm[(k + j) % len(bpm)][1])
        
        for j in range(0,10):
            psiijave[j] = phase.calc_phase_mean(psiijall[j][k],1)
            psiijstd[j] = phase.calc_phase_std(psiijall[j][k],1)
            result_["".join([plane, bnj[0], bnj[j + 1]])] = [psiijave[j],psiijstd[j],psiijmdl[j][k]]
            
        muave_ += psiijave[0]
        try:    result_[bpm[k][1]]=[psiijave[0],psiijstd[0],psiijave[1],psiijstd[1],psiijmdl[0][k],psiijmdl[1][k],bpm[k+1][1]]
        except: result_[bpm[k][1]]=[psiijave[0],psiijstd[0],psiijave[1],psiijstd[1],psiijmdl[0][k],psiijmdl[1][k],bpm[0][1]]    #-- The last BPM
        
        
        # ========== OLD 
        psi12ave = phase.calc_phase_mean(psi12all[k],1)
        psi12std = phase.calc_phase_std(psi12all[k],1)
        psi13ave = phase.calc_phase_mean(psi13all[k],1)
        psi13std = phase.calc_phase_std(psi13all[k],1)
        psi14ave = phase.calc_phase_mean(psi14all[k],1)
        psi14std = phase.calc_phase_std(psi14all[k],1)
        psi15ave = phase.calc_phase_mean(psi15all[k],1)
        psi15std = phase.calc_phase_std(psi15all[k],1)
        psi16ave = phase.calc_phase_mean(psi16all[k],1)
        psi16std = phase.calc_phase_std(psi16all[k],1)
        psi17ave = phase.calc_phase_mean(psi17all[k],1)
        psi17std = phase.calc_phase_std(psi17all[k],1)
        psi18ave = phase.calc_phase_mean(psi18all[k],1)
        psi18std = phase.calc_phase_std(psi18all[k],1)
        psi19ave = phase.calc_phase_mean(psi19all[k],1)
        psi19std = phase.calc_phase_std(psi19all[k],1)
        psi110ave = phase.calc_phase_mean(psi110all[k],1)
        psi110std = phase.calc_phase_std(psi110all[k],1)
        psi111ave = phase.calc_phase_mean(psi111all[k],1)
        psi111std = phase.calc_phase_std(psi111all[k],1)
        
        
        
        
        muave=muave+psi12ave
        try:    result[bpm[k][1]]=[psi12ave,psi12std,psi13ave,psi13std,psi12mdl[k],psi13mdl[k],bpm[k+1][1]]
        except: result[bpm[k][1]]=[psi12ave,psi12std,psi13ave,psi13std,psi12mdl[k],psi13mdl[k],bpm[0][1]]    #-- The last BPM

        bn1 = str.upper(bpm[k%len(bpm)][1])
        bn2 = str.upper(bpm[(k+1)%len(bpm)][1])
        bn3 = str.upper(bpm[(k+2)%len(bpm)][1])
        bn4 = str.upper(bpm[(k+3)%len(bpm)][1])
        bn5 = str.upper(bpm[(k+4)%len(bpm)][1])
        bn6 = str.upper(bpm[(k+5)%len(bpm)][1])
        bn7 = str.upper(bpm[(k+6)%len(bpm)][1])
        bn8 = str.upper(bpm[(k+7)%len(bpm)][1])
        bn9 = str.upper(bpm[(k+8)%len(bpm)][1])
        bn10 = str.upper(bpm[(k+9)%len(bpm)][1])
        bn11 = str.upper(bpm[(k+10)%len(bpm)][1])

        if plane=='H':
            result["".join(['H',bn1,bn2])] = [psi12ave,psi12std,psi12mdl[k]]
            result["".join(['H',bn1,bn3])] = [psi13ave,psi13std,psi13mdl[k]]
            result["".join(['H',bn1,bn4])] = [psi14ave,psi14std,psi14mdl[k]]
            result["".join(['H',bn1,bn5])] = [psi15ave,psi15std,psi15mdl[k]]
            result["".join(['H',bn1,bn6])] = [psi16ave,psi16std,psi16mdl[k]]
            result["".join(['H',bn1,bn7])] = [psi17ave,psi17std,psi17mdl[k]]
            result["".join(['H',bn1,bn8])] = [psi18ave,psi18std,psi18mdl[k]]
            result["".join(['H',bn1,bn9])] = [psi19ave,psi19std,psi19mdl[k]]
            result["".join(['H',bn1,bn10])] = [psi110ave,psi110std,psi110mdl[k]]
            result["".join(['H',bn1,bn11])] = [psi111ave,psi111std,psi111mdl[k]]
        elif plane=='V':
            result["".join(['V',bn1,bn2])] = [psi12ave,psi12std,psi12mdl[k]]
            result["".join(['V',bn1,bn3])] = [psi13ave,psi13std,psi13mdl[k]]
            result["".join(['V',bn1,bn4])] = [psi14ave,psi14std,psi14mdl[k]]
            result["".join(['V',bn1,bn5])] = [psi15ave,psi15std,psi15mdl[k]]
            result["".join(['V',bn1,bn6])] = [psi16ave,psi16std,psi16mdl[k]]
            result["".join(['V',bn1,bn7])] = [psi17ave,psi17std,psi17mdl[k]]
            result["".join(['V',bn1,bn8])] = [psi18ave,psi18std,psi18mdl[k]]
            result["".join(['V',bn1,bn9])] = [psi19ave,psi19std,psi19mdl[k]]
            result["".join(['V',bn1,bn10])] = [psi110ave,psi110std,psi110mdl[k]]
            result["".join(['V',bn1,bn11])] = [psi111ave,psi111std,psi111mdl[k]]
            
             # >>>>>>>>>>
    return [result_,muave_,bpm]


def get_free_beta_from_amp_eq(MADTwiss_ac, Files, Qd, Q, psid_ac2bpmac, plane, bd, op):
    #-- Select common BPMs
    all_bpms = Utilities.bpm.model_intersect(
        Utilities.bpm.intersect(Files),
        MADTwiss_ac,
    )
    all_bpms = [(b[0], str.upper(b[1])) for b in all_bpms]

    good_bpms_for_kick = intersect_bpm_list_with_arc_bpms(
        intersect_bpms_list_with_bad_known_bpms(all_bpms)
    )

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op == "1" and bd == 1:
        s_lastbpm = MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L2.B1']]
    if op == "1" and bd == -1:
        s_lastbpm = MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L8.B2']]

    #-- Determine the BPM closest to the AC dipole and its position
    bpm_ac1 = ""
    bpm_ac2 = ""
    for b in psid_ac2bpmac.keys():
        if '5L4' in b:
            bpm_ac1 = b
        if '6L4' in b:
            bpm_ac2 = b

    bpm_names = list(zip(*all_bpms)[1])
    if bpm_ac1 in bpm_names:
        k_bpmac = bpm_names.index(bpm_ac1)
        bpmac = bpm_ac1
    elif bpm_ac2 in bpm_names:
        k_bpmac = bpm_names.index(bpm_ac2)
        bpmac = bpm_ac2
    else:
        print >> sys.stderr, 'WARN: BPMs next to AC dipoles missing. Was looking for: bpm_ac1="{0}" and bpm_ac2="{1}"'.format(bpm_ac1, bpm_ac2)
        return [{}, 0.0, [], [float('nan'), float('nan')]]

    #-- Model beta and phase advance
    if plane == 'H':
        betmdl = np.array(
            [MADTwiss_ac.BETX[MADTwiss_ac.indx[b[1]]] for b in all_bpms]
        )
    if plane == 'V':
        betmdl = np.array(
            [MADTwiss_ac.BETY[MADTwiss_ac.indx[b[1]]] for b in all_bpms]
        )

    #-- Global parameters of the driven motion
    r = sin(np.pi * (Qd - Q)) / sin(np.pi * (Qd + Q))

    # TODO: Use std to compute errorbars.
    sqrt2j, sqrt2j_std = get_kick_from_bpm_list_w_ACdipole(
        MADTwiss_ac, good_bpms_for_kick, Files, plane
    )

    #-- Loop for files
    betall = np.zeros((len(all_bpms), len(Files)))
    for i in range(len(Files)):
        if plane == 'H':
            amp = np.array(
                [2 * Files[i].AMPX[Files[i].indx[b[1]]] for b in all_bpms]
            )
            psid = bd * 2 * np.pi * np.array(
                [Files[i].MUX[Files[i].indx[b[1]]] for b in all_bpms]
            )  # bd flips B2 phase to B1 direction
        if plane == 'V':
            amp = np.array(
                [2 * Files[i].AMPY[Files[i].indx[b[1]]] for b in all_bpms]
            )
            psid = bd * 2 * np.pi * np.array(
                [Files[i].MUY[Files[i].indx[b[1]]] for b in all_bpms]
            )  # bd flips B2 phase to B1 direction

        # This loop is just to fix the phase jump at the beginning of the ring.
        for k in range(len(all_bpms)):
            try:
                if all_bpms[k][0] > s_lastbpm:
                    psid[k] += 2 * np.pi * Qd
            except:
                pass

        psid = psid - (psid[k_bpmac] - psid_ac2bpmac[bpmac])
        Psid = psid + np.pi * Qd
        Psid[k_bpmac:] = Psid[k_bpmac:] - 2 * np.pi * Qd
        bet = ((amp / sqrt2j[i]) ** 2 *
               (1 + r ** 2 + 2 * r * np.cos(2 * Psid)) / (1 - r ** 2))
        for bpm_index in range(len(all_bpms)):
            betall[bpm_index][i] = bet[bpm_index]

    #-- Output
    result = {}
    bb = []
    for k in range(len(all_bpms)):
        betave = np.mean(betall[k])
        betstd = np.std(betall[k])
        bb.append((betave - betmdl[k]) / betmdl[k])
        result[all_bpms[k][1]] = [betave, betstd, all_bpms[k][0]]
    bb = math.sqrt(np.mean(np.array(bb) ** 2))

    return [result, bb, all_bpms]


def intersect_bpm_list_with_arc_bpms(bpms_list):
    bpm_arcs = []
    # Selecting ARC BPMs
    for b in bpms_list:
        if ((b[1][4]) == '1' and (b[1][5]) >= '4') or (b[1][4]) == '2':
            bpm_arcs.append(b)
    return bpm_arcs


BAD_BPM_LIST = ['BPM.15R8.B1', 'BPM.16R3.B1', 'BPM.31L5.B1', 'BPM.23L6.B1',
                'BPM.22R8.B1', 'BPM.11R6.B1', 'BPM.18R6.B1', 'BPM.34R5.B1']


def intersect_bpms_list_with_bad_known_bpms(bpms_list):
    bpm_arcs_clean = []
    for b in bpms_list:
        if b[1] not in BAD_BPM_LIST:
            bpm_arcs_clean.append(b)
    return bpm_arcs_clean


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
        betmdl = np.array(
            [MADTwiss_ac.BETX[MADTwiss_ac.indx[bpm[1]]] for bpm in bpm_list]
        )
    if plane == 'V':
        betmdl = np.array(
            [MADTwiss_ac.BETY[MADTwiss_ac.indx[bpm[1]]] for bpm in bpm_list]
        )

    actions_sqrt = []
    actions_sqrt_err = []

    for meas in measurements:
        if plane == 'H':
            amp = np.array(
                [2 * meas.AMPX[meas.indx[bpm[1]]] for bpm in bpm_list]
            )
        if plane == 'V':
            amp = np.array(
                [2 * meas.AMPY[meas.indx[bpm[1]]] for bpm in bpm_list]
            )
        actions_sqrt.append(np.average(amp / np.sqrt(betmdl)))
        actions_sqrt_err.append(np.std(amp / np.sqrt(betmdl)))

    return actions_sqrt, actions_sqrt_err


def GetFreeCoupling_Eq(MADTwiss,FilesX,FilesY,Qh,Qv,Qx,Qy,psih_ac2bpmac,psiv_ac2bpmac,bd):

    #-- Details of this algorithms is in http://www.agsrhichome.bnl.gov/AP/ap_notes/ap_note_410.pdf

    #-- Check linx/liny files, may be redundant
    if len(FilesX)!=len(FilesY): return [{},[]]

    #-- Select common BPMs
    bpm=Utilities.bpm.model_intersect(Utilities.bpm.intersect(FilesX+FilesY),MADTwiss)
    bpm=[(b[0],str.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    #if op=="1" and bd== 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    #if op=="1" and bd==-1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    #-- Determine the BPM closest to the AC dipole and its position
    bpmac1 = psih_ac2bpmac.keys()[0]
    bpmac2 = psih_ac2bpmac.keys()[1]

    try:
        k_bpmac = list(zip(*bpm)[1]).index(bpmac1)
        bpmac = bpmac1
    except:
        try:
            k_bpmac = list(zip(*bpm)[1]).index(bpmac2)
            bpmac = bpmac2
        except:
            print >> sys.stderr,'WARN: BPMs next to AC dipoles missing. AC dipole effects not calculated with analytic eqs for coupling'
            return [{},[]]
    
    print "WARN: GetFreeCoupling ---- maybe problem with different kicker positions for H and V kicks"
    
    #-- Global parameters of the driven motion
    dh =Qh-Qx
    dv =Qv-Qy
    rh =sin(np.pi*(Qh-Qx))/sin(np.pi*(Qh+Qx))
    rv =sin(np.pi*(Qv-Qy))/sin(np.pi*(Qv+Qy))
    rch=sin(np.pi*(Qh-Qy))/sin(np.pi*(Qh+Qy))
    rcv=sin(np.pi*(Qx-Qv))/sin(np.pi*(Qx+Qv))

    #-- Loop for files
    f1001Abs =np.zeros((len(bpm),len(FilesX)))
    f1010Abs =np.zeros((len(bpm),len(FilesX)))
    f1001xArg=np.zeros((len(bpm),len(FilesX)))
    f1001yArg=np.zeros((len(bpm),len(FilesX)))
    f1010xArg=np.zeros((len(bpm),len(FilesX)))
    f1010yArg=np.zeros((len(bpm),len(FilesX)))
    for i in range(len(FilesX)):

        #-- Read amplitudes and phases
        amph  =     np.array([FilesX[i].AMPX[FilesX[i].indx[b[1]]]    for b in bpm])
        ampv  =     np.array([FilesY[i].AMPY[FilesY[i].indx[b[1]]]    for b in bpm])
        amph01=     np.array([FilesX[i].AMP01[FilesX[i].indx[b[1]]]   for b in bpm])
        ampv10=     np.array([FilesY[i].AMP10[FilesY[i].indx[b[1]]]   for b in bpm])
        psih  =2*np.pi*np.array([FilesX[i].MUX[FilesX[i].indx[b[1]]]     for b in bpm])
        psiv  =2*np.pi*np.array([FilesY[i].MUY[FilesY[i].indx[b[1]]]     for b in bpm])
        psih01=2*np.pi*np.array([FilesX[i].PHASE01[FilesX[i].indx[b[1]]] for b in bpm])
        psiv10=2*np.pi*np.array([FilesY[i].PHASE10[FilesY[i].indx[b[1]]] for b in bpm])
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
        psih=psih-(psih[k_bpmac]-psih_ac2bpmac[bpmac])
        psiv=psiv-(psiv[k_bpmac]-psiv_ac2bpmac[bpmac])

        Psih=psih-np.pi*Qh
        Psih[:k_bpmac]=Psih[:k_bpmac]+2*np.pi*Qh
        Psiv=psiv-np.pi*Qv
        Psiv[:k_bpmac]=Psiv[:k_bpmac]+2*np.pi*Qv

        Psix=np.arctan((1-rh)/(1+rh)*np.tan(Psih))%np.pi
        Psiy=np.arctan((1-rv)/(1+rv)*np.tan(Psiv))%np.pi
        for k in range(len(bpm)):
            if Psih[k]%(2*np.pi)>np.pi: Psix[k]=Psix[k]+np.pi
            if Psiv[k]%(2*np.pi)>np.pi: Psiy[k]=Psiy[k]+np.pi

        psix=Psix-np.pi*Qx
        psix[k_bpmac:]=psix[k_bpmac:]+2*np.pi*Qx
        psiy=Psiy-np.pi*Qy
        psiy[k_bpmac:]=psiy[k_bpmac:]+2*np.pi*Qy

        #-- Construct f1001h, f1001v, f1010h, f1010v (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
        f1001h=1/math.sqrt(1-rv**2)*(np.exp(-1j*(Psiv-Psiy))*f1001hv+rv*np.exp( 1j*(Psiv+Psiy))*f1010hv)
        f1010h=1/math.sqrt(1-rv**2)*(np.exp( 1j*(Psiv-Psiy))*f1010hv+rv*np.exp(-1j*(Psiv+Psiy))*f1001hv)
        f1001v=1/math.sqrt(1-rh**2)*(np.exp( 1j*(Psih-Psix))*f1001vh+rh*np.exp(-1j*(Psih+Psix))*np.conjugate(f1010vh))
        f1010v=1/math.sqrt(1-rh**2)*(np.exp( 1j*(Psih-Psix))*f1010vh+rh*np.exp(-1j*(Psih+Psix))*np.conjugate(f1001vh))

        #-- Construct f1001 and f1010 from h and v BPMs (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
        g1001h          =np.exp(-1j*((psih-psih[k_bpmac])-(psiy-psiy[k_bpmac])))*(ampv/amph*amph[k_bpmac]/ampv[k_bpmac])*f1001h[k_bpmac]
        g1001h[:k_bpmac]=1/(np.exp(2*np.pi*1j*(Qh-Qy))-1)*(f1001h-g1001h)[:k_bpmac]
        g1001h[k_bpmac:]=1/(1-np.exp(-2*np.pi*1j*(Qh-Qy)))*(f1001h-g1001h)[k_bpmac:]

        g1010h          =np.exp(-1j*((psih-psih[k_bpmac])+(psiy-psiy[k_bpmac])))*(ampv/amph*amph[k_bpmac]/ampv[k_bpmac])*f1010h[k_bpmac]
        g1010h[:k_bpmac]=1/(np.exp(2*np.pi*1j*(Qh+Qy))-1)*(f1010h-g1010h)[:k_bpmac]
        g1010h[k_bpmac:]=1/(1-np.exp(-2*np.pi*1j*(Qh+Qy)))*(f1010h-g1010h)[k_bpmac:]

        g1001v          =np.exp(-1j*((psix-psix[k_bpmac])-(psiv-psiv[k_bpmac])))*(amph/ampv*ampv[k_bpmac]/amph[k_bpmac])*f1001v[k_bpmac]
        g1001v[:k_bpmac]=1/(np.exp(2*np.pi*1j*(Qx-Qv))-1)*(f1001v-g1001v)[:k_bpmac]
        g1001v[k_bpmac:]=1/(1-np.exp(-2*np.pi*1j*(Qx-Qv)))*(f1001v-g1001v)[k_bpmac:]

        g1010v          =np.exp(-1j*((psix-psix[k_bpmac])+(psiv-psiv[k_bpmac])))*(amph/ampv*ampv[k_bpmac]/amph[k_bpmac])*f1010v[k_bpmac]
        g1010v[:k_bpmac]=1/(np.exp(2*np.pi*1j*(Qx+Qv))-1)*(f1010v-g1010v)[:k_bpmac]
        g1010v[k_bpmac:]=1/(1-np.exp(-2*np.pi*1j*(Qx+Qv)))*(f1010v-g1010v)[k_bpmac:]

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
        for k in range(len(bpm)):
            f1001Abs[k][i] =math.sqrt(abs(f1001x[k]*f1001y[k]))
            f1010Abs[k][i] =math.sqrt(abs(f1010x[k]*f1010y[k]))
            f1001xArg[k][i]=np.angle(f1001x[k])%(2*np.pi)
            f1001yArg[k][i]=np.angle(f1001y[k])%(2*np.pi)
            f1010xArg[k][i]=np.angle(f1010x[k])%(2*np.pi)
            f1010yArg[k][i]=np.angle(f1010y[k])%(2*np.pi)

    #-- Output
    fwqw={}
    goodbpm=[]
    for k in range(len(bpm)):

        #-- Bad BPM flag based on phase
        badbpm=0
        f1001xArgAve = phase.calc_phase_mean(f1001xArg[k],2*np.pi)
        f1001yArgAve = phase.calc_phase_mean(f1001yArg[k],2*np.pi)
        f1010xArgAve = phase.calc_phase_mean(f1010xArg[k],2*np.pi)
        f1010yArgAve = phase.calc_phase_mean(f1010yArg[k],2*np.pi)
        if min(abs(f1001xArgAve-f1001yArgAve),2*np.pi-abs(f1001xArgAve-f1001yArgAve))>np.pi/2: badbpm=1
        if min(abs(f1010xArgAve-f1010yArgAve),2*np.pi-abs(f1010xArgAve-f1010yArgAve))>np.pi/2: badbpm=1

        #-- Output
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
            fwqw[bpm[k][1]] = [[f1001Ave          ,f1001AbsStd       ,f1010Ave          ,f1010AbsStd       ],
                             [f1001ArgAve/(2*np.pi),f1001ArgStd/(2*np.pi),f1010ArgAve/(2*np.pi),f1010ArgStd/(2*np.pi)]]  #-- Phases renormalized to [0,1)
            goodbpm.append(bpm[k])

    #-- Global parameters not implemented yet
    fwqw['Global']=['"null"','"null"']

    return [fwqw,goodbpm]

def GetFreeIP2_Eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,oa,op):

    #-- Common BPMs
    bpm = Utilities.bpm.model_intersect(Utilities.bpm.intersect(Files),MADTwiss)
    bpm=[(b[0],str.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and bd== 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    if op=="1" and bd==-1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    #-- Determine the BPM closest to the AC dipole and its position
    for b in psid_ac2bpmac.keys():
        if '5L4' in b: bpmac1=b
        if '6L4' in b: bpmac2=b
    try:
        k_bpmac=list(zip(*bpm)[1]).index(bpmac1)
        bpmac=bpmac1
    except:
        try:
            k_bpmac=list(zip(*bpm)[1]).index(bpmac2)
            bpmac=bpmac2
        except:
            return [{},[]]

    #-- Global parameters of the driven motion
    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))

    #-- Determine Psid (w.r.t the AC dipole) for each file
    Psidall=[]
    for i in range(len(Files)):
        if plane=='H': psid=bd*2*np.pi*np.array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid=bd*2*np.pi*np.array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*np.pi*Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+np.pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*np.pi*Qd
        Psidall.append(Psid)

    #-- Loop for IPs
    result={}
    for ip in ('1','2','5','8'):

        bpml='BPMSW.1L'+ip+'.'+oa[3:]
        bpmr='BPMSW.1R'+ip+'.'+oa[3:]
        if (bpml in zip(*bpm)[1]) and (bpmr in zip(*bpm)[1]):

            #-- Model values
            L=0.5*(MADTwiss.S[MADTwiss.indx[bpmr]]-MADTwiss.S[MADTwiss.indx[bpml]])
            if L<0: L+=0.5*MADTwiss.LENGTH
            if plane=='H':
                betlmdl=MADTwiss.BETX[MADTwiss.indx[bpml]]
                alflmdl=MADTwiss.ALFX[MADTwiss.indx[bpml]]
            if plane=='V':
                betlmdl=MADTwiss.BETY[MADTwiss.indx[bpml]]
                alflmdl=MADTwiss.ALFY[MADTwiss.indx[bpml]]
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
                        al=Files[i].AMPX[Files[i].indx[bpml]]
                        ar=Files[i].AMPX[Files[i].indx[bpmr]]
                    if plane=='V':
                        al=Files[i].AMPY[Files[i].indx[bpml]]
                        ar=Files[i].AMPY[Files[i].indx[bpmr]]
                    Psidl=Psidall[i][list(zip(*bpm)[1]).index(bpml)]
                    Psidr=Psidall[i][list(zip(*bpm)[1]).index(bpmr)]
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

    all_bpms_x = Utilities.bpm.model_intersect(Utilities.bpm.intersect(files[0]), MADTwiss_ac )
    all_bpms_y = Utilities.bpm.model_intersect(Utilities.bpm.intersect(files[1]), MADTwiss_ac )
    all_bpms_x = [(b[0], str.upper(b[1])) for b in all_bpms_x]
    all_bpms_y = [(b[0], str.upper(b[1])) for b in all_bpms_y]

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
