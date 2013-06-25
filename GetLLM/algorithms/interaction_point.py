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
from numpy import sin, cos, tan

import Utilities.bpm


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


#===================================================================================================
# helper-functions
#===================================================================================================

def get_ip(ip_num, measured, model):
    bpm_left, bpm_right = _find_bpm(ip_num, model, measured)

    if DEBUG:
        print "getIP", ip_num

    if bpm_left is None or bpm_right is None:
        print "skipping ip%1s calculation, no BPM found" % ip_num
        betahor = [ip_num, 0, 0, 0, 0, 0, 0]
        betaver = [ip_num, 0, 0, 0, 0, 0, 0]
        return [betahor, betaver]

    # model
    sxl = model.S[model.indx[bpm_left]]
    sip = model.S[model.indx["ip"+ip_num]]
    sxr = model.S[model.indx[bpm_right]]
    betaipx = model.BETX[model.indx["ip"+ip_num]]
    betaipy = model.BETY[model.indx["ip"+ip_num]]

    # measured value
    betaxl = measured[0][bpm_left][0]
    betayl = measured[1][bpm_left][0]

    betaxr = measured[0][bpm_right][0]
    betayr = measured[1][bpm_right][0]

    deltaphimodel = abs(model.MUX[model.indx[bpm_right]]-model.MUX[model.indx[bpm_left]])

    L = ((sip-sxl)+(sxr-sip))/2
    betastar = (2*math.sqrt(betaxl)*math.sqrt(betaxr)*sin(deltaphimodel*2*np.pi))/(betayl+betayr-2*math.sqrt(betaxl)*math.sqrt(betaxr)*cos(2*np.pi*deltaphimodel))*L
    location = ((betaxl-betaxr)/(betaxl+betaxr-2*math.sqrt(betaxl)*math.sqrt(betaxr)*cos(2*np.pi*deltaphimodel)))*L

    deltaphi = (math.atan((L-location)/betastar)+math.atan((L+location)/betastar))/(2*np.pi)

    betahor = [ip_num, betastar, location, deltaphi, betaipx, deltaphimodel, 0]

    if DEBUG:
        print "horizontal betastar for ", ip_num, " is ", str(betastar), " at location ", str(location), " of ip_num center with phase advance ", str(deltaphi)

    #vertical
    deltaphimodel = abs(model.MUY[model.indx[bpm_right]]-model.MUY[model.indx[bpm_left]])


    betastar = (2*math.sqrt(betayl)*math.sqrt(betayr)*sin(deltaphimodel*2*np.pi))/(betayl+betayr-2*math.sqrt(betayl)*math.sqrt(betayr)*cos(2*np.pi*deltaphimodel))*L
    location = ((betayl-betayr)/(betayl+betayr-2*math.sqrt(betayl)*math.sqrt(betayr)*cos(2*np.pi*deltaphimodel)))*L

    deltaphi = (math.atan((L-location)/betastar)+math.atan((L+location)/betastar))/(2*np.pi)

    betaver = [ip_num, betastar, location, deltaphi, betaipy, deltaphimodel, 0]

    if DEBUG:
        print "vertical betastar for ", ip_num, " is ", str(betastar), " at location ", str(location), " of ip_num center with phase advance ", str(deltaphi)

    return [betahor, betaver]


def _find_bpm(ip_num, model, measured):
    bpm_left = None
    bpm_right = None
    
    for bpm_name in model.NAME:
        if "BPMSW.1L"+ip_num in bpm_name:
            bpm_left = bpm_name
            try:
                measured[0][bpm_left][0]  # Test if bpm_left exists
            except KeyError:
                traceback.print_exc()
                bpm_left = None
        if "BPMSW.1R"+ip_num in bpm_name:
            bpm_right = bpm_name
            try:
                test = measured[0][bpm_right][0]  # @UnusedVariable
            except KeyError:
                bpm_right = None
                
    return [bpm_left, bpm_right]


def get_ip_2(mad_twiss, files, Q, plane, beam_direction, accel, lhc_phase):
    #-- Common BPMs
    bpm = Utilities.bpm.model_intersect(Utilities.bpm.intersect(files), mad_twiss)
    bpm = [(b[0], str.upper(b[1])) for b in bpm]
    
    bpm_names = [ b[1] for b in bpm]

    #-- Loop for IPs
    result = {}
    for ip in ('1', '2', '5', '8'):

        bpml = 'BPMSW.1L'+ip+'.'+accel[3:]
        bpmr = 'BPMSW.1R'+ip+'.'+accel[3:]

        if (bpml in bpm_names) and (bpmr in bpm_names):

            #-- Model values
            L = 0.5*(mad_twiss.S[mad_twiss.indx[bpmr]] - mad_twiss.S[mad_twiss.indx[bpml]])
            if L < 0: 
                L += 0.5*mad_twiss.LENGTH  #-- For sim starting in the middle of an IP
            if plane == 'H':
                betlmdl = mad_twiss.BETX[mad_twiss.indx[bpml]]
                alflmdl = mad_twiss.ALFX[mad_twiss.indx[bpml]]
            if plane == 'V':
                betlmdl = mad_twiss.BETY[mad_twiss.indx[bpml]]
                alflmdl = mad_twiss.ALFY[mad_twiss.indx[bpml]]
            betsmdl = betlmdl/(1+alflmdl**2)
            betmdl = betlmdl-2*alflmdl*L+L**2/betsmdl
            alfmdl = alflmdl-L/betsmdl
            dsmdl = alfmdl*betsmdl

            #-- Measurement for each file
            betall = []
            alfall = []
            betsall = []
            dsall = []
            rt2j_all = []
            for t_f in files: # t_f := twiss_file
                try:
                    if plane == 'H':
                        a_l = t_f.AMPX[t_f.indx[bpml]]
                        a_r = t_f.AMPX[t_f.indx[bpmr]]
                        if bpm_names.index(bpmr) > bpm_names.index(bpml):
                            dpsi = 2*np.pi*beam_direction*(t_f.MUX[t_f.indx[bpmr]]-t_f.MUX[t_f.indx[bpml]])
                        else:
                            dpsi = 2*np.pi*(Q+beam_direction*(t_f.MUX[t_f.indx[bpmr]]-t_f.MUX[t_f.indx[bpml]]))
                    elif plane == 'V':
                        a_l = t_f.AMPY[t_f.indx[bpml]]
                        a_r = t_f.AMPY[t_f.indx[bpmr]]
                        if bpm_names.index(bpmr) > bpm_names.index(bpml):
                            dpsi = 2*np.pi*beam_direction*(t_f.MUY[t_f.indx[bpmr]]-t_f.MUY[t_f.indx[bpml]])
                        else:
                            dpsi = 2*np.pi*(Q+beam_direction*(t_f.MUY[t_f.indx[bpmr]]-t_f.MUY[t_f.indx[bpml]]))
                    else:
                        raise ValueError("plane is neither 'H' nor 'V'.")

                    #-- To compensate the phase shift by tune
                    if lhc_phase == '1':
                        if (beam_direction==1 and ip=='2') or (beam_direction==-1 and ip=='8'): 
                            dpsi += 2*np.pi*Q

                    #-- bet, alf, and math.sqrt(2J) from amp and phase advance
                    bet = L*(a_l**2+a_r**2+2*a_l*a_r*cos(dpsi))/(2*a_l*a_r*sin(dpsi))
                    alf = (a_l**2-a_r**2)/(2*a_l*a_r*sin(dpsi))
                    bets = bet/(1+alf**2)
                    d_s = alf*bets
                    rt2j = math.sqrt(a_l*a_r*sin(dpsi)/(2*L))
                    betall.append(bet)
                    alfall.append(alf)
                    betsall.append(bets)
                    dsall.append(d_s)
                    rt2j_all.append(rt2j)
                except ValueError:
                    traceback.print_exc() # math domain error
                except ZeroDivisionError:
                    traceback.print_exc()
                except Exception:
                    traceback.print_exc()

            #-- Ave and Std
            betall = np.array(betall)
            betave = np.mean(betall)
            betstd = math.sqrt(np.mean((betall-betave)**2))
            
            alfall = np.array(alfall)
            alfave = np.mean(alfall)
            alfstd = math.sqrt(np.mean((alfall-alfave)**2))
            
            betsall = np.array(betsall)
            betsave = np.mean(betsall)
            betsstd = math.sqrt(np.mean((betsall-betsave)**2))
            
            dsall = np.array(dsall)
            dsave = np.mean(dsall)
            dsstd = math.sqrt(np.mean((dsall-dsave)**2))
            
            rt2j_all = np.array(rt2j_all)
            rt2j_ave = np.mean(rt2j_all)
            rt2j_std = math.sqrt(np.mean((rt2j_all-rt2j_ave)**2))
            
            result['IP'+ip] = [betave, betstd, betmdl, alfave, alfstd, alfmdl, betsave, betsstd, betsmdl, dsave, dsstd, dsmdl, rt2j_ave, rt2j_std]

    return result


def get_ip_from_phase(MADTwiss,psix,psiy,oa):

    IP=('1','2','5','8')
    result={}
    for i in IP:
        bpml='BPMSW.1L'+i+'.'+oa[3:]
        bpmr=bpml.replace('L','R')
        try:
            if psix[bpml][-1]==bpmr:
                #-- Model
                L       =0.5*(MADTwiss.S[MADTwiss.indx[bpmr]]-MADTwiss.S[MADTwiss.indx[bpml]])
                dpsixmdl=MADTwiss.MUX[MADTwiss.indx[bpmr]]-MADTwiss.MUX[MADTwiss.indx[bpml]]
                dpsiymdl=MADTwiss.MUY[MADTwiss.indx[bpmr]]-MADTwiss.MUY[MADTwiss.indx[bpml]]
                betxmdl =MADTwiss.BETX[MADTwiss.indx[bpml]]/(1+MADTwiss.ALFX[MADTwiss.indx[bpml]]**2)
                betymdl =MADTwiss.BETY[MADTwiss.indx[bpml]]/(1+MADTwiss.ALFY[MADTwiss.indx[bpml]]**2)
                #-- For sim starting in the middle of an IP
                if L < 0:
                    L += 0.5 * MADTwiss.LENGTH
                    dpsixmdl += MADTwiss.Q1
                    dpsiymdl += MADTwiss.Q2
                #-- Measurement
                dpsix   =psix['BPMSW.1L'+i+'.'+oa[3:]][0]
                dpsiy   =psiy['BPMSW.1L'+i+'.'+oa[3:]][0]
                dpsixstd=psix['BPMSW.1L'+i+'.'+oa[3:]][1]
                dpsiystd=psiy['BPMSW.1L'+i+'.'+oa[3:]][1]
                betx    =L/tan(np.pi*dpsix)
                bety    =L/tan(np.pi*dpsiy)
                betxstd =L*np.pi*dpsixstd/(2*sin(np.pi*dpsix)**2)
                betystd =L*np.pi*dpsiystd/(2*sin(np.pi*dpsiy)**2)
                result['IP'+i]=[2*L,betx,betxstd,betxmdl,bety,betystd,betymdl,dpsix,dpsixstd,dpsixmdl,dpsiy,dpsiystd,dpsiymdl]
        except: pass
        #-- This part due to the format difference of phasef2 (from the model)
        try:
            if psix[bpml][-2]==bpmr:
                #-- Model
                L       =0.5*(MADTwiss.S[MADTwiss.indx[bpmr]]-MADTwiss.S[MADTwiss.indx[bpml]])
                dpsixmdl=MADTwiss.MUX[MADTwiss.indx[bpmr]]-MADTwiss.MUX[MADTwiss.indx[bpml]]
                dpsiymdl=MADTwiss.MUY[MADTwiss.indx[bpmr]]-MADTwiss.MUY[MADTwiss.indx[bpml]]
                betxmdl =MADTwiss.BETX[MADTwiss.indx[bpml]]/(1+MADTwiss.ALFX[MADTwiss.indx[bpml]]**2)
                betymdl =MADTwiss.BETY[MADTwiss.indx[bpml]]/(1+MADTwiss.ALFY[MADTwiss.indx[bpml]]**2)
                #-- For sim starting in the middle of an IP
                if L < 0:
                    L += 0.5 * MADTwiss.LENGTH
                    dpsixmdl += MADTwiss.Q1
                    dpsiymdl += MADTwiss.Q2
                #-- Measurement
                dpsix   =psix['BPMSW.1L'+i+'.'+oa[3:]][0]
                dpsiy   =psiy['BPMSW.1L'+i+'.'+oa[3:]][0]
                dpsixstd=psix['BPMSW.1L'+i+'.'+oa[3:]][1]
                dpsiystd=psiy['BPMSW.1L'+i+'.'+oa[3:]][1]
                betx    =L/tan(np.pi*dpsix)
                bety    =L/tan(np.pi*dpsiy)
                betxstd =L*np.pi*dpsixstd/(2*sin(np.pi*dpsix)**2)
                betystd =L*np.pi*dpsiystd/(2*sin(np.pi*dpsiy)**2)
                result['IP'+i]=[2*L,betx,betxstd,betxmdl,bety,betystd,betymdl,dpsix,dpsixstd,dpsixmdl,dpsiy,dpsiystd,dpsiymdl]
        except: pass

    return result


#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def get_free_ip_2(MADTwiss,MADTwiss_ac,IP,plane,oa):

    for i in ('1','2','5','8'):

        bpml='BPMSW.1L'+i+'.'+oa[3:]
        bpmr='BPMSW.1R'+i+'.'+oa[3:]
        if 'IP'+i in IP:

            L=0.5*(MADTwiss.S[MADTwiss.indx[bpmr]]-MADTwiss.S[MADTwiss.indx[bpml]])
            if L<0: L+=0.5*MADTwiss.LENGTH
            #-- bet and alf at the left BPM
            if plane=='H':
                betl =MADTwiss.BETX[MADTwiss.indx[bpml]]; betdl=MADTwiss_ac.BETX[MADTwiss_ac.indx[bpml]]
                alfl =MADTwiss.ALFX[MADTwiss.indx[bpml]]; alfdl=MADTwiss_ac.ALFX[MADTwiss_ac.indx[bpml]]
            if plane=='V':
                betl =MADTwiss.BETY[MADTwiss.indx[bpml]]; betdl=MADTwiss_ac.BETY[MADTwiss_ac.indx[bpml]]
                alfl =MADTwiss.ALFY[MADTwiss.indx[bpml]]; alfdl=MADTwiss_ac.ALFY[MADTwiss_ac.indx[bpml]]
            #-- IP parameters propagated from the left BPM
            bets =betl/(1+alfl**2)       ; betds=betdl/(1+alfdl**2)
            bet  =betl-2*alfl*L+L**2/bets; betd =betdl-2*alfdl*L+L**2/betds
            alf  =alfl-L/bets            ; alfd =alfdl-L/betds
            ds   =alf*bets               ; dsd  =alfd*betds
            #-- Apply corrections
            IP['IP'+i][0]=IP['IP'+i][0]+bet-betd  ; IP['IP'+i][2] =bet
            IP['IP'+i][3]=IP['IP'+i][3]+alf-alfd  ; IP['IP'+i][5] =alf
            IP['IP'+i][6]=IP['IP'+i][6]+bets-betds; IP['IP'+i][8] =bets
            IP['IP'+i][9]=IP['IP'+i][9]+ds-dsd    ; IP['IP'+i][11]=ds

    return IP


