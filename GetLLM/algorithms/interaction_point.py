'''
Created on 27 May 2013

@author: ?, vimaier

@version: 0.0.1

GetLLM.algorithms.interaction_point.py stores helper functions for phase calculations for GetLLM.
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
from model.accelerators.accelerator import AccExcitationMode
import utils.bpm
from utils import logging_tools
import compensate_excitation


LOGGER = logging_tools.get_logger(__name__)


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# main part
#===================================================================================================

def calculate_ip(getllm_d, twiss_d, tune_d, phase_d, beta_d, mad_twiss, mad_ac, files_dict):
    '''
    Calculates ip and fills the following TfsFiles:
        getIP.out        
        getIPx.out        getIPx_free.out        getIPx_free2.out        
        getIPy.out        getIPy_free.out        getIPy_free2.out        
        getIPfromphase.out        getIPfromphase_free.out        getIPfromphase_free2.out
        
    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            lhc_phase, accel and beam_direction are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances
        'phase_d': _PhaseData (In-param, values will only be read)
            Holds results from get_phases
        'beta_d': _BetaData (In-param, values will only be read)
            Holds results from get_beta. Beta from amp and ratios will be set.
    '''
    LOGGER.info('Calculating IP')
    tfs_file = files_dict['getIP.out']
    tfs_file.add_column_names(["NAME", "BETASTARH", "BETASTARHMDL", "H", "PHIH", "PHIXH", "PHIHMDL", "BETASTARV", "BETASTARVMDL", "V", "PHIV", "PHIYV", "PHIVMDL"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    ips = ["1", "2", "3", "4", "5", "6", "7", "8"]
    measured = [beta_d.x_amp, beta_d.y_amp]
    for num_ip in ips:
        try:
            betahor, betaver = _get_ip(num_ip, measured, mad_twiss)
            list_row_entries = ['"IP' + num_ip + '"', betahor[1], betahor[4], betahor[2], betahor[3], betahor[6], betahor[5], betaver[1], betaver[4], betaver[2], betaver[3], betaver[6], betaver[5]]
            tfs_file.add_table_row(list_row_entries)
        except (KeyError, ValueError):
            pass
    
    #-- Parameters at IP1, IP2, IP5, and IP8
    ip_x = _get_ip_2(mad_ac, twiss_d.zero_dpp_x, tune_d.q1, 'H', getllm_d.accelerator, twiss_d.zero_dpp_commonbpms_x, getllm_d.lhc_phase)
    ip_y = _get_ip_2(mad_ac, twiss_d.zero_dpp_y, tune_d.q2, 'V', getllm_d.accelerator, twiss_d.zero_dpp_commonbpms_y, getllm_d.lhc_phase)

    fill_ip_tfs_file(tfs_file=files_dict['getIPx.out'],
                     column_names=["NAME", "BETX", "BETXSTD", "BETXMDL", "ALFX", "ALFXSTD", "ALFXMDL", "BETX*", "BETX*STD", "BETX*MDL", "SX*", "SX*STD", "SX*MDL", "rt(2JX)", "rt(2JX)STD"],
                     column_datatypes=["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"],
                     results=ip_x)
    fill_ip_tfs_file(tfs_file=files_dict['getIPy.out'],
                     column_names=["NAME", "BETY", "BETYSTD", "BETYMDL", "ALFY", "ALFYSTD", "ALFYMDL", "BETY*", "BETY*STD", "BETY*MDL", "SY*", "SY*STD", "SY*MDL", "rt(2JY)", "rt(2JY)STD"],
                     column_datatypes=["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"],
                     results=ip_y)

    #-- ac to free parameters at IP1, IP2, IP5, and IP8
    if getllm_d.accelerator.excitation is not AccExcitationMode.FREE:
        #-- From Eq
        ip_x_f = compensate_excitation.GetFreeIP2_Eq(mad_twiss, twiss_d.zero_dpp_x, tune_d.q1, tune_d.q1f,
                                                     phase_d.ac2bpmac_x, 'H', getllm_d.accelerator,
                                                     twiss_d.zero_dpp_commonbpms_x, getllm_d.lhc_phase)
        ip_y_f = compensate_excitation.GetFreeIP2_Eq(mad_twiss, twiss_d.zero_dpp_y, tune_d.q2, tune_d.q2f,
                                                     phase_d.ac2bpmac_y, 'V', getllm_d.accelerator,
                                                     twiss_d.zero_dpp_commonbpms_y, getllm_d.lhc_phase)
        
        col_names_x = ["NAME", "BETX", "BETXSTD", "BETXMDL", "ALFX", "ALFXSTD", "ALFXMDL", "BETX*", "BETX*STD", "BETX*MDL", "SX*", "SX*STD", "SX*MDL", "rt(2JXD)", "rt(2JXD)STD"]
        col_datatypes_x = ["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"]
        fill_ip_tfs_file(files_dict['getIPx_free.out'], col_names_x, col_datatypes_x, ip_x_f)
        
        col_names_y = ["NAME", "BETY", "BETYSTD", "BETYMDL", "ALFY", "ALFYSTD", "ALFYMDL", "BETY*", "BETY*STD", "BETY*MDL", "SY*", "SY*STD", "SY*MDL", "rt(2JYD)", "rt(2JYD)STD"]
        col_datatypes_y = ["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"]
        fill_ip_tfs_file(files_dict['getIPy_free.out'], col_names_y, col_datatypes_y, ip_y_f)
        
        
        #-- From model
        ip_x_f_2 = _get_free_ip_2(mad_twiss, mad_ac, ip_x, 'H', getllm_d.accelerator)
        ip_y_f_2 = _get_free_ip_2(mad_twiss, mad_ac, ip_y, 'V', getllm_d.accelerator)
        
        fill_ip_tfs_file(files_dict['getIPx_free2.out'], col_names_x,col_datatypes_x, results=ip_x_f_2)
        fill_ip_tfs_file(files_dict['getIPy_free2.out'], col_names_y, col_datatypes_y, results=ip_y_f_2)

    #-- IP beta* and phase from phase only
    if not phase_d.phase_advances_free_x is None and not phase_d.phase_advances_free_y is None: # Designed to run with both 
        try:
            ip_from_phase = _get_ip_from_phase(mad_twiss, phase_d.phase_advances_free_x, phase_d.phase_advances_free_y,
                                               getllm_d.accelerator)
            col_names = ["NAME", "2L", "BETX*", "BETX*STD", "BETX*MDL", "BETY*", "BETY*STD", "BETY*MDL", "PHX", "PHXSTD", "PHXMDL", "PHY", "PHYSTD", "PHYMDL"]
            col_datatypes = ["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"]
            fill_ip_tfs_file(files_dict['getIPfromphase_free.out'], col_names, col_datatypes, ip_from_phase)
        except:
            traceback.print_exc()
            print 'No output from IP from phase. H or V file missing?'

    #-- ac to free beta*
    if getllm_d.accelerator.excitation is not AccExcitationMode.FREE:
        #-- from eqs
        try:
            ip_from_phase_f =_get_ip_from_phase(mad_ac, phase_d.phase_advances_x, phase_d.phase_advances_y,
                                                getllm_d.accelerator)
            fill_ip_tfs_file(files_dict['getIPfromphase.out'], col_names, col_datatypes, ip_from_phase_f)
        except:
            traceback.print_exc()

        try:
            ip_from_phase_f2 = _get_ip_from_phase(mad_twiss, phase_d.phase_advances_free2_x,
                                                  phase_d.phase_advances_free2_y, getllm_d.accelerator)
            fill_ip_tfs_file(files_dict['getIPfromphase_free2.out'], col_names, col_datatypes, ip_from_phase_f2)
        except:
            traceback.print_exc()

# END calculate_ip ---------------------------------------------------------------------------------




#===================================================================================================
# helper-functions
#===================================================================================================

def fill_ip_tfs_file(tfs_file, column_names, column_datatypes, results):
    tfs_file.add_column_names(column_names)
    tfs_file.add_column_datatypes(column_datatypes)
    for ip_name in 'IP1', 'IP5', 'IP8', 'IP2':
        list_row_entries = ['"' + ip_name + '"']
        if ip_name in results:
            for ip_value in results[ip_name]:
                list_row_entries.append(ip_value)
            
            tfs_file.add_table_row(list_row_entries)


def _get_ip(ip_num, measured, model):
    bpm_left, bpm_right = _get_closest_bpm_names_to_ip(ip_num, model, measured)

    if DEBUG:
        print "_get_ip", ip_num

    if bpm_left is None or bpm_right is None:
        print "skipping ip%1s calculation, no BPM found" % ip_num
        raise ValueError
        
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


def _get_closest_bpm_names_to_ip(ip_num, model, measured):
    bpm_left = None
    bpm_right = None

    for bpm_name in model.NAME:
        if "BPMSW.1L" + ip_num in bpm_name:
            bpm_left = bpm_name
            if _bpm_is_not_in_measured_data(bpm_name, measured):
                bpm_left = None
        if "BPMSW.1R" + ip_num in bpm_name:
            bpm_right = bpm_name
            if _bpm_is_not_in_measured_data(bpm_name, measured):
                bpm_left = None

    return [bpm_left, bpm_right]


def _bpm_is_not_in_measured_data(bpm_name, measured):
    return bpm_name not in measured[0] or 0 >= len([bpm_name][0])


def _get_ip_2(mad_twiss, files, Q, plane, accel, bpms, lhc_phase):
    #-- Common BPMs
    beam_direction = accel.get_beam_direction()

    #-- Loop for IPs
    result = {}
    for ip in ('1', '2', '5', '8'):

        # TODO: the following, maybe this whole function, has to be moved to the accelerator class since it is
        # accelerator dependent calculations

        bpml = 'BPMSW.1L' + ip + '.' + accel.__class__.__name__[-2:]
        bpmr = 'BPMSW.1R' + ip + '.' + accel.__class__.__name__[-2:]
        print bpml

        if (bpml in bpms.index) and (bpmr in bpms.index):

            #-- Model values
            L = 0.5*(mad_twiss.loc[bpmr, "S"] - mad_twiss.loc[bpml, "S"])
            if L < 0: 
                L += 0.5*mad_twiss.LENGTH  #-- For sim starting in the middle of an IP
            if plane == 'H':
                betlmdl = mad_twiss.loc[bpml, "BETX"]
                alflmdl = mad_twiss.loc[bpml, "ALFX"]
            if plane == 'V':
                betlmdl = mad_twiss.loc[bpml, "BETY"]
                alflmdl = mad_twiss.loc[bpml, "ALFY"]
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
                    if bpms.index.get_loc(bpmr) > bpms.index.get_loc(bpml):
                        tune = 0
                    else:
                        tune = Q
                    if plane == 'H':
                        amp_l = t_f.loc[bpml, "AMPX"]
                        amp_r = t_f.loc[bpmr, "AMPX"]
                        dpsi = 2*np.pi*(tune+beam_direction*(t_f.loc[bpmr, "MUX"]-t_f.loc[bpml, "MUX"]))
                    elif plane == 'V':
                        amp_l = t_f.loc[bpml, "AMPY"]
                        amp_r = t_f.loc[bpmr, "AMPY"]
                        dpsi = 2*np.pi*(tune+beam_direction*(t_f.loc[bpmr, "MUY"]-t_f.loc[bpml, "MUY"]))
                    else:
                        raise ValueError("plane is neither 'H' nor 'V'.")

                    #-- To compensate the phase shift by tune
                    if lhc_phase == '1':
                        if (beam_direction==1 and ip=='2') or (beam_direction==-1 and ip=='8'): 
                            dpsi += 2*np.pi*Q

                    #-- bet, alf, and math.sqrt(2J) from amp and phase advance
                    bet = L*(amp_l**2+amp_r**2+2*amp_l*amp_r*cos(dpsi))/(2*amp_l*amp_r*sin(dpsi))
                    alf = (amp_l**2-amp_r**2)/(2*amp_l*amp_r*sin(dpsi))
                    bets = bet/(1+alf**2)
                    d_s = alf*bets
                    sin_dpsi = sin(dpsi)
                    if 0 >sin_dpsi:
                        if DEBUG:
                            print "_get_ip_2: Negative sin_dpsi("+str(sin_dpsi)+") for IP:"+str(ip)
                        continue
                    rt2j = math.sqrt(amp_l*amp_r*sin_dpsi/(2*L))
                    betall.append(bet)
                    alfall.append(alf)
                    betsall.append(bets)
                    dsall.append(d_s)
                    rt2j_all.append(rt2j)
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


def _get_ip_from_phase(MADTwiss, psix, psiy, accelerator):

    raise NotImplementedError("this function hasn't been upgraded to accelerator class")
    IP=('1','2','5','8')
    result={}
    for i in IP:
        bpml = 'BPMSW.1L' + ip + '.' + accelerator.__class__.__name__[-2:]
        bpmr = 'BPMSW.1R' + ip + '.' + accelerator.__class__.__name__[-2:]
        try:
            if psix[bpml][-1]==bpmr:
                #-- Model
                L       =0.5*(MADTwiss.loc[bpmr, "S"] - MADTwiss.loc[bpml, "S"])
                dpsixmdl = MADTwiss.loc[bpmr, "MUX"] - MADTwiss.loc[bpml, "MUX"]
                dpsiymdl = MADTwiss.loc[bpmr, "MUY"] - MADTwiss.loc[bpml, "MUY"]
                betxmdl =MADTwiss.loc[bpml, "BETX"] / (1 + MADTwiss.loc[bpml, "ALFX"] ** 2)
                betymdl =MADTwiss.loc[bpml, "BETY"] / (1 + MADTwiss.loc[bpml, "ALFY"] ** 2)
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
        except KeyError: 
            pass
        except: 
            traceback.print_exc()
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
        
        except KeyError: 
            pass
        except: 
            traceback.print_exc()

    return result


#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def _get_free_ip_2(mad_twiss, mad_ac, ip, plane, accel):

    for i in ('1', '2', '5', '8'):
        bpml = 'BPMSW.1L' + i + '.' + accel.__class__.__name__[-2:]
        bpmr = 'BPMSW.1R' + i + '.' + accel.__class__.__name__[-2:]
        if 'IP'+i in ip:

            L = 0.5*(mad_twiss.S[mad_twiss.indx[bpmr]]-mad_twiss.S[mad_twiss.indx[bpml]])
            if L < 0:
                L += 0.5*mad_twiss.LENGTH
            #-- bet and alf at the left BPM
            if plane == 'H':
                betl = mad_twiss.BETX[mad_twiss.indx[bpml]]
                betdl = mad_ac.BETX[mad_ac.indx[bpml]]
                alfl = mad_twiss.ALFX[mad_twiss.indx[bpml]]
                alfdl = mad_ac.ALFX[mad_ac.indx[bpml]]
            if plane == 'V':
                betl = mad_twiss.BETY[mad_twiss.indx[bpml]]
                betdl = mad_ac.BETY[mad_ac.indx[bpml]]
                alfl = mad_twiss.ALFY[mad_twiss.indx[bpml]]
                alfdl = mad_ac.ALFY[mad_ac.indx[bpml]]
            #-- IP parameters propagated from the left BPM
            bets = betl/(1+alfl**2)       
            betds = betdl/(1+alfdl**2)
            bet = betl-2*alfl*L+L**2/bets
            betd = betdl-2*alfdl*L+L**2/betds
            alf = alfl-L/bets
            alfd = alfdl-L/betds
            ds = alf*bets
            dsd = alfd*betds
            #-- Apply corrections
            ip['IP'+i][0] = ip['IP'+i][0]+bet-betd
            ip['IP'+i][2] = bet
            ip['IP'+i][3] = ip['IP'+i][3]+alf-alfd
            ip['IP'+i][5] = alf
            ip['IP'+i][6] = ip['IP'+i][6]+bets-betds
            ip['IP'+i][8] = bets
            ip['IP'+i][9] = ip['IP'+i][9]+ds-dsd
            ip['IP'+i][11] = ds

    return ip


