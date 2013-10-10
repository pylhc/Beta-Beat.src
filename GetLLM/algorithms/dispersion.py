'''
Created on 27 May 2013

@author: ?, vimaier

@version: 0.0.1

GetLLM.algorithms.dispersion.py stores helper functions for dispersion calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import math

import numpy as np
from numpy import sin, cos

import Utilities.bpm


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


#===================================================================================================
# main part
#===================================================================================================
def calculate_dispersion(getllm_d, twiss_d, tune_d, mad_twiss, files_dict, beta_x_from_amp, list_of_co_x, list_of_co_y):
    '''
    Calculates dispersion and fills the following TfsFiles:
        getNDx.out        getDx.out        getDy.out

    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            accel is used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances
    '''
    print 'Calculating dispersion'
    if twiss_d.has_zero_dpp_x() and twiss_d.has_non_zero_dpp_x():
        [nda, d_x, dpx, bpms] = _norm_disp_x(mad_twiss, twiss_d.zero_dpp_x, twiss_d.non_zero_dpp_x, list_of_co_x, beta_x_from_amp, getllm_d.cut_for_closed_orbit)
        tfs_file = files_dict['getNDx.out']
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "NDX", "STDNDX", "DX", "DPX", "NDXMDL", "DXMDL", "DPXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            ndmdl = mad_twiss.DX[mad_twiss.indx[bn1]] / math.sqrt(mad_twiss.BETX[mad_twiss.indx[bn1]])
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.non_zero_dpp_x), nda[bn1][0], nda[bn1][1], d_x[bn1][0], dpx[bn1], ndmdl, mad_twiss.DX[mad_twiss.indx[bn1]], mad_twiss.DPX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

        [dxo, bpms] = _dispersion_from_orbit(twiss_d.zero_dpp_x, twiss_d.non_zero_dpp_x, list_of_co_x, getllm_d.cut_for_closed_orbit, getllm_d.bpm_unit)
        dpx = _get_dpx(mad_twiss, dxo, bpms)
        tfs_file = files_dict['getDx.out']
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "DX", "STDDX", "DPX", "DXMDL", "DPXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.non_zero_dpp_x), dxo[bn1][0], dxo[bn1][1], dpx[bn1], mad_twiss.DX[mad_twiss.indx[bn1]], mad_twiss.DPX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

    if twiss_d.has_zero_dpp_y() and twiss_d.has_non_zero_dpp_y():
        [dyo, bpms] = _dispersion_from_orbit(twiss_d.zero_dpp_y, twiss_d.non_zero_dpp_y, list_of_co_y, getllm_d.cut_for_closed_orbit, getllm_d.bpm_unit)
        dpy = _get_dpy(mad_twiss, dyo, bpms)
        tfs_file = files_dict['getDy.out']
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "DY", "STDDY", "DPY", "DYMDL", "DPYMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.non_zero_dpp_y), dyo[bn1][0], dyo[bn1][1], dpy[bn1], mad_twiss.DY[mad_twiss.indx[bn1]], mad_twiss.DPY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
# END calculate_dispersion -------------------------------------------------------------------------


#===================================================================================================
# helper-functions
#===================================================================================================

def _norm_disp_x(mad_twiss, list_zero_dpp_x, list_non_zero_dpp_x, list_co_x, betax, cut_co):

    nzdpp=len(list_non_zero_dpp_x) # How many non zero dpp
    zdpp=len(list_zero_dpp_x)  # How many zero dpp
    if zdpp==0 or nzdpp ==0:
        print >> sys.stderr, 'Warning: No data for dp/p!=0 or even for dp/p=0.'
        print >> sys.stderr, 'No output for the  dispersion.'
        dum0={}
        dum1=[]
        return [dum0, dum1]


    coac=list_co_x[0] # COX dictionary after cut bad BPMs
    coact={}
    for i in coac:
        if (coac[i][1] < cut_co):
            coact[i]=coac[i]
    coac=coact


    nda={} # Dictionary for the output containing [average Disp, rms error]

    ALL=list_zero_dpp_x+list_non_zero_dpp_x
    commonbpmsALL=Utilities.bpm.intersect(ALL)
    commonbpmsALL=Utilities.bpm.model_intersect(commonbpmsALL, mad_twiss)

    mydp=[]
    gf=[]
    for j in list_non_zero_dpp_x:
        mydp.append(float(j.DPP))
        gf.append(0.0) # just to initialize the array, may not be a clever way...
    mydp=np.array(mydp)
    wf=np.array(abs(mydp))/sum(abs(mydp))*len(mydp) #Weight for the average depending on DPP


    # Find the global factor
    nd=[]
    ndmdl=[]
    badco=0
    for i in range(0,len(commonbpmsALL)):
        bn1=str.upper(commonbpmsALL[i][1])
        bns1=commonbpmsALL[i][0]
        ndmdli=mad_twiss.DX[mad_twiss.indx[bn1]]/math.sqrt(mad_twiss.BETX[mad_twiss.indx[bn1]])
        ndmdl.append(ndmdli)

        try:
            coi=coac[bn1]

            ampi=0.0
            for j in list_zero_dpp_x:
                ampi+=j.AMPX[j.indx[bn1]]
            ampi=ampi/zdpp

            ndi=[]
            for j in range(0,nzdpp): # the range(0,nzdpp) instead of list_non_zero_dpp_x is used because the index j is used in the loop
                codpp=list_co_x[j+1]
                orbit=codpp[bn1][0]-coi[0]
                ndm=orbit/ampi
                gf[j]+=ndm
                ndi.append(ndm)
            nd.append(ndi)
        except:
            badco+=1
            coi=0

    ndmdl=np.array(ndmdl)
    avemdl=np.average(ndmdl)

    gf=np.array(gf)
    gf=gf/avemdl/(len(commonbpmsALL)-badco)


    # Find normalized dispersion and Dx construction
    nd=np.array(nd)
    Dx={}
    dummy=0.0 # dummy for DXSTD
    bpms=[]
    badco=0
    for i in range(0,len(commonbpmsALL)):
        ndi=[]
        bn1=str.upper(commonbpmsALL[i][1])
        bns1=commonbpmsALL[i][0]
        try:
            coac[bn1]
            for j in range(0,nzdpp): # the range(0,nzdpp) instead of list_zero_dpp_x is used because the index j is used in the loop
                ndi.append(nd[i-badco][j]/gf[j])
            ndi=np.array(ndi)
            ndstd=math.sqrt(np.average(ndi*ndi)-(np.average(ndi))**2.0+2.2e-16)
            ndas=np.average(wf*ndi)
            nda[bn1]=[ndas,ndstd]
            Dx[bn1]=[nda[bn1][0]*math.sqrt(betax[bn1][0]),dummy]
            bpms.append([bns1,bn1])
        except:
            badco+=1
    dpx=_get_dpx(mad_twiss, Dx, bpms)

    return [nda,Dx,dpx,bpms]


def _get_dpx(mad_twiss, d_x, commonbpms):

    dpx = {}
    for i in range(0,len(commonbpms)):
        bn1 = str.upper(commonbpms[i][1])
        if i == len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2 = str.upper(commonbpms[0][1])
            # phase model between BPM1 and BPM2
            p_mdl_12 = 2.*np.pi*(mad_twiss.Q1+mad_twiss.MUX[mad_twiss.indx[bn2]]-mad_twiss.MUX[mad_twiss.indx[bn1]])
        else:
            bn2 = str.upper(commonbpms[i+1][1])
            # phase model between BPM1 and BPM2
            p_mdl_12 = 2.*np.pi*(mad_twiss.MUX[mad_twiss.indx[bn2]]-mad_twiss.MUX[mad_twiss.indx[bn1]])
        betmdl1 = mad_twiss.BETX[mad_twiss.indx[bn1]]
        betmdl2 = mad_twiss.BETX[mad_twiss.indx[bn2]]
        alpmdl1 = mad_twiss.ALFX[mad_twiss.indx[bn1]]
        dxmdl1 = mad_twiss.DX[mad_twiss.indx[bn1]]
        dpxmdl1 = mad_twiss.DPX[mad_twiss.indx[bn1]]
        dxmdl2 = mad_twiss.DX[mad_twiss.indx[bn2]]

        M11 = math.sqrt(betmdl2/betmdl1)*(cos(p_mdl_12)+alpmdl1*sin(p_mdl_12))
        M12 = math.sqrt(betmdl1*betmdl2)*sin(p_mdl_12)
        M13 = dxmdl2-M11*dxmdl1-M12*dpxmdl1
        # use the beta from amplitude
        dpx[bn1] = (-M13+d_x[bn2][0]-M11*d_x[bn1][0])/M12

    return dpx


def _get_dpy(mad_twiss, d_y, commonbpms):
    dpy={}
    for i in range(0,len(commonbpms)):
        bn1 = str.upper(commonbpms[i][1])
        if i == len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2 = str.upper(commonbpms[0][1])
            # phase model between BPM1 and BPM2
            phmdl12 = 2.*np.pi*(mad_twiss.Q2+mad_twiss.MUY[mad_twiss.indx[bn2]]-mad_twiss.MUY[mad_twiss.indx[bn1]])
        else:
            bn2 = str.upper(commonbpms[i+1][1])
            phmdl12 = 2.*np.pi*(mad_twiss.MUY[mad_twiss.indx[bn2]]-mad_twiss.MUY[mad_twiss.indx[bn1]])
        betmdl1 = mad_twiss.BETY[mad_twiss.indx[bn1]]
        betmdl2 = mad_twiss.BETY[mad_twiss.indx[bn2]]
        alpmdl1 = mad_twiss.ALFY[mad_twiss.indx[bn1]]
        dymdl1 = mad_twiss.DY[mad_twiss.indx[bn1]]
        dpymdl1 = mad_twiss.DPY[mad_twiss.indx[bn1]]
        dymdl2 = mad_twiss.DY[mad_twiss.indx[bn2]]

        M11 = math.sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12 = math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        M13 = dymdl2-M11*dymdl1-M12*dpymdl1
        #M13 = -M11*dymdl1-M12*dpymdl1
        # use the beta from amplitude
        dpy[bn1] = (-M13+d_y[bn2][0]-M11*d_y[bn1][0])/M12
        #dpy[bn1] = (-M13-M11*d_y[bn1][0])/M12

    return dpy


def _dispersion_from_orbit(ListOfZeroDPP,ListOfNonZeroDPP,ListOfCO,COcut,BPMU):

    if BPMU=='um': scalef=1.0e-6
    elif BPMU=='mm': scalef=1.0e-3
    elif BPMU=='cm': scalef=1.0e-2
    elif BPMU=='m': scalef=1.0


    coac=ListOfCO[0]
    coact={}
    for i in coac:
        if (coac[i][1] < COcut):
            coact[i]=coac[i]

    coac=coact # COY dictionary after cut bad BPMs

    ALL=ListOfZeroDPP+ListOfNonZeroDPP
    commonbpmsALL=Utilities.bpm.intersect(ALL)



    mydp=[]
    for j in ListOfNonZeroDPP:
        mydp.append(float(j.DPP))
    mydp=np.array(mydp)
    wf=np.array(abs(mydp))/sum(abs(mydp))*len(mydp) #Weitghs for the average


    dco={} # Dictionary for the output containing [(averaged)Disp, rms error]
    bpms=[]
    for i in range(0,len(commonbpmsALL)):
        bn1=str.upper(commonbpmsALL[i][1])
        bns1=commonbpmsALL[i][0]

        try:
            coi=coac[bn1]
            dcoi=[]
            for j in ListOfNonZeroDPP:
                dcoi.append((j.CO[j.indx[bn1]]-coi[0])*scalef/float(j.DPP))
            dcoi=np.array(dcoi)
            dcostd=math.sqrt(np.average(dcoi*dcoi)-(np.average(dcoi))**2.0+2.2e-16)
            dcos=np.average(wf*dcoi)
            dco[bn1]=[dcos,dcostd]
            bpms.append([bns1,bn1])
        except:
            coi=0
    return [dco,bpms]

