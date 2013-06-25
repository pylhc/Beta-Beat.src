'''
Created on 27 May 2013

@author: vimaier

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
# helper-functions
#===================================================================================================

def NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, betax, cut_co):

    nzdpp=len(ListOfNonZeroDPPX) # How many non zero dpp
    zdpp=len(ListOfZeroDPPX)  # How many zero dpp
    if zdpp==0 or nzdpp ==0:
        print >> sys.stderr, 'Warning: No data for dp/p!=0 or even for dp/p=0.'
        print >> sys.stderr, 'No output for the  dispersion.'
        dum0={}
        dum1=[]
        return [dum0, dum1]


    coac=ListOfCOX[0] # COX dictionary after cut bad BPMs
    coact={}
    for i in coac:
        if (coac[i][1] < cut_co):
            coact[i]=coac[i]
    coac=coact


    nda={} # Dictionary for the output containing [average Disp, rms error]

    ALL=ListOfZeroDPPX+ListOfNonZeroDPPX
    commonbpmsALL=Utilities.bpm.intersect(ALL)
    commonbpmsALL=Utilities.bpm.model_intersect(commonbpmsALL, MADTwiss)

    mydp=[]
    gf=[]
    for j in ListOfNonZeroDPPX:
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
        ndmdli=MADTwiss.DX[MADTwiss.indx[bn1]]/math.sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
        ndmdl.append(ndmdli)

        try:
            coi=coac[bn1]

            ampi=0.0
            for j in ListOfZeroDPPX:
                ampi+=j.AMPX[j.indx[bn1]]
            ampi=ampi/zdpp

            ndi=[]
            for j in range(0,nzdpp): # the range(0,nzdpp) instead of ListOfNonZeroDPPX is used because the index j is used in the loop
                codpp=ListOfCOX[j+1]
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
            for j in range(0,nzdpp): # the range(0,nzdpp) instead of ListOfZeroDPPX is used because the index j is used in the loop
                ndi.append(nd[i-badco][j]/gf[j])
            ndi=np.array(ndi)
            ndstd=math.sqrt(np.average(ndi*ndi)-(np.average(ndi))**2.0+2.2e-16)
            ndas=np.average(wf*ndi)
            nda[bn1]=[ndas,ndstd]
            Dx[bn1]=[nda[bn1][0]*math.sqrt(betax[bn1][0]),dummy]
            bpms.append([bns1,bn1])
        except:
            badco+=1
    DPX=GetDPX(MADTwiss, Dx, bpms)

    return [nda,Dx,DPX,bpms]


#---------- DPX, New!
def GetDPX(MADTwiss,Dx,commonbpms):

    DPX={}
    for i in range(0,len(commonbpms)):
        bn1=str.upper(commonbpms[i][1])
        if i==len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2=str.upper(commonbpms[0][1])
            phmdl12=2.*np.pi*(MADTwiss.Q1+MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
        else:
            bn2=str.upper(commonbpms[i+1][1])
            phmdl12=2.*np.pi*(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
        betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
        betmdl2=MADTwiss.BETX[MADTwiss.indx[bn2]]
        alpmdl1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
        dxmdl1=MADTwiss.DX[MADTwiss.indx[bn1]]
        dpxmdl1=MADTwiss.DPX[MADTwiss.indx[bn1]]
        dxmdl2=MADTwiss.DX[MADTwiss.indx[bn2]]

        M11=math.sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        M13=dxmdl2-M11*dxmdl1-M12*dpxmdl1
        # use the beta from amplitude
        DPX[bn1]=(-M13+Dx[bn2][0]-M11*Dx[bn1][0])/M12

    return DPX



#----------- DPY, New!
def GetDPY(MADTwiss,Dy,commonbpms):
    DPY={}
    for i in range(0,len(commonbpms)):
        bn1=str.upper(commonbpms[i][1])
        if i==len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2=str.upper(commonbpms[0][1])
            phmdl12=2.*np.pi*(MADTwiss.Q2+MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
        else:
            bn2=str.upper(commonbpms[i+1][1])
            phmdl12=2.*np.pi*(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
        betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
        betmdl2=MADTwiss.BETY[MADTwiss.indx[bn2]]
        alpmdl1=MADTwiss.ALFY[MADTwiss.indx[bn1]]
        dymdl1=MADTwiss.DY[MADTwiss.indx[bn1]]
        dpymdl1=MADTwiss.DPY[MADTwiss.indx[bn1]]
        dymdl2=MADTwiss.DY[MADTwiss.indx[bn2]]

        M11=math.sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        M13=dymdl2-M11*dymdl1-M12*dpymdl1
        #M13=-M11*dymdl1-M12*dpymdl1
        # use the beta from amplitude
        DPY[bn1]=(-M13+Dy[bn2][0]-M11*Dy[bn1][0])/M12
        #DPY[bn1]=(-M13-M11*Dy[bn1][0])/M12

    return DPY


def DispersionfromOrbit(ListOfZeroDPP,ListOfNonZeroDPP,ListOfCO,COcut,BPMU):

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

#===================================================================================================
# ac-dipole stuff
#===================================================================================================

