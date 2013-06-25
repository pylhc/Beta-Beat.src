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
import math

import numpy as np
from numpy import sin, cos, tan

import Utilities.bpm


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


#===================================================================================================
# helper-functions
#===================================================================================================

def BetaFromPhase(MADTwiss,ListOfFiles,phase,plane):

    alfa={}
    beta={}
    commonbpms = Utilities.bpm.intersect(ListOfFiles)
    commonbpms = Utilities.bpm.model_intersect(commonbpms,MADTwiss)
    alfii=[]
    betii=[]
    delbeta=[]
    for i in range(0,len(commonbpms)):
        if i==len(commonbpms)-2: # The last monitor but one
            bn1=str.upper(commonbpms[i][1])
            bn2=str.upper(commonbpms[i+1][1])
            bn3=str.upper(commonbpms[0][1])
        elif i==len(commonbpms)-1: # The last monitor
            bn1=str.upper(commonbpms[i][1])
            bn2=str.upper(commonbpms[0][1])
            bn3=str.upper(commonbpms[1][1])
        else : # Others
            bn1=str.upper(commonbpms[i][1])
            bn2=str.upper(commonbpms[i+1][1])
            bn3=str.upper(commonbpms[i+2][1])
        ph2pi12=2.*np.pi*phase[bn1][0]
        ph2pi23=2.*np.pi*phase[bn2][0]
        ph2pi13=2.*np.pi*phase[bn1][2]
        # Find the model transfer matrices for beta1
        phmdl12=2.*np.pi*phase[bn1][4]
        phmdl13=2.*np.pi*phase[bn1][5]
        phmdl23=2.*np.pi*phase[bn2][4]
        if plane=='H':
            betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
            betmdl2=MADTwiss.BETX[MADTwiss.indx[bn2]]
            betmdl3=MADTwiss.BETX[MADTwiss.indx[bn3]]
            alpmdl1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
            alpmdl2=MADTwiss.ALFX[MADTwiss.indx[bn2]]
            alpmdl3=MADTwiss.ALFX[MADTwiss.indx[bn3]]
        elif plane=='V':
            betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
            betmdl2=MADTwiss.BETY[MADTwiss.indx[bn2]]
            betmdl3=MADTwiss.BETY[MADTwiss.indx[bn3]]
            alpmdl1=MADTwiss.ALFY[MADTwiss.indx[bn1]]
            alpmdl2=MADTwiss.ALFY[MADTwiss.indx[bn2]]
            alpmdl3=MADTwiss.ALFY[MADTwiss.indx[bn3]]
        if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
            print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
            sys.exit(1)

        # Find beta1 and alpha1 from phases assuming model transfer matrix
        # Matrix M: BPM1-> BPM2
        # Matrix N: BPM1-> BPM3
        M11=math.sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        N11=math.sqrt(betmdl3/betmdl1)*(cos(phmdl13)+alpmdl1*sin(phmdl13))
        N12=math.sqrt(betmdl1*betmdl3)*sin(phmdl13)

        denom=M11/M12-N11/N12+1e-16
        numer=1/tan(ph2pi12)-1/tan(ph2pi13)
        beti1=numer/denom
        # beterr1=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # beterr1=beterr1+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # beterr1=beterr1/abs(denom)
        betstd1=        (2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        betstd1=betstd1+(2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        betstd1=math.sqrt(betstd1)/abs(denom)

        denom=M12/M11-N12/N11+1e-16
        numer=-M12/M11/tan(ph2pi12)+N12/N11/tan(ph2pi13)
        alfi1=numer/denom
        # alferr1=abs(M12/M11*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # alferr1=alferr1+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # alferr1=alferr1/abs(denom)
        alfstd1=        (M12/M11*2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        alfstd1=alfstd1+(N12/N11*2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        alfstd1=math.sqrt(alfstd1)/denom

        # Find beta2 and alpha2 from phases assuming model transfer matrix
        # Matrix M: BPM1-> BPM2
        # Matrix N: BPM2-> BPM3
        M22=math.sqrt(betmdl1/betmdl2)*(cos(phmdl12)-alpmdl2*sin(phmdl12))
        M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        N11=math.sqrt(betmdl3/betmdl2)*(cos(phmdl23)+alpmdl2*sin(phmdl23))
        N12=math.sqrt(betmdl2*betmdl3)*sin(phmdl23)

        denom=M22/M12+N11/N12+1e-16
        numer=1/tan(ph2pi12)+1/tan(ph2pi23)
        beti2=numer/denom
        # beterr2=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # beterr2=beterr2+abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # beterr2=beterr2/abs(denom)
        betstd2=        (2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        betstd2=betstd2+(2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        betstd2=math.sqrt(betstd2)/abs(denom)

        denom=M12/M22+N12/N11+1e-16
        numer=M12/M22/tan(ph2pi12)-N12/N11/tan(ph2pi23)
        alfi2=numer/denom
        # alferr2=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # alferr2=alferr2+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # alferr2=alferr2/abs(denom)
        alfstd2=        (M12/M22*2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        alfstd2=alfstd2+(N12/N11*2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        alfstd2=math.sqrt(alfstd2)/abs(denom)

        # Find beta3 and alpha3 from phases assuming model transfer matrix
        # Matrix M: BPM2-> BPM3
        # Matrix N: BPM1-> BPM3
        M22=math.sqrt(betmdl2/betmdl3)*(cos(phmdl23)-alpmdl3*sin(phmdl23))
        M12=math.sqrt(betmdl2*betmdl3)*sin(phmdl23)
        N22=math.sqrt(betmdl1/betmdl3)*(cos(phmdl13)-alpmdl3*sin(phmdl13))
        N12=math.sqrt(betmdl1*betmdl3)*sin(phmdl13)

        denom=M22/M12-N22/N12+1e-16
        numer=1/tan(ph2pi23)-1/tan(ph2pi13)
        beti3=numer/denom
        # beterr3=abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # beterr3=beterr3+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # beterr3=beterr3/abs(denom)
        betstd3=        (2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        betstd3=betstd3+(2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        betstd3=math.sqrt(betstd3)/abs(denom)

        denom=M12/M22-N12/N22+1e-16
        numer=M12/M22/tan(ph2pi23)-N12/N22/tan(ph2pi13)
        alfi3=numer/denom
        # alferr3=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # alferr3=alferr3+abs(N22/N12*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # alferr3=alferr3/abs(denom)
        alfstd3=        (M12/M22*2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        alfstd3=alfstd3+(N12/N22*2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        alfstd3=math.sqrt(alfstd3)/abs(denom)

        betii.append([beti1,betstd1,beti2,betstd2,beti3,betstd3])
        alfii.append([alfi1,alfstd1,alfi2,alfstd2,alfi3,alfstd3])

    # Output the beta at all monitors even for the first, second, last and last but one using beta1,2 and 3!
    for i in range(0,len(commonbpms)):
        bn1=str.upper(commonbpms[i][1])
        if i==0: # The first monitor
            ib1=0
            ib2=len(commonbpms)-1
            ib3=len(commonbpms)-2

        elif i==1: # The second monitor
            ib1=i
            ib2=0
            ib3=len(commonbpms)-1
        else: # Others
            ib1=i
            ib2=i-1
            ib3=i-2
        # Averaging over beta/alpha1,2 and 3. Find beta/alpha errors
        beti=(betii[ib1][0]+betii[ib2][2]+betii[ib3][4])/3.
        betstd=math.sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/math.sqrt(3.)
        # betstd=math.sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/math.sqrt(3.*len(ListOfFiles))  # If we want to plot std "over files" !!!!!
        try:
            beterr=math.sqrt((betii[ib1][0]**2.+betii[ib2][2]**2.+betii[ib3][4]**2.)/3.-beti**2.)
        except:
            beterr=0
        alfi=(alfii[ib1][0]+alfii[ib2][2]+alfii[ib3][4])/3.
        alfstd=math.sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/math.sqrt(3.)
        # alfstd=math.sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/math.sqrt(3.*len(ListOfFiles))  # If we want to plot std "over files" !!!!!
        try:
            alferr=math.sqrt((alfii[ib1][0]**2.+alfii[ib2][2]**2.+alfii[ib3][4]**2.)/3.-alfi**2.)
        except:
            alferr=0

        beta[bn1]=(beti,beterr,betstd)
        alfa[bn1]=(alfi,alferr,alfstd)
        if plane=='H':
            betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
        delbeta.append((beti-betmdl1)/betmdl1)


    delbeta=np.array(delbeta)
    rmsbb=math.sqrt(np.average(delbeta*delbeta))

    return [beta,rmsbb,alfa,commonbpms]


def BetaFromAmplitude(MADTwiss,ListOfFiles,plane):

    beta={}
    root2J=[]
    commonbpms=Utilities.bpm.intersect(ListOfFiles)
    commonbpms=Utilities.bpm.model_intersect(commonbpms,MADTwiss)
    SumA=0.0
    Amp=[]
    Amp2=[]
    Kick2=[]
    for i in range(0,len(commonbpms)): # this loop have become complicated after modifications... anybody simplify?
        bn1=str.upper(commonbpms[i][1])
        if plane=='H':
            tembeta=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            tembeta=MADTwiss.BETY[MADTwiss.indx[bn1]]
        Ampi=0.0
        Ampj2=[]
        root2Ji=0.0
        jj=0
        for j in ListOfFiles:
            if i==0:
                Kick2.append(0)
            if plane=='H':
                Ampi+=j.AMPX[j.indx[bn1]]
                Ampj2.append(j.AMPX[j.indx[bn1]]**2)
                root2Ji+=j.PK2PK[j.indx[bn1]]/2.

            elif plane=='V':
                Ampi+=j.AMPY[j.indx[bn1]]
                Ampj2.append(j.AMPY[j.indx[bn1]]**2)
                root2Ji+=j.PK2PK[j.indx[bn1]]/2.
            Kick2[jj]+=Ampj2[jj]/tembeta
            jj+=1
        Ampi=Ampi/len(ListOfFiles)
        root2Ji=root2Ji/len(ListOfFiles)
        Amp.append(Ampi)
        Amp2.append(Ampj2)



        SumA+=Ampi**2/tembeta
        root2J.append(root2Ji/math.sqrt(tembeta))


    Kick=SumA/len(commonbpms) # Assuming the average of beta is constant
    Kick2=np.array(Kick2)
    Kick2=Kick2/len(commonbpms)
    Amp2=np.array(Amp2)
    root2J=np.array(root2J)
    root2Jave=np.average(root2J)
    root2Jrms=math.sqrt(np.average(root2J*root2J)-root2Jave**2+2.2e-16)

    delbeta=[]
    for i in range(0,len(commonbpms)):
        bn1=str.upper(commonbpms[i][1])
        location=commonbpms[i][0]
        for j in range(0,len(ListOfFiles)):
            Amp2[i][j]=Amp2[i][j]/Kick2[j]
        #print np.average(Amp2[i]*Amp2[i]),np.average(Amp2[i])**2
        try:
            betstd=math.sqrt(np.average(Amp2[i]*Amp2[i])-np.average(Amp2[i])**2+2.2e-16)
        except:
            betstd=0

        beta[bn1]=[Amp[i]**2/Kick,betstd,location]
        if plane=='H':
            betmdl=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            betmdl=MADTwiss.BETY[MADTwiss.indx[bn1]]
        delbeta.append((beta[bn1][0]-betmdl)/betmdl)

    invariantJ=[root2Jave,root2Jrms]

    delbeta=np.array(delbeta)
    rmsbb=math.sqrt(np.average(delbeta*delbeta))
    return [beta,rmsbb,commonbpms,invariantJ]

#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def getFreeBeta(modelfree,modelac,betal,rmsbb,alfal,bpms,plane): # to check "+"
    if DEBUG:
        print "Calculating free beta using model"
    bpms=Utilities.bpm.model_intersect(bpms, modelfree)
    bpms=Utilities.bpm.model_intersect(bpms, modelac)
    betan={}
    alfan={}
    for bpma in bpms:

        bpm=bpma[1].upper()
        beta,beterr,betstd=betal[bpm]
        alfa,alferr,alfstd=alfal[bpm]

        if plane=="H":
            betmf=modelfree.BETX[modelfree.indx[bpm]]
            betma=modelac.BETX[modelac.indx[bpm]]
            bb=(betma-betmf)/betmf
            alfmf=modelfree.ALFX[modelfree.indx[bpm]]
            alfma=modelac.ALFX[modelac.indx[bpm]]
            aa=(alfma-alfmf)/alfmf
        else:
            betmf=modelfree.BETY[modelfree.indx[bpm]]
            betma=modelac.BETY[modelac.indx[bpm]]
            alfmf=modelfree.ALFY[modelfree.indx[bpm]]
            alfma=modelac.ALFY[modelac.indx[bpm]]
            bb=(betma-betmf)/betmf
            aa=(alfma-alfmf)/alfmf

        betan[bpm]=beta*(1+bb),beterr,betstd # has to be plus!
        alfan[bpm]=alfa*(1+aa),alferr,alfstd

    return betan,rmsbb,alfan,bpms


def getFreeAmpBeta(betai,rmsbb,bpms,invJ,modelac,modelfree,plane): # "-"

    #
    # Why difference in betabeta calculation ??
    #
    #

    betas={}

    if DEBUG:
        print "Calculating free beta from amplitude using model"

    for bpm in bpms:

        bpmm=bpm[1].upper()
        beta=betai[bpmm][0]

        if plane=="H":
            betmf=modelfree.BETX[modelfree.indx[bpmm]]
            betma=modelac.BETX[modelac.indx[bpmm]]
            bb=(betmf-betma)/betma

        else:
            betmf=modelfree.BETY[modelfree.indx[bpmm]]
            betma=modelac.BETY[modelac.indx[bpmm]]
            bb=(betmf-betma)/betma

        betas[bpmm]=[beta*(1+bb),betai[bpmm][1],betai[bpmm][2]]

    return betas,rmsbb,bpms,invJ

