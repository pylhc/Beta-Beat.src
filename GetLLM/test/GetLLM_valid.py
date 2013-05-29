# Python script to obtain Linear Lattice functions and More -> GetLLM
# Version-up history:
# V1.0, 11/Feb/2008 by Masa. Aiba
# V1.1, 18/Feb/2008:
# - Debugging, add model phase and tunes to output
# - add function to obtain DY
# - add chromatic parameter (phase for non zero DPP)
# V1.2, 22/Feb/2008:
# - test version for beta with all BPM
# V1.3, 29/Feb/2008:
# - beta from phases is improved, averaging beta1, 2 and 3
# V1.31, 12/Mar/2008:
# - debugged on alpha3
# V1.4, 12/Mar/2008:
# - modify output to fit latest TSF format and to meet requests from Rogelio
# - fix buggs in r.m.s. beta-beat and added to the output of getbetax/y.out
# V1.5 Rogelio, 13 March 2008:
# - Update to option parser, include BPMdictionary to filter BPMs not in Model
# V1.51, 13/Mar/2008:
# - Modify output to fit latest TSF format again. Add STD to beta.
# V1.6, 15/Jul/2008:
# - Add the integer part of tunes - assuming that the phase advance is always
#   less than 1.0.
# V1.71 27/Jul/2008:
# - Add GetCO. Filter in dispersion calculation to exclude bad bpms.
# V1.8, 13/Aug/2008 Ref. note by A. Franchi, R. T. Garcia, G. Vanbavinckhove:
# - Add GetCoupling
# - "Computation of the Coupling Resonance Driving term f1001 and the coupling
#   coefficient C from turn-by-turn single-BPM data", 28/May/2008
# - The GetCoupling.py is initiated by Glenn V. and imported into GetLLM /
#   finalized by Masa. Aiba
# - Some bugs are fixed - you can run even without liny file and find the
#   results for horizontal plane only.
# V1.81, 1/Sep/2008:
# - For an accelerator in which the beam goes opposite direction to the model
#   as in LHCB2, the beam direction parameter (bd) is added to GetPhases.
# - Bug in the phi13 for the last monitor and the last but one is fixed.
# V1.9, 21/Oct/2008:
# - Add the beta from spectrum height.
# V1.91, 08/Dec/2008:
# - Add option - SUSSIX or SVD for input file
# V1.99, 13/Feb/09:
# - Add DPX! The MEMORIAL version for Version 1.** since all the linear lattice
#   parameters become available!
# - Add the option for the Harmonic analysis.
# - Refine coding, add several lines of comment
# V2.0, 17/Feb/2009:
# - Major version up to output "More", that is, other than the linear lattice
#   parameters.
# - Add off-momentum lattice (dbeta/beta)/(dp/p)
# - Add SPS coupling with "Pseudo-double plane BPM"-it need at this moment a
#   list containing
# - pre-paired H-V monitor. This should be replaced by a clever algorithm to
#   find good pairs automatically.
# V2.01, 10/Mar/2009:
# - Fix bug on SPS double plane BPM monitor, in which missing BPM could cause
#   an error.
# - Modify BetaFromAmplitude to output invariant J (Rogelio / finalised by MA)
# V2.02, 10/Mar/2009:
# - Fix bug in getcoupling, bad input for phasex and phasey
# V2.10, 13/Mar/2009, Glenn Vanbavinckhove:
# - Added function for finding sextupole lines (amp and phases) + chiterms amp
# V2.11. 26/Mar/2009
# - Fix bug in Getphase (model phase advance from the last to first monitor).
# V2.12, 28/03/2009
# - Following bugs fixed : avoid negative square, change option -h to -l
# - Add -r option as in correct.py
# - Change the way to import BPM pair file for SPS. (import -> execfile)
# V2.13, 06/Apl/2009
# - Fix bug in weight function to accept negative dpp
# V2.14, 08/Apl/2009
# - Fix bug in Normalized dispersion to treat COcut correctly.
# V2.15:
# - Enable coupling in RHIC as in SPS
# V2.16, 28/May/2009:
# - Add STDBET for the beta from amplitude.
# - Add option for off momentum beta-beating to choose algorithm, that is, beta
#   from phase or amp
# - Add a routine to detect wrong data having two lines in linx/y file with the
#   same BPM name.
# - Add a routine to avoid zero division due to exactly n*pi phase advance in
#   beta from phase (see the last part of GetPhases).
# V2.21, 23/June/2009:
# - Add STDBET Model for off momentum beta beat phase.
# V2.25:
# - Fixed bd flag (must be -1 for beam2)
# V2.25:
# - Adding VERSION variable to be always output and modified in subsequent
#   versions, do not forget!!!
# V2.26:
# - Adding F2000 (two different methods linear and non-linear)
# V2.27:
# - Adding new method for IP calculation
# V2.28:
# - Changing the rejection of bad BPM for the coupling phase - averaging the
#   phase over the sets of data first , then cut if the q1 and q2 are very
#   different.
# - 24/Feb/2010 the change is not yet checked.
# - Hyphens in the @ field of tfs files is not allowed:  The previous label
#   "RMS-beta-beat" has been  moved to "RMSbetabeat"
# V2.29 5/Mar/2010
# - Change the default value for COcut from 1000 to 4000 as it was too small
# V2.30 13/Sept/2010:
# - Updating for AC-Dipole, implementing chromatic coupling (not done for RHIC
#   or SPS)
# V2.31 8/November/2010:
# - Update for AC-Dipole: gives free beta,phase,coupling(global factor)
# V2.32 15/January/2010:
# - Taking models
# V2.33 7/February/2011:
# - implementing coupling correction for AC-dipole.
# V2.34 7/04/2011:
# - Updating to deal with chromatic twiss files.
# V2.35 9/06/2011:
# - Functions to cancel the AC dipole effect for beta, phase, and total phase,
#   based on equations, are added.
# - A function to calculate beta* from the phase advance between Q1s is added.
# - Phase shift by tune is compensated for the LHC experiment data.
# V2.36 30/09/2011:
# - Rescaling algorithm for BetaFromAmp and action (by Andy) is implemented.
# - 2nd function to calculate IP parameters from Q1s is added.
# - The compensation of the phase shift by tune for LHC exp data is modified.
# - A function to calculate action of the AC dipole excitation is added.
# V2.38 08/03/2012:
# - added main() function
# - using raise ValueError() instead of sys.exit() some places
# 13/09/2012:
# - merged in patch 2.35-2.37
# V2.38b 03/dec/2012, tbach:
# - reformatted comments
# - changed all parts of code, where the program exits inside an exception to
#   display the exception, because the messages are not always helpful
# - removed ";" from all over the code (still some left)
# - 207 errors, 983 warning, 574 infos left from static code analysis...


# Usage1 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./
# Usage2 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./




# Some rules for variable name:
# - Dictionary is used to contain the output of function
# - Variable containing 'm' is a value directly obtained from measurment data
# - Variable containing 'mdl' is a value related to model
#


####
#######
#########
VERSION='V2.38b PRO'
#########
#######
####

import sys
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')

import traceback
from metaclass import *
from numpy import *
import numpy

from math import *
import cmath
import pickle,os

from string import *

from math import *
from linreg import *

# tentative solution for SPS pseudo double plane BPM
# from SPSBPMpair import *




#######################################################
#                 Functions
#######################################################

#------------

def modelIntersect(expbpms, model):
    bpmsin=[]
    #print "start Intersect, expbpms #:", len(expbpms)
    if len(expbpms)==0:
        print >> sys.stderr, "Zero exp BPMs sent to modelIntersect"
        sys.exit(1)
    for bpm in expbpms:
        try:
            check=model.indx[bpm[1].upper()]
            bpmsin.append(bpm)
        except:
            print bpm, "Not in Model"

    if len(bpmsin)==0:
        print >> sys.stderr, "Zero intersection of Exp and Model"
        print >> sys.stderr, "Please, provide a good Dictionary or correct data"
        print >> sys.stderr, "Now we better leave!"
        sys.exit(1)
    return bpmsin


def intersect(ListOfFile):
    '''Pure intersection of all bpm names in all files '''
    if len(ListOfFile)==0:
        print >> sys.stderr, "Nothing to intersect!!!!"
        sys.exit(1)
    z=ListOfFile[0].NAME
    if len(z)==0:
        print >> sys.stderr, "No exp BPMs..."
        sys.exit(1)
    for b in ListOfFile:
        z=filter(lambda x: x in z   , b.NAME)
    #SORT by S
    result=[]
    x0=ListOfFile[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))

    result.sort()
    return result


#------------ Get phases

def phiLastAndLastButOne(phi,ftune):
    if ftune > 0.0:
        phit=phi+ftune
        if phit>1.0: phit=phit-1.0
    elif ftune <= 0.0:
        phit=phi+(1.0+ftune)
        if phit>1.0: phit=phit-1.0
    return phit

def PhaseMean(phase0,norm):  #-- phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi

    phase0   =array(phase0)%norm
    phase1   =(phase0+0.5*norm)%norm-0.5*norm
    phase0ave=mean(phase0)
    phase1ave=mean(phase1)
    phase0std=sqrt(mean((phase0-phase0ave)**2))
    phase1std=sqrt(mean((phase1-phase1ave)**2))
    if phase0std<phase1std: return phase0ave
    else                  : return phase1ave%norm

def PhaseStd(phase0,norm):  #-- phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi

    phase0   =array(phase0)%norm
    phase1   =(phase0+0.5*norm)%norm-0.5*norm
    phase0ave=mean(phase0)
    phase1ave=mean(phase1)
    phase0std=sqrt(mean((phase0-phase0ave)**2))
    phase1std=sqrt(mean((phase1-phase1ave)**2))
    return min(phase0std,phase1std)

def GetPhasesTotal(MADTwiss,ListOfFiles,Q,plane,bd,oa,op):

    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms, MADTwiss)
    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and oa=="LHCB1": s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    if op=="1" and oa=="LHCB2": s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    bn1=upper(commonbpms[0][1])
    phaseT={}
    print "Reference BPM:", bn1, "Plane:", plane
    for i in range(0,len(commonbpms)):
        #bn2=upper(commonbpms[i+1][1]) ?
        bn2=upper(commonbpms[i][1])
        if plane=='H':
            phmdl12=(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]) % 1
        if plane=='V':
            phmdl12=(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]) % 1

        phi12=[]
        for j in ListOfFiles:
            # Phase is in units of 2pi
            if plane=='H':
                phm12=(j.MUX[j.indx[bn2]]-j.MUX[j.indx[bn1]]) % 1
            if plane=='V':
                phm12=(j.MUY[j.indx[bn2]]-j.MUY[j.indx[bn1]]) % 1
            #-- To fix the phase shift by Q in LHC
            try:
                if commonbpms[i][0]>s_lastbpm: phm12+=bd*Q
            except: pass
            phi12.append(phm12)
        phi12=array(phi12)
        # for the beam circulating reversely to the model
        if bd==-1: phi12=1.0-phi12

        #phstd12=sqrt(average(phi12*phi12)-(average(phi12))**2.0+2.2e-15)
        #phi12=average(phi12)
        phstd12=PhaseStd(phi12,1.0)
        phi12  =PhaseMean(phi12,1.0)
        phaseT[bn2]=[phi12,phstd12,phmdl12,bn1]

    return [phaseT,commonbpms]


def GetPhases(MADTwiss,ListOfFiles,Q,plane,outputpath,bd,oa,op):
    #print "Hi Get", len(ListOfFiles)

    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms, MADTwiss)
    #print len(commonbpms)
    #sys.exit()

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and oa=="LHCB1": s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    if op=="1" and oa=="LHCB2": s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    mu=0.0
    tunem=[]
    phase={} # Dictionary for the output containing [average phase, rms error]
    for i in range(0,len(commonbpms)): # To find the integer part of tune as well, the loop is up to the last monitor
        if i==len(commonbpms)-2: # The last monitor but one
            bn1=upper(commonbpms[i][1])
            bn2=upper(commonbpms[i+1][1])
            bn3=upper(commonbpms[0][1])
            bns1=commonbpms[i][0]
            bns2=commonbpms[i+1][0]
        elif i==len(commonbpms)-1: # The last monitor
            bn1=upper(commonbpms[i][1])
            bn2=upper(commonbpms[0][1])
            bn3=upper(commonbpms[1][1])
            bns1=commonbpms[i][0]
            bns2=commonbpms[0][0]
        else : # Others
            bn1=upper(commonbpms[i][1])
            bn2=upper(commonbpms[i+1][1])
            bn3=upper(commonbpms[i+2][1])
            bns1=commonbpms[i][0]
            bns2=commonbpms[i+1][0]
        if (bn1==bn2):
            print >> sys.stderr, "There seem two lines with the same BPM name "+bn1+" in linx/y file."
            print >> sys.stderr, "Please check your input data....leaving GetLLM."
            sys.exit(1)
        if plane=='H':
            phmdl12=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
            phmdl13=MADTwiss.MUX[MADTwiss.indx[bn3]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
        elif plane=='V':
            phmdl12=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
            phmdl13=MADTwiss.MUY[MADTwiss.indx[bn3]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
        if i==len(commonbpms)-2:
            if plane=='H':
                madtune=MADTwiss.Q1 % 1.0
            elif plane=='V':
                madtune=MADTwiss.Q2 % 1.0
            if madtune>0.5:
                madtune=madtune-1.0
            phmdl13=phmdl13 % 1.0
            phmdl13=phiLastAndLastButOne(phmdl13,madtune)
        elif i==len(commonbpms)-1:
            if plane=='H':
                madtune=MADTwiss.Q1 % 1.0
            elif plane=='V':
                madtune=MADTwiss.Q2 % 1.0
            if madtune>0.5:
                madtune=madtune-1.0
            phmdl12=phmdl12 % 1.0
            phmdl13=phmdl13 % 1.0
            phmdl12=phiLastAndLastButOne(phmdl12,madtune)
            phmdl13=phiLastAndLastButOne(phmdl13,madtune)




        phi12=[]
        phi13=[]
        tunemi=[]
        for j in ListOfFiles:
            # Phase is in units of 2pi
            if plane=='H':
                phm12=(j.MUX[j.indx[bn2]]-j.MUX[j.indx[bn1]]) # the phase advance between BPM1 and BPM2
                phm13=(j.MUX[j.indx[bn3]]-j.MUX[j.indx[bn1]]) # the phase advance between BPM1 and BPM3
                tunemi.append(j.TUNEX[j.indx[bn1]])
            elif plane=='V':
                phm12=(j.MUY[j.indx[bn2]]-j.MUY[j.indx[bn1]]) # the phase advance between BPM1 and BPM2
                phm13=(j.MUY[j.indx[bn3]]-j.MUY[j.indx[bn1]]) # the phase advance between BPM1 and BPM3
                tunemi.append(j.TUNEY[j.indx[bn1]])
            #-- To fix the phase shift by Q in LHC
            try:
                if MADTwiss.S[MADTwiss.indx[bn1]]<=s_lastbpm and MADTwiss.S[MADTwiss.indx[bn2]] >s_lastbpm: phm12+= bd*Q
                if MADTwiss.S[MADTwiss.indx[bn1]]<=s_lastbpm and MADTwiss.S[MADTwiss.indx[bn3]] >s_lastbpm: phm13+= bd*Q
                if MADTwiss.S[MADTwiss.indx[bn1]] >s_lastbpm and MADTwiss.S[MADTwiss.indx[bn2]]<=s_lastbpm: phm12+=-bd*Q
                if MADTwiss.S[MADTwiss.indx[bn1]] >s_lastbpm and MADTwiss.S[MADTwiss.indx[bn3]]<=s_lastbpm: phm13+=-bd*Q
            except: pass
            if phm12<0: phm12+=1
            if phm13<0: phm13+=1
            phi12.append(phm12)
            phi13.append(phm13)

        phi12=array(phi12)
        phi13=array(phi13)
        if bd==-1: # for the beam circulating reversely to the model
            phi12=1.0-phi12
            phi13=1.0-phi13

        #if any(phi12)>0.9 and i !=len(commonbpms): # Very small phase advance could result in larger than 0.9 due to measurement error
        #       print 'Warning: there seems too large phase advance! '+bn1+' to '+bn2+' = '+str(phi12)+'in plane '+plane+', recommended to check.'
        phstd12=PhaseStd(phi12,1.0)
        phstd13=PhaseStd(phi13,1.0)
        phi12=PhaseMean(phi12,1.0)
        phi13=PhaseMean(phi13,1.0)
        #phstd12=sqrt(average(phi12*phi12)-(average(phi12))**2.0+2.2e-15)
        #phstd13=sqrt(average(phi13*phi13)-(average(phi13))**2.0+2.2e-15)
        #phi12=average(phi12)
        #phi13=average(phi13)
        tunemi=array(tunemi)
        if i<len(commonbpms)-1 :
            tunem.append(average(tunemi))

        # Note that the phase advance between the last monitor and the first monitor should be find by taking into account the fractional part of tune.
        if i==len(commonbpms)-2:
            tunem=array(tunem)
            tune=average(tunem)
            phi13=phiLastAndLastButOne(phi13,tune)
        elif i==len(commonbpms)-1:
            phi12=phiLastAndLastButOne(phi12,tune)
            phi13=phiLastAndLastButOne(phi13,tune)
        mu=mu+phi12

        small=0.0000001
        if (abs(phmdl12) < small):
            phmdl12=small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn2+" in MAD model is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        if (abs(phmdl13) < small):
            phmdl13=small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn3+" in MAD model is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        if (abs(phi12) < small ):
            phi12 = small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn2+" in measurement is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        if (abs(phi13) < small):
            phi13 = small
            print "Note: Phase advance (Plane"+plane+") between "+bn1+" and "+bn3+" in measurement is EXACTLY n*pi. GetLLM slightly differ the phase advance here, artificially."
            print "Beta from amplitude around this monitor will be slightly varied."
        phase[bn1]=[phi12,phstd12,phi13,phstd13,phmdl12,phmdl13,bn2]

    #print len(phase)
    #sys.exit()

    return [phase,tune,mu,commonbpms]

#-------- Beta from pahses

def BetaFromPhase(MADTwiss,ListOfFiles,phase,plane):

    alfa={}
    beta={}
    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms,MADTwiss)
    alfii=[]
    betii=[]
    delbeta=[]
    for i in range(0,len(commonbpms)):
        if i==len(commonbpms)-2: # The last monitor but one
            bn1=upper(commonbpms[i][1])
            bn2=upper(commonbpms[i+1][1])
            bn3=upper(commonbpms[0][1])
        elif i==len(commonbpms)-1: # The last monitor
            bn1=upper(commonbpms[i][1])
            bn2=upper(commonbpms[0][1])
            bn3=upper(commonbpms[1][1])
        else : # Others
            bn1=upper(commonbpms[i][1])
            bn2=upper(commonbpms[i+1][1])
            bn3=upper(commonbpms[i+2][1])
        ph2pi12=2.*pi*phase[bn1][0]
        ph2pi23=2.*pi*phase[bn2][0]
        ph2pi13=2.*pi*phase[bn1][2]
        # Find the model transfer matrices for beta1
        phmdl12=2.*pi*phase[bn1][4]
        phmdl13=2.*pi*phase[bn1][5]
        phmdl23=2.*pi*phase[bn2][4]
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
        M11=sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
        N11=sqrt(betmdl3/betmdl1)*(cos(phmdl13)+alpmdl1*sin(phmdl13))
        N12=sqrt(betmdl1*betmdl3)*sin(phmdl13)

        denom=M11/M12-N11/N12+1e-16
        numer=1/tan(ph2pi12)-1/tan(ph2pi13)
        beti1=numer/denom
        # beterr1=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # beterr1=beterr1+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # beterr1=beterr1/abs(denom)
        betstd1=        (2*pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        betstd1=betstd1+(2*pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        betstd1=sqrt(betstd1)/abs(denom)

        denom=M12/M11-N12/N11+1e-16
        numer=-M12/M11/tan(ph2pi12)+N12/N11/tan(ph2pi13)
        alfi1=numer/denom
        # alferr1=abs(M12/M11*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # alferr1=alferr1+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # alferr1=alferr1/abs(denom)
        alfstd1=        (M12/M11*2*pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        alfstd1=alfstd1+(N12/N11*2*pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        alfstd1=sqrt(alfstd1)/denom

        # Find beta2 and alpha2 from phases assuming model transfer matrix
        # Matrix M: BPM1-> BPM2
        # Matrix N: BPM2-> BPM3
        M22=sqrt(betmdl1/betmdl2)*(cos(phmdl12)-alpmdl2*sin(phmdl12))
        M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
        N11=sqrt(betmdl3/betmdl2)*(cos(phmdl23)+alpmdl2*sin(phmdl23))
        N12=sqrt(betmdl2*betmdl3)*sin(phmdl23)

        denom=M22/M12+N11/N12+1e-16
        numer=1/tan(ph2pi12)+1/tan(ph2pi23)
        beti2=numer/denom
        # beterr2=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # beterr2=beterr2+abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # beterr2=beterr2/abs(denom)
        betstd2=        (2*pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        betstd2=betstd2+(2*pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        betstd2=sqrt(betstd2)/abs(denom)

        denom=M12/M22+N12/N11+1e-16
        numer=M12/M22/tan(ph2pi12)-N12/N11/tan(ph2pi23)
        alfi2=numer/denom
        # alferr2=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # alferr2=alferr2+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # alferr2=alferr2/abs(denom)
        alfstd2=        (M12/M22*2*pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        alfstd2=alfstd2+(N12/N11*2*pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        alfstd2=sqrt(alfstd2)/abs(denom)

        # Find beta3 and alpha3 from phases assuming model transfer matrix
        # Matrix M: BPM2-> BPM3
        # Matrix N: BPM1-> BPM3
        M22=sqrt(betmdl2/betmdl3)*(cos(phmdl23)-alpmdl3*sin(phmdl23))
        M12=sqrt(betmdl2*betmdl3)*sin(phmdl23)
        N22=sqrt(betmdl1/betmdl3)*(cos(phmdl13)-alpmdl3*sin(phmdl13))
        N12=sqrt(betmdl1*betmdl3)*sin(phmdl13)

        denom=M22/M12-N22/N12+1e-16
        numer=1/tan(ph2pi23)-1/tan(ph2pi13)
        beti3=numer/denom
        # beterr3=abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # beterr3=beterr3+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # beterr3=beterr3/abs(denom)
        betstd3=        (2*pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        betstd3=betstd3+(2*pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        betstd3=sqrt(betstd3)/abs(denom)

        denom=M12/M22-N12/N22+1e-16
        numer=M12/M22/tan(ph2pi23)-N12/N22/tan(ph2pi13)
        alfi3=numer/denom
        # alferr3=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # alferr3=alferr3+abs(N22/N12*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # alferr3=alferr3/abs(denom)
        alfstd3=        (M12/M22*2*pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        alfstd3=alfstd3+(N12/N22*2*pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        alfstd3=sqrt(alfstd3)/abs(denom)

        betii.append([beti1,betstd1,beti2,betstd2,beti3,betstd3])
        alfii.append([alfi1,alfstd1,alfi2,alfstd2,alfi3,alfstd3])

    # Output the beta at all monitors even for the first, second, last and last but one using beta1,2 and 3!
    for i in range(0,len(commonbpms)):
        bn1=upper(commonbpms[i][1])
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
        betstd=sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/sqrt(3.)
        # betstd=sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/sqrt(3.*len(ListOfFiles))  # If we want to plot std "over files" !!!!!
        try:
            beterr=sqrt((betii[ib1][0]**2.+betii[ib2][2]**2.+betii[ib3][4]**2.)/3.-beti**2.)
        except:
            beterr=0
        alfi=(alfii[ib1][0]+alfii[ib2][2]+alfii[ib3][4])/3.
        alfstd=sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/sqrt(3.)
        # alfstd=sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/sqrt(3.*len(ListOfFiles))  # If we want to plot std "over files" !!!!!
        try:
            alferr=sqrt((alfii[ib1][0]**2.+alfii[ib2][2]**2.+alfii[ib3][4]**2.)/3.-alfi**2.)
        except:
            alferr=0

        beta[bn1]=(beti,beterr,betstd)
        alfa[bn1]=(alfi,alferr,alfstd)
        if plane=='H':
            betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
        delbeta.append((beti-betmdl1)/betmdl1)


    delbeta=array(delbeta)
    rmsbb=sqrt(average(delbeta*delbeta))

    return [beta,rmsbb,alfa,commonbpms]


#------------- Beta from amplitude

def BetaFromAmplitude(MADTwiss,ListOfFiles,plane):

    beta={}
    root2J=[]
    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms,MADTwiss)
    SumA=0.0
    Amp=[]
    Amp2=[]
    Kick2=[]
    for i in range(0,len(commonbpms)): # this loop have become complicated after modifications... anybody simplify?
        bn1=upper(commonbpms[i][1])
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
        root2J.append(root2Ji/sqrt(tembeta))


    Kick=SumA/len(commonbpms) # Assuming the average of beta is constant
    Kick2=array(Kick2)
    Kick2=Kick2/len(commonbpms)
    Amp2=array(Amp2)
    root2J=array(root2J)
    root2Jave=average(root2J)
    root2Jrms=sqrt(average(root2J*root2J)-root2Jave**2+2.2e-16)

    #print Amp2/Kick2


    delbeta=[]
    for i in range(0,len(commonbpms)):
        bn1=upper(commonbpms[i][1])
        location=commonbpms[i][0]
        for j in range(0,len(ListOfFiles)):
            Amp2[i][j]=Amp2[i][j]/Kick2[j]
        #print average(Amp2[i]*Amp2[i]),average(Amp2[i])**2
        try:
            betstd=sqrt(average(Amp2[i]*Amp2[i])-average(Amp2[i])**2+2.2e-16)
        except:
            betstd=0

        beta[bn1]=[Amp[i]**2/Kick,betstd,location]
        if plane=='H':
            betmdl=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            betmdl=MADTwiss.BETY[MADTwiss.indx[bn1]]
        delbeta.append((beta[bn1][0]-betmdl)/betmdl)

    invariantJ=[root2Jave,root2Jrms]

    delbeta=array(delbeta)
    rmsbb=sqrt(average(delbeta*delbeta))
    return [beta,rmsbb,commonbpms,invariantJ]

#-------------------------

def GetCO(MADTwiss, ListOfFiles):

    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms, MADTwiss)
    co={} # Disctionary for output
    for i in range(0,len(commonbpms)):
        bn1=upper(commonbpms[i][1])
        bns1=commonbpms[i][0]
        coi=0.0
        coi2=0.0
        for j in ListOfFiles:
            coi=coi + j.CO[j.indx[bn1]]
            coi2=coi2 + j.CO[j.indx[bn1]]**2
        coi=coi/len(ListOfFiles)
        corms=sqrt(coi2/len(ListOfFiles)-coi**2+2.2e-16)
        co[bn1]=[coi,corms]
    return [co, commonbpms]


def NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, betax, COcut):

    nzdpp=len(ListOfNonZeroDPPX) # How many non zero dpp
    zdpp=len(ListOfZeroDPPX)  # How many zero dpp
    if zdpp==0 or nzdpp ==0:
        print 'Warning: No data for dp/p!=0 or even for dp/p=0.'
        print 'No output for the  dispersion.'
        dum0={}
        dum1=[]
        return [dum0, dum1]


    coac=ListOfCOX[0] # COX dictionary after cut bad BPMs
    coact={}
    for i in coac:
        if (coac[i][1] < COcut):
            coact[i]=coac[i]
    coac=coact


    nda={} # Dictionary for the output containing [average Disp, rms error]

    ALL=ListOfZeroDPPX+ListOfNonZeroDPPX
    commonbpmsALL=intersect(ALL)
    commonbpmsALL=modelIntersect(commonbpmsALL, MADTwiss)

    mydp=[]
    gf=[]
    for j in ListOfNonZeroDPPX:
        mydp.append(float(j.DPP))
        gf.append(0.0) # just to initialize the array, may not be a clever way...
    mydp=array(mydp)
    wf=array(abs(mydp))/sum(abs(mydp))*len(mydp) #Weight for the average depending on DPP


    # Find the global factor
    nd=[]
    ndmdl=[]
    badco=0
    for i in range(0,len(commonbpmsALL)):
        bn1=upper(commonbpmsALL[i][1])
        bns1=commonbpmsALL[i][0]
        ndmdli=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
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

    ndmdl=array(ndmdl)
    avemdl=average(ndmdl)

    gf=array(gf)
    gf=gf/avemdl/(len(commonbpmsALL)-badco)


    # Find normalized dispersion and Dx construction
    nd=array(nd)
    Dx={}
    dummy=0.0 # dummy for DXSTD
    bpms=[]
    badco=0
    for i in range(0,len(commonbpmsALL)):
        ndi=[]
        bn1=upper(commonbpmsALL[i][1])
        bns1=commonbpmsALL[i][0]
        try:
            coac[bn1]
            for j in range(0,nzdpp): # the range(0,nzdpp) instead of ListOfZeroDPPX is used because the index j is used in the loop
                ndi.append(nd[i-badco][j]/gf[j])
            ndi=array(ndi)
            ndstd=sqrt(average(ndi*ndi)-(average(ndi))**2.0+2.2e-16)
            ndas=average(wf*ndi)
            nda[bn1]=[ndas,ndstd]
            Dx[bn1]=[nda[bn1][0]*sqrt(betax[bn1][0]),dummy]
            bpms.append([bns1,bn1])
        except:
            badco+=1
    DPX=GetDPX(MADTwiss, Dx, bpms)

    return [nda,Dx,DPX,bpms]


#---------- DPX, New!
def GetDPX(MADTwiss,Dx,commonbpms):

    DPX={}
    for i in range(0,len(commonbpms)):
        bn1=upper(commonbpms[i][1])
        if i==len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2=upper(commonbpms[0][1])
            phmdl12=2.*pi*(MADTwiss.Q1+MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
        else:
            bn2=upper(commonbpms[i+1][1])
            phmdl12=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
        betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
        betmdl2=MADTwiss.BETX[MADTwiss.indx[bn2]]
        alpmdl1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
        alpmdl2=MADTwiss.ALFX[MADTwiss.indx[bn2]]
        dxmdl1=MADTwiss.DX[MADTwiss.indx[bn1]]
        dpxmdl1=MADTwiss.DPX[MADTwiss.indx[bn1]]
        dxmdl2=MADTwiss.DX[MADTwiss.indx[bn2]]
        dpxmdl2=MADTwiss.DPX[MADTwiss.indx[bn2]]

        M11=sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
        M13=dxmdl2-M11*dxmdl1-M12*dpxmdl1
        # use the beta from amplitude
        DPX[bn1]=(-M13+Dx[bn2][0]-M11*Dx[bn1][0])/M12

    return DPX



#----------- DPY, New!
def GetDPY(MADTwiss,Dy,commonbpms):
    DPY={}
    for i in range(0,len(commonbpms)):
        bn1=upper(commonbpms[i][1])
        if i==len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2=upper(commonbpms[0][1])
            phmdl12=2.*pi*(MADTwiss.Q2+MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
        else:
            bn2=upper(commonbpms[i+1][1])
            phmdl12=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
        betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
        betmdl2=MADTwiss.BETY[MADTwiss.indx[bn2]]
        alpmdl1=MADTwiss.ALFY[MADTwiss.indx[bn1]]
        alpmdl2=MADTwiss.ALFY[MADTwiss.indx[bn2]]
        dymdl1=MADTwiss.DY[MADTwiss.indx[bn1]]
        dpymdl1=MADTwiss.DPY[MADTwiss.indx[bn1]]
        dymdl2=MADTwiss.DY[MADTwiss.indx[bn2]]
        dpymdl2=MADTwiss.DPY[MADTwiss.indx[bn2]]

        M11=sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
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
    commonbpmsALL=intersect(ALL)

    nzdpp=len(ListOfNonZeroDPP) # How many non zero dpp
    zdpp=len(ListOfZeroDPP)  # How many zero dpp


    mydp=[]
    for j in ListOfNonZeroDPP:
        mydp.append(float(j.DPP))
    mydp=array(mydp)
    wf=array(abs(mydp))/sum(abs(mydp))*len(mydp) #Weitghs for the average


    dco={} # Dictionary for the output containing [(averaged)Disp, rms error]
    bpms=[]
    for i in range(0,len(commonbpmsALL)):
        bn1=upper(commonbpmsALL[i][1])
        bns1=commonbpmsALL[i][0]

        try:
            coi=coac[bn1]
            dcoi=[]
            for j in ListOfNonZeroDPP:
                dcoi.append((j.CO[j.indx[bn1]]-coi[0])*scalef/float(j.DPP))
            dcoi=array(dcoi)
            dcostd=sqrt(average(dcoi*dcoi)-(average(dcoi))**2.0+2.2e-16)
            dcos=average(wf*dcoi)
            dco[bn1]=[dcos,dcostd]
            bpms.append([bns1,bn1])
        except:
            coi=0
    return [dco,bpms]


#-----------

def GetCoupling1(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2):

    # not applicable to db=-1 for the time being...

    tp=2.0*pi

    # find operation point
    try:
        fdi=open(outputpath+'Drive.inp','r')  # Drive.inp file is normally in the outputpath directory in GUI operation
        for line in fdi:
            if "TUNE X" in line:
                fracxinp=line.split("=")
                fracx=fracxinp[1]
            if "TUNE Y" in line:
                fracyinp=line.split("=")
                fracy=fracyinp[1]
        fdi.close()
    except:
        fracx=Q1 # Otherwise, the fractional parts are assumed to be below 0.5
        fracy=Q2

    if fracx<0.0 :
        fracx=1.0-Q1
    else:
        fracx=Q1
    if fracy<0.0 :
        fracx=1.0-Q2
    else:
        fracy=Q2

    if fracx>fracy:
        sign_QxmQy=1.0
    else:
        sign_QxmQy=-1.0

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement.
    if len(ListOfZeroDPPX)!=len(ListOfZeroDPPY):
        print 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]


    XplusY=ListOfZeroDPPX+ListOfZeroDPPY
    dbpms=intersect(XplusY)
    dbpms=modelIntersect(dbpms, MADTwiss)


    # caculate fw and qw, exclude bpms having wrong phases

    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)):
        bn1=upper(dbpms[i][1])

        fij=[]
        qij=[]
        q1j=[]
        q2j=[]
        badbpm=0
        for j in range(0,len(ListOfZeroDPPX)):
            jx=ListOfZeroDPPX[j]
            jy=ListOfZeroDPPY[j]
            C01ij=jx.AMP01[jx.indx[bn1]]
            C10ij=jy.AMP10[jy.indx[bn1]]
            fij.append(0.5*atan(sqrt(C01ij*C10ij)))

            #q1=(jx.MUX[jx.indx[bn1]]-jy.PHASE10[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
            #q2=(jx.PHASE01[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]]-0.25)%1.0
            #q1=(0.5-q1)%1.0 # This sign change in the real part is to comply with MAD output
            #q2=(0.5-q2)%1.0
            q1j.append((jx.MUX[jx.indx[bn1]]-jy.PHASE10[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
            q2j.append((jx.PHASE01[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]]-0.25)%1.0)
            q1j[j]=(0.5-q1j[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2j[j]=(0.5-q2j[j])%1.0

            #if abs(q1-q2)<0.25:
            #       qij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            #       qij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
            #       badbpm=1
            #       countBadPhase += 1
            #       #print "Bad Phases in BPM ",bn1, "total so far", countBadPhase
        q1j=array(q1j)
        q2j=array(q2j)
        q1=average(q1j)
        q2=average(q2j)

        if abs(q1-q2)<0.25:  # Very rough cut !!!!!!!!!!!!!!!!!!!
            qi=(q1+q2)/2.0
        elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            qi=q1 # Note that q1 and q2 are confined 0. to 1.
        else:
            badbpm=1
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1, "total so far", countBadPhase



        if badbpm==0:
            fij=array(fij)
            fi=average(fij)
            fistd=sqrt(average(fij*fij)-(average(fij))**2.0+2.2e-16)
            #qij=array(qij)
            #qi=average(qij)
            #qistd=sqrt(average(qij*qij)-(average(qij))**2.0+2.2e-16)
            qistd=sqrt(average(q1j*q1j)-(average(q1j))**2.0+2.2e-16) # Not very exact...
            fi=fi*complex(cos(tp*qi),sin(tp*qi))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])
            fwqw[bn1]=[[fi,fistd],[qi,qistd]]


    dbpms=dbpmt


    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)):
        jx=ListOfZeroDPPX[j]
        jy=ListOfZeroDPPY[j]
        bn1=upper(dbpms[i][1])
        CG=CG+sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(jx.MUX[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]])


    CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
    QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]


    return [fwqw,dbpms]

#-----------

def ComplexSecondaryLine(delta, cw, cw1, pw, pw1):
    tp=2.0*pi
    a1=complex(1.0,-tan(tp*delta))
    a2=cw*complex(cos(tp*pw),sin(tp*pw))
    a3=-1.0/cos(tp*delta)*complex(0.0,1.0)
    a4=cw1*complex(cos(tp*pw1),sin(tp*pw1))
    SL=a1*a2+a3*a4
    sizeSL=sqrt(SL.real**2+SL.imag**2)
    phiSL=(arctan2(SL.imag , SL.real)/tp) %1.0
    #SL=complex(-SL.real,SL.imag)    # This sign change in the real part is to comply with MAD output
    return [sizeSL,phiSL]


def ComplexSecondaryLineExtended(delta,edelta, amp1,amp2, phase1,phase2):
    '''
     Input : - delta: phase advance between two BPMs
             - edelta: error on the phase advance between two BPMs
             - amp1: amplitude of secondary line at ith BPM
             - amp2: amplitude of secondary line at i+1th BPM
             - phase1: phase of secondary line at ith BPM
             - phase2: phase of secondary line at i+1th BPM
     Return: - amp: amplitude of the complex signal
             - phase: phase of the complex signal
             - eamp: error on amplitude of the complex signal
             - ephase: error on phase of the complex signal
    '''


    # functions
    tp=2.0*pi
    C=cos(delta*tp)
    S=sin(delta*tp)
    T=tan(delta*tp)
    SC=sin(delta*tp)/((cos(delta*tp*2)+1)/2)

    # signal
    cs1=cos(tp*phase1)
    ss1=sin(tp*phase1)
    cs2=cos(tp*phase2)
    ss2=sin(tp*phase2)

    sig1=amp1*complex(cs1,ss1)
    sig2=amp2*complex(cs2,ss2)

    # computing complex secondary line (h-)
    sig=sig1*complex(1,-T)-sig2*complex(0,1)*(1/C)

    amp=abs(sig)
    phase=(arctan2(sig.imag,sig.real)/tp) %1.0

    # computing error secondary line (h-)
    esig=(sig1*complex(1,-(1/C))-sig2*complex(0,1)*(SC))*edelta

    eamp=abs(esig)
    ephase=(arctan2(esig.imag,esig.real)/tp) %1.0

    return [amp,phase,eamp,ephase]


def GetCoupling2(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, phasex, phasey, bd, oa):


    # find operation point
    try:
        fdi=open(outputpath+'Drive.inp','r')  # Drive.inp file is normally in the outputpath directory in GUI operation
        for line in fdi:
            if "TUNE X" in line:
                fracxinp=line.split("=")
                fracx=fracxinp[1]
            if "TUNE Y" in line:
                fracyinp=line.split("=")
                fracy=fracyinp[1]
        fdi.close()
    except:
        fracx=Q1 # Otherwise, the fractional parts are assumed to be below 0.5
        fracy=Q2

    if fracx<0.0 :
        fracx=1.0-Q1
    else:
        fracx=Q1
    if fracy<0.0 :
        fracx=1.0-Q2
    else:
        fracy=Q2

    if fracx>fracy:
        sign_QxmQy=1.0
    else:
        sign_QxmQy=-1.0

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(ListOfZeroDPPX)!=len(ListOfZeroDPPY):
        print 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]


    XplusY=ListOfZeroDPPX+ListOfZeroDPPY
    dbpms=intersect(XplusY)
    dbpms=modelIntersect(dbpms, MADTwiss)


    # caculate fw and qw, exclude bpms having wrong phases

    tp=2.0*pi
    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)-1):
        bn1=upper(dbpms[i][1])
        bn2=upper(dbpms[i+1][1])

        delx= phasex[bn1][0] - 0.25  # Missprint in the coupling note
        dely= phasey[bn1][0] - 0.25

        f1001ij=[]
        #q1001ij=[]
        f1010ij=[]
        #q1010ij=[]
        q1js=[]
        q2js=[]
        q1jd=[]
        q2jd=[]
        badbpm=0
        for j in range(0,len(ListOfZeroDPPX)):
            jx=ListOfZeroDPPX[j]
            jy=ListOfZeroDPPY[j]
            [SA0p1ij,phi0p1ij]=ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]], jx.PHASE01[jx.indx[bn1]], jx.PHASE01[jx.indx[bn2]])
            [SA0m1ij,phi0m1ij]=ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]], -jx.PHASE01[jx.indx[bn1]], -jx.PHASE01[jx.indx[bn2]])
            [TBp10ij,phip10ij]=ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]], jy.PHASE10[jy.indx[bn1]], jy.PHASE10[jy.indx[bn2]])
            [TBm10ij,phim10ij]=ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]], -jy.PHASE10[jy.indx[bn1]], -jy.PHASE10[jy.indx[bn2]])


            #print SA0p1ij,phi0p1ij,SA0m1ij,phi0m1ij,TBp10ij,phip10ij,TBm10ij,phim10ij
            f1001ij.append(0.5*sqrt(TBp10ij*SA0p1ij/2.0/2.0))
            f1010ij.append(0.5*sqrt(TBm10ij*SA0m1ij/2.0/2.0))

            if bd==1:
                q1jd.append((phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append((-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0)
            elif bd==-1:
                q1jd.append((phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append(-(-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0)
            #print q1,q2
            q1jd[j]=(0.5-q1jd[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2jd[j]=(0.5-q2jd[j])%1.0


            #if abs(q1-q2)<0.25:
                #q1001ij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
                #q1001ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
                #badbpm=1
                #q1001ij.append(q1)
                #countBadPhase += 1
                #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

            if bd==1:
                q1js.append((phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append((phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0)
            if bd==-1:
                q1js.append((phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append(-(phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0)
            #print q1,q2
            q1js[j]=(0.5-q1js[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2js[j]=(0.5-q2js[j])%1.0

            #if abs(q1-q2)<0.25:
                #q1010ij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
                #q1010ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
                #badbpm=1
                #if (oa=="SPS" or oa=="RHIC"):
                #       badbpm=0
                #q1010ij.append(q1)
                #countBadPhase += 1
                #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

        q1jd=array(q1jd)
        q2jd=array(q2jd)
        q1d=average(q1jd)
        q2d=average(q2jd)

        q1js=array(q1js)
        q2js=array(q2js)
        q1s=average(q1js)
        q2s=average(q2js)


        if abs(q1d-q2d)<0.25:
            q1001i=(q1d+q2d)/2.0
        elif abs(q1d-q2d)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            q1001i=q1d # Note that q1 and q2 are confined 0. to 1.
        else:
            badbpm=1
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase
        if abs(q1s-q2s)<0.25:
            q1010i=(q1s+q2s)/2.0
        elif abs(q1s-q2s)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            q1010i=q1s # Note that q1 and q2 are confined 0. to 1.
        else:
            badbpm=1
        if (oa=="SPS" or "RHIC" in oa):
            # No check for the SPS or RHIC
            badbpm=0
            q1001i=q1d
            #q1001r=q1s
            #q1010i=q1d
            q1010i=q1s
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase



        if badbpm==0:
            f1001ij=array(f1001ij)
            f1001i=average(f1001ij)
            f1001istd=sqrt(average(f1001ij*f1001ij)-(average(f1001ij))**2.0+2.2e-16)
            f1010ij=array(f1010ij)
            f1010i=average(f1010ij)
            try:
                f1010istd=sqrt(average(f1010ij*f1010ij)-(average(f1010ij))**2.0+2.2e-16)
            except:
                f1010istd=0
            #q1001ij=array(q1001ij)
            #q1001i=average(q1001ij)
            try:
                q1001istd=sqrt(average(q1jd*q1jd)-(q1001i)**2.0+2.2e-16) # Not very correct
            except:
                q1001istd=0
            #q1010ij=array(q1010ij)
            #q1010i=average(q1010ij)
            try:
                q1010istd=sqrt(average(q1js*q1js)-(q1010i)**2.0+2.2e-16)
            except:
                q1010istd=0
            f1001i=f1001i*complex(cos(tp*q1001i),sin(tp*q1001i))
            f1010i=f1010i*complex(cos(tp*q1010i),sin(tp*q1010i))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])
            if bd==1:
                fwqw[bn1]=[[f1001i,f1001istd,f1010i,f1010istd],[q1001i,q1001istd,q1010i,q1010istd]]
            elif bd==-1:
                fwqw[bn1]=[[f1010i,f1010istd,f1001i,f1001istd],[q1010i,q1010istd,q1001i,q1001istd]]



    dbpms=dbpmt

    # possible correction ??
    #bn0=upper(dbpms[0][1])
    #up1=fwqw[bn0][0][0]
    #up2=fwqw[bn0][0][2]
    #for i in range(1,len(dbpms)):
        #bn0=upper(dbpms[i-1][1])
        #bn1=upper(dbpms[i][1])
        #df1001=sqrt(fwqw[bn0][0][0].real**2+fwqw[bn0][0][0].imag**2)/sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        #df1010=sqrt(fwqw[bn0][0][2].real**2+fwqw[bn0][0][2].imag**2)/sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2)
        #fwqw[bn0][0][0]=up1
        #fwqw[bn0][0][2]=up2
        #up1=complex(df1001*fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag)
        #up2=complex(df1010*fwqw[bn1][0][2].real,fwqw[bn1][0][2].imag)

    #fwqw[bn1][0][0]=up1
    #fwqw[bn1][0][2]=up2
    # end of possible correction



    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)-1):
        jx=ListOfZeroDPPX[0]
        jy=ListOfZeroDPPY[0]
        bn1=upper(dbpms[i][1])
        CG=CG+sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(jx.MUX[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]])

    if len(dbpms)==0:
        print 'Warning: There is no BPM to output linear coupling properly... leaving Getcoupling.'
        fwqw['Global']=[-99,-99] #Quick fix Evian 2012
        return [fwqw,dbpms]
    else:
        CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
        QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]


    return [fwqw,dbpms]

#-----------

def GetCoupling2b(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, phasex, phasey, bd, oa):

    # This fixes the averaging of multiple files in GetCoupling2.
    # Everything else is the same. By Ryoichi

    # find operation point
    try:
        fdi=open(outputpath+'Drive.inp','r')  # Drive.inp file is normally in the outputpath directory in GUI operation
        for line in fdi:
            if "TUNE X" in line:
                fracxinp=line.split("=")
                fracx=fracxinp[1]
            if "TUNE Y" in line:
                fracyinp=line.split("=")
                fracy=fracyinp[1]
        fdi.close()
    except:
        fracx=Q1 # Otherwise, the fractional parts are assumed to be below 0.5
        fracy=Q2

    if fracx<0.0 :
        fracx=1.0-Q1
    else:
        fracx=Q1

    if fracy<0.0 :
        fracx=1.0-Q2
    else:
        fracy=Q2

    if fracx>fracy:
        sign_QxmQy=1.0
    else:
        sign_QxmQy=-1.0

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(ListOfZeroDPPX)!=len(ListOfZeroDPPY):
        print 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]


    XplusY=ListOfZeroDPPX+ListOfZeroDPPY
    dbpms=intersect(XplusY)
    dbpms=modelIntersect(dbpms, MADTwiss)

    # caculate fw and qw, exclude bpms having wrong phases

    tp=2.0*pi
    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)-1):
        bn1=upper(dbpms[i][1])
        bn2=upper(dbpms[i+1][1])

        delx= phasex[bn1][0] - 0.25  # Missprint in the coupling note
        dely= phasey[bn1][0] - 0.25

        f1001ij=[]
        #q1001ij=[]
        f1010ij=[]
        #q1010ij=[]
        q1js=[]
        q2js=[]
        q1jd=[]
        q2jd=[]
        badbpm=0
        for j in range(0,len(ListOfZeroDPPX)):
            jx=ListOfZeroDPPX[j]
            jy=ListOfZeroDPPY[j]
            [SA0p1ij,phi0p1ij]=ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]], jx.PHASE01[jx.indx[bn1]], jx.PHASE01[jx.indx[bn2]])
            [SA0m1ij,phi0m1ij]=ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]], -jx.PHASE01[jx.indx[bn1]], -jx.PHASE01[jx.indx[bn2]])
            [TBp10ij,phip10ij]=ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]], jy.PHASE10[jy.indx[bn1]], jy.PHASE10[jy.indx[bn2]])
            [TBm10ij,phim10ij]=ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]], -jy.PHASE10[jy.indx[bn1]], -jy.PHASE10[jy.indx[bn2]])


            #print SA0p1ij,phi0p1ij,SA0m1ij,phi0m1ij,TBp10ij,phip10ij,TBm10ij,phim10ij
            f1001ij.append(0.5*sqrt(TBp10ij*SA0p1ij/2.0/2.0))
            f1010ij.append(0.5*sqrt(TBm10ij*SA0m1ij/2.0/2.0))

            if bd==1:
                q1jd.append((phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append((-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0)
            elif bd==-1:
                q1jd.append((phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append(-(-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0)
            #print q1,q2
            q1jd[j]=(0.5-q1jd[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2jd[j]=(0.5-q2jd[j])%1.0


            #if abs(q1-q2)<0.25:
                #q1001ij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
                #q1001ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
                #badbpm=1
                #q1001ij.append(q1)
                #countBadPhase += 1
                #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

            if bd==1:
                q1js.append((phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append((phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0)
            if bd==-1:
                q1js.append((phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append(-(phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0)
            #print q1,q2
            q1js[j]=(0.5-q1js[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2js[j]=(0.5-q2js[j])%1.0

            #if abs(q1-q2)<0.25:
                #q1010ij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
                #q1010ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
                #badbpm=1
                #if (oa=="SPS" or oa=="RHIC"):
                #       badbpm=0
                #q1010ij.append(q1)
                #countBadPhase += 1
                #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

        #--------- Different from getcoupling2 ->

        q1jd=array(q1jd)
        q2jd=array(q2jd)
        q1d=PhaseMean(q1jd,1.0)
        q2d=PhaseMean(q2jd,1.0)

        q1js=array(q1js)
        q2js=array(q2js)
        q1s=PhaseMean(q1js,1.0)
        q2s=PhaseMean(q2js,1.0)

        if min(abs(q1d-q2d),1.0-abs(q1d-q2d))>0.25 or min(abs(q1s-q2s),1.0-abs(q1s-q2s))>0.25:
            badbpm=1
            countBadPhase += 1

        if (oa=="SPS" or oa=="RHIC"):
            # No check for the SPS or RHIC
            badbpm=0
            q1010i=q1d
            q1010i=q1s
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

        if badbpm==0:

            f1001ij=array(f1001ij)
            f1001i=average(f1001ij)
            f1001istd=sqrt(average(f1001ij*f1001ij)-(average(f1001ij))**2.0+2.2e-16)
            f1010ij=array(f1010ij)
            f1010i=average(f1010ij)
            f1010istd=sqrt(average(f1010ij*f1010ij)-(average(f1010ij))**2.0+2.2e-16)

            q1001i=PhaseMean(array([q1d,q2d]),1.0)
            q1010i=PhaseMean(array([q1s,q2s]),1.0)
            q1001istd=PhaseStd(numpy.append(q1jd,q2jd),1.0)
            q1010istd=PhaseStd(numpy.append(q1js,q2js),1.0)

            f1001i=f1001i*complex(cos(tp*q1001i),sin(tp*q1001i))
            f1010i=f1010i*complex(cos(tp*q1010i),sin(tp*q1010i))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])

            if bd==1:
                fwqw[bn1]=[[f1001i,f1001istd,f1010i,f1010istd],[q1001i,q1001istd,q1010i,q1010istd]]
            elif bd==-1:
                fwqw[bn1]=[[f1010i,f1010istd,f1001i,f1001istd],[q1010i,q1010istd,q1001i,q1001istd]]

        #--------- <- Different from getcoupling2

    dbpms=dbpmt

    # possible correction ??
    #bn0=upper(dbpms[0][1])
    #up1=fwqw[bn0][0][0]
    #up2=fwqw[bn0][0][2]
    #for i in range(1,len(dbpms)):
        #bn0=upper(dbpms[i-1][1])
        #bn1=upper(dbpms[i][1])
        #df1001=sqrt(fwqw[bn0][0][0].real**2+fwqw[bn0][0][0].imag**2)/sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        #df1010=sqrt(fwqw[bn0][0][2].real**2+fwqw[bn0][0][2].imag**2)/sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2)
        #fwqw[bn0][0][0]=up1
        #fwqw[bn0][0][2]=up2
        #up1=complex(df1001*fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag)
        #up2=complex(df1010*fwqw[bn1][0][2].real,fwqw[bn1][0][2].imag)

    #fwqw[bn1][0][0]=up1
    #fwqw[bn1][0][2]=up2
    # end of possible correction

    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)-1):
        jx=ListOfZeroDPPX[0]
        jy=ListOfZeroDPPY[0]
        bn1=upper(dbpms[i][1])
        CG=CG+sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(jx.MUX[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]])

    if len(dbpms)==0:
        print 'Warning: There is no BPM to output linear coupling properly... leaving Getcoupling.'
        fwqw['Global']=[CG,QG] #Quick fix Evian 2012
        return [fwqw,dbpms]
    else:
        CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
        QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]

    return [fwqw,dbpms]

#---------------------------
def PseudoDoublePlaneMonitors(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, BPMdictionary):



    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(ListOfZeroDPPX)!=len(ListOfZeroDPPY):
        print 'Leaving PseudoDoublePlaneMonitors as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]

    bpmh=intersect(ListOfZeroDPPX)
    bpmv=intersect(ListOfZeroDPPY)
    bpmh=modelIntersect(bpmh, MADTwiss)
    bpmv=modelIntersect(bpmv, MADTwiss)


    fbpmx=[]
    fbpmy=[]
    for i in range(0,len(ListOfZeroDPPX)):
        filex='temp'+str(i)+'_linx'
        filey='temp'+str(i)+'_liny'
        fbpmxi=open(filex,'w')
        fbpmyi=open(filey,'w')
        fbpmx.append(fbpmxi)
        fbpmy.append(fbpmyi)
        fbpmx[i].write('* NAME   S      TUNEX  MUX    AMPX   AMP01  PHASE01\n')
        fbpmy[i].write('* NAME   S      TUNEY  MUY    AMPY   AMP10  PHASE10\n')
        fbpmx[i].write('$ %s     %le    %le    %le    %le    %le    %le\n')
        fbpmy[i].write('$ %s     %le    %le    %le    %le    %le    %le\n')

    bpmhp=[]
    for i in range(0,len(bpmh)):
        smin=1.0e10
        jsave=0
        for j in range (0,len(bpmv),10):
            sdiff=abs(bpmh[i][0]-bpmv[j][0])
            if sdiff<smin:
                smin=sdiff
                jsave=j
        jlower=jsave-9
        jupper=jsave+9
        if jupper > len(bpmv):
            jupper=len(bpmv)
        for j in range (jlower,jupper):
            sdiff=abs(bpmh[i][0]-bpmv[j][0])
            if sdiff<smin:
                smin=sdiff
                jsave=j

        bpmhp.append([bpmh[i][0],bpmh[i][1],bpmv[jsave][1],0])

    #bpmvp=[]
    #for i in range(0,len(bpmv)):
        #smin=1.0e10
        #jsave=0
        #for j in range (0,len(bpmh),10):
            #sdiff=abs(bpmv[i][0]-bpmh[j][0])
            #if sdiff<smin:
                #smin=sdiff
                #jsave=j
        #jlower=jsave-9
        #jupper=jsave+9
        #if jupper > len(bpmh):
            #jupper=len(bpmh)
        #for j in range (jlower,jupper):
            #sdiff=abs(bpmv[i][0]-bpmh[j][0])
            #if sdiff<smin:
                #smin=sdiff
                #jsave=j

        #bpmvp.append([bpmv[i][0],bpmv[i][1],bpmh[jsave][1],1])


    #dbpms=combinebpms(bpmhp,bpmvp)
    dbpms=bpmhp

    # tentative solution
    dbpms=bpmpair() # model BPM name
    countofmissingBPMs=0
    for i in range(0,len(dbpms)):
        wname=upper(dbpms[i][1]) # horizontal BPM basis of the pairing (model name)
        pname=upper(dbpms[i][2]) # vertical pair of the horizontal as in SPSBPMpairs (model name)
        #print wname
        ws=dbpms[i][0]  # Location
        #print "name ",wname, pname
        #Check whether the inputs (linx/y) have BPM name of model or experiment
        try:
            exwname=BPMdictionary[wname][0] #Experimental BPM name of horizontal To be paired
            #print exwname

            expname=BPMdictionary[pname][1] #Experimental BPM name of vertical  (one of them does not exist!) to be paired

            #print expname

        except:
            if len(BPMdictionary)!=0:
                countofmissingBPMs = countofmissingBPMs + 1
                print wname, "or", pname, "not found in the BPMdictionary. Total so far = ",countofmissingBPMs
        try:
            for j in range(0,len(ListOfZeroDPPX)):
                jx=ListOfZeroDPPX[j]
                jy=ListOfZeroDPPY[j]
                #if dbpms[i][3]==0:
                dphix=MADTwiss.MUX[MADTwiss.indx[upper(pname)]]-MADTwiss.MUX[MADTwiss.indx[upper(wname)]] # dphix is not used anyway
                dphiy=MADTwiss.MUY[MADTwiss.indx[upper(pname)]]-MADTwiss.MUY[MADTwiss.indx[upper(wname)]]
                # Going to try using model names, to be able to use simulation data
                try:
                    wampx=jx.AMPX[jx.indx[wname]]
                    wampy=jy.AMPY[jy.indx[pname]]
                    wamp01=jx.AMP01[jx.indx[wname]]
                    wamp10=jy.AMP10[jy.indx[pname]]
                    wtunex=jx.TUNEX[jx.indx[wname]]
                    wtuney=jy.TUNEY[jy.indx[pname]]
                    wmux=jx.MUX[jx.indx[wname]]
                    wmuy=(jy.MUY[jy.indx[pname]]-dphiy)%1.0
                    if (wmuy > 0.5): wmuy=wmuy-1.0
                    wphase01=jx.PHASE01[jx.indx[wname]]
                    wphase10=(jy.PHASE10[jy.indx[pname]]-dphiy)%1.0
                    if (wphase10 > 0.5): wphase10=wphase10-1.0
                # This seems to be experiment data, going to try with experimental names
                except:
                    wampx=jx.AMPX[jx.indx[exwname]]
                    wampy=jy.AMPY[jy.indx[expname]]
                    wamp01=jx.AMP01[jx.indx[exwname]]
                    wamp10=jy.AMP10[jy.indx[expname]]
                    wtunex=jx.TUNEX[jx.indx[exwname]]
                    wtuney=jy.TUNEY[jy.indx[expname]]
                    wmux=jx.MUX[jx.indx[exwname]]
                    wmuy=(jy.MUY[jy.indx[expname]]-dphiy)%1.0
                    if (wmuy > 0.5): wmuy=wmuy-1.0
                    wphase01=jx.PHASE01[jx.indx[exwname]]
                    wphase10=(jy.PHASE10[jy.indx[expname]]-dphiy)%1.0
                    if (wphase10 > 0.5): wphase10=wphase10-1.0
                #elif dbpms[i][3]==1:
                    #wampx=jx.AMPX[jx.indx[pname]]
                    #wampy=jy.AMPY[jy.indx[wname]]
                    #wamp01=jx.AMP01[jx.indx[pname]]
                    #wamp10=jy.AMP10[jy.indx[wname]]
                    #wtunex=jx.TUNEX[jx.indx[pname]]
                    #wtuney=jy.TUNEY[jy.indx[wname]]
                    #dphix=MADTwiss.MUX[MADTwiss.indx[upper(pname)]]-MADTwiss.MUX[MADTwiss.indx[upper(wname)]]
                    #dphiy=MADTwiss.MUY[MADTwiss.indx[upper(pname)]]-MADTwiss.MUY[MADTwiss.indx[upper(wname)]]
                    #wmux=(jx.MUX[jx.indx[pname]]-dphix)%1.0
                    #if (wmux > 0.5): wmux=wmux-1
                    #wmuy=jy.MUY[jy.indx[wname]]
                    #wphase01=(jx.PHASE01[jx.indx[pname]]-dphix)%1.0
                    #wphase10=jy.PHASE10[jy.indx[wname]]
                    #if (wphase01 > 0.5): wphase01=wphase01-1
                #elif dbpms[i][3]==2:
                    #wampx=jx.AMPX[jx.indx[wname]]
                    #wampy=jy.AMPY[jy.indx[wname]]
                    #wamp01=jx.AMP01[jx.indx[wname]]
                    #wamp10=jy.AMP10[jy.indx[wname]]
                    #wtunex=jx.TUNEX[jx.indx[wname]]
                    #wtuney=jy.TUNEY[jy.indx[wname]]
                    #wmux=jx.MUX[jx.indx[wname]]
                    #wmuy=jy.MUY[jy.indx[wname]]
                    #wphase01=jx.PHASE01[jx.indx[wname]]
                    #wphase10=jy.PHASE10[jy.indx[wname]]
                fbpmx[j].write('"'+wname+'" '+str(ws)+' '+str(wtunex)+' '+str(wmux)+' '+str(wampx)+' '+str(wamp01)+' '+str(wphase01)+'\n')
                fbpmy[j].write('"'+wname+'" '+str(ws)+' '+str(wtuney)+' '+str(wmuy)+' '+str(wampy)+' '+str(wamp10)+' '+str(wphase10)+'\n')
        except:
            if len(BPMdictionary)!=0:
                countofmissingBPMs = countofmissingBPMs + 1
                print wname, "or", pname, "not found in the DATA. Total so far = ",countofmissingBPMs


    PseudoListX=[]
    PseudoListY=[]
    for j in range(0,len(ListOfZeroDPPX)):
        fbpmx[j].close()
        fbpmy[j].close()
        filex='temp'+str(j)+'_linx'
        filey='temp'+str(j)+'_liny'
        PseudoListX.append(twiss(filex))
        PseudoListY.append(twiss(filey))


    return [PseudoListX,PseudoListY]

#----------------------- for finding the lines of the sextupoles (@ Glenn Vanbavinckhove)
def f2h(amp,ampphase,termj,factor,term,M2M):    # converts from f-term to h-term

    # conversion to include f2h
    tp=2.0*pi
    H=(amp/(termj*factor))*(1-e**complex(0,term*tp))

    Ampi=H.imag
    Ampr=H.real
    Amp=abs(H)/M2M
    phase=(atan2(Ampi,Ampr))/tp


    fh=[Ampi,Ampr,Amp,phase]

    return fh

def function(x):


    fterm=4*abs(x[0])

    main=2*sinh(fterm)

    single=sinh(fterm)

    first=cosh(fterm)*sin(x[1])

    second=cosh(fterm)*sin(x[1]+2*phi12_nlinear)

    fun1=float((main*(single+first))-bb[0])

    fun2=float((main*(single+second))-bb[1])

    fun=[fun1,fun2]

    return fun

phi12_nlinear=0
bb=[0,0]

def Getquadrupole(MADTwiss,names,plane,phase,amp):

    f2000_tot={}

    #print MADTwiss

    for bcount in range(len(names)):

        #print bcount

        if bcount==len(names)-1:
            bpm=names[bcount]
            bpm2=names[0]
        else:
            bpm=names[bcount]
            bpm2=names[bcount+1]

#               print bpm,bpm2

        if plane=='H':
            BETM=MADTwiss.BETX[MADTwiss.indx[bpm]]
            BET=amp[bpm][0]
            BETM1=MADTwiss.BETX[MADTwiss.indx[bpm2]]
            BET1=amp[bpm2][0]
            phi12_nlinear=(MADTwiss.PHIX[MADTwiss.indx[bpm]])*2*pi
            tphi01=2*(MADTwiss.MUX[MADTwiss.indx[bpm2]]-MADTwiss.MUX[MADTwiss.indx[bpm]])*2*pi
        else:
            BETM=MADTwiss.BETY[MADTwiss.indx[bpm]]
            BET=amp[bpm][0]
            BETM1=MADTwiss.BETY[MADTwiss.indx[bpm2]]
            BET1=amp[bpm2][0]
            phi12_nlinear=(MADTwiss.PHIY[MADTwiss.indx[bpm]])*2*pi
            tphi01=2*(MADTwiss.MUY[MADTwiss.indx[bpm2]]-MADTwiss.MUY[MADTwiss.indx[bpm]])*2*pi

        #linear

        bb=(BET-BETM)/BETM
        bb1=(BET1-BETM1)/BETM1

        r=bb/bb1
        P2000_linear=atan(r*sin(tphi01)/(1-r*cos(tphi01)))


        F2000_linear=abs(bb/8/sin(P2000_linear))


        # non-linear
        #resul=optimize.fsolve(function,[0.1,0.1])

        #resul[1]=resul[1]%1

        f2000_tot[bpm]=[F2000_linear,P2000_linear,0,0,0,0,0,0,names[bcount]]



    return f2000_tot



def Getsextupole(MADTwiss,amp20list,phase,tune,j,k):
    '''
    function written to calculate resonance driving terms
    '''

    # constructing complex amplitude and phase using two BPM method

    bpms=intersect(amp20list)
    bpms=modelIntersect(bpms,MADTwiss)

    [beta,rmsbb,bpms,invariantJx]=BetaFromAmplitude(MADTwiss,amp20list,'H')
    sqrt2jx=invariantJx[0]

    Q=tune+float(str(MADTwiss.Q1).split(".")[0])

    afactor=(1-cos(2*(j-k)*pi*Q))#(2*sin(pi*(j-k)*Q))
    #print (2*sin(pi*(j-k)*Q)),(1-cos(6*pi*Q))
    #sys.exit()
    pfactor=(pi*(j-k)*Q)

    htot={}

    for i in range(len(bpms)):

        if i<(len(bpms)-1):
            bpm=bpms[i][1]
            bpm1=bpms[i+1][1]
            s=bpms[i][0]
        else:
            bpm=bpms[i][1]
            bpm1=bpms[0][1]
            s=bpms[i][0]

        amp_i_list=[]
        phase_i_list=[]

        hlist=[]
        hplist=[]

        flist=[]
        fplist=[]


        for fileamp in amp20list:

            amp_201=fileamp.AMP_20[fileamp.indx[bpm]]*fileamp.AMPX[fileamp.indx[bpm]]
            amp_202=fileamp.AMP_20[fileamp.indx[bpm1]]*fileamp.AMPX[fileamp.indx[bpm1]]

            phase_201=fileamp.PHASE_20[fileamp.indx[bpm]]
            phase_202=fileamp.PHASE_20[fileamp.indx[bpm1]]

            delta=phase[bpm.upper()][0]-0.25
            edelta=phase[bpm.upper()][1]

            #computing complex line
            ampi,phasei,eampi,ephasei=ComplexSecondaryLineExtended(delta,edelta,amp_201,amp_202,phase_201,phase_202)


            if ampi!=0.0:

                amp_i_list.append(ampi)
                phase_i_list.append(phasei)

                if (j==3 and k==0):
                    factor=sqrt(2)### factor
                    fterm=ampi/(factor*2*j*sqrt2jx**2)
                    pterm=(phasei-phase[bpm.upper()][0]+0.25)%1

                    hterm=fterm/afactor

                    hpterm=(pterm-pfactor)%1


                elif (j==2 and k==1):
                    factor=sqrt(2)### factor
                    fterm=ampi/(factor*2*j*sqrt2jx**2)
                    pterm=(phasei-phase[bpm][0]+0.25)%1

                    hterm=fterm/afactor

                    hpterm=(pterm-pfactor)%1

                flist.append(fterm)
                fplist.append(pterm)
                hlist.append(hterm)
                hplist.append(hpterm)


        if len(amp_i_list)!=0.0:
            al=mean(amp_i_list)
            alstd=std(amp_i_list)

            pl=mean(phase_i_list)
            plstd=mean(phasei)

            fl=mean(flist)
            fstd=std(flist)

            fpl=mean(fplist)
            fpstd=std(fplist)

            hl=mean(hlist)
            hstd=std(hlist)

            hpl=mean(hplist)
            hpstd=std(hplist)


            htot[bpm]=[bpm,s,al,alstd,pl,plstd,fl,fstd,fpl,fpstd,hl,hstd,hpl,hpstd]


    return htot,afactor,pfactor


def Getoctopole(MADTwiss,plane,listF,phaseI,Q,fname,fM,NAMES):
    '''
    for finding secondary lines of the octuple (@ Glenn Vanbavinckhove)
    '''

        # intersects BPMs
    dbpms=intersect(listF[0])
    dbpms=modelIntersect(dbpms,MADTwiss)



    # value definition
    tp=2.0*pi

    hMODELT=[]
    hMODELTi=[]
    hMODELTr=[]
    h_phase_MODELT=[]

    AT=[]
    A_RMST=[]

    phaseT=[]
    phase_RMST=[]

    hT=[]
    hTi=[]
    hTr=[]
    h_RMST=[]

    h_phaseT=[]
    h_phase_RMST=[]

    invarianceJx=[]
    invarianceJy=[]

    # finding the invariances
    for j in range(0,len(listF[0])):
        singleFilex=[listF[0][j]]
        singleFiley=[listF[1][j]]

        [beta,rmsbb,bpms,invariantJx]=BetaFromAmplitude(MADTwiss,singleFilex,'H')

        [beta,rmsbb,bpms,invariantJy]=BetaFromAmplitude(MADTwiss,singleFiley,'V')

        invarianceJx.append(invariantJx)
        invarianceJy.append(invariantJy)


    # for the model
    for i in range(0,len(dbpms)):

        bpm=upper(dbpms[i][1])

        bpmC=MADTwiss.NAME[MADTwiss.indx[bpm]]


        for j in range(0,len(NAMES)):
            try:
                name=NAMES[j]

                if name==bpmC:

                    amp=abs(fM[j])
                    ampr=fM[i].real
                    ampi=fM[j].imag
                    phase=arctan2(ampi,ampr)%1

                    hMODELT.append(amp)
                    hMODELTr.append(ampr)
                    hMODELTi.append(ampi)
                    h_phase_MODELT.append(phase)



            except:
                print 'name '+str(NAMES[j])+' is not found in dictionary'
            hMODEL=[hMODELT,hMODELTi,hMODELTr,h_phase_MODELT]

    #calculation of f,q,h,qh
    for i in range(0,len(dbpms)-1):

        bn1=upper(dbpms[i][1])
        bn2=upper(dbpms[i+1][1])

        #print bn1
        #print phaseT


        dell= phaseI[bn1][0] - 0.25



        # internal value definition
        AS=[]
        A_SRMS=[]
        phaseS=[]
        phase_RMSS=[]

        hS=[]
        hSi=[]
        hSr=[]
        h_RMSS=[]
        h_phaseS=[]
        h_phase_RMSS=[]

        phaseMM=h_phase_MODELT[i]

        for j in range(0,len(listF[0])):

            file=listF[0][j]

            # for f4000
            if fname=='f4000':

                [A,phi]=ComplexSecondaryLine(dell, file.AMP_30[file.indx[bn1]], file.AMP_30[file.indx[bn2]], file.PHASE_30[file.indx[bn1]], file.PHASE_30[file.indx[bn2]])

                factor=float(8*invarianceJx[j][0]**1.5)   # 1 to fit with model
                term=float(4*Q[0])
                termj=4
                M2M=0.5

            #------ converting
            phase0=file.MUX[file.indx[bn1]]
            h=f2h(A,phi,termj,factor,term,M2M)

            #----- adding the terms
            AS.append(A)
            phaseS.append(phi)
            hSi.append(h[0])
            hSr.append(h[1])
            hS.append(h[2])
            h_phaseS.append(h[3])

        # array and taking average for all the input files for one BPM
        AS=array(AS)
        A_SRMS=sqrt(average(AS*AS)-(average(AS))**2+2.2e-16)

        phaseS=array(phaseS)
        try:
            phase_RMSS=sqrt(average(phaseS*phaseS)-(average(phaseS))**2+2.2e-16)
        except:
            phase_RMSS=0

        hS=array(hS)
        hSi=array(hSi)
        hSr=array(hSr)
        try:
            h_RMSS=sqrt(average(hS*hS)-(average(hS))**2+2.2e-16)
        except:
            h_RMSS=0

        h_phaseS=array(h_phaseS)
        try:
            phase_rms=average(h_phaseS*h_phaseS)-(average(h_phaseS))**2+2.2e-16
        except:
            phase_rms=0
        h_phase_RMSS=sqrt(phase_rms)

        # real output
        AT.append(average(AS))
        A_RMST.append(A_SRMS)

        phaseT.append(average(phaseS))
        phase_RMST.append(phase_RMSS)

        hT.append(average(hS))
        hTi.append(average(hSi))
        hTr.append(average(hSr))
        h_RMST.append(h_RMSS)

        h_phaseT.append(average(h_phaseS))
        h_phase_RMST.append(h_phase_RMSS)

        A=[AT,A_RMST,phaseT,phase_RMST]
        h=[hT,hTi,hTr,h_RMST,h_phaseT,h_phase_RMST]




    return [A,h,hMODEL,dbpms]

#------------------------- for finding the chi terms
def computeChiTerms(amp,phase_20,phase,terms,J,plane,ima,rea):

    #computes the chiterms for different inputs
    twoPi=2*pi

    delta1=((phase[1]-phase[0]-0.25)*twoPi)
    delta2=((phase[2]-phase[1]-0.25)*twoPi)

    inp=0.13 # ????
    #term1=((amp[0]*e**complex(0,phase_20[0]*twoPi)))/cos(delta1)
    #term2=((amp[1]*e**complex(0,phase_20[1]*twoPi)))*(tan(delta1)+tan(delta2))
    #term3=((amp[2]*e**complex(0,phase_20[2]*twoPi)))/cos(delta2)
    term1=((amp[0]*e**complex(0,(phase_20[0]+inp)*twoPi)))/cos(delta1)
    term2=((amp[1]*e**complex(0,(phase_20[1]+inp)*twoPi)))*(tan(delta1)+tan(delta2))
    term3=((amp[2]*e**complex(0,(phase_20[2]+inp)*twoPi)))/cos(delta2)
    chiTOT=(term1+term2+term3)

    chiAMP=abs(chiTOT)

    chiAMPi=chiTOT.imag
    chiAMPr=chiTOT.real
    #print chiTOT.imag,chiTOT.real
    chiPHASE=(((arctan2(chiTOT.imag,chiTOT.real)))/twoPi)%1
    #chiPHASE=(0.5-chiPHASE)%1



    JX=J[0]**(2.*(terms[0]+terms[1]-2.)/2.)
    JY=J[1]**(2.*(terms[2]+terms[3])/2.)


    Invariance=JX*JY
    Facot4AMP=Invariance*4/2 # to for conversion complex, invariance = ((2*JX)^(j+k-2)/2)*((2*JY)^(l+m)/2)


    chiAMP=chiAMP/Facot4AMP
    chiAMPi=chiAMPi/Facot4AMP
    chiAMPr=chiAMPr/Facot4AMP

    #print 'measured ima + real '+ str(chiAMPi)+' '+str(chiAMPr) + ' model ima + real '+str(ima)+' '+str(rea)

    return [chiAMP,chiAMPi,chiAMPr,chiPHASE]

def getChiTerms(MADTwiss,filesF,plane,name,ListOfZeroDPPX,ListOfZeroDPPY):

    # bmps
    files=filesF[0]

    dbpms=intersect(files)
    dbpms=modelIntersect(dbpms, MADTwiss)


    # initiliasing variables
    twoPi=2*pi

    XIT=[]
    XITi=[]
    XITr=[]
    XIrmsT=[]
    XI_phase_T=[]
    XI_phaseRMS_T=[]

    POS1=[]
    POS2=[]
    POS3=[]

    XITMODEL=[]
    XITMODELi=[]
    XITMODELr=[]
    XITMODEL_phase=[]

    BPMS=[]
    invarianceJx=[]
    invarianceJy=[]

    for i in range(0,len(dbpms)): # ask rogelio

        bn1=upper(dbpms[i][1])

        BPMS.append(bn1)

    #### invariance
    for j in range(0,len(ListOfZeroDPPX)):
        [betax,rmsbbx,bpms,invariantJX]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPX,'H')
        [betay,rmsbby,bpms,invariantJY]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPY,'V')
        invarianceJx.append(invariantJX[0])
        invarianceJy.append(invariantJY[0])

    #print invarianceJx
    #### model chi
    MADTwiss.chiterms(BPMS)
    if name=='chi3000':
        MODEL=MADTwiss.chi
    elif name=='chi4000':
        MODEL=MADTwiss.chi4000


    for i in range(0,len(MODEL)):

        MODEL[i]=MODEL[i]
        amp=abs(MODEL[i])
        ampi=MODEL[i].imag
        ampr=MODEL[i].real

        if(MODEL[i].real==0. ):

            phase=0

        else:

            phase=arctan2(MODEL[i].imag,MODEL[i].real)%1


        XITMODEL.append(amp)
        XITMODELi.append(ampi)
        XITMODELr.append(ampr)
        XITMODEL_phase.append(phase)

    XIMODEl=[XITMODEL,XITMODELi,XITMODELr,XITMODEL_phase]

    for i in range(0,len(dbpms)-2):

        XI=[]
        XIi=[]
        XIr=[]
        XIrms=[]
        XI_phase=[]
        XI_phaseRMS=[]

        bn1=upper(dbpms[i][1])
        bn2=upper(dbpms[i+1][1])
        bn3=upper(dbpms[i+2][1])

        filej=ListOfZeroDPPX[0]

        pos1=filej.S[filej.indx[bn1]]
        pos2=filej.S[filej.indx[bn2]]
        pos3=filej.S[filej.indx[bn3]]


        POS1.append(pos1)
        POS2.append(pos2)
        POS3.append(pos3)

        imaM=XITMODELi[i]
        realM=XITMODELr[i]

        for j in range(0,len(files)):
            jx=files[j]
            listJX=[jx]




            # for chi3000
            if name=='chi3000':
                phase1=jx.PHASE_20[jx.indx[bn1]]
                phase2=jx.PHASE_20[jx.indx[bn2]]
                phase3=jx.PHASE_20[jx.indx[bn3]]
                phase_SL=[phase1,phase2,phase3]

                terms=[3,0,0,0]
                amp1=jx.AMP_20[jx.indx[bn1]]
                amp2=jx.AMP_20[jx.indx[bn2]]
                amp3=jx.AMP_20[jx.indx[bn3]]
                amp=[amp1,amp2,amp3]

            # for chi4000
            elif name=='chi4000':
                phase1=jx.PHASE_30[jx.indx[bn1]]
                phase2=jx.PHASE_30[jx.indx[bn2]]
                phase3=jx.PHASE_30[jx.indx[bn3]]
                phase_SL=[phase1,phase2,phase3]

                terms=[4,0,0,0]
                amp1=jx.AMP_30[jx.indx[bn1]]
                amp2=jx.AMP_30[jx.indx[bn2]]
                amp3=jx.AMP_30[jx.indx[bn3]]
                amp=[amp1,amp2,amp3]

            phase11=jx.MUX[jx.indx[bn1]]
            phase12=jx.MUX[jx.indx[bn2]]
            phase13=jx.MUX[jx.indx[bn3]]
            phase=[phase11,phase12,phase13]


            J=[ invarianceJx[j],invarianceJy[j]]



            chi=computeChiTerms(amp,phase_SL,phase,terms,J,'H',imaM,realM)



            XI.append(chi[0])
            XIi.append(chi[1])
            XIr.append(chi[2])
            XI_phase.append(chi[3])



        XI=array(XI)
        XIi=array(XIi)
        XIr=array(XIr)
        try:
            XIrms=sqrt(average(XI*XI)-average(XI)**2+2.2e-16)
        except:
            XIrms=0
        XI_phase=array(XI_phase)
        try:
            XI_phaseRMS=sqrt(average(XI_phase*XI_phase)-average(XI_phase)**2+2.2e-16)
        except:
            XI_phaseRMS=0


        XIT.append(average(XI))
        XITi.append(average(XIi))
        XITr.append(average(XIr))
        XIrmsT.append(XIrms)
        XI_phase_T.append(average(XI_phase))
        XI_phaseRMS_T.append(XI_phaseRMS)

        POS=[POS1,POS2,POS3]

        XItot=[XIT,XITi,XITr,XIrmsT,XI_phase_T,XI_phaseRMS_T]

    return [dbpms,POS,XItot,XIMODEl]

def getchi1010(MADTwiss,filesF,plane,name,bn1,ListOfZeroDPPX,ListOfZeroDPPY):

        # bmps
    files=filesF[0]
    filesy=filesF[1]

    dbpms=intersect(files+filesy)
    dbpms=modelIntersect(dbpms, MADTwiss)

    dbpmsy=intersect(filesy+files)
    dbpmsy=modelIntersect(dbpmsy, MADTwiss)


    # initiliasing variables
    twoPi=2*pi

    XIT=[]
    XITi=[]
    XITr=[]
    XIrmsT=[]
    XI_phase_T=[]
    XI_phaseRMS_T=[]

    POS1=[]
    POS2=[]
    POS3=[]

    XITMODEL=[]
    XITMODELi=[]
    XITMODELr=[]
    XITMODEL_phase=[]

    BPMS=[]
    invarianceJx=[]
    invarianceJy=[]

    for i in range(0,len(dbpms)): # ask rogelio


        BPMS.append(bn1)

    #### invariance
    for j in range(0,len(files)):
        [betax,rmsbbx,bpms,invariantJX]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPX,'H')
        [betay,rmsbby,bpms,invariantJY]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPY,'V')
        invarianceJx.append(invariantJX[0])
        invarianceJy.append(invariantJY[0])


    for i in range(0,len(dbpms)):

        XI=[]
        XIrms=[]
        XI_phase=[]
        XI_phaseRMS=[]

        bn=upper(dbpms[i][1])
        bny=upper(dbpmsy[i][1])

        for j in range(0,len(files)):

            jx=files[j]
            listJX=[jx]

            jy=filesy[j]
            listJy=[jy]


            amp10x=jx.AMP01[jx.indx[bn]]
            amp10y=jy.AMP10[jy.indx[bny]]
            phase10x=jx.PHASE01[jx.indx[bn]]
            phasex=jx.MUX[jx.indx[bn]]

            XI1010=0.25*sqrt(amp10x*amp10y)
            phase1010=phase10x+phasex

            XI.append(XI1010)
            XI_phase.append(phase1010)

        XI=array(XI)
        XIrms=sqrt(average(XI*XI)-average(XI)**2+2.2e-16)
        XI_phase=array(XI_phase)
        XI_phaseRMS=sqrt(average(XI_phase*XI_phase)-average(XI_phase)**2+2.2e-16)


        XIT.append(average(XI))
        XIrmsT.append(XIrms)
        XI_phase_T.append(average(XI_phase))
        XI_phaseRMS_T.append(XI_phaseRMS)

        XItot=[XIT,XIrmsT,XI_phase_T,XI_phaseRMS_T]

    return [dbpms,XItot]

#---- construct offmomentum
def ConstructOffMomentumModel(MADTwiss,dpp, dictionary):

    j=MADTwiss
    bpms=intersect([MADTwiss])

    Qx=j.Q1+dpp*j.DQ1
    Qy=j.Q2+dpp*j.DQ2

    ftemp=open("./TempTwiss.dat","w")
    ftemp.write("@ Q1 %le "+str(Qx)+"\n")
    ftemp.write("@ Q2 %le "+str(Qy)+"\n")
    ftemp.write("@ DPP %le "+str(dpp)+"\n")
    ftemp.write("* NAME S BETX BETY ALFX ALFY MUX MUY\n")
    ftemp.write("$ %s %le %le  %le  %le  %le  %le %le\n")


    for i in range(0,len(bpms)):
        bn=upper(bpms[i][1])
        bns=bpms[i][0]

        # dbeta and dalpha will be extract via metaclass. As it is for the time being.
        ax=j.WX[j.indx[bn]]*cos(2.0*pi*j.PHIX[j.indx[bn]])
        bx=j.WX[j.indx[bn]]*sin(2.0*pi*j.PHIX[j.indx[bn]])
        bx1=bx+j.ALFX[j.indx[bn]]*ax
        NBETX=j.BETX[j.indx[bn]]*(1.0+ax*dpp)
        NALFX=j.ALFX[j.indx[bn]]+bx1*dpp
        NMUX=j.MUX[j.indx[bn]]+j.DMUX[j.indx[bn]]*dpp

        ay=j.WY[j.indx[bn]]*cos(2.0*pi*j.PHIY[j.indx[bn]])
        by=j.WY[j.indx[bn]]*sin(2.0*pi*j.PHIY[j.indx[bn]])
        by1=by+j.ALFY[j.indx[bn]]*ay
        NBETY=j.BETY[j.indx[bn]]*(1.0+ay*dpp)
        NALFY=j.ALFY[j.indx[bn]]+by1*dpp
        NMUY=j.MUY[j.indx[bn]]+j.DMUY[j.indx[bn]]*dpp

        ftemp.write('"'+bn+'" '+str(bns)+" "+str(NBETX)+" "+str(NBETY)+" "+str(NALFX)+" "+str(NALFY)+" "+str(NMUX)+" "+str(NMUY)+"\n")

    ftemp.close()
    dpptwiss=twiss("./TempTwiss.dat",dictionary)


    return dpptwiss


#---- finding kick
def getkick(files,MADTwiss):

    invarianceJx=[]
    invarianceJy=[]

    tunex=[]
    tuney=[]

    tunexRMS=[]
    tuneyRMS=[]

    dpp=[]



        # finding the invariances
    for j in range(0,len(files[0])):
        x=files[0][j]
        y=files[1][j]


        [beta,rmsbb,bpms,invariantJx]=BetaFromAmplitude(MADTwiss,[x],'H')

        [beta,rmsbb,bpms,invariantJy]=BetaFromAmplitude(MADTwiss,[y],'V')

        invarianceJx.append(invariantJx)
        invarianceJy.append(invariantJy)

        try:
            dpp.append(x.DPP)
        except:
            dpp.append(0.0)
        tunex.append(x.Q1)
        tuney.append(y.Q2)
        tunexRMS.append(x.Q1RMS)
        tuneyRMS.append(y.Q2RMS)




    tune=[tunex,tuney]
    tuneRMS=[tunexRMS,tuneyRMS]

    return [invarianceJx,invarianceJy,tune,tuneRMS,dpp]

def BPMfinder(IP,model,measured):

    # last index
    indxes=model.S
    #print len(indxes)
    indxlast=len(indxes)-1

    S=model.S[model.indx["IP"+IP]]
    indx=model.indx["IP"+IP]
    bpml="null"
    bpmh="null"
    for ind in range(len(indxes)):
        #print ind
        name=model.NAME[ind]
        if "BPMSW.1L"+IP in name:
            bpml=name
            try:
                test =measured[0][bpml][0]
            except:
                bpml="null"
        if "BPMSW.1R"+IP in name:
            bpmh=name
            try:
                test =measured[0][bpmh][0]
            except:
                bpmh="null"
    return [bpml,bpmh]



def getIP(IP,measured,model,phase,bpms):

    BPMleft,BPMright=BPMfinder(IP,model,measured)

    #print "IN ip"

    if "null" in BPMleft or "null" in BPMright:

        print "skipping IP"+IP+" calculation, no BPM found"
        betahor=[IP,0,0,0,0,0,0]
        betaver=[IP,0,0,0,0,0,0]
        #sys.exit()

    else:
        # model
        sxl=model.S[model.indx[BPMleft]]
        sip=model.S[model.indx["IP"+IP]]
        sxr=model.S[model.indx[BPMright]]
        betaipx=model.BETX[model.indx["IP"+IP]]
        #alxl=model.ALFX[[model.indx[BPMleft]]
        #alyl=model.ALFY[[model.indx[BPMleft]]
        #alxr=model.ALFX[[model.indx[BPMright]]
        #alyr=model.ALFY[[model.indx[BPMright]]
        #print betaipx
        betaipy=model.BETY[model.indx["IP"+IP]]
        #print betaipy

        # measured value
        betaxl=measured[0][BPMleft][0]
        betayl=measured[1][BPMleft][0]

        betaxr=measured[0][BPMright][0]
        betayr=measured[1][BPMright][0]

        deltaphimodel=abs(model.MUX[model.indx[BPMright]]-model.MUX[model.indx[BPMleft]])


        L=((sip-sxl)+(sxr-sip))/2
        betastar=(2*sqrt(betaxl)*sqrt(betaxr)*sin(deltaphimodel*2*pi))/(betayl+betayr-2*sqrt(betaxl)*sqrt(betaxr)*cos(2*pi*deltaphimodel))*L
        location=((betaxl-betaxr)/(betaxl+betaxr-2*sqrt(betaxl)*sqrt(betaxr)*cos(2*pi*deltaphimodel)))*L

        deltaphi=(atan((L-location)/betastar)+atan((L+location)/betastar))/(2*pi)

        betahor=[IP,betastar,location,deltaphi,betaipx,deltaphimodel,0]


        print "horizontal betastar for ",IP," is ",str(betastar)," at location ",str(location), " of IP center with phase advance ",str(deltaphi)

        #vertical
        deltaphimodel=abs(model.MUY[model.indx[BPMright]]-model.MUY[model.indx[BPMleft]])


        betastar=(2*sqrt(betayl)*sqrt(betayr)*sin(deltaphimodel*2*pi))/(betayl+betayr-2*sqrt(betayl)*sqrt(betayr)*cos(2*pi*deltaphimodel))*L
        location=((betayl-betayr)/(betayl+betayr-2*sqrt(betayl)*sqrt(betayr)*cos(2*pi*deltaphimodel)))*L

        deltaphi=(atan((L-location)/betastar)+atan((L+location)/betastar))/(2*pi)

        betaver=[IP,betastar,location,deltaphi,betaipy,deltaphimodel,0]

        print "vertical betastar for ",IP," is ",str(betastar)," at location ",str(location), " of IP center with phase advance ",str(deltaphi)


    return [betahor,betaver]

def GetIP2(MADTwiss,Files,Q,plane,bd,oa,op):

    #-- Common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],upper(b[1])) for b in bpm]

    #-- Loop for IPs
    result={}
    for ip in ('1','2','5','8'):

        bpml='BPMSW.1L'+ip+'.'+oa[3:]
        bpmr='BPMSW.1R'+ip+'.'+oa[3:]

        if (bpml in zip(*bpm)[1]) and (bpmr in zip(*bpm)[1]):

            #-- Model values
            L=0.5*(MADTwiss.S[MADTwiss.indx[bpmr]]-MADTwiss.S[MADTwiss.indx[bpml]])
            if L<0: L+=0.5*MADTwiss.LENGTH  #-- For sim starting in the middle of an IP
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
                try:
                    if plane=='H':
                        al=Files[i].AMPX[Files[i].indx[bpml]]
                        ar=Files[i].AMPX[Files[i].indx[bpmr]]
                        if list(zip(*bpm)[1]).index(bpmr)>list(zip(*bpm)[1]).index(bpml):
                            dpsi=2*pi*bd*(Files[i].MUX[Files[i].indx[bpmr]]-Files[i].MUX[Files[i].indx[bpml]])
                        else:
                            dpsi=2*pi*(Q+bd*(Files[i].MUX[Files[i].indx[bpmr]]-Files[i].MUX[Files[i].indx[bpml]]))
                        #-- To compensate the phase shift by tune
                        if op=='1':
                            if (bd==1 and ip=='2') or (bd==-1 and ip=='8'): dpsi+=2*pi*Q
                    if plane=='V':
                        al=Files[i].AMPY[Files[i].indx[bpml]]
                        ar=Files[i].AMPY[Files[i].indx[bpmr]]
                        if list(zip(*bpm)[1]).index(bpmr)>list(zip(*bpm)[1]).index(bpml):
                            dpsi=2*pi*bd*(Files[i].MUY[Files[i].indx[bpmr]]-Files[i].MUY[Files[i].indx[bpml]])
                        else:
                            dpsi=2*pi*(Q+bd*(Files[i].MUY[Files[i].indx[bpmr]]-Files[i].MUY[Files[i].indx[bpml]]))
                        #-- To compensate the phase shift by tune
                        if op=='1':
                            if (bd==1 and ip=='2') or (bd==-1 and ip=='8'): dpsi+=2*pi*Q

                    #-- bet, alf, and sqrt(2J) from amp and phase advance
                    bet =L*(al**2+ar**2+2*al*ar*cos(dpsi))/(2*al*ar*sin(dpsi))
                    alf =(al**2-ar**2)/(2*al*ar*sin(dpsi))
                    bets=bet/(1+alf**2)
                    ds  =alf*bets
                    rt2J=sqrt(al*ar*sin(dpsi)/(2*L))
                    betall.append(bet)
                    alfall.append(alf)
                    betsall.append(bets)
                    dsall.append(ds)
                    rt2Jall.append(rt2J)
                except:
                    pass

            #-- Ave and Std
            betall =array(betall) ; betave =mean(betall) ; betstd =sqrt(mean((betall-betave)**2))
            alfall =array(alfall) ; alfave =mean(alfall) ; alfstd =sqrt(mean((alfall-alfave)**2))
            betsall=array(betsall); betsave=mean(betsall); betsstd=sqrt(mean((betsall-betsave)**2))
            dsall  =array(dsall)  ; dsave  =mean(dsall)  ; dsstd  =sqrt(mean((dsall-dsave)**2))
            rt2Jall=array(rt2Jall); rt2Jave=mean(rt2Jall); rt2Jstd=sqrt(mean((rt2Jall-rt2Jave)**2))
            result['IP'+ip]=[betave,betstd,betmdl,alfave,alfstd,alfmdl,betsave,betsstd,betsmdl,dsave,dsstd,dsmdl,rt2Jave,rt2Jstd]

    return result

def GetIPFromPhase(MADTwiss,psix,psiy,oa):

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
                betx    =L/tan(pi*dpsix)
                bety    =L/tan(pi*dpsiy)
                betxstd =L*pi*dpsixstd/(2*sin(pi*dpsix)**2)
                betystd =L*pi*dpsiystd/(2*sin(pi*dpsiy)**2)
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
                betx    =L/tan(pi*dpsix)
                bety    =L/tan(pi*dpsiy)
                betxstd =L*pi*dpsixstd/(2*sin(pi*dpsix)**2)
                betystd =L*pi*dpsiystd/(2*sin(pi*dpsiy)**2)
                result['IP'+i]=[2*L,betx,betxstd,betxmdl,bety,betystd,betymdl,dpsix,dpsixstd,dpsixmdl,dpsiy,dpsiystd,dpsiymdl]
        except: pass

    return result

def getCandGammaQmin(fqwq,bpms,tunex,tuney,twiss):

    QQ1=float(str(twiss.Q1).split('.')[0])
    QQ2=float(str(twiss.Q2).split('.')[0])

    tunex=float(tunex)+QQ1
    tuney=float(tuney)+QQ2

    tunefactor=(cos(2*pi*tunex)-cos(2*pi*tuney))/(pi*(sin(2*pi*tunex)+sin(2*pi*tuney)))

    coupleterms={}
    Qmin=[]


    if len(bpms)==0:
        print "No bpms in getCandGammaQmin. Returning emtpy stuff"
        return coupleterms,0,0,bpms

    for bpm in bpms:

        bpmm=bpm[1].upper()

        detC=1-(1/(1+4*(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2)))


        check2=0.25+abs(fqwq[bpmm][0][0])**2


        if check2>abs(fqwq[bpmm][0][2])**2: # checking if sum or difference resonance is dominant!
            gamma=sqrt(1/(1/(1+4*(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2))))
            #print detC
            ffactor= 2*gamma*tunefactor*sqrt(abs(detC)) # cannot take abs
            C11=-(fqwq[bpmm][0][0].imag-fqwq[bpmm][0][2].imag)*2*gamma
            C12=-(fqwq[bpmm][0][0].real+fqwq[bpmm][0][2].real)*2*gamma
            C21=(fqwq[bpmm][0][0].real+fqwq[bpmm][0][2].real)*2*gamma
            C22=(fqwq[bpmm][0][0].imag-fqwq[bpmm][0][2].imag)*2*gamma
        else: # negative gamma
            gamma=-1
            ffactor=-1
            C11=C12=C21=C22=-1


        Qmin.append(ffactor)

        if (abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2)>0.0:

            err=(2*((abs(fqwq[bpmm][0][1])*abs(fqwq[bpmm][0][0]))+(abs(fqwq[bpmm][0][3])*abs(fqwq[bpmm][0][2]))))/(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2)

        else:
            err=-1


        coupleterms[bpmm]=[detC,err,gamma,err,C11,C12,C21,C22]

    if gamma==-1:
        print "WARN: Sum resonance is dominant! "

    Qmin=array(Qmin)

    Qminerr=sqrt(average(Qmin*Qmin)-(average(Qmin))**2+2.2e-16)
    Qminav=average(Qmin)




    return coupleterms,Qminav,Qminerr,bpms






###### ac-dipole stuff

def getFreeBeta(modelfree,modelac,betal,rmsbb,alfal,bpms,plane): # to check "+"

    # ! how to deal with error !

    #print "Calculating free beta using model"

    bpms=modelIntersect(bpms, modelfree)
    bpms=modelIntersect(bpms, modelac)
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

    #sys.exit()
    return betan,rmsbb,alfan,bpms

def getFreeCoupling(tunefreex,tunefreey,tunedrivenx,tunedriveny,fterm,twiss,bpms):

    print "Calculating free fterms"
    couple={}
    couple['Global']=[fterm['Global'][0],fterm['Global'][1]]

    QQ1=float(str(twiss.Q1).split('.')[0])
    QQ2=float(str(twiss.Q2).split('.')[0])

    if(tunefreey>0.50):
        tunefreey=1-tunefreey
        tunefreey=abs(QQ2+tunefreey)
    else:
        tunefreey=abs(QQ2+abs(tunefreey))
    if(tunefreex>0.50):
        tunefreex=1-float(tunefreex)
        tunefreex=abs(QQ1+tunefreex)
    else:
        tunefreex=abs(QQ1+abs(tunefreex))

    if(tunedrivenx>0.50):
        tunedrivenx=1-tunedrivenx
    if(tunedriveny>0.50):
        tunedriveny=1-tunedriveny

    tunedrivenx=abs(QQ1+abs(tunedrivenx))
    tunedriveny=abs(QQ2+abs(tunedriveny))


    # diff f1001
    factor_top_diff=sqrt(sin(pi*(tunedrivenx-tunefreey))*sin(pi*(tunefreex-tunedriveny)))
    factor_bottom_diff=sin(pi*(tunefreex-tunefreey))

    factor_diff=abs((factor_top_diff/factor_bottom_diff))

    print "Factor for coupling diff ",factor_diff

    # sum f1010
    factor_top_sum=sqrt(sin(pi*(tunedrivenx+tunefreey))*sin(pi*(tunefreex+tunedriveny)))
    factor_bottom_sum=sin(pi*(tunefreex+tunefreey))

    factor_sum=abs((factor_top_sum/factor_bottom_sum))

    print "Factor for coupling sum ",factor_sum

    for bpm in bpms:

        bpmm=bpm[1].upper()
        [amp,phase]=fterm[bpmm]

        #print amp[2]

        ampp=[amp[0]*factor_diff,amp[1],amp[2]*factor_sum,amp[3]]
        pphase=[phase[0]*factor_diff,phase[1],phase[2]*factor_sum,phase[3]]

        couple[bpmm]=[ampp,pphase]

    return couple,bpms

def getFreeAmpBeta(betai,rmsbb,bpms,invJ,modelac,modelfree,plane): # "-"

    #
    # Why difference in betabeta calculation ??
    #
    #

    betas={}

    #print "Calculating free beta from amplitude using model"

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

        #print beta,beta*(1+bb)

        betas[bpmm]=[beta*(1+bb),betai[bpmm][1],betai[bpmm][2]]

    return betas,rmsbb,bpms,invJ

def getfreephase(phase,Qac,Q,bpms,MADTwiss_ac,MADTwiss,plane):

    #print "Calculating free phase using model"

    phasef={}
    phi=[]

    for bpm in bpms:

        s=bpm[0]
        bn1=bpm[1].upper()

        phi12,phstd12,phi13,phstd13,phmdl12,phmdl13,bn2=phase[bn1]
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

def getfreephaseTotal(phase,bpms,plane,MADTwiss,MADTwiss_ac):

    #print "Calculating free total phase using model"

    first=bpms[0][1]

    phasef={}

    for bpm in bpms:
        s=bpm[0]
        bn2=bpm[1].upper()

        if plane=="H":

            ph_ac_m=(MADTwiss_ac.MUX[MADTwiss_ac.indx[bn2]]-MADTwiss_ac.MUX[MADTwiss_ac.indx[first]])%1
            ph_m=(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[first]])%1

        else:
            ph_ac_m=(MADTwiss_ac.MUY[MADTwiss_ac.indx[bn2]]-MADTwiss_ac.MUY[MADTwiss_ac.indx[first]])%1
            ph_m=(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[first]])%1

        phi12,phstd12,phmdl12,bn1=phase[bn2]


        phi12=phi12-(ph_ac_m-ph_m)
        phstd12=phstd12

        phasef[bn2]=phi12,phstd12,ph_m

    return phasef,bpms

# free coupling from equations

def GetFreeIP2(MADTwiss,MADTwiss_ac,IP,plane,oa):

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

#---------  The following is functions to compensate the AC dipole effect based on analytic formulae (by R. Miyamoto)

def GetACPhase_AC2BPMAC(MADTwiss,Qd,Q,plane,oa):
    if   oa=='LHCB1':
        bpmac1='BPMYA.5L4.B1'
        bpmac2='BPMYB.6L4.B1'
    elif oa=='LHCB2':
        bpmac1='BPMYB.5L4.B2'
        bpmac2='BPMYA.6L4.B2'
    else:
        return {}

    if plane=='H':
        psi_ac2bpmac1=MADTwiss.MUX[MADTwiss.indx[bpmac1]]-MADTwiss.MUX[MADTwiss.indx['MKQA.6L4.'+oa[3:]]]  #-- B1 direction for B2
        psi_ac2bpmac2=MADTwiss.MUX[MADTwiss.indx[bpmac2]]-MADTwiss.MUX[MADTwiss.indx['MKQA.6L4.'+oa[3:]]]  #-- B1 direction for B2
    if plane=='V':
        psi_ac2bpmac1=MADTwiss.MUY[MADTwiss.indx[bpmac1]]-MADTwiss.MUY[MADTwiss.indx['MKQA.6L4.'+oa[3:]]]  #-- B1 direction for B2
        psi_ac2bpmac2=MADTwiss.MUY[MADTwiss.indx[bpmac2]]-MADTwiss.MUY[MADTwiss.indx['MKQA.6L4.'+oa[3:]]]  #-- B1 direction for B2

    r=sin(pi*(Qd-Q))/sin(pi*(Qd+Q))
    psid_ac2bpmac1=arctan((1+r)/(1-r)*tan(2*pi*psi_ac2bpmac1-pi*Q))%pi-pi+pi*Qd
    psid_ac2bpmac2=arctan((1+r)/(1-r)*tan(2*pi*psi_ac2bpmac2+pi*Q))%pi-pi*Qd

    return {bpmac1:psid_ac2bpmac1,bpmac2:psid_ac2bpmac2}


def GetFreePhaseTotal_Eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,op):

    #-- Select common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],upper(b[1])) for b in bpm]

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

    #-- Model phase advances
    if plane=='H': psimdl=array([(MADTwiss.MUX[MADTwiss.indx[b[1]]]-MADTwiss.MUX[MADTwiss.indx[bpm[0][1]]])%1 for b in bpm])
    if plane=='V': psimdl=array([(MADTwiss.MUY[MADTwiss.indx[b[1]]]-MADTwiss.MUY[MADTwiss.indx[bpm[0][1]]])%1 for b in bpm])

    #-- Global parameters of the driven motion
    r=sin(pi*(Qd-Q))/sin(pi*(Qd+Q))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    psiall=zeros((len(bpm),len(Files)))
    for i in range(len(Files)):
        if plane=='H': psid=bd*2*pi*array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid=bd*2*pi*array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*pi*Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*pi*Qd
        Psi=numpy.arctan((1-r)/(1+r)*numpy.tan(Psid))%pi
        for k in range(len(bpm)):
            if Psid[k]%(2*pi)>pi: Psi[k]=Psi[k]+pi
        psi=Psi-Psi[0]
        psi[k_bpmac:]=psi[k_bpmac:]+2*pi*Q
        for k in range(len(bpm)): psiall[k][i]=psi[k]/(2*pi)  #-- phase range back to [0,1)

    #-- Output
    result={}
    for k in range(len(bpm)):
        psiave=PhaseMean(psiall[k],1)
        psistd=PhaseStd(psiall[k],1)
        result[bpm[k][1]]=[psiave,psistd,psimdl[k],bpm[0][1]]

    return [result,bpm]


def GetFreePhase_Eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,op):

    #-- Select common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],upper(b[1])) for b in bpm]

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
            print 'WARN: BPMs next to AC dipoles missing. AC dipole effects not calculated for '+plane+' with eqs !'
            return [{},'',[]]

    #-- Model phase advances
    if plane=='H': psimdl=array([MADTwiss.MUX[MADTwiss.indx[b[1]]] for b in bpm])
    if plane=='V': psimdl=array([MADTwiss.MUY[MADTwiss.indx[b[1]]] for b in bpm])
    psi12mdl=(append(psimdl[1:],psimdl[0] +Q)-psimdl)%1
    psi13mdl=(append(psimdl[2:],psimdl[:2]+Q)-psimdl)%1

    #-- Global parameters of the driven motion
    r=sin(pi*(Qd-Q))/sin(pi*(Qd+Q))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    psi12all=zeros((len(bpm),len(Files)))
    psi13all=zeros((len(bpm),len(Files)))
    for i in range(len(Files)):
        if plane=='H': psid=bd*2*pi*array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid=bd*2*pi*array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*pi*Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*pi*Qd
        Psi=numpy.arctan((1-r)/(1+r)*numpy.tan(Psid))%pi
        for k in range(len(bpm)):
            if Psid[k]%(2*pi)>pi: Psi[k]=Psi[k]+pi
        psi=Psi-Psi[0]
        psi[k_bpmac:]=psi[k_bpmac:]+2*pi*Q
        psi12=(append(psi[1:],psi[0] +2*pi*Q)-psi)/(2*pi)  #-- phase range back to [0,1)
        #psi12=(append(psi[1:],psi[0])-psi)/(2*pi)  #-- phase range back to [0,1)
        psi13=(append(psi[2:],psi[:2]+2*pi*Q)-psi)/(2*pi)  #-- phase range back to [0,1)
        for k in range(len(bpm)):
            psi12all[k][i]=psi12[k]
            psi13all[k][i]=psi13[k]

    #-- Output
    result={}
    muave=0.0  #-- mu is the same as psi but w/o mod
    for k in range(len(bpm)):
        psi12ave=PhaseMean(psi12all[k],1)
        psi12std=PhaseStd(psi12all[k],1)
        psi13ave=PhaseMean(psi13all[k],1)
        psi13std=PhaseStd(psi13all[k],1)
        muave=muave+psi12ave
        try:    result[bpm[k][1]]=[psi12ave,psi12std,psi13ave,psi13std,psi12mdl[k],psi13mdl[k],bpm[k+1][1]]
        except: result[bpm[k][1]]=[psi12ave,psi12std,psi13ave,psi13std,psi12mdl[k],psi13mdl[k],bpm[0][1]]    #-- The last BPM

    return [result,muave,bpm]


def GetFreeBetaFromAmp_Eq(MADTwiss_ac,Files,Qd,Q,psid_ac2bpmac,plane,bd,op):

    #-- Select common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss_ac)
    bpm=[(b[0],upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and bd== 1: s_lastbpm=MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L2.B1']]
    if op=="1" and bd==-1: s_lastbpm=MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L8.B2']]

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
            return [{},'',[],[]]

    #-- Model beta and phase advance
    if plane=='H': betmdl=array([MADTwiss_ac.BETX[MADTwiss_ac.indx[b[1]]] for b in bpm])
    if plane=='V': betmdl=array([MADTwiss_ac.BETY[MADTwiss_ac.indx[b[1]]] for b in bpm])

    #-- Global parameters of the driven motion
    r=sin(pi*(Qd-Q))/sin(pi*(Qd+Q))

    #-- Loop for files
    betall=zeros((len(bpm),len(Files)))
    Adall=zeros((len(bpm),len(Files)))
    for i in range(len(Files)):
        if plane=='H':
            amp =array([2*Files[i].AMPX[Files[i].indx[b[1]]] for b in bpm])
            psid=bd*2*pi*array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V':
            amp =array([2*Files[i].AMPY[Files[i].indx[b[1]]] for b in bpm])
            psid=bd*2*pi*array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*pi*Qd  #-- To fix the phase shift by Q
            except: pass
        Ad  =amp/map(sqrt,betmdl)
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*pi*Qd
        bet =(amp/mean(Ad))**2*(1+r**2+2*r*numpy.cos(2*Psid))/(1-r**2)
        for k in range(len(bpm)):
            betall[k][i]=bet[k]
            Adall[k][i]=Ad[k]

    #-- Output
    result={}
    bb=[]
    Adave=[]
    for k in range(len(bpm)):
        betave=mean(betall[k])
        betstd=sqrt(mean((betall[k]-betave)**2))
        bb.append((betave-betmdl[k])/betmdl[k])
        Adave.append(mean(Adall[k]))
        result[bpm[k][1]]=[betave,betstd,bpm[k][0]]
    bb=sqrt(mean(array(bb)**2))
    Ad=[mean(Adave),sqrt(mean((Adave-mean(Adave))**2))]

    return [result,bb,bpm,Ad]


def GetFreeCoupling_Eq(MADTwiss,FilesX,FilesY,Qh,Qv,Qx,Qy,psih_ac2bpmac,psiv_ac2bpmac,bd):

    #-- Details of this algorithms is in http://www.agsrhichome.bnl.gov/AP/ap_notes/ap_note_410.pdf

    #-- Check linx/liny files, may be redundant
    if len(FilesX)!=len(FilesY): return [{},[]]

    #-- Select common BPMs
    bpm=modelIntersect(intersect(FilesX+FilesY),MADTwiss)
    bpm=[(b[0],upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    #if op=="1" and bd== 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    #if op=="1" and bd==-1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    #-- Determine the BPM closest to the AC dipole and its position
    for b in psih_ac2bpmac.keys():
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
            print 'WARN: BPMs next to AC dipoles missing. AC dipole effects not calculated with analytic eqs for coupling'
            return [{},[]]

    #-- Global parameters of the driven motion
    dh =Qh-Qx
    dv =Qv-Qy
    rh =sin(pi*(Qh-Qx))/sin(pi*(Qh+Qx))
    rv =sin(pi*(Qv-Qy))/sin(pi*(Qv+Qy))
    rch=sin(pi*(Qh-Qy))/sin(pi*(Qh+Qy))
    rcv=sin(pi*(Qx-Qv))/sin(pi*(Qx+Qv))

    #-- Loop for files
    f1001Abs =zeros((len(bpm),len(FilesX)))
    f1010Abs =zeros((len(bpm),len(FilesX)))
    f1001xArg=zeros((len(bpm),len(FilesX)))
    f1001yArg=zeros((len(bpm),len(FilesX)))
    f1010xArg=zeros((len(bpm),len(FilesX)))
    f1010yArg=zeros((len(bpm),len(FilesX)))
    for i in range(len(FilesX)):

        #-- Read amplitudes and phases
        amph  =     array([FilesX[i].AMPX[FilesX[i].indx[b[1]]]    for b in bpm])
        ampv  =     array([FilesY[i].AMPY[FilesY[i].indx[b[1]]]    for b in bpm])
        amph01=     array([FilesX[i].AMP01[FilesX[i].indx[b[1]]]   for b in bpm])
        ampv10=     array([FilesY[i].AMP10[FilesY[i].indx[b[1]]]   for b in bpm])
        psih  =2*pi*array([FilesX[i].MUX[FilesX[i].indx[b[1]]]     for b in bpm])
        psiv  =2*pi*array([FilesY[i].MUY[FilesY[i].indx[b[1]]]     for b in bpm])
        psih01=2*pi*array([FilesX[i].PHASE01[FilesX[i].indx[b[1]]] for b in bpm])
        psiv10=2*pi*array([FilesY[i].PHASE10[FilesY[i].indx[b[1]]] for b in bpm])
        #-- I'm not sure this is correct for the coupling so I comment out this part for now (by RM 9/30/11).
        #for k in range(len(bpm)):
        #       try:
        #               if bpm[k][0]>s_lastbpm:
        #                       psih[k]  +=bd*2*pi*Qh  #-- To fix the phase shift by Qh
        #                       psiv[k]  +=bd*2*pi*Qv  #-- To fix the phase shift by Qv
        #                       psih01[k]+=bd*2*pi*Qv  #-- To fix the phase shift by Qv
        #                       psiv10[k]+=bd*2*pi*Qh  #-- To fix the phase shift by Qh
        #       except: pass

        #-- Construct Fourier components
        #   * be careful for that the note is based on x+i(alf*x*bet*x')).
        #   * Calculating Eqs (87)-(92) by using Eqs (47) & (48) (but in the Fourier space) in the note.
        #   * Note that amph(v)01 is normalized by amph(v) and it is un-normalized in the following.
        dpsih  =append(psih[1:]  ,2*pi*Qh+psih[0]  )-psih
        dpsiv  =append(psiv[1:]  ,2*pi*Qv+psiv[0]  )-psiv
        dpsih01=append(psih01[1:],2*pi*Qv+psih01[0])-psih01
        dpsiv10=append(psiv10[1:],2*pi*Qh+psiv10[0])-psiv10

        X_m10=2*amph*numpy.exp(-1j*psih)
        Y_0m1=2*ampv*numpy.exp(-1j*psiv)
        X_0m1=amph*numpy.exp(-1j*psih01)/(1j*numpy.sin(dpsih))*(amph01*numpy.exp(1j*dpsih)-append(amph01[1:],amph01[0])*numpy.exp(-1j*dpsih01))
        X_0p1=amph*numpy.exp( 1j*psih01)/(1j*numpy.sin(dpsih))*(amph01*numpy.exp(1j*dpsih)-append(amph01[1:],amph01[0])*numpy.exp( 1j*dpsih01))
        Y_m10=ampv*numpy.exp(-1j*psiv10)/(1j*numpy.sin(dpsiv))*(ampv10*numpy.exp(1j*dpsiv)-append(ampv10[1:],ampv10[0])*numpy.exp(-1j*dpsiv10))
        Y_p10=ampv*numpy.exp( 1j*psiv10)/(1j*numpy.sin(dpsiv))*(ampv10*numpy.exp(1j*dpsiv)-append(ampv10[1:],ampv10[0])*numpy.exp( 1j*dpsiv10))

        #-- Construct f1001hv, f1001vh, f1010hv (these include sqrt(betv/beth) or sqrt(beth/betv))
        f1001hv=-conjugate(1/(2j)*Y_m10/X_m10)  #-- - sign from the different def
        f1001vh=-1/(2j)*X_0m1/Y_0m1             #-- - sign from the different def
        f1010hv=-1/(2j)*Y_p10/conjugate(X_m10)  #-- - sign from the different def
        f1010vh=-1/(2j)*X_0p1/conjugate(Y_0m1)  #-- - sign from the different def
##              f1001hv=conjugate(1/(2j)*Y_m10/X_m10)
##              f1001vh=1/(2j)*X_0m1/Y_0m1
##              f1010hv=1/(2j)*Y_p10/conjugate(X_m10)
##              f1010vh=1/(2j)*X_0p1/conjugate(Y_0m1)

        #-- Construct phases psih, psiv, Psih, Psiv w.r.t. the AC dipole
        psih=psih-(psih[k_bpmac]-psih_ac2bpmac[bpmac])
        psiv=psiv-(psiv[k_bpmac]-psiv_ac2bpmac[bpmac])

        Psih=psih-pi*Qh
        Psih[:k_bpmac]=Psih[:k_bpmac]+2*pi*Qh
        Psiv=psiv-pi*Qv
        Psiv[:k_bpmac]=Psiv[:k_bpmac]+2*pi*Qv

        Psix=numpy.arctan((1-rh)/(1+rh)*numpy.tan(Psih))%pi
        Psiy=numpy.arctan((1-rv)/(1+rv)*numpy.tan(Psiv))%pi
        for k in range(len(bpm)):
            if Psih[k]%(2*pi)>pi: Psix[k]=Psix[k]+pi
            if Psiv[k]%(2*pi)>pi: Psiy[k]=Psiy[k]+pi

        psix=Psix-pi*Qx
        psix[k_bpmac:]=psix[k_bpmac:]+2*pi*Qx
        psiy=Psiy-pi*Qy
        psiy[k_bpmac:]=psiy[k_bpmac:]+2*pi*Qy

        #-- Construct f1001h, f1001v, f1010h, f1010v (these include sqrt(betv/beth) or sqrt(beth/betv))
        f1001h=1/numpy.sqrt(1-rv**2)*(numpy.exp(-1j*(Psiv-Psiy))*f1001hv+rv*numpy.exp( 1j*(Psiv+Psiy))*f1010hv)
        f1010h=1/numpy.sqrt(1-rv**2)*(numpy.exp( 1j*(Psiv-Psiy))*f1010hv+rv*numpy.exp(-1j*(Psiv+Psiy))*f1001hv)
        f1001v=1/numpy.sqrt(1-rh**2)*(numpy.exp( 1j*(Psih-Psix))*f1001vh+rh*numpy.exp(-1j*(Psih+Psix))*conjugate(f1010vh))
        f1010v=1/numpy.sqrt(1-rh**2)*(numpy.exp( 1j*(Psih-Psix))*f1010vh+rh*numpy.exp(-1j*(Psih+Psix))*conjugate(f1001vh))

        #-- Construct f1001 and f1010 from h and v BPMs (these include sqrt(betv/beth) or sqrt(beth/betv))
        g1001h          =numpy.exp(-1j*((psih-psih[k_bpmac])-(psiy-psiy[k_bpmac])))*(ampv/amph*amph[k_bpmac]/ampv[k_bpmac])*f1001h[k_bpmac]
        g1001h[:k_bpmac]=1/(numpy.exp(2*pi*1j*(Qh-Qy))-1)*(f1001h-g1001h)[:k_bpmac]
        g1001h[k_bpmac:]=1/(1-numpy.exp(-2*pi*1j*(Qh-Qy)))*(f1001h-g1001h)[k_bpmac:]

        g1010h          =numpy.exp(-1j*((psih-psih[k_bpmac])+(psiy-psiy[k_bpmac])))*(ampv/amph*amph[k_bpmac]/ampv[k_bpmac])*f1010h[k_bpmac]
        g1010h[:k_bpmac]=1/(numpy.exp(2*pi*1j*(Qh+Qy))-1)*(f1010h-g1010h)[:k_bpmac]
        g1010h[k_bpmac:]=1/(1-numpy.exp(-2*pi*1j*(Qh+Qy)))*(f1010h-g1010h)[k_bpmac:]

        g1001v          =numpy.exp(-1j*((psix-psix[k_bpmac])-(psiv-psiv[k_bpmac])))*(amph/ampv*ampv[k_bpmac]/amph[k_bpmac])*f1001v[k_bpmac]
        g1001v[:k_bpmac]=1/(numpy.exp(2*pi*1j*(Qx-Qv))-1)*(f1001v-g1001v)[:k_bpmac]
        g1001v[k_bpmac:]=1/(1-numpy.exp(-2*pi*1j*(Qx-Qv)))*(f1001v-g1001v)[k_bpmac:]

        g1010v          =numpy.exp(-1j*((psix-psix[k_bpmac])+(psiv-psiv[k_bpmac])))*(amph/ampv*ampv[k_bpmac]/amph[k_bpmac])*f1010v[k_bpmac]
        g1010v[:k_bpmac]=1/(numpy.exp(2*pi*1j*(Qx+Qv))-1)*(f1010v-g1010v)[:k_bpmac]
        g1010v[k_bpmac:]=1/(1-numpy.exp(-2*pi*1j*(Qx+Qv)))*(f1010v-g1010v)[k_bpmac:]

        f1001x=numpy.exp(1j*(psih-psix))*f1001h
        f1001x=f1001x-rh*numpy.exp(-1j*(psih+psix))/rch*conjugate(f1010h)
        f1001x=f1001x-2j*numpy.sin(pi*dh)*numpy.exp(1j*(Psih-Psix))*g1001h
        f1001x=f1001x-2j*numpy.sin(pi*dh)*numpy.exp(-1j*(Psih+Psix))/rch*conjugate(g1010h)
        f1001x=1/numpy.sqrt(1-rh**2)*numpy.sin(pi*(Qh-Qy))/numpy.sin(pi*(Qx-Qy))*f1001x

        f1010x=numpy.exp(1j*(psih-psix))*f1010h
        f1010x=f1010x-rh*numpy.exp(-1j*(psih+psix))*rch*conjugate(f1001h)
        f1010x=f1010x-2j*numpy.sin(pi*dh)*numpy.exp(1j*(Psih-Psix))*g1010h
        f1010x=f1010x-2j*numpy.sin(pi*dh)*numpy.exp(-1j*(Psih+Psix))*rch*conjugate(g1001h)
        f1010x=1/numpy.sqrt(1-rh**2)*numpy.sin(pi*(Qh+Qy))/numpy.sin(pi*(Qx+Qy))*f1010x

        f1001y=numpy.exp(-1j*(psiv-psiy))*f1001v
        f1001y=f1001y+rv*numpy.exp(1j*(psiv+psiy))/rcv*f1010v
        f1001y=f1001y+2j*numpy.sin(pi*dv)*numpy.exp(-1j*(Psiv-Psiy))*g1001v
        f1001y=f1001y-2j*numpy.sin(pi*dv)*numpy.exp(1j*(Psiv+Psiy))/rcv*g1010v
        f1001y=1/numpy.sqrt(1-rv**2)*numpy.sin(pi*(Qx-Qv))/numpy.sin(pi*(Qx-Qy))*f1001y

        f1010y=numpy.exp(1j*(psiv-psiy))*f1010v
        f1010y=f1010y+rv*numpy.exp(-1j*(psiv+psiy))*rcv*f1001v
        f1010y=f1010y-2j*numpy.sin(pi*dv)*numpy.exp(1j*(Psiv-Psiy))*g1010v
        f1010y=f1010y+2j*numpy.sin(pi*dv)*numpy.exp(-1j*(Psiv+Psiy))*rcv*g1001v
        f1010y=1/numpy.sqrt(1-rv**2)*numpy.sin(pi*(Qx+Qv))/numpy.sin(pi*(Qx+Qy))*f1010y

        #-- For B2, must be double checked
        if bd == -1:
            f1001x=-conjugate(f1001x)
            f1001y=-conjugate(f1001y)
            f1010x=-conjugate(f1010x)
            f1010y=-conjugate(f1010y)

        #-- Separate to amplitudes and phases, amplitudes averaged to cancel sqrt(betv/beth) and sqrt(beth/betv)
        for k in range(len(bpm)):
            f1001Abs[k][i] =numpy.sqrt(abs(f1001x[k]*f1001y[k]))
            f1010Abs[k][i] =numpy.sqrt(abs(f1010x[k]*f1010y[k]))
            f1001xArg[k][i]=angle(f1001x[k])%(2*pi)
            f1001yArg[k][i]=angle(f1001y[k])%(2*pi)
            f1010xArg[k][i]=angle(f1010x[k])%(2*pi)
            f1010yArg[k][i]=angle(f1010y[k])%(2*pi)

    #-- Output
    fwqw={}
    goodbpm=[]
    for k in range(len(bpm)):

        #-- Bad BPM flag based on phase
        badbpm=0
        f1001xArgAve=PhaseMean(f1001xArg[k],2*pi)
        f1001yArgAve=PhaseMean(f1001yArg[k],2*pi)
        f1010xArgAve=PhaseMean(f1010xArg[k],2*pi)
        f1010yArgAve=PhaseMean(f1010yArg[k],2*pi)
        if min(abs(f1001xArgAve-f1001yArgAve),2*pi-abs(f1001xArgAve-f1001yArgAve))>pi/2: badbpm=1
        if min(abs(f1010xArgAve-f1010yArgAve),2*pi-abs(f1010xArgAve-f1010yArgAve))>pi/2: badbpm=1

        #-- Output
        if badbpm==0:
            f1001AbsAve=mean(f1001Abs[k])
            f1010AbsAve=mean(f1010Abs[k])
            f1001ArgAve=PhaseMean(append(f1001xArg[k],f1001yArg[k]),2*pi)
            f1010ArgAve=PhaseMean(append(f1010xArg[k],f1010yArg[k]),2*pi)
            f1001Ave   =f1001AbsAve*numpy.exp(1j*f1001ArgAve)
            f1010Ave   =f1010AbsAve*numpy.exp(1j*f1010ArgAve)
            f1001AbsStd=sqrt(mean((f1001Abs[k]-f1001AbsAve)**2))
            f1010AbsStd=sqrt(mean((f1010Abs[k]-f1010AbsAve)**2))
            f1001ArgStd=PhaseStd(append(f1001xArg[k],f1001yArg[k]),2*pi)
            f1010ArgStd=PhaseStd(append(f1010xArg[k],f1010yArg[k]),2*pi)
            fwqw[bpm[k][1]]=[[f1001Ave          ,f1001AbsStd       ,f1010Ave          ,f1010AbsStd       ],
                             [f1001ArgAve/(2*pi),f1001ArgStd/(2*pi),f1010ArgAve/(2*pi),f1010ArgStd/(2*pi)]]  #-- Phases renormalized to [0,1)
            goodbpm.append(bpm[k])

    #-- Global parameters not implemented yet
    fwqw['Global']=['"null"','"null"']

    return [fwqw,goodbpm]

def GetFreeIP2_Eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,oa,op):

    #-- Common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],upper(b[1])) for b in bpm]

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
    r=sin(pi*(Qd-Q))/sin(pi*(Qd+Q))

    #-- Determine Psid (w.r.t the AC dipole) for each file
    Psidall=[]
    for i in range(len(Files)):
        if plane=='H': psid=bd*2*pi*array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V': psid=bd*2*pi*array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*pi*Qd  #-- To fix the phase shift by Q
            except: pass
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*pi*Qd
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
                try:    #-- Maybe not needed, to avoid like sqrt(-...)
                    if plane=='H':
                        al=Files[i].AMPX[Files[i].indx[bpml]]
                        ar=Files[i].AMPX[Files[i].indx[bpmr]]
                    if plane=='V':
                        al=Files[i].AMPY[Files[i].indx[bpml]]
                        ar=Files[i].AMPY[Files[i].indx[bpmr]]
                    Psidl=Psidall[i][list(zip(*bpm)[1]).index(bpml)]
                    Psidr=Psidall[i][list(zip(*bpm)[1]).index(bpmr)]
                    dpsid=Psidr-Psidl

                    #-- betd, alfd, and sqrt(2Jd) at BPM_left from amp and phase advance
                    betdl=2*L*al/(ar*sin(dpsid))
                    alfdl=(al-ar*cos(dpsid))/(ar*sin(dpsid))
                    rt2J =sqrt(al*ar*sin(dpsid)/(2*L))
                    #-- Convert to free bet and alf
                    betl=(1+r**2+2*r*numpy.cos(2*Psidl))/(1-r**2)*betdl
                    alfl=((1+r**2+2*r*numpy.cos(2*Psidl))*alfdl+2*r*numpy.sin(2*Psidl))/(1-r**2)
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
            betall =array(betall) ; betave =mean(betall) ; betstd =sqrt(mean((betall-betave)**2))
            alfall =array(alfall) ; alfave =mean(alfall) ; alfstd =sqrt(mean((alfall-alfave)**2))
            betsall=array(betsall); betsave=mean(betsall); betsstd=sqrt(mean((betsall-betsave)**2))
            dsall  =array(dsall)  ; dsave  =mean(dsall)  ; dsstd  =sqrt(mean((dsall-dsave)**2))
            rt2Jall=array(rt2Jall); rt2Jave=mean(rt2Jall); rt2Jstd=sqrt(mean((rt2Jall-rt2Jave)**2))
            result['IP'+ip]=[betave,betstd,betmdl,alfave,alfstd,alfmdl,betsave,betsstd,betsmdl,dsave,dsstd,dsmdl,rt2Jave,rt2Jstd]

    return result

def getkickac(MADTwiss_ac,files,Qh,Qv,Qx,Qy,psih_ac2bpmac,psiv_ac2bpmac,bd,op):

    invarianceJx=[]
    invarianceJy=[]
    tunex       =[]
    tuney       =[]
    tunexRMS    =[]
    tuneyRMS    =[]
    dpp=[]

    for j in range(len(files[0])):

        x=files[0][j]
        y=files[1][j]
        [beta,rmsbb,bpms,invariantJx]=GetFreeBetaFromAmp_Eq(MADTwiss_ac,[x],Qh,Qx,psih_ac2bpmac,'H',bd,op)
        [beta,rmsbb,bpms,invariantJy]=GetFreeBetaFromAmp_Eq(MADTwiss_ac,[y],Qv,Qy,psiv_ac2bpmac,'V',bd,op)
        invarianceJx.append(invariantJx)
        invarianceJy.append(invariantJy)
        try:
            dpp.append(x.DPP)
        except:
            dpp.append(0.0)
        tunex.append(x.Q1)
        tuney.append(y.Q2)
        tunexRMS.append(x.Q1RMS)
        tuneyRMS.append(y.Q2RMS)

    tune   =[tunex,tuney]
    tuneRMS=[tunexRMS,tuneyRMS]

    return [invarianceJx,invarianceJy,tune,tuneRMS,dpp]

######### end ac-dipole stuff



#----------------- end glenn part

#---- Functions for Andy's BetaFromAmp re-scaling

def filterbpm(ListOfBPM):
    '''Filter non-arc BPM'''
    if len(ListOfBPM)==0:
        print >> sys.stderr, "Nothing to filter!!!!"
        sys.exit(1)
    result=[]
    for b in ListOfBPM:
        if ('BPM.' in b[1] or 'bpm.' in b[1]):
            result.append(b)
    return result

def union(a, b):
    ''' return the union of two lists '''
    return list(set(a) | set(b))

def _fix_output(outputpath):
    if not os.path.isdir(outputpath):
        os.makedirs(outputpath)
    if '/'!=outputpath[-1]:
        outputpath+='/'
    return outputpath

def _write_llm_tfs_header(filename,mad_files):
    fout=file(filename,'w')
    fout.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    for mad_file in mad_files:
        fout.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fout.write('@ FILES %s "')
    return fout

#######################################################
#                   Main part                         #
#######################################################


#-- Reading sys.argv
def parse_args():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-a", "--accel",
                    help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
                    metavar="ACCEL", default="LHCB1",dest="ACCEL")
    parser.add_option("-d", "--dictionary",
                    help="File with the BPM dictionary",
                    metavar="DICT", default="0", dest="dict")
    parser.add_option("-m", "--model",
                    help="Twiss free model file *.dat. For free oscillations, ONLY the file *.dat should be present in the path. For AC dipole, the path should also contain a file *_ac.dat (BOTH files a needed in this case).",
                    metavar="TwissFile", default="0", dest="Twiss")
    parser.add_option("-f", "--files",
                    help="Files from analysis, separated by comma",
                    metavar="TwissFile", default="0", dest="files")
    parser.add_option("-o", "--output",
                    help="Output Path",
                    metavar="OUT", default="./", dest="output")
    parser.add_option("-c", "--cocut",
                    help="Cut for closed orbit measurement [um]",
                    metavar="COCUT", default=4000, dest="COcut")
    parser.add_option("-n", "--nbcpl",
                    help="Analysis option for coupling, 1 bpm or 2 bpms",
                    metavar="NBCPL", default=2, dest="NBcpl")
    parser.add_option("-t", "--tbtana",
                    help="Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA",
                    metavar="TBTANA", default="SUSSIX", dest="TBTana")
    parser.add_option("-b", "--bpmu",
                    help="BPMunit: um, mm, cm, m (default um)",
                    metavar="BPMUNIT", default="um", dest="BPMUNIT")
    parser.add_option("-l", "--nonlinear",
                    help="Switch to output higher oerder resonance stuffs, on=1(default)/off=0",
                    metavar="HIGHER", default="1" , dest="higher")
    parser.add_option("-p", "--lhcphase",
                    help="Compensate phase shifts by tunes for the LHC experiment data, off=0(default)/on=1",
                    metavar="LHCPHASE", default="0" , dest="lhcphase")

    (options, args) = parser.parse_args()
    return options,args

def main(outputpath,files_to_analyse,twiss_model_file,dict_file="0",accel="LHCB1",lhcphase="0",BPMU="um",COcut=4000,NBcpl=2,TBTana="SUSSIX",higher_order=1):
    '''
     GetLLM main function.

     :param outputpath: The output path to store results
     :param files_to_analyse: List of files, comma separated string
     :param NBcpl: For selecting the coupling measurement method
     :param higher_order: output higher order resonance stuff
    '''

    print "Starting GetLLM ", VERSION

    outputpath=_fix_output(outputpath)

    if dict_file=="0":
        BPMdictionary={}
    else:
        execfile(dict_file)
        BPMdictionary=dictionary   # temporaryly since presently name is not BPMdictionary

    listOfInputFiles=files_to_analyse.split(",")

    # Beam direction
    bd=1
    if accel=="LHCB2":
        bd=-1 # THIS IS CORRECT, be careful with tune sign in SUSSIX and eigenmode order in SVD
    elif accel=="LHCB4":
        bd=1  # IS THIS CORRECT? I (rogelio) try for Simon...
        accel="LHCB2" #This is because elements later are named B2 anyway, not B4

    #-- finding base model
    try:
        MADTwiss=twiss(twiss_model_file,BPMdictionary) # MODEL from MAD
        print "Base model found!"
    except:
        print >> sys.stderr, "twiss file loading failed for:",twiss_model_file
        print >> sys.stderr, traceback.format_exc()
        sys.exit(1)

    #-- finding the ac dipole model
    try:
        MADTwiss_ac=twiss(twiss_model_file.replace(".dat","_ac.dat"))
        acswitch="1"
        print "Driven Twiss file found. AC dipole effects calculated with the effective model (get***_free2.out)"
    except:
        MADTwiss_ac=MADTwiss
        acswitch="0"
        print "WARN: AC dipole effects not calculated. Driven twiss file does not exsist !"

    #-- Test if the AC dipole (MKQA) is in the model of LHC
    if acswitch=='1':
        if 'LHC' in accel:
            if 'MKQA.6L4.'+accel[3:] in MADTwiss.NAME:
                print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
            else:
                try:
                    MADTwissElem=twiss(twiss_model_file.replace(".dat","_elements.dat"))
                    print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
                except:
                    print 'WARN: AC dipoles not in the model. AC dipole effects not calculated with analytic equations !'
        else: print 'WARN: AC dipole effects calculated with analytic equations only for LHC for now'

    if BPMU=='um': COcut=COcut
    elif BPMU=='mm': COcut=COcut/1.0e3
    elif BPMU=='cm': COcut=COcut/1.0e4
    elif BPMU=='m': COcut=COcut/1.0e6




    if TBTana=="SUSSIX":
        Suffix1='_linx'
        Suffix2='_liny'
    elif TBTana=='SVD':
        Suffix1='_svdx'
        Suffix2='_svdy'
    elif TBTana=='HA':
        Suffix1='_hax'
        Suffix2='_hay'


    fphasex=open(outputpath+'getphasex.out','w')
    fphasey=open(outputpath+'getphasey.out','w')
    fphasex.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fphasey.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fphasex.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fphasey.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fphasex.write('@ FILES %s "')
    fphasey.write('@ FILES %s "')

    fphasexT=open(outputpath+'getphasetotx.out','w')
    fphaseyT=open(outputpath+'getphasetoty.out','w')
    fphasexT.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fphaseyT.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fphasexT.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fphaseyT.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fphasexT.write('@ FILES %s "')
    fphaseyT.write('@ FILES %s "')

    if acswitch=="1":

        fphasexf=open(outputpath+'getphasex_free.out','w')
        fphaseyf=open(outputpath+'getphasey_free.out','w')
        fphasexf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphaseyf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphasexf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphaseyf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphasexf.write('@ FILES %s "')
        fphaseyf.write('@ FILES %s "')
        fphasexf2=open(outputpath+'getphasex_free2.out','w')
        fphaseyf2=open(outputpath+'getphasey_free2.out','w')
        fphasexf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphaseyf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphasexf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphaseyf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphasexf2.write('@ FILES %s "')
        fphaseyf2.write('@ FILES %s "')

        fphasexTf=open(outputpath+'getphasetotx_free.out','w')
        fphaseyTf=open(outputpath+'getphasetoty_free.out','w')
        fphasexTf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphaseyTf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphasexTf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphaseyTf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphasexTf.write('@ FILES %s "')
        fphaseyTf.write('@ FILES %s "')
        fphasexTf2=open(outputpath+'getphasetotx_free2.out','w')
        fphaseyTf2=open(outputpath+'getphasetoty_free2.out','w')
        fphasexTf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphaseyTf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fphasexTf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphaseyTf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fphasexTf2.write('@ FILES %s "')
        fphaseyTf2.write('@ FILES %s "')

    fbetax=open(outputpath+'getbetax.out','w')
    fbetay=open(outputpath+'getbetay.out','w')
    fbetax.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fbetay.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fbetax.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fbetay.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fbetax.write('@ FILES %s "')
    fbetay.write('@ FILES %s "')

    if acswitch=="1":
        fbetaxf=open(outputpath+'getbetax_free.out','w')
        fbetayf=open(outputpath+'getbetay_free.out','w')
        fbetaxf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fbetayf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fbetaxf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fbetayf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fbetaxf.write('@ FILES %s "')
        fbetayf.write('@ FILES %s "')
        fbetaxf2=open(outputpath+'getbetax_free2.out','w')
        fbetayf2=open(outputpath+'getbetay_free2.out','w')
        fbetaxf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fbetayf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fbetaxf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fbetayf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fbetaxf2.write('@ FILES %s "')
        fbetayf2.write('@ FILES %s "')

    fabetax=open(outputpath+'getampbetax.out','w')
    fabetay=open(outputpath+'getampbetay.out','w')
    fabetax.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fabetay.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fabetax.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fabetay.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fabetax.write('@ FILES %s "')
    fabetay.write('@ FILES %s "')

    if acswitch=="1":
        fabetaxf=open(outputpath+'getampbetax_free.out','w')
        fabetayf=open(outputpath+'getampbetay_free.out','w')
        fabetaxf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fabetayf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fabetaxf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fabetayf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fabetaxf.write('@ FILES %s "')
        fabetayf.write('@ FILES %s "')
        fabetaxf2=open(outputpath+'getampbetax_free2.out','w')
        fabetayf2=open(outputpath+'getampbetay_free2.out','w')
        fabetaxf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fabetayf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fabetaxf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fabetayf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fabetaxf2.write('@ FILES %s "')
        fabetayf2.write('@ FILES %s "')

    fcox=open(outputpath+'getCOx.out','w')
    fcoy=open(outputpath+'getCOy.out','w')
    fcox.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fcoy.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fcox.write('@ FILES %s "')
    fcoy.write('@ FILES %s "')

    fNDx=open(outputpath+'getNDx.out','w')
    fDx=open(outputpath+'getDx.out','w')
    fDy=open(outputpath+'getDy.out','w')
    fNDx.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fDx.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fDy.write('@ GetLLMVersion %s "'+VERSION+'"\n')
    fNDx.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fDx.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fDy.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fNDx.write('@ FILES %s "')
    fDx.write('@ FILES %s "')
    fDy.write('@ FILES %s "')

    fcouple=open(outputpath+'getcouple.out','w')
    fcouple.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
    fcouple.write('@ FILES %s "')
    if acswitch=="1":
        fcouplef=open(outputpath+'getcouple_free.out','w')
        fcouplef.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fcouplef.write('@ FILES %s "')
        fcouplef2=open(outputpath+'getcouple_free2.out','w')
        fcouplef2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fcouplef2.write('@ FILES %s "')


    fcoupleterms=open(outputpath+'getcoupleterms.out','w')
    fcoupleterms.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    if "LHC" in accel:
        fIP=open(outputpath+'getIP.out','w')
        fIP.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fIP.write('* NAME  BETASTARH  BETASTARHMDL   H   PHIH PHIXH   PHIHMDL  BETASTARV  BETASTARVMDL  V   PHIV  PHIYV  PHIVMDL\n')
        fIP.write('$  %s  %le    %le   %le  %le    %le    %le    %le   %le   %le   %le   %le   %le  \n')
        fIPx=open(outputpath+'getIPx.out','w')
        fIPx.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fIPx.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fIPx.write('@ FILES %s "')
        fIPy=open(outputpath+'getIPy.out','w')
        fIPy.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fIPy.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fIPy.write('@ FILES %s "')
        fIPfromphase=open(outputpath+'getIPfromphase.out','w')
        fIPfromphase.write('@ GetLLMVersion %s "'+VERSION+'"\n')
        fIPfromphase.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
        fIPfromphase.write('@ FILES %s "')
        if acswitch=='1':
            fIPxf=open(outputpath+'getIPx_free.out','w')
            fIPxf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
            fIPxf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fIPxf.write('@ FILES %s "')
            fIPyf=open(outputpath+'getIPy_free.out','w')
            fIPyf.write('@ GetLLMVersion %s "'+VERSION+'"\n')
            fIPyf.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fIPyf.write('@ FILES %s "')
            fIPxf2=open(outputpath+'getIPx_free2.out','w')
            fIPxf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
            fIPxf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fIPxf2.write('@ FILES %s "')
            fIPyf2=open(outputpath+'getIPy_free2.out','w')
            fIPyf2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
            fIPyf2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fIPyf2.write('@ FILES %s "')
            fIPfromphasef=open(outputpath+'getIPfromphase_free.out','w')
            fIPfromphasef.write('@ GetLLMVersion %s "'+VERSION+'"\n')
            fIPfromphasef.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fIPfromphasef.write('@ FILES %s "')
            fIPfromphasef2=open(outputpath+'getIPfromphase_free2.out','w')
            fIPfromphasef2.write('@ GetLLMVersion %s "'+VERSION+'"\n')
            fIPfromphasef2.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fIPfromphasef2.write('@ FILES %s "')

    FileOfZeroDPPX=[]
    FileOfZeroDPPY=[]
    FileOfNonZeroDPPX=[]
    FileOfNonZeroDPPY=[]
    ListOfZeroDPPX=[]
    ListOfNonZeroDPPX=[]
    ListOfZeroDPPY=[]
    ListOfNonZeroDPPY=[]
    woliny=0  # Let's assume there is liny for the moment
    woliny2=0
    for filein in listOfInputFiles:
        if filein[-3:]=='.gz':
            file1=filein[:-3]+Suffix1+'.gz'
        else:
            file1=filein+Suffix1
        x1=twiss(file1)
        try:
            dppi=x1.DPP
        except:
            dppi=0.0

        if type(dppi)!=float:
            print 'Warning: DPP may not be given as a number in ',file1 ,'...trying to forcibly cast it as a number'
            try:
                dppi=float(dppi)
                print 'dppi= ',dppi
            except:
                print >> sys.stderr, 'but failing. DPP in ',file1 , ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)

        if dppi==0.0:
            ListOfZeroDPPX.append(twiss(file1))
            FileOfZeroDPPX.append(file1)
            fphasex.write(file1+'')
            fphasexT.write(file1+'')
            fbetax.write(file1+'')
            fabetax.write(file1+'')
            fcox.write(file1+'')
            fNDx.write(file1+'')
            fDx.write(file1+'')
            fcouple.write(filein+'')
            if "LHC" in accel:
                fIPx.write(filein+'')
                fIPy.write(filein+'')
                fIPfromphase.write(filein+'')
                if acswitch=='1':
                    fIPxf.write(filein+'')
                    fIPyf.write(filein+'')
                    fIPxf2.write(filein+'')
                    fIPyf2.write(filein+'')
                    fIPfromphasef.write(filein+'')
                    fIPfromphasef2.write(filein+'')
            if acswitch=="1":
                fphasexf.write(file1+'')
                fphasexf2.write(file1+'')
                fphasexTf.write(file1+'')
                fphasexTf2.write(file1+'')
                fbetaxf.write(file1+'')
                fbetaxf2.write(file1+'')
                fabetaxf.write(file1+'')
                fabetaxf2.write(file1+'')
                fcouplef.write(filein+'')
                fcouplef2.write(filein+'')
        else:
            ListOfNonZeroDPPX.append(twiss(file1))
            FileOfNonZeroDPPX.append(file1)
            fNDx.write(file1+' ')
            fDx.write(file1+' ')

        try:
            if filein[-3:]=='.gz':
                file1=filein[:-3]+Suffix2+'.gz'
            else:
                file1=filein+Suffix2
            y1=twiss(file1)
            try:
                dppi=y1.DPP
            except:
                dppi=0.0

            if type(dppi)!=float:
                print 'Warning: DPP may not be given as a number in ',file1 ,'...trying to forcibly cast it as a number'
                try:
                    dppi=float(dppi)
                    print 'dppi= ',dppi
                except:
                    print 'but failing. DPP in ',file1 , ' is something wrong. String?'
                    raise ValueError('leaving GetLLM')

            if dppi==0.0:
                ListOfZeroDPPY.append(twiss(file1))
                FileOfZeroDPPY.append(file1)
                fphasey.write(file1+'')
                fphaseyT.write(file1+'')
                fbetay.write(file1+'')
                fabetay.write(file1+'')
                fcoy.write(file1+' ')
                fDy.write(file1+' ')
                if acswitch=="1":
                    fphaseyf.write(file1+'')
                    fphaseyf2.write(file1+'')
                    fphaseyTf.write(file1+'')
                    fphaseyTf2.write(file1+'')
                    fbetayf.write(file1+'')
                    fbetayf2.write(file1+'')
                    fabetayf.write(file1+'')
                    fabetayf2.write(file1+'')
            else:
                ListOfNonZeroDPPY.append(twiss(file1))
                FileOfNonZeroDPPY.append(file1)
                fDy.write(file1+' ')
        except:
            print 'Warning: There seems no '+str(file1)+' file in the specified directory.'




    woliny=0
    woliny2=0
    wolinx=0
    wolinx2=0


    if len(ListOfZeroDPPY)==0 :
        woliny=1  #FLAG meaning there is no _liny file for zero DPPY!
    if len(ListOfNonZeroDPPY)==0 :
        woliny2=1  #FLAG meaning there is no _liny file for non-zero DPPY!
    if len(ListOfNonZeroDPPX)==0 :
        wolinx2=1

    if len(ListOfZeroDPPX)==0 :
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'
        wolinx=1

    if (len(ListOfNonZeroDPPX)!=0) and (len(ListOfZeroDPPX)==0):
        ListOfZeroDPPX=ListOfNonZeroDPPX
        ListOfZeroDPPY=ListOfNonZeroDPPY
        wolinx=0
        woliny=0
        woliny2=1
        wolinx2=1
        print "Previous warning suppressed, running in chromatic mode"
        fphasex.write('chrommode')
        fbetax.write('chrommode')
        fabetax.write('chrommode')
        fcox.write('chrommode')
        fNDx.write('chrommode')
        fDx.write('chrommode')
        fcouple.write('chrommode')
        if acswitch=="1":
            fcouplef.write('chrommode')
            fcouplef2.write('chrommode')
        fphasey.write('chrommode')
        fbetay.write('chrommode')
        fabetay.write('chrommode')
        fcoy.write('chrommode')
        fDy.write('chrommode')


    if acswitch=="1":
        Q1f=abs(float(str(MADTwiss.Q1).split('.')[0])-MADTwiss.Q1)        #-- Free Q1 (tempolarlly, overwritten later)
        Q2f=abs(float(str(MADTwiss.Q2).split('.')[0])-MADTwiss.Q2)        #-- Free Q2 (tempolarlly, overwritten later)
        Q1 =abs(float(str(MADTwiss_ac.Q1).split('.')[0])-MADTwiss_ac.Q1)  #-- Drive Q1 (tempolarlly, overwritten later)
        Q2 =abs(float(str(MADTwiss_ac.Q2).split('.')[0])-MADTwiss_ac.Q2)  #-- Drive Q2 (tempolarlly, overwritten later)
        d1 =Q1-Q1f                                                        #-- Used later to calculate free Q1
        d2 =Q2-Q2f                                                        #-- Used later to calculate free Q2
    else:
        Q1f=ListOfZeroDPPX[0].Q1
        Q2f=ListOfZeroDPPY[0].Q2


    fphasex.write('"'+'\n')
    fphasey.write('"'+'\n')
    fphasexT.write('"'+'\n')
    fphaseyT.write('"'+'\n')
    fbetax.write('"'+'\n')
    fbetay.write('"'+'\n')
    fabetax.write('"'+'\n')
    fabetay.write('"'+'\n')
    fcox.write('"'+'\n')
    fcoy.write('"'+'\n')
    fNDx.write('"'+'\n')
    fDx.write('"'+'\n')
    fDy.write('"'+'\n')
    fcouple.write('"'+'\n')
    if "LHC" in accel:
        fIPx.write('"'+'\n')
        fIPy.write('"'+'\n')
        fIPfromphase.write('"'+'\n')
        if acswitch=='1':
            fIPxf.write('"'+'\n')
            fIPyf.write('"'+'\n')
            fIPxf2.write('"'+'\n')
            fIPyf2.write('"'+'\n')
            fIPfromphasef.write('"'+'\n')
            fIPfromphasef2.write('"'+'\n')
    if acswitch=="1":
        fphasexf.write('"'+'\n')
        fphaseyf.write('"'+'\n')
        fphasexf2.write('"'+'\n')
        fphaseyf2.write('"'+'\n')
        fphasexTf.write('"'+'\n')
        fphaseyTf.write('"'+'\n')
        fphasexTf2.write('"'+'\n')
        fphaseyTf2.write('"'+'\n')
        fbetaxf.write('"'+'\n')
        fbetayf.write('"'+'\n')
        fbetaxf2.write('"'+'\n')
        fbetayf2.write('"'+'\n')
        fabetaxf.write('"'+'\n')
        fabetayf.write('"'+'\n')
        fabetaxf2.write('"'+'\n')
        fabetayf2.write('"'+'\n')
        fcouplef.write('"'+'\n')
        fcouplef2.write('"'+'\n')

    # Construct pseudo-double plane BPMs
    if (accel=="SPS" or "RHIC" in accel) and wolinx!=1 and woliny!=1 :
        execfile(twiss_model_file.replace("twiss.dat","BPMpair.py"))
        [PseudoListX,PseudoListY]=PseudoDoublePlaneMonitors(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, BPMdictionary)


    #-------- Check monitor compatibility between data and model

    ALL=ListOfNonZeroDPPX+ListOfZeroDPPX+ListOfNonZeroDPPY+ListOfZeroDPPY
    for j in range(0,len(ALL)) :
        z=ALL[j].NAME
        for bpm in z:
            try:
                check=MADTwiss.NAME[MADTwiss.indx[bpm]]
            except:
                try:
                    check=MADTwiss.NAME[MADTwiss.indx[upper(bpm)]]
                except:
                    print 'Monitor '+bpm+' cannot be found in the model!'
                    #exit()


    #-------- START Phase
    print 'Calculating phase'

    #---- Calling GetPhases first to save tunes
    if wolinx!=1 and woliny!=1:
        #-- Calculate temp value of tune
        Q1temp=[]
        Q2temp=[]
        for i in ListOfZeroDPPX: Q1temp.append(mean(i.TUNEX))
        for i in ListOfZeroDPPY: Q2temp.append(mean(i.TUNEY))
        Q1temp=mean(Q1temp)
        Q2temp=mean(Q2temp)
        if len(ListOfZeroDPPX[0].NAME)==0:
            print "No BPMs in linx file"
            sys.exit(1)
        if  len(ListOfZeroDPPY[0].NAME)==0:
            print "No BPMs in liny file"
            sys.exit(1)
        [phasex,Q1,MUX,bpmsx]=GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1temp,'H',outputpath,bd,accel,lhcphase)
        [phasey,Q2,MUY,bpmsy]=GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2temp,'V',outputpath,bd,accel,lhcphase)
        print "KK end"
    elif wolinx!=1:
        #-- Calculate temp value of tune
        Q1temp=[]
        for i in ListOfZeroDPPX: Q1temp.append(mean(i.TUNEX))
        Q1temp=mean(Q1temp)

        [phasex,Q1,MUX,bpmsx]=GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1temp,'H',outputpath,bd,accel,lhcphase)
        print 'liny missing and output x only ...'
    elif woliny!=1:
        #-- Calculate temp value of tune
        Q2temp=[]
        for i in ListOfZeroDPPY: Q2temp.append(mean(i.TUNEY))
        Q2temp=mean(Q2temp)

        [phasey,Q2,MUY,bpmsy]=GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2temp,'V',outputpath,bd,accel,lhcphase)
        print 'linx missing and output y only ...'

    #---- Re-run GetPhase to fix the phase shift by Q for exp data of LHC
    if lhcphase=="1":
        if wolinx!=1 and woliny!=1:
            [phasex,Q1,MUX,bpmsx]=GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',outputpath,bd,accel,lhcphase)
            [phasey,Q2,MUY,bpmsy]=GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',outputpath,bd,accel,lhcphase)
        elif wolinx!=1:
            [phasex,Q1,MUX,bpmsx]=GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',outputpath,bd,accel,lhcphase)
        elif woliny!=1:
            [phasey,Q2,MUY,bpmsy]=GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',outputpath,bd,accel,lhcphase)


    #---- ac to free phase from eq and the model
    if acswitch=='1':
        if wolinx!=1:
            Q1f=Q1-d1  #-- Free H-tune
            try:    acphasex_ac2bpmac=GetACPhase_AC2BPMAC(MADTwissElem,Q1,Q1f,'H',accel)
            except: acphasex_ac2bpmac=GetACPhase_AC2BPMAC(MADTwiss,Q1,Q1f,'H',accel)
            [phasexf,muxf,bpmsxf]=GetFreePhase_Eq(MADTwiss,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',bd,lhcphase)
            [phasexf2,muxf2,bpmsxf2]=getfreephase(phasex,Q1,Q1f,bpmsx,MADTwiss_ac,MADTwiss,"H")
        if woliny!=1:
            Q2f=Q2-d2  #-- Free V-tune
            try:    acphasey_ac2bpmac=GetACPhase_AC2BPMAC(MADTwissElem,Q2,Q2f,'V',accel)
            except: acphasey_ac2bpmac=GetACPhase_AC2BPMAC(MADTwiss,Q2,Q2f,'V',accel)
            [phaseyf,muyf,bpmsyf]=GetFreePhase_Eq(MADTwiss,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',bd,lhcphase)
            [phaseyf2,muyf2,bpmsyf2]=getfreephase(phasey,Q2,Q2f,bpmsy,MADTwiss_ac,MADTwiss,"V")

    #---- H plane result
    if wolinx!=1:

        phasexlist=[]
        phasex['DPP']=0.0
        phasexlist.append(phasex)
        fphasex.write('@ Q1 %le '+str(Q1)+'\n')
        fphasex.write('@ MUX %le '+str(MUX)+'\n')
        try:
            fphasex.write('@ Q2 %le '+str(Q2)+'\n')
            fphasex.write('@ MUY %le '+str(MUY)+'\n')
        except:
            fphasex.write('@ Q2 %le '+'0.0'+'\n')
            fphasex.write('@ MUY %le '+'0.0'+'\n')
        fphasex.write('* NAME   NAME2  S   S1   COUNT  PHASEX  STDPHX  PHXMDL MUXMDL\n')
        fphasex.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
        for i in range(len(bpmsx)):
            bn1=upper(bpmsx[i][1])
            bns1=bpmsx[i][0]
            phmdl=phasex[bn1][4]
            if i==len(bpmsx)-1:
                bn2=upper(bpmsx[0][1])
                bns2=bpmsx[0][0]
            else:
                bn2=upper(bpmsx[i+1][1])
                bns2=bpmsx[i+1][0]
            fphasex.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasex[bn1][0])+' '+str(phasex[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]])+'\n' )
        fphasex.close()

        #-- ac to free phase
        if acswitch=='1':

            #-- from eq
            try:
                fphasexf.write('@ Q1 %le '+str(Q1f)+'\n')
                fphasexf.write('@ MUX %le '+str(muxf)+'\n')
                try:
                    fphasexf.write('@ Q2 %le '+str(Q2f)+'\n')
                    fphasexf.write('@ MUY %le '+str(muyf)+'\n')
                except:
                    fphasexf.write('@ Q2 %le '+'0.0'+'\n')
                    fphasexf.write('@ MUY %le '+'0.0'+'\n')
                fphasexf.write('* NAME   NAME2  S   S1   COUNT  PHASEX  STDPHX  PHXMDL MUXMDL\n')
                fphasexf.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
                for i in range(len(bpmsxf)):
                    bn1=upper(bpmsxf[i][1])
                    bns1=bpmsxf[i][0]
                    phmdlf=phasexf[bn1][4]
                    if i==len(bpmsxf)-1:
                        bn2=upper(bpmsxf[0][1])
                        bns2=bpmsxf[0][0]
                    else:
                        bn2=upper(bpmsxf[i+1][1])
                        bns2=bpmsxf[i+1][0]
                    fphasexf.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasexf[bn1][0])+' '+str(phasexf[bn1][1])+' '+str(phmdlf)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            except: pass
            fphasexf.close()

            #-- from the model
            fphasexf2.write('@ Q1 %le '+str(Q1f)+'\n')
            fphasexf2.write('@ MUX %le '+str(muxf2)+'\n')
            try:
                fphasexf2.write('@ Q2 %le '+str(Q2f)+'\n')
                fphasexf2.write('@ MUY %le '+str(muyf2)+'\n')
            except:
                fphasexf2.write('@ Q2 %le '+'0.0'+'\n')
                fphasexf2.write('@ MUY %le '+'0.0'+'\n')
            fphasexf2.write('* NAME   NAME2  S   S1   COUNT  PHASEX  STDPHX  PHXMDL MUXMDL\n')
            fphasexf2.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpmsxf2)):
                bn1=upper(bpmsxf2[i][1])
                bns1=bpmsxf2[i][0]
                phmdlf2=phasexf2[bn1][2]
                bn2=phasexf2[bn1][3]
                bns2=phasexf2[bn1][4]
                fphasexf2.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasexf2[bn1][0])+' '+str(phasexf2[bn1][1])+' '+str(phmdlf2)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fphasexf2.close()

    #---- V plane result
    if woliny!=1:

        phaseylist=[]
        phasey['DPP']=0.0
        phaseylist.append(phasey)
        try:
            fphasey.write('@ Q1 %le '+str(Q1)+'\n')
            fphasey.write('@ MUX %le '+str(MUX)+'\n')
        except:
            fphasey.write('@ Q1 %le '+'0.0'+'\n')
            fphasey.write('@ MUX %le '+'0.0'+'\n')
        fphasey.write('@ Q2 %le '+str(Q2)+'\n')
        fphasey.write('@ MUY %le '+str(MUY)+'\n')
        fphasey.write('* NAME   NAME2  S   S1   COUNT  PHASEY  STDPHY  PHYMDL MUYMDL\n')
        fphasey.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
        for i in range(len(bpmsy)):
            bn1=upper(bpmsy[i][1])
            bns1=bpmsy[i][0]
            phmdl=phasey[bn1][4]
            if i==len(bpmsy)-1:
                bn2=upper(bpmsy[0][1])
                bns2=bpmsy[0][0]
            else:
                bn2=upper(bpmsy[i+1][1])
                bns2=bpmsy[i+1][0]
            fphasey.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phasey[bn1][0])+' '+str(phasey[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]])+'\n' )
        fphasey.close()

        #-- ac to free phase
        if acswitch=='1':

            #-- from eq
            try:
                try:
                    fphaseyf.write('@ Q1 %le '+str(Q1f)+'\n')
                    fphaseyf.write('@ MUX %le '+str(muxf)+'\n')
                except:
                    fphaseyf.write('@ Q1 %le '+'0.0'+'\n')
                    fphaseyf.write('@ MUX %le '+'0.0'+'\n')
                fphaseyf.write('@ Q2 %le '+str(Q2f)+'\n')
                fphaseyf.write('@ MUY %le '+str(muyf)+'\n')
                fphaseyf.write('* NAME   NAME2  S   S1   COUNT  PHASEY  STDPHY  PHYMDL MUYMDL\n')
                fphaseyf.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
                for i in range(len(bpmsyf)):
                    bn1=upper(bpmsyf[i][1])
                    bns1=bpmsyf[i][0]
                    phmdlf=phaseyf[bn1][4]
                    if i==len(bpmsyf)-1:
                        bn2=upper(bpmsyf[0][1])
                        bns2=bpmsyf[0][0]
                    else:
                        bn2=upper(bpmsyf[i+1][1])
                        bns2=bpmsyf[i+1][0]
                    fphaseyf.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phaseyf[bn1][0])+' '+str(phaseyf[bn1][1])+' '+str(phmdlf)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            except: pass
            fphaseyf.close()

            #-- from the model
            try:
                fphaseyf2.write('@ Q1 %le '+str(Q1f)+'\n')
                fphaseyf2.write('@ MUX %le '+str(muxf2)+'\n')
            except:
                fphaseyf2.write('@ Q1 %le '+'0.0'+'\n')
                fphaseyf2.write('@ MUX %le '+'0.0'+'\n')
            fphaseyf2.write('@ Q2 %le '+str(Q2f)+'\n')
            fphaseyf2.write('@ MUY %le '+str(muyf2)+'\n')
            fphaseyf2.write('* NAME   NAME2  S   S1   COUNT  PHASEY  STDPHY  PHYMDL MUYMDL\n')
            fphaseyf2.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpmsyf2)):
                bn1=upper(bpmsyf2[i][1])
                bns1=bpmsyf2[i][0]
                phmdlf2=phaseyf2[bn1][2]
                bn2=phaseyf2[bn1][3]
                bns2=phaseyf2[bn1][4]
                fphaseyf2.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phaseyf2[bn1][0])+' '+str(phaseyf2[bn1][1])+' '+str(phmdlf2)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            fphaseyf2.close()


    #-------- START Total Phase
    print 'Calculating total phase'

    #---- H plane result
    if wolinx!=1:

        [phasexT,bpmsxT]=GetPhasesTotal(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',bd,accel,lhcphase)
        fphasexT.write('@ Q1 %le '+str(Q1)+'\n')
        fphasexT.write('@ MUX %le '+str(MUX)+'\n')
        try:
            fphasexT.write('@ Q2 %le '+str(Q2)+'\n')
            fphasexT.write('@ MUY %le '+str(MUY)+'\n')
        except:
            fphasexT.write('@ Q2 %le '+'0.0'+'\n')
            fphasexT.write('@ MUY %le '+'0.0'+'\n')
        fphasexT.write('* NAME   NAME2  S   S1   COUNT  PHASEX  STDPHX  PHXMDL MUXMDL\n')
        fphasexT.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
        for i in range(0,len(bpmsxT)):
            bn1=upper(bpmsxT[i][1])
            bns1=bpmsxT[i][0]
            phmdl=phasexT[bn1][2]
            bn2=upper(bpmsxT[0][1])
            bns2=bpmsxT[0][0]
            fphasexT.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasexT[bn1][0])+' '+str(phasexT[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]])+'\n' )
        fphasexT.close()

        #-- ac to free total phase
        if acswitch=='1':

            #-- from eq
            try:
                [phasexTf,bpmsxTf]=GetFreePhaseTotal_Eq(MADTwiss,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',bd,lhcphase)
                fphasexTf.write('@ Q1 %le '+str(Q1f)+'\n')
                fphasexTf.write('@ MUX %le '+str(muxf)+'\n')
                try:
                    fphasexTf.write('@ Q2 %le '+str(Q2f)+'\n')
                    fphasexTf.write('@ MUY %le '+str(muyf)+'\n')
                except:
                    fphasexTf.write('@ Q2 %le '+'0.0'+'\n')
                    fphasexTf.write('@ MUY %le '+'0.0'+'\n')
                fphasexTf.write('* NAME   NAME2  S   S1   COUNT  PHASEX  STDPHX  PHXMDL MUXMDL\n')
                fphasexTf.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
                for i in range(0,len(bpmsxTf)):
                    bn1=upper(bpmsxTf[i][1])
                    bns1=bpmsxTf[i][0]
                    phmdlf=phasexTf[bn1][2]
                    bn2=upper(bpmsxTf[0][1])
                    bns2=bpmsxTf[0][0]
                    fphasexTf.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasexTf[bn1][0])+' '+str(phasexTf[bn1][1])+' '+str(phmdlf)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            except: pass
            fphasexTf.close()

            #-- from the model
            [phasexTf2,bpmsxTf2]=getfreephaseTotal(phasexT,bpmsxT,"H",MADTwiss,MADTwiss_ac)
            fphasexTf2.write('@ Q1 %le '+str(Q1f)+'\n')
            fphasexTf2.write('@ MUX %le '+str(muxf2)+'\n')
            try:
                fphasexTf2.write('@ Q2 %le '+str(Q2f)+'\n')
                fphasexTf2.write('@ MUY %le '+str(muyf2)+'\n')
            except:
                fphasexTf2.write('@ Q2 %le '+'0.0'+'\n')
                fphasexTf2.write('@ MUY %le '+'0.0'+'\n')
            fphasexTf2.write('* NAME   NAME2  S   S1   COUNT  PHASEX  STDPHX  PHXMDL MUXMDL\n')
            fphasexTf2.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpmsxTf2)):
                bn1=upper(bpmsxTf2[i][1])
                bns1=bpmsxTf2[i][0]
                phmdlf2=phasexTf2[bn1][2]
                bn2=upper(bpmsxTf2[0][1])
                bns2=bpmsxTf2[0][0]
                fphasexTf2.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasexTf2[bn1][0])+' '+str(phasexTf2[bn1][1])+' '+str(phmdlf2)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fphasexTf2.close()


    #---- V plane result
    if woliny!=1:

        [phaseyT,bpmsyT]=GetPhasesTotal(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',bd,accel,lhcphase)
        try:
            fphaseyT.write('@ Q1 %le '+str(Q1)+'\n')
            fphaseyT.write('@ MUX %le '+str(MUX)+'\n')
        except:
            fphaseyT.write('@ Q1 %le '+'0.0'+'\n')
            fphaseyT.write('@ MUX %le '+'0.0'+'\n')
        fphaseyT.write('@ Q2 %le '+str(Q2)+'\n')
        fphaseyT.write('@ MUY %le '+str(MUY)+'\n')
        fphaseyT.write('* NAME   NAME2  S   S1   COUNT  PHASEY  STDPHY  PHYMDL MUYMDL\n')
        fphaseyT.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
        for i in range(0,len(bpmsyT)):
            bn1=upper(bpmsyT[i][1])
            bns1=bpmsyT[i][0]
            phmdl=phaseyT[bn1][2]
            bn2=upper(bpmsyT[0][1])
            bns2=bpmsyT[0][0]
            fphaseyT.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+ str(phaseyT[bn1][0])+' '+str(phaseyT[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]])+'\n' )
        fphaseyT.close()

        #-- ac to free total phase
        if acswitch=='1':

            #-- from eq
            try:
                [phaseyTf,bpmsyTf]=GetFreePhaseTotal_Eq(MADTwiss,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',bd,lhcphase)
                try:
                    fphaseyTf.write('@ Q1 %le '+str(Q1f)+'\n')
                    fphaseyTf.write('@ MUX %le '+str(muxf)+'\n')
                except:
                    fphaseyTf.write('@ Q1 %le '+'0.0'+'\n')
                    fphaseyTf.write('@ MUX %le '+'0.0'+'\n')
                fphaseyTf.write('@ Q2 %le '+str(Q2f)+'\n')
                fphaseyTf.write('@ MUY %le '+str(muyf)+'\n')
                fphaseyTf.write('* NAME   NAME2  S   S1   COUNT  PHASEY  STDPHY  PHYMDL MUYMDL\n')
                fphaseyTf.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
                for i in range(0,len(bpmsyTf)):
                    bn1=upper(bpmsyTf[i][1])
                    bns1=bpmsyTf[i][0]
                    phmdlf=phaseyTf[bn1][2]
                    bn2=upper(bpmsyTf[0][1])
                    bns2=bpmsyTf[0][0]
                    fphaseyTf.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phaseyTf[bn1][0])+' '+str(phaseyTf[bn1][1])+' '+str(phmdlf)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            except: pass
            fphaseyTf.close()

            #-- from the model
            [phaseyTf2,bpmsyTf2]=getfreephaseTotal(phaseyT,bpmsyT,"V",MADTwiss,MADTwiss_ac)
            try:
                fphaseyTf2.write('@ Q1 %le '+str(Q1f)+'\n')
                fphaseyTf2.write('@ MUX %le '+str(muxf2)+'\n')
            except:
                fphaseyTf2.write('@ Q1 %le '+'0.0'+'\n')
                fphaseyTf2.write('@ MUX %le '+'0.0'+'\n')
            fphaseyTf2.write('@ Q2 %le '+str(Q2f)+'\n')
            fphaseyTf2.write('@ MUY %le '+str(muyf2)+'\n')
            fphaseyTf2.write('* NAME   NAME2  S   S1   COUNT  PHASEY  STDPHY  PHYMDL MUYMDL\n')
            fphaseyTf2.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpmsyTf2)):
                bn1=upper(bpmsyTf2[i][1])
                bns1=bpmsyTf2[i][0]
                phmdlf2=phaseyTf2[bn1][2]
                bn2=upper(bpmsyTf2[0][1])
                bns2=bpmsyTf2[0][0]
                fphaseyTf2.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phaseyTf2[bn1][0])+' '+str(phaseyTf2[bn1][1])+' '+str(phmdlf2)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            fphaseyTf2.close()


    #-------- START Beta
    print 'Calculating beta'
    betaxlist=[]
    betaylist=[]

    #---- H plane
    if wolinx!=1:

        [betax,rmsbbx,alfax,bpms]=BetaFromPhase(MADTwiss_ac,ListOfZeroDPPX,phasex,'H')
        betax['DPP']=0
        betaxlist.append(betax)
        betaxPhaseCopy=betax  #-- For Andy's BetaFromAmp re-scaling
        fbetax.write('@ Q1 %le '+str(Q1)+'\n')
        try:    fbetax.write('@ Q2 %le '+str(Q2)+'\n')
        except: fbetax.write('@ Q2 %le '+'0.0'+'\n')
        fbetax.write('@ RMSbetabeat %le '+str(rmsbbx)+'\n')
        fbetax.write('* NAME   S    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
        fbetax.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fbetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1][0])+' '+str(betax[bn1][1])+' '+str(betax[bn1][2])+' '+str(alfax[bn1][0])+' '+str(alfax[bn1][1])+' '+str(alfax[bn1][2])+' '+str(MADTwiss_ac.BETX[MADTwiss_ac.indx[bn1]])+' '+str(MADTwiss_ac.ALFX[MADTwiss_ac.indx[bn1]])+' '+str(MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]])+'\n' )
        fbetax.close()

        #-- ac to free beta
        if acswitch=='1':

            #-- from eq
            try:
                [betaxf,rmsbbxf,alfaxf,bpmsf]=BetaFromPhase(MADTwiss,ListOfZeroDPPX,phasexf,'H')
                betaxPhaseCopyf=betaxf  #-- For Andy's BetaFromAmp re-scaling
                fbetaxf.write('@ Q1 %le '+str(Q1f)+'\n')
                try:    fbetaxf.write('@ Q2 %le '+str(Q2f)+'\n')
                except: fbetaxf.write('@ Q2 %le '+'0.0'+'\n')
                fbetaxf.write('@ RMSbetabeat %le '+str(rmsbbxf)+'\n')
                fbetaxf.write('* NAME   S    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
                fbetaxf.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
                for i in range(0,len(bpmsf)):
                    bn1=upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    fbetaxf.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betaxf[bn1][0])+' '+str(betaxf[bn1][1])+' '+str(betaxf[bn1][2])+' '+str(alfaxf[bn1][0])+' '+str(alfaxf[bn1][1])+' '+str(alfaxf[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            except: pass
            fbetaxf.close()

            #-- from the model
            [betaxf2,rmsbbxf2,alfaxf2,bpmsf2]=getFreeBeta(MADTwiss_ac,MADTwiss,betax,rmsbbx,alfax,bpms,'H')
            fbetaxf2.write('@ Q1 %le '+str(Q1f)+'\n')
            try:    fbetaxf2.write('@ Q2 %le '+str(Q2f)+'\n')
            except: fbetaxf2.write('@ Q2 %le '+'0.0'+'\n')
            fbetaxf2.write('@ RMSbetabeat %le '+str(rmsbbxf2)+'\n')
            fbetaxf2.write('* NAME   S    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
            fbetaxf2.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
            for i in range(0,len(bpmsf2)):
                bn1=upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                fbetaxf2.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betaxf2[bn1][0])+' '+str(betaxf2[bn1][1])+' '+str(betaxf2[bn1][2])+' '+str(alfaxf2[bn1][0])+' '+str(alfaxf2[bn1][1])+' '+str(alfaxf2[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fbetaxf2.close()

    #---- V plane
    if woliny!=1:

        [betay,rmsbby,alfay,bpms]=BetaFromPhase(MADTwiss_ac,ListOfZeroDPPY,phasey,'V')
        betay['DPP']=0
        betaylist.append(betay)
        betayPhaseCopy=betay  #-- For Andy's BetaFromAmp re-scaling
        try:    fbetay.write('@ Q1 %le '+str(Q1)+'\n')
        except: fbetay.write('@ Q1 %le '+'0.0'+'\n')
        fbetay.write('@ Q2 %le '+str(Q2)+'\n')
        fbetay.write('@ RMSbetabeat %le '+str(rmsbby)+'\n')
        fbetay.write('* NAME   S    COUNT  BETY   ERRBETY STDBETY ALFY   ERRALFY STDALFY BETYMDL ALFYMDL MUYMDL\n')
        fbetay.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fbetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(betay[bn1][2])+' '+str(alfay[bn1][0])+' '+str(alfay[bn1][1])+' '+str(alfay[bn1][2])+' '+str(MADTwiss_ac.BETY[MADTwiss_ac.indx[bn1]])+' '+str(MADTwiss_ac.ALFY[MADTwiss_ac.indx[bn1]])+' '+str(MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]])+'\n' )
        fbetay.close()

        #-- ac to free beta
        if acswitch=='1':

            #-- from eq
            try:
                [betayf,rmsbbyf,alfayf,bpmsf]=BetaFromPhase(MADTwiss,ListOfZeroDPPY,phaseyf,'V')
                betayPhaseCopyf=betayf  #-- For Andy's BetaFromAmp re-scaling
                try:    fbetayf.write('@ Q1 %le '+str(Q1f)+'\n')
                except: fbetayf.write('@ Q1 %le '+'0.0'+'\n')
                fbetayf.write('@ Q2 %le '+str(Q2f)+'\n')
                fbetayf.write('@ RMSbetabeat %le '+str(rmsbbyf)+'\n')
                fbetayf.write('* NAME   S    COUNT  BETY   ERRBETY STDBETY ALFY   ERRALFY STDALFY BETYMDL ALFYMDL MUYMDL\n')
                fbetayf.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
                for i in range(0,len(bpmsf)):
                    bn1=upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    fbetayf.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betayf[bn1][0])+' '+str(betayf[bn1][1])+' '+str(betayf[bn1][2])+' '+str(alfayf[bn1][0])+' '+str(alfayf[bn1][1])+' '+str(alfayf[bn1][2])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            except: pass
            fbetayf.close()

            #-- from the model
            [betayf2,rmsbbyf2,alfayf2,bpmsf2]=getFreeBeta(MADTwiss_ac,MADTwiss,betay,rmsbby,alfay,bpms,'V')
            try:    fbetayf2.write('@ Q1 %le '+str(Q1f)+'\n')
            except: fbetayf2.write('@ Q1 %le '+'0.0'+'\n')
            fbetayf2.write('@ Q2 %le '+str(Q2f)+'\n')
            fbetayf2.write('@ RMSbetabeat %le '+str(rmsbbyf2)+'\n')
            fbetayf2.write('* NAME   S    COUNT  BETY   ERRBETY STDBETY ALFY   ERRALFY STDALFY BETYMDL ALFYMDL MUYMDL\n')
            fbetayf2.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
            for i in range(0,len(bpmsf2)):
                bn1=upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                fbetayf2.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betayf2[bn1][0])+' '+str(betayf2[bn1][1])+' '+str(betayf2[bn1][2])+' '+str(alfayf2[bn1][0])+' '+str(alfayf2[bn1][1])+' '+str(alfayf2[bn1][2])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            fbetayf2.close()

    #------- Start beta from amplitude
    print 'Calculating beta from amplitude'
    betaxalist=[]
    betayalist=[]

    #---- H plane
    if wolinx!=1:

        [betax,rmsbbx,bpms,invJx]=BetaFromAmplitude(MADTwiss_ac,ListOfZeroDPPX,'H')
        betax['DPP']=0
        beta2_save=betax
        betaxalist.append(betax)

        #-- Rescaling
        betax_ratio=0
        skipped_bpmx=[]
        arcbpms=filterbpm(bpms)
        for bpm in arcbpms:
            s=bpm[0] # first entry is position
            name=upper(bpm[1]) # second entry is the name
            #Skip BPM with strange data
            if abs(betaxPhaseCopy[name][0]/betax[name][0])>100: skipped_bpmx.append(name)
            elif (betax[name][0]<0 or betaxPhaseCopy[name][0]<0): skipped_bpmx.append(name)
            else: betax_ratio=betax_ratio+(betaxPhaseCopy[name][0]/betax[name][0])
        try: betax_ratio=betax_ratio/(len(arcbpms)-len(skipped_bpmx))
        except: betax_ratio=1
        betax_rescale={}
        for bpm in map(upper,zip(*bpms)[1]): betax_rescale[bpm]=[betax_ratio*betax[bpm][0],betax_ratio*betax[bpm][1],betax[bpm][2]]

        fabetax.write('@ Q1 %le '+str(Q1)+'\n')
        try:    fabetax.write('@ Q2 %le '+str(Q2)+'\n')
        except: fabetax.write('@ Q2 %le '+'0.0'+'\n')
        fabetax.write('@ RMSbetabeat %le '+str(rmsbbx)+'\n')
        fabetax.write('@ RescalingFactor %le '+str(betax_ratio)+'\n')
        fabetax.write('* NAME   S    COUNT  BETX   BETXSTD BETXMDL MUXMDL BETXRES BETXSTDRES\n')
        fabetax.write('$ %s     %le    %le    %le    %le     %le     %le  %le  %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fabetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1][0])+' '+str(betax[bn1][1])+' '+str(MADTwiss_ac.BETX[MADTwiss_ac.indx[bn1]])+' '+str(MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]])+' '+str(betax_rescale[bn1][0])+' '+str(betax_rescale[bn1][1])+'\n')
        fabetax.close()

        #-- ac to free amp beta
        if acswitch=='1':

            #-- from eq
            try:
                [betaxf,rmsbbxf,bpmsf,invJxf]=GetFreeBetaFromAmp_Eq(MADTwiss_ac,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',bd,lhcphase)

                #-- Rescaling
                betaxf_ratio=0
                skipped_bpmxf=[]
                arcbpms=filterbpm(bpmsf)
                for bpm in arcbpms:
                    s=bpm[0] # first entry is position
                    name=upper(bpm[1]) # second entry is the name
                    #Skip BPM with strange data
                    if abs(betaxPhaseCopyf[name][0]/betaxf[name][0])>10: skipped_bpmxf.append(name)
                    elif abs(betaxPhaseCopyf[name][0]/betaxf[name][0])<0.1: skipped_bpmxf.append(name)
                    elif (betaxf[name][0]<0 or betaxPhaseCopyf[name][0]<0): skipped_bpmxf.append(name)
                    else: betaxf_ratio=betaxf_ratio+(betaxPhaseCopyf[name][0]/betaxf[name][0])
                try: betaxf_ratio=betaxf_ratio/(len(arcbpms)-len(skipped_bpmxf))
                except: betaxf_ratio=1

                fabetaxf.write('@ Q1 %le '+str(Q1f)+'\n')
                try:    fabetaxf.write('@ Q2 %le '+str(Q2f)+'\n')
                except: fabetaxf.write('@ Q2 %le '+'0.0'+'\n')
                fabetaxf.write('@ RMSbetabeat %le '+str(rmsbbxf)+'\n')
                fabetaxf.write('@ RescalingFactor %le '+str(betaxf_ratio)+'\n')
                fabetaxf.write('* NAME   S    COUNT  BETX   BETXSTD BETXMDL MUXMDL BETXRES BETXSTDRES\n')
                fabetaxf.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le\n')
                for i in range(0,len(bpmsf)):
                    bn1=upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    fabetaxf.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betaxf[bn1][0])+' '+str(betaxf[bn1][1])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+' '+str(betaxf_ratio*betaxf[bn1][0])+' '+str(betaxf_ratio*betaxf[bn1][1])+'\n')
            except: pass
            fabetaxf.close()

            #-- from the model
            [betaxf2,rmsbbxf2,bpmsf2,invJxf2]=getFreeAmpBeta(betax,rmsbbx,bpms,invJx,MADTwiss_ac,MADTwiss,'H')
            betaxf2_rescale=getFreeAmpBeta(betax_rescale,rmsbbx,bpms,invJx,MADTwiss_ac,MADTwiss,'H')[0]
            fabetaxf2.write('@ Q1 %le '+str(Q1f)+'\n')
            try:    fabetaxf2.write('@ Q2 %le '+str(Q2f)+'\n')
            except: fabetaxf2.write('@ Q2 %le '+'0.0'+'\n')
            fabetaxf2.write('@ RMSbetabeat %le '+str(rmsbbxf2)+'\n')
            fabetaxf2.write('@ RescalingFactor %le '+str(betax_ratio)+'\n')
            fabetaxf2.write('* NAME   S    COUNT  BETX   BETXSTD BETXMDL MUXMDL BETXRES BETXSTDRES\n')
            fabetaxf2.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le\n')
            for i in range(0,len(bpmsf2)):
                bn1=upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                fabetaxf2.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betaxf2[bn1][0])+' '+str(betaxf2[bn1][1])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+' '+str(betaxf2_rescale[bn1][0])+' '+str(betaxf2_rescale[bn1][1])+'\n')
            fabetaxf2.close()

    #---- V plane
    if woliny!=1:

        [betay,rmsbby,bpms,invJy]=BetaFromAmplitude(MADTwiss_ac,ListOfZeroDPPY,'V')
        betay['DPP']=0
        betayalist.append(betay)

        #-- Rescaling
        betay_ratio=0
        skipped_bpmy=[]
        arcbpms=filterbpm(bpms)
        for bpm in arcbpms:
            s=bpm[0] # first entry is position
            name=upper(bpm[1]) # second entry is the name
            #Skip BPM with strange data
            if abs(betayPhaseCopy[name][0]/betay[name][0])>100: skipped_bpmy.append(name)
            elif (betay[name][0]<0 or betayPhaseCopy[name][0]<0): skipped_bpmy.append(name)
            else: betay_ratio=betay_ratio+(betayPhaseCopy[name][0]/betay[name][0])
        try: betay_ratio=betay_ratio/(len(arcbpms)-len(skipped_bpmy))
        except: betay_ratio=1
        betay_rescale={}
        for bpm in map(upper,zip(*bpms)[1]): betay_rescale[bpm]=[betay_ratio*betay[bpm][0],betay_ratio*betay[bpm][1],betay[bpm][2]]

        try:    fabetay.write('@ Q1 %le '+str(Q1)+'\n')
        except: fabetay.write('@ Q1 %le '+'0.0'+'\n')
        fabetay.write('@ Q2 %le '+str(Q2)+'\n')
        fabetay.write('@ RMSbetabeat %le '+str(rmsbby)+'\n')
        fabetay.write('@ RescalingFactor %le '+str(betay_ratio)+'\n')
        fabetay.write('* NAME   S    COUNT  BETY   BETYSTD BETYMDL MUYMDL BETYRES BETYSTDRES\n')
        fabetay.write('$ %s     %le    %le    %le    %le     %le     %le   %le   %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fabetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(MADTwiss_ac.BETY[MADTwiss_ac.indx[bn1]])+' '+str(MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]])+' '+str(betay_rescale[bn1][0])+' '+str(betay_rescale[bn1][1])+'\n')
        fabetay.close()

        #-- ac to free amp beta
        if acswitch=='1':

            #-- from eq
            try:
                [betayf,rmsbbyf,bpmsf,invJyf]=GetFreeBetaFromAmp_Eq(MADTwiss_ac,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',bd,accel)

                #-- Rescaling
                betayf_ratio=0
                skipped_bpmyf=[]
                arcbpms=filterbpm(bpmsf)
                for bpm in arcbpms:
                    s=bpm[0] # first entry is position
                    name=upper(bpm[1]) # second entry is the name
                    #Skip BPM with strange data
                    if abs(betayPhaseCopyf[name][0]/betayf[name][0])>10: skipped_bpmyf.append(name)
                    elif (betayf[name][0]<0 or betayPhaseCopyf[name][0]<0): skipped_bpmyf.append(name)
                    elif abs(betayPhaseCopyf[name][0]/betayf[name][0])<0.1: skipped_bpmyf.append(name)
                    else: betayf_ratio=betayf_ratio+(betayPhaseCopyf[name][0]/betayf[name][0])
                try: betayf_ratio=betayf_ratio/(len(arcbpms)-len(skipped_bpmyf))
                except: betayf_ratio=1

                try:    fabetayf.write('@ Q1 %le '+str(Q1f)+'\n')
                except: fabetayf.write('@ Q1 %le '+'0.0'+'\n')
                fabetayf.write('@ Q2 %le '+str(Q2f)+'\n')
                fabetayf.write('@ RMSbetabeat %le '+str(rmsbbyf)+'\n')
                fabetayf.write('@ RescalingFactor %le '+str(betayf_ratio)+'\n')
                fabetayf.write('* NAME   S    COUNT  BETY   BETYSTD BETYMDL MUYMDL BETYRES BETYSTDRES\n')
                fabetayf.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le\n')
                for i in range(0,len(bpmsf)):
                    bn1=upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    fabetayf.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betayf[bn1][0])+' '+str(betayf[bn1][1])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+' '+str(betayf_ratio*betayf[bn1][0])+' '+str(betayf_ratio*betayf[bn1][1])+'\n')
            except: pass
            fabetayf.close()

            #-- from the model
            [betayf2,rmsbbyf2,bpmsf2,invJyf2]=getFreeAmpBeta(betay,rmsbby,bpms,invJy,MADTwiss_ac,MADTwiss,'V')
            betayf2_rescale=getFreeAmpBeta(betay_rescale,rmsbby,bpms,invJy,MADTwiss_ac,MADTwiss,'V')[0]
            try:    fabetayf2.write('@ Q1 %le '+str(Q1f)+'\n')
            except: fabetayf2.write('@ Q1 %le '+'0.0'+'\n')
            fabetayf2.write('@ Q2 %le '+str(Q2f)+'\n')
            fabetayf2.write('@ RMSbetabeat %le '+str(rmsbbyf2)+'\n')
            fabetayf2.write('@ RescalingFactor %le '+str(betay_ratio)+'\n')
            fabetayf2.write('* NAME   S    COUNT  BETY   BETYSTD BETYMDL MUYMDL BETYRES BETYSTDRES\n')
            fabetayf2.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le\n')
            for i in range(0,len(bpmsf2)):
                bn1=upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                fabetayf2.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betayf2[bn1][0])+' '+str(betayf2[bn1][1])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+' '+str(betayf2_rescale[bn1][0])+' '+str(betayf2_rescale[bn1][1])+'\n')
            fabetayf2.close()

    #-------- START IP
    print 'Calculating IP'

    if "LHC" in accel:
        ips=["1","2","3","4","5","6","7","8"]
        try:
            measured=[betax,betay]
            phases=[phasex,phasey]
            bpmss=[bpmsx,bpmsy]
        except:
            pass
        for ip in ips:
            try:
                #print "entering"
                betahor,betaver=getIP(ip,measured,MADTwiss,phases,bpmss)
            except:
                betahor=[0,0,0,0,0,0,0]
                betaver=[0,0,0,0,0,0,0]
            #print str(betahor[6])
            fIP.write("\"IP"+ip+"\" "+str(betahor[1])+" "+str(betahor[4])+" "+str(betahor[2])+" "+str(betahor[3])+" "+str(betahor[6])+" "+str(betahor[5])+" "+str(betaver[1])+" "+str(betaver[4])+" "+str(betaver[2])+" "+str(betaver[3])+" "+str(betaver[6])+" "+str(betaver[5])+"\n")

        fIP.close()

        #-- Parameters at IP1, IP2, IP5, and IP8
        IPx=GetIP2(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',bd,accel,lhcphase)
        IPy=GetIP2(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',bd,accel,lhcphase)
        fIPx.write('* NAME  BETX  BETXSTD  BETXMDL  ALFX  ALFXSTD  ALFXMDL  BETX*  BETX*STD  BETX*MDL  SX*  SX*STD  SX*MDL  rt(2JX)  rt(2JX)STD\n')
        fIPx.write(" ".join(["$","%s","%le","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le",  "%le"]))
        fIPx.write('\n')
        fIPy.write('* NAME  BETY  BETYSTD  BETYMDL  ALFY  ALFYSTD  ALFYMDL  BETY*  BETY*STD  BETY*MDL  SY*  SY*STD  SY*MDL  rt(2JY)  rt(2JY)STD\n')
        fIPy.write(" ".join(["$","%s","%le","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le",  "%le"]))
        fIPy.write('\n')
        #FIXME: weiter bei getipfromphase
        for i in ('IP1','IP5','IP8','IP2'):
            try:
                fIPx.write('"'+i+'"'+' ')
                for k in IPx[i]: fIPx.write(str(k)+' ')
                fIPx.write('\n')
            except: pass
            try:
                fIPy.write('"'+i+'"'+' ')
                for k in IPy[i]: fIPy.write(str(k)+' ')
                fIPy.write('\n')
            except: pass
        fIPx.close()
        fIPy.close()

        #-- ac to free parameters at IP1, IP2, IP5, and IP8
        if acswitch=='1':

            #-- From Eq
            IPxf=GetFreeIP2_Eq(MADTwiss,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',bd,accel,lhcphase)
            IPyf=GetFreeIP2_Eq(MADTwiss,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',bd,accel,lhcphase)
            fIPxf.write('* NAME  BETX  BETXSTD  BETXMDL  ALFX  ALFXSTD  ALFXMDL  BETX*  BETX*STD  BETX*MDL  SX*  SX*STD  SX*MDL  rt(2JXD)  rt(2JXD)STD\n')
            fIPxf.write(" ".join(["$","%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le","\n"]))
            fIPyf.write('* NAME  BETY  BETYSTD  BETYMDL  ALFY  ALFYSTD  ALFYMDL  BETY*  BETY*STD  BETY*MDL  SY*  SY*STD  SY*MDL  rt(2JYD)  rt(2JYD)STD\n')
            fIPyf.write(" ".join(["$","%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le","\n"]))
            for i in ('IP1','IP5','IP8','IP2'):
                try:
                    fIPxf.write('"'+i+'"'+' ')
                    for k in IPxf[i]: fIPxf.write(str(k)+' ')
                    fIPxf.write('\n')
                except: pass
                try:
                    fIPyf.write('"'+i+'"'+' ')
                    for k in IPyf[i]: fIPyf.write(str(k)+' ')
                    fIPyf.write('\n')
                except: pass
            fIPxf.close()
            fIPyf.close()

            #-- From model
            IPxf2=GetFreeIP2(MADTwiss,MADTwiss_ac,IPx,'H',accel)
            IPyf2=GetFreeIP2(MADTwiss,MADTwiss_ac,IPy,'V',accel)
            fIPxf2.write('* NAME  BETX  BETXSTD  BETXMDL  ALFX  ALFXSTD  ALFXMDL  BETX*  BETX*STD  BETX*MDL  SX*  SX*STD  SX*MDL  rt(2JXD)  rt(2JXD)STD\n')
            fIPxf2.write(" ".join(["$","%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le","\n"]))
            fIPyf2.write('* NAME  BETY  BETYSTD  BETYMDL  ALFY  ALFYSTD  ALFYMDL  BETY*  BETY*STD  BETY*MDL  SY*  SY*STD  SY*MDL  rt(2JYD)  rt(2JYD)STD\n')
            fIPyf2.write(" ".join(["$","%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le","\n"]))
            for i in ('IP1','IP5','IP8','IP2'):
                try:
                    fIPxf2.write('"'+i+'"'+' ')
                    for k in IPxf2[i]: fIPxf2.write(str(k)+' ')
                    fIPxf2.write('\n')
                except: pass
                try:
                    fIPyf2.write('"'+i+'"'+' ')
                    for k in IPyf2[i]: fIPyf2.write(str(k)+' ')
                    fIPyf2.write('\n')
                except: pass

            fIPxf2.close()
            fIPyf2.close()

            GetFreeIP2(MADTwiss,MADTwiss_ac,IPx,'H',accel)

        #-- IP beta* and phase from phase only
        try:    IPfromphase=GetIPFromPhase(MADTwiss_ac,phasex,phasey,accel)
        except: print 'No output from IP from phase. H or V file missing?'
        fIPfromphase.write('* NAME  2L  BETX*  BETX*STD  BETX*MDL  BETY*  BETY*STD  BETY*MDL  PHX  PHXSTD  PHXMDL  PHY  PHYSTD  PHYMDL\n')
        fIPfromphase.write(" ".join(["$","%s","%le","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le","\n"]))
        for i in ('IP1','IP5','IP8','IP2'):
            fIPfromphase.write('"'+i+'"'+' ')
            try:
                for k in IPfromphase[i]: fIPfromphase.write(str(k)+' ')
                fIPfromphase.write('\n')
            except: fIPfromphase.write('\n')
        fIPfromphase.close()

        #-- ac to free beta*
        if acswitch=='1':

            #-- from eqs
            try:    IPfromphasef=GetIPFromPhase(MADTwiss,phasexf,phaseyf,accel)
            except: pass
            fIPfromphasef.write('* NAME  2L  BETX*  BETX*STD  BETX*MDL  BETY*  BETY*STD  BETY*MDL  PHX  PHXSTD  PHXMDL  PHY  PHYSTD  PHYMDL\n')
            fIPfromphasef.write(" ".join(["$","%s","%le","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le","\n"]))
            for i in ('IP1','IP5','IP8','IP2'):
                fIPfromphasef.write('"'+i+'"'+' ')
                try:
                    for k in IPfromphasef[i]: fIPfromphasef.write(str(k)+' ')
                    fIPfromphasef.write('\n')
                except: fIPfromphasef.write('\n')
            fIPfromphasef.close()

            #-- from the model
            try:    IPfromphasef2=GetIPFromPhase(MADTwiss,phasexf2,phaseyf2,accel)
            except: pass
            fIPfromphasef2.write('* NAME  2L  BETX*  BETX*STD  BETX*MDL  BETY*  BETY*STD  BETY*MDL  PHX  PHXSTD  PHXMDL  PHY  PHYSTD  PHYMDL\n')
            fIPfromphasef2.write(" ".join(["$","%s","%le","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le","\n"]))
            for i in ('IP1','IP5','IP8','IP2'):
                fIPfromphasef2.write('"'+i+'"'+' ')
                try:
                    for k in IPfromphasef2[i]: fIPfromphasef2.write(str(k)+' ')
                    fIPfromphasef2.write('\n')
                except: fIPfromphasef2.write('\n')
            fIPfromphasef2.close()

    #-------- START Orbit
    ListOfCOX=[]
    if wolinx!=1:

        [cox,bpms]=GetCO(MADTwiss, ListOfZeroDPPX)
        # The output file can be directly used for orbit correction with MADX
        fcox.write('@ TABLE %05s "ORBIT"\n')
        fcox.write('@ TYPE %05s "ORBIT"\n')
        fcox.write('@ SEQUENCE %05s "'+accel+'"\n')
        fcox.write('@ Q1 %le '+str(Q1)+'\n')
        try:    fcox.write('@ Q2 %le '+str(Q2)+'\n')
        except: fcox.write('@ Q2 %le '+'0.0'+'\n')
        fcox.write('* NAME   S   COUNT  X      STDX   XMDL   MUXMDL\n')
        fcox.write('$ %s     %le    %le    %le    %le    %le    %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fcox.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(cox[bn1][0])+' '+str(cox[bn1][1])+' '+str(MADTwiss.X[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )

        ListOfCOX.append(cox)
    fcox.close()




    ListOfCOY=[]
    if woliny!=1:
        [coy,bpms]=GetCO(MADTwiss, ListOfZeroDPPY)
        # The output file can be directly used for orbit correction with MADX
        fcoy.write('@ TABLE %05s "ORBIT"\n')
        fcoy.write('@ TYPE %05s "ORBIT"\n')
        fcoy.write('@ SEQUENCE %05s "'+accel+'"\n')
        try:    fcoy.write('@ Q1 %le '+str(Q1)+'\n')
        except: fcoy.write('@ Q1 %le '+'0.0'+'\n')
        fcoy.write('@ Q2 %le '+str(Q2)+'\n')
        fcoy.write('* NAME   S   COUNT  Y      STDY   YMDL   MUYMDL\n')
        fcoy.write('$ %s     %le    %le    %le    %le    %le    %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fcoy.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(coy[bn1][0])+' '+str(coy[bn1][1])+' '+str(MADTwiss.Y[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )


        ListOfCOY.append(coy)
    fcoy.close()



    #-------- Orbit for non-zero DPP
    if wolinx2!=1:

        k=0
        for j in ListOfNonZeroDPPX:
            SingleFile=[]
            SingleFile.append(j)
            file1=outputpath+'getCOx_dpp_'+str(k+1)+'.out'
            fcoDPP=open(file1,'w')
            fcoDPP.write('@ MAD_FILE: %s "'+twiss_model_file+'"'+'\n')
            fcoDPP.write('@ FILE %s "')
            fcoDPP.write(FileOfNonZeroDPPX[k]+' "'+'\n')
            fcoDPP.write('@ DPP %le '+str(float(j.DPP))+'\n')
            try:
                fcoDPP.write('@ Q1 %le '+str(Q1)+'\n')
            except:
                fcoDPP.write('@ Q1 %le '+'0.0'+'\n')
            try:
                fcoDPP.write('@ Q2 %le '+str(Q2)+'\n')
            except:
                fcoDPP.write('@ Q2 %le '+'0.0'+'\n')
            [codpp,bpms]=GetCO(MADTwiss, SingleFile)
            fcoDPP.write('* NAME   S   COUNT  X      STDX   XMDL   MUXMDL\n')
            fcoDPP.write('$ %s     %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                fcoDPP.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(codpp[bn1][0])+' '+str(codpp[bn1][1])+' '+str(MADTwiss.X[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fcoDPP.close()
            ListOfCOX.append(codpp)
            k+=1

    if woliny2!=1:
        k=0
        for j in ListOfNonZeroDPPY:
            SingleFile=[]
            SingleFile.append(j)
            file1=outputpath+'getCOy_dpp_'+str(k+1)+'.out'
            fcoDPP=open(file1,'w')
            fcoDPP.write('@ MAD_FILE: %s "'+twiss_model_file+'"'+'\n')
            fcoDPP.write('@ FILE %s "')
            fcoDPP.write(FileOfNonZeroDPPY[k]+' "'+'\n')
            fcoDPP.write('@ DPP %le '+str(float(j.DPP))+'\n')
            try:
                fcoDPP.write('@ Q1 %le '+str(Q1)+'\n')
            except:
                fcoDPP.write('@ Q1 %le '+'0.0'+'\n')
            try:
                fcoDPP.write('@ Q2 %le '+str(Q2)+'\n')
            except:
                fcoDPP.write('@ Q2 %le '+'0.0'+'\n')
            [codpp,bpms]=GetCO(MADTwiss, SingleFile)
            fcoDPP.write('* NAME   S   COUNT  Y      STDY   YMDL   MUYMDL\n')
            fcoDPP.write('$ %s     %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                fcoDPP.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(codpp[bn1][0])+' '+str(codpp[bn1][1])+' '+str(MADTwiss.Y[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            fcoDPP.close()
            ListOfCOY.append(codpp)
            k+=1



    #-------- START Dispersion

    if wolinx!=1 and wolinx2!=1:


        [nda,Dx,DPX,bpms]=NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, beta2_save, COcut)
        fNDx.write('@ Q1 %le '+str(Q1)+'\n')
        try:
            fNDx.write('@ Q2 %le '+str(Q2)+'\n')
        except:
            fNDx.write('@ Q2 %le '+'0.0'+'\n')
        fNDx.write('* NAME   S    COUNT  NDX    STDNDX DX     DPX    NDXMDL DXMDL  DPXMDL MUXMDL\n')
        fNDx.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            ndmdl=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
            fNDx.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPX))+' '+str(nda[bn1][0])+' '+str(nda[bn1][1])+' '+str(Dx[bn1][0])+' '+str(DPX[bn1])+' '+str(ndmdl)+' '+str(MADTwiss.DX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.DPX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )




        [dxo,bpms]=DispersionfromOrbit(ListOfZeroDPPX,ListOfNonZeroDPPX,ListOfCOX,COcut,BPMU)

        DPX=GetDPX(MADTwiss,dxo,bpms)
        fDx.write('@ Q1 %le '+str(Q1)+'\n')
        fDx.write('@ Q2 %le '+str(Q2)+'\n')
        fDx.write('* NAME   S    COUNT  DX     STDDX  DPX    DXMDL  DPXMDL MUXMDL\n')
        fDx.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le\n')
        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fDx.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPX))+' '+str(dxo[bn1][0])+' '+str(dxo[bn1][1])+' '+str(DPX[bn1])+' '+str(MADTwiss.DX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.DPX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )


    fNDx.close()
    fDx.close()



    if woliny!=1 and woliny2!=1:
        [dyo,bpms]=DispersionfromOrbit(ListOfZeroDPPY,ListOfNonZeroDPPY,ListOfCOY,COcut,BPMU)
        DPY=GetDPY(MADTwiss,dyo,bpms)
        fDy.write('@ Q1 %le '+str(Q1)+'\n')
        fDy.write('@ Q2 %le '+str(Q2)+'\n')
        fDy.write('* NAME   S    COUNT  DY     STDDY  DPY  DYMDL  DPYMDL  MUYMDL\n')
        fDy.write('$ %s     %le    %le    %le    %le    %le  %le   %le    %le\n')

        for i in range(0,len(bpms)):
            bn1=upper(bpms[i][1])
            bns1=bpms[i][0]
            fDy.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPY))+' '+str(dyo[bn1][0])+' '+str(dyo[bn1][1])+' '+str(DPY[bn1])+' '+str(MADTwiss.DY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.DPY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

    fDy.close()

    #-------- START coupling.
    print "Calculating coupling"

    if wolinx!=1 and woliny!=1:

        #-- Coupling in the model
        try:    MADTwiss.Cmatrix()
        except: pass

        #-- Main part
        if   NBcpl==1:

            [fwqw,bpms]=GetCoupling1(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2)
            fcouple.write('@ CG %le '+str(fwqw['Global'][0])+'\n')
            fcouple.write('@ QG %le '+str(fwqw['Global'][1])+'\n')
            fcouple.write('* NAME   S    COUNT  F1001W FWSTD  Q1001W QWSTD MDLF1001R MDLF1001I\n')
            fcouple.write('$ %s     %le    %le    %le    %le    %le    %le   %le       %le\n')
            for i in range(len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                try:    fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(MADTwiss.f1001[MADTwiss.indx(bn1)].real)+' '+str(MADTwiss.f1001[MADTwiss.indx(bn1)].imag)+' '+str(MADTwiss_ac.f1010[MADTwiss_ac.indx(bn1)].real)+' '+str(MADTwiss_ac.f1010[MADTwiss_ac.indx(bn1)].imag)+'\n')
                #-- Output zero if the model does not have couping parameters
                except: fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' 0.0 0.0'+'\n')
            fcouple.close()

        elif NBcpl==2:

            if accel=="SPS" or "RHIC" in accel:
                [phasexp,Q1,MUX,bpmsx]=GetPhases(MADTwiss,PseudoListX,'H',outputpath,bd,accel,lhcphase)
                [phaseyp,Q2,MUY,bpmsy]=GetPhases(MADTwiss,PseudoListY,'V',outputpath,bd,accel,lhcphase)
                #[fwqw,bpms]=GetCoupling2(MADTwiss,PseudoListX,PseudoListY,Q1,Q2,phasexp,phaseyp,bd,accel)
                [fwqw,bpms]=GetCoupling2b(MADTwiss,PseudoListX,PseudoListY,Q1,Q2,phasexp,phaseyp,bd,accel)
            else:
                #[fwqw,bpms]=GetCoupling2(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2,phasexlist[0],phaseylist[0],bd,accel)
                [fwqw,bpms]=GetCoupling2b(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2,phasexlist[0],phaseylist[0],bd,accel)
            fcouple.write('@ CG %le '+str(fwqw['Global'][0])+'\n')
            fcouple.write('@ QG %le '+str(fwqw['Global'][1])+'\n')
            fcouple.write('* NAME   S    COUNT  F1001W FWSTD1 F1001R F1001I F1010W FWSTD2 F1010R F1010I MDLF1001R MDLF1001I MDLF1010R MDLF1010I\n')
            fcouple.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le       %le       %le       %le\n')
            for i in range(len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                try:    fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2))+' '+str(fwqw[bn1][0][3])+' '+str(fwqw[bn1][0][2].real)+' '+str(fwqw[bn1][0][2].imag)+' '+str(MADTwiss_ac.f1001[MADTwiss_ac.indx[bn1]].real)+' '+str(MADTwiss_ac.f1001[MADTwiss_ac.indx[bn1]].imag)+' '+str(MADTwiss_ac.f1010[MADTwiss_ac.indx[bn1]].real)+' '+str(MADTwiss_ac.f1010[MADTwiss_ac.indx[bn1]].imag)+'\n')
                #-- Output zero if the model does not have couping parameters
                except: fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2))+' '+str(fwqw[bn1][0][3])+' '+str(fwqw[bn1][0][2].real)+' '+str(fwqw[bn1][0][2].imag)+' 0.0 0.0 0.0 0.0'+'\n')
            fcouple.close()

            #-- ac to free coupling
            if acswitch=="1":

                #-- analytic eqs
                try:
                    [fwqwf,bpmsf]=GetFreeCoupling_Eq(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2,Q1f,Q2f,acphasex_ac2bpmac,acphasey_ac2bpmac,bd)
                    fcouplef.write('@ CG %le '+str(fwqw['Global'][0])+'\n')
                    fcouplef.write('@ QG %le '+str(fwqw['Global'][1])+'\n')
                    fcouplef.write('* NAME   S    COUNT  F1001W FWSTD1 F1001R F1001I F1010W FWSTD2 F1010R F1010I MDLF1001R MDLF1001I MDLF1010R MDLF1010I\n')
                    fcouplef.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le       %le       %le       %le\n')
                    for i in range(len(bpmsf)):
                        bn1=upper(bpmsf[i][1])
                        bns1=bpmsf[i][0]
                        try:    fcouplef.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqwf[bn1][0][0].real**2+fwqwf[bn1][0][0].imag**2))+' '+str(fwqwf[bn1][0][1])+' '+str(fwqwf[bn1][0][0].real)+' '+str(fwqwf[bn1][0][0].imag)+' '+str(sqrt(fwqwf[bn1][0][2].real**2+fwqwf[bn1][0][2].imag**2))+' '+str(fwqwf[bn1][0][3])+' '+str(fwqwf[bn1][0][2].real)+' '+str(fwqwf[bn1][0][2].imag)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].imag)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].imag)+'\n')
                        #-- Output zero if the model does not have couping parameters
                        except: fcouplef.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqwf[bn1][0][0].real**2+fwqwf[bn1][0][0].imag**2))+' '+str(fwqwf[bn1][0][1])+' '+str(fwqwf[bn1][0][0].real)+' '+str(fwqwf[bn1][0][0].imag)+' '+str(sqrt(fwqwf[bn1][0][2].real**2+fwqwf[bn1][0][2].imag**2))+' '+str(fwqwf[bn1][0][3])+' '+str(fwqwf[bn1][0][2].real)+' '+str(fwqwf[bn1][0][2].imag)+' 0.0 0.0 0.0 0.0'+'\n')
                except: pass
                fcouplef.close()

                #-- global factor
                [fwqwf2,bpmsf2]=getFreeCoupling(Q1f,Q2f,Q1,Q2,fwqw,MADTwiss,bpms)
                fcouplef2.write('@ CG %le '+str(fwqw['Global'][0])+'\n')
                fcouplef2.write('@ QG %le '+str(fwqw['Global'][1])+'\n')
                fcouplef2.write('* NAME   S    COUNT  F1001W FWSTD1 F1001R F1001I F1010W FWSTD2 F1010R F1010I MDLF1001R MDLF1001I MDLF1010R MDLF1010I\n')
                fcouplef2.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le       %le       %le       %le\n')
                for i in range(len(bpmsf2)):
                    bn1=upper(bpmsf2[i][1])
                    bns1=bpmsf2[i][0]
                    try:    fcouplef2.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqwf2[bn1][0][0].real**2+fwqwf2[bn1][0][0].imag**2))+' '+str(fwqwf2[bn1][0][1])+' '+str(fwqwf2[bn1][0][0].real)+' '+str(fwqwf2[bn1][0][0].imag)+' '+str(sqrt(fwqwf2[bn1][0][2].real**2+fwqwf2[bn1][0][2].imag**2))+' '+str(fwqwf2[bn1][0][3])+' '+str(fwqwf2[bn1][0][2].real)+' '+str(fwqwf2[bn1][0][2].imag)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].imag)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].imag)+'\n')
                    #-- Output zero if the model does not have couping parameters
                    except: fcouplef2.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqwf2[bn1][0][0].real**2+fwqwf2[bn1][0][0].imag**2))+' '+str(fwqwf2[bn1][0][1])+' '+str(fwqwf2[bn1][0][0].real)+' '+str(fwqwf2[bn1][0][0].imag)+' '+str(sqrt(fwqwf2[bn1][0][2].real**2+fwqwf2[bn1][0][2].imag**2))+' '+str(fwqwf2[bn1][0][3])+' '+str(fwqwf2[bn1][0][2].real)+' '+str(fwqwf2[bn1][0][2].imag)+' 0.0 0.0 0.0 0.0'+'\n')
                fcouplef2.close()

        else:
            raise ValueError('Number of monitors for coupling analysis should be 1 or 2 (option -n)')

        #-- Convert to C-matrix:
        if acswitch=='1':
            try:    [coupleterms,Qminav,Qminerr,bpms]=getCandGammaQmin(fwqwf,bpmsf,Q1f,Q2f,MADTwiss)
            except: [coupleterms,Qminav,Qminerr,bpms]=getCandGammaQmin(fwqwf2,bpmsf2,Q1f,Q2f,MADTwiss)
        else: [coupleterms,Qminav,Qminerr,bpms]=getCandGammaQmin(fwqw,bpms,Q1f,Q2f,MADTwiss)
        print >> fcoupleterms,"@ DQMIN %le ",Qminav
        print >> fcoupleterms,"@ DQMINE %le ",Qminerr
        print >> fcoupleterms,"* NAME S DETC DETCE  GAMMA GAMMAE C11 C12 C21 C22"
        print >> fcoupleterms,"$ %s %le %le %le %le %le %le %le %le %le"
        for bpm in bpms:
            bps=bpm[0]
            bpmm=bpm[1].upper()
            print >> fcoupleterms,bpmm,bps,coupleterms[bpmm][0],coupleterms[bpmm][1],coupleterms[bpmm][2],coupleterms[bpmm][3],coupleterms[bpmm][4],coupleterms[bpmm][5],coupleterms[bpmm][6],coupleterms[bpmm][7]
        fcoupleterms.close()

        #-- For chromatic coupling
        fwqw['DPP']=0
        couplelist=[fwqw]

    #-------- Phase, Beta and coupling for non-zero DPP

    print "Phase and Beta for non-zero DPP"
    print "lenght of zerothingie "+ str(len(ListOfNonZeroDPPX))
    print "lenght of zerothingie "+ str(len(ListOfNonZeroDPPY))


    if wolinx2!=1:
        plane='H'
        k=0
        for j in ListOfNonZeroDPPX:
            dpop=float(j.DPP)
            SingleFile=[]
            SingleFile.append(j)
            file1=outputpath+'getphasex_dpp_'+str(k+1)+'.out'
            fphDPP=open(file1,'w')
            fphDPP.write('@ MAD_FILE: %s "'+twiss_model_file+'"'+'\n')
            fphDPP.write('@ FILE %s "'+FileOfNonZeroDPPX[k]+'"'+'\n')
            fphDPP.write('@ DPP %le '+str(dpop)+'\n')
            try:
                fphDPP.write('@ Q1 %le '+str(Q1)+'\n')
            except:
                fphDPP.write('@ Q1 %le '+'0.0'+'\n')
            try:
                fphDPP.write('@ Q2 %le '+str(Q2)+'\n')
            except:
                fphDPP.write('@ Q2 %le '+'0.0'+'\n')
            DPPTwiss=ConstructOffMomentumModel(MADTwiss,dpop,BPMdictionary)
            [phasex,Q1DPP,MUX,bpms]=GetPhases(DPPTwiss,SingleFile,Q1,plane,outputpath,bd,accel,lhcphase)
            phasex['DPP']=dpop
            phasexlist.append(phasex)
            fphDPP.write('@ Q1DPP %le '+str(Q1DPP)+'\n')
            fphDPP.write('* NAME   NAME2  S   S1   COUNT  PHASE  STDPH  PHXMDL MUXMDL\n')
            fphDPP.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                if i==len(bpms)-1:
                    bn2=upper(bpms[0][1])
                    bns2=bpms[0][0]

                else:
                    bn2=upper(bpms[i+1][1])
                    bns2=bpms[i+1][0]
                try:
                    phmdl=phasexlist[0][bn1][4]
                except:
                    phmdl=0.0
                #phmdl=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
                fphDPP.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(1)+' '+str(phasex[bn1][0])+' '+str(phasex[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fphDPP.close()

            betax={}
            alfax={}
            rmsbbx=0.
            [betax,rmsbbx,alfax,bpms]=BetaFromPhase(MADTwiss,SingleFile,phasex,plane)
            betax['DPP']=dpop
            betaxlist.append(betax)
            betaxa={}
            [betaxa,rmsbbx,bpms,invJx]=BetaFromAmplitude(MADTwiss,SingleFile,plane)
            betaxa['DPP']=dpop
            betaxalist.append(betaxa)
            file2=outputpath+'getbetax_dpp_'+str(k+1)+'.out'
            fbetaxDPP=open(file2,'w')
            fbetaxDPP.write('@ MAD_FILE: %s "'+twiss_model_file+'"'+'\n')
            fbetaxDPP.write('@ FILE %s "'+FileOfNonZeroDPPX[k]+'"'+'\n')
            fbetaxDPP.write('@ DPP %le '+str(dpop)+'\n')
            try:
                fbetaxDPP.write('@ Q1 %le '+str(Q1)+'\n')
            except:
                fbetaxDPP.write('@ Q1 %le '+'0.0'+'\n')
            try:
                fbetaxDPP.write('@ Q2 %le '+str(Q2)+'\n')
            except:
                fbetaxDPP.write('@ Q2 %le '+'0.0'+'\n')
            #fbetaxDPP.write('@ RMSbetabeat %le '+str(rmsbbx)+'\n')
            fbetaxDPP.write('* NAME   S    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
            fbetaxDPP.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
            for i in range(0,len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                fbetaxDPP.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1][0])+' '+str(betax[bn1][1])+' '+str(betax[bn1][2])+' '+str(alfax[bn1][0])+' '+str(alfax[bn1][1])+' '+str(alfax[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fbetaxDPP.close()
            k+=1


    if woliny2!=1:
        plane='V'
        k=0

        for j in ListOfNonZeroDPPY:
            dpop=float(j.DPP)
            SingleFile=[]
            SingleFile.append(j)
            file1=outputpath+'getphasey_dpp_'+str(k+1)+'.out'
            fphDPP=open(file1,'w')
            fphDPP.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
            fphDPP.write('@ FILE %s "')
            fphDPP.write(FileOfNonZeroDPPY[k]+' ')
            fphDPP.write('"'+'\n')
            try:
                fphDPP.write('@ Q1 %le '+str(Q1)+'\n')
            except:
                fphDPP.write('@ Q1 %le '+'0.0'+'\n')
            try:
                fphDPP.write('@ Q2 %le '+str(Q2)+'\n')
            except:
                fphDPP.write('@ Q2 %le '+'0.0'+'\n')
            DPPTwiss=ConstructOffMomentumModel(MADTwiss,dpop,BPMdictionary)

            [phasey,Q2DPP,MUY,bpms]=GetPhases(DPPTwiss,SingleFile,Q2,plane,outputpath,bd,accel,lhcphase)
            phasey['DPP']=dpop
            phaseylist.append(phasey)
            fphDPP.write('@ Q2DPP %le '+str(Q2DPP)+'\n')
            fphDPP.write('* NAME   NAME2  S   S1  COUNT  PHASE  STDPH  PHYMDL MUYMDL\n')
            fphDPP.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
            for i in range(0,len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                if i==len(bpms)-1:
                    bn2=upper(bpms[0][1])
                    bns2=bpms[0][0]

                else:
                    bn2=upper(bpms[i+1][1])
                    bns2=bpms[i+1][0]
                try:
                    phmdl=phaseylist[0][bn1][4]
                except:
                    phmdl=0.0
                #phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
                fphDPP.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(1)+' '+str(phasey[bn1][0])+' '+str(phasey[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
            fphDPP.close()

            betay={}
            alfay={}
            rmsbby=0.
            [betay,rmsbby,alfay,bpms]=BetaFromPhase(DPPTwiss,SingleFile,phasey,plane)
            betay['DPP']=dpop
            betaylist.append(betay)
            betaya={}
            [betaya,rmsbby,bpms,invJy]=BetaFromAmplitude(DPPTwiss,SingleFile,plane)
            betaya['DPP']=dpop
            betayalist.append(betaya)
            file2=outputpath+'getbetay_dpp_'+str(k+1)+'.out'
            fbetayDPP=open(file2,'w')
            fbetayDPP.write('@ MAD_FILE: %s "'+twiss_model_file+'"'+'\n')
            fbetayDPP.write('@ FILE %s "'+FileOfNonZeroDPPY[k]+'"'+'\n')
            fbetayDPP.write('@ DPP %le '+str(dpop)+'\n')
            try:
                fbetayDPP.write('@ Q1 %le '+str(Q1)+'\n')
            except:
                fbetayDPP.write('@ Q1 %le '+'0.0'+'\n')
            try:
                fbetayDPP.write('@ Q2 %le '+str(Q2)+'\n')
            except:
                fbetayDPP.write('@ Q2 %le '+'0.0'+'\n')
            #fbetayDPP.write('@ RMSbetabeat %le '+str(rmsbbx)+'\n')
            fbetayDPP.write('* NAME   S    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
            fbetayDPP.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
            for i in range(0,len(bpms)):
                bn1=upper(bpms[i][1])
                bns1=bpms[i][0]
                fbetayDPP.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(betay[bn1][2])+' '+str(alfay[bn1][0])+' '+str(alfay[bn1][1])+' '+str(alfay[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
            fbetayDPP.close()

            k+=1


    if woliny2!=1 and wolinx2!=1:

        if len(ListOfNonZeroDPPX)!=len(ListOfNonZeroDPPY):

            raise ValueError("list of dppx is not equal list of dppy")


        for j in range(len(ListOfNonZeroDPPX)):
            dpop=float(ListOfNonZeroDPPX[j].DPP)
            SingleFilex=[]
            SingleFiley=[]
            SingleFilex.append(ListOfNonZeroDPPX[j])
            SingleFiley.append(ListOfNonZeroDPPY[j])
            ### coupling
            try:
                MADTwiss.Cmatrix()
            except:
                0.0

            if accel=="SPS" or "RHIC" in accel:
                plane='H'
                [phasexp,Q1,MUX,bpmsx]=GetPhases(MADTwiss,PseudoListX,plane,outputpath,bd,accel,lhcphase)
                plane='V'
                [phaseyp,Q2,MUY,bpmsy]=GetPhases(MADTwiss,PseudoListY,plane,outputpath,bd,accel,lhcphase)
                # [fwqw,bpms]=GetCoupling2(MADTwiss, PseudoListX, PseudoListY, Q1, Q2, phasexp, phaseyp, bd, accel)
                [fwqw,bpms]=GetCoupling2b(MADTwiss, PseudoListX, PseudoListY, Q1, Q2, phasexp, phaseyp, bd, accel)
            elif NBcpl==1:
                [fwqw,bpms]=GetCoupling1(MADTwiss, SingleFilex, SingleFiley, Q1, Q2)
            elif NBcpl==2:
                print phasexlist[j+1]['DPP'],dpop
                # [fwqw,bpms]=GetCoupling2(MADTwiss, SingleFilex, SingleFiley, Q1, Q2, phasexlist[j+1], phaseylist[j+1], bd, accel)
                [fwqw,bpms]=GetCoupling2b(MADTwiss, SingleFilex, SingleFiley, Q1, Q2, phasexlist[j+1], phaseylist[j+1], bd, accel)
                if acswitch=="1":
                    [fwqw,bpms]=getFreeCoupling(Q1f,Q2f,Q1,Q2,fwqw,MADTwiss,bpms)

            else:
                raise ValueError('Number of monitors for coupling analysis (option -n) should be 1 or 2.')

            fwqw['DPP']=dpop
            couplelist.append(fwqw)

            ####
    #---------------------------------------- Start getsextupoles @ Glenn Vanbavinckhove

    if not higher_order:
        print "Not analysing higher order..."
        return 0

    fsex3000=open(outputpath+'getsex3000.out','w')
    fsex3000.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    # Removed since it is not used in GetLLM (vimaier)
#     fsex1200=open(outputpath+'getsex1200.out','w')
#     fsex1200.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
#     fsex2100=open(outputpath+'getsex1200.out','w')
#     fsex2100.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    foct4000=open(outputpath+'getoct4000.out','w')
    foct4000.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    fchi3000=open(outputpath+'getchi3000.out','w')
    fchi3000.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    fchi1010=open(outputpath+'getchi1010.out','w')
    fchi1010.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    fchi4000=open(outputpath+'getchi4000.out','w')
    fchi4000.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')

    fkick=open(outputpath+'getkick.out','w')
    if acswitch=='1':
        fkickac=open(outputpath+'getkickac.out','w')

    if TBTana=="SUSSIX":
    #-> 1) f3000 line (-2,0)
    #-> 2) f1200 line  (2,0)
    #-> 3) f2100 line  (0,0)

    # global stuff

    # 1)

        htot,afactor,pfactor=Getsextupole(MADTwiss,ListOfZeroDPPX,phasexlist[0],Q1f,3,0)

        print >> fsex3000,"@","f2h_factor","%le",afactor
        print >> fsex3000,"@","p_f2h_factor","%le",pfactor

        print >> fsex3000,"*","NAME","S","AMP_20","AMP_20std","PHASE_20","PHASE_20std","f3000","f3000std","phase_f_3000","phase_f_3000std","h3000","h3000_std","phase_h_3000","phase_h_3000_std"
        print >> fsex3000,"$","%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"

        for bpm in htot:
            li=htot[bpm]
            print >>fsex3000,li[0],li[1],li[2],li[3],li[4],li[5],li[6],li[7],li[8],li[9],li[10],li[11],li[12],li[13]

        fsex3000.close()



        # --------------------------------------- end getsextupoles
        #---------------------------------------- begin getchiterms @ Glenn Vanbavinckhove
        #-> 1) chi3000
        #-> 2) chi1010
        #-> 2) chi4000

        # 1) chi3000

        fchi3000.write('* NAME    S    S1    S2    X3000    X3000i    X3000r    X3000RMS   X3000PHASE   X3000PHASERMS   X3000M    X3000Mi   X3000Mr    X3000MPHASE \n')
        fchi3000.write('$ %s   %le    %le   %le   %le   %le   %le   %le   %le %le   %le   %le   %le   %le \n')

        files=[ListOfZeroDPPX,ListOfZeroDPPY]
        name='chi3000'
        plane='H'

        [dbpms,POS,XItot,XIMODEL]=getChiTerms(MADTwiss,files,plane,name,ListOfZeroDPPX,ListOfZeroDPPY)

        for i in range(0,len(dbpms)-2):

            bn=upper(dbpms[i][1])

            fchi3000.write('"'+bn+'" '+str(POS[0][i])+' '+str(POS[1][i])+' '+str(POS[2][i])+' '+str(XItot[0][i])+' '+' '+str(XItot[1][i])+' '+str(XItot[2][i])+' '+str(XItot[3][i])+' '+str(XItot[4][i])+' '+str(XItot[5][i])+' '+str(XIMODEL[0][i])+' '+str(XIMODEL[1][i])+' '+str(XIMODEL[2][i])+' '+str(XIMODEL[3][i])+'\n')

        fchi3000.close()

    # 2) chi1010

        if  accel!='SPS':

            fchi1010.write('* NAME  S    X1010   X1010RMS   X1010PHASE   X1010PHASERMS   X1010M   X1010MPHASE \n')
            fchi1010.write('$ %s   %le  %le    %le   %le   %le   %le   %le  \n')

            files=[ListOfZeroDPPX,ListOfZeroDPPY]
            name='chi1010'
            plane='H'

            [dbpms,XItot]=getchi1010(MADTwiss,files,plane,name,bn1,ListOfZeroDPPX,ListOfZeroDPPY)

            for i in range(0,len(dbpms)-2):
                bn=upper(dbpms[i][1])
                bns=dbpms[i][0]
                fchi1010.write('"'+bn+'" '+str(bns)+' '+str(XItot[0][i])+' '+str(XItot[1][i])+' '+str(XItot[2][i])+' '+str(XItot[3][i])+' '+str('0')+' '+str('0')+' '+'\n')

            fchi1010.close()
    # 1) chi4000

        fchi4000.write('* NAME    S    S1    S2    X4000    X4000i    X4000r    X4000RMS   X4000PHASE   X4000PHASERMS   X4000M    X4000Mi   X4000Mr    X4000MPHASE \n')
        fchi4000.write('$ %s   %le    %le   %le   %le   %le   %le   %le   %le %le   %le   %le   %le   %le \n')

        #files=[ListOfZeroDPPX,ListOfZeroDPPY]
        #name='chi4000'
        #plane='H'

        #[dbpms,POS,XItot,XIMODEL]=getChiTerms(MADTwiss,files,plane,name,ListOfZeroDPPX,ListOfZeroDPPY)

        #for i in range(0,len(dbpms)-2):

    #               bn=upper(dbpms[i][1])

        #       fchi4000.write('"'+bn+'" '+str(POS[0][i])+' '+str(POS[1][i])+' '+str(POS[2][i])+' '+str(XItot[0][i])+' '+' '+str(XItot[1][i])+' '+str(XItot[2][i])+' '+str(XItot[3][i])+' '+str(XItot[4][i])+' '+str(XItot[5][i])+' '+str(XIMODEL[0][i])+' '+str(XIMODEL[1][i])+' '+str(XIMODEL[2][i])+' '+str(XIMODEL[3][i])+'\n')

        fchi4000.close()




    #---------------------------------------- end chiterms
    #-----------------------------------------begin octupole
    #->  1) f4000 (-3,0)
    #

        foct4000.write('* NAME    S    AMP_30    AMP_30RMS   PHASE_30   PHASE_30RMS   H4000   H4000I   H4000R   H4000RMS  H4000PHASE  H4000PHASERMS    H4000M    H4000MI    H4000MR    HMPHASE4000  \n')
        foct4000.write('$ %s   %le   %le   %le   %le   %le   %le   %le   %le   %le   %le   %le   %le   %le   %le   %le  \n')

        files=[ListOfZeroDPPX,ListOfZeroDPPY]
        Q=[Q1,Q2]
        plane='H'
        name='f4000'

    #       [A,h,hMODEL,dbpms]=Getoctopole(MADTwiss,plane,files,phasexlist[0],Q,name,f2100M,NAMES)

    #       for i in range(0,len(dbpms)-1):
    #
    #               bn=upper(dbpms[i][1])
    #               bns=dbpms[i][0]
    #               foct4000.write('"'+bn+'" '+str(bns)+' '+str(A[0][i])+' '+str(A[1][i])+' '+str(A[2][i])+' '+str(A[3][i])+' '+str(h[0][i])+' '+str(h[1][i])+' '+str(h[2][i])+' '+str(h[3][i])+' '+str(h[4][i])+' '+str(h[5][i])+' '+str(hMODEL[0][i])+' '+str(hMODEL[1][i])+' '+str(hMODEL[2][i])+' '+str(hMODEL[3][i])+' \n')

        foct4000.close()

    #-----------------------------------------end octupole

    #----------------------------- begin get Q,JX,delta

    files=[ListOfZeroDPPX+ListOfNonZeroDPPX,ListOfZeroDPPY+ListOfNonZeroDPPY]


    fkick.write('@ RescalingFactor_for_X %le '+str(betax_ratio)+'\n')
    fkick.write('@ RescalingFactor_for_Y %le '+str(betay_ratio)+'\n')
    fkick.write('*  DPP  QX  QXRMS  QY  QYRMS  sqrt2JX  sqrt2JXSTD  sqrt2JY  sqrt2JYSTD  2JX  2JXSTD  2JY  2JYSTD  sqrt2JXRES  sqrt2JXSTDRES  sqrt2JYRES  sqrt2JYSTDRES  2JXRES  2JXSTDRES  2JYRES  2JYSTDRES\n')
    fkick.write('$  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n')

    [invarianceJx,invarianceJy,tune,tuneRMS,dpp]=getkick(files,MADTwiss)

    for i in range(0,len(dpp)):
        fkick.write(str(dpp[i])+' '+str(tune[0][i])+' '+str(tuneRMS[0][i])+' '+str(tune[1][i])+' '+str(tuneRMS[1][i])+' '+str(invarianceJx[i][0])+' '+str(invarianceJx[i][1])+' '+str(invarianceJy[i][0])+' '+str(invarianceJy[i][1])+' '+str(invarianceJx[i][0]**2)+' '+str(2*invarianceJx[i][0]*invarianceJx[i][1])+' '+str(invarianceJy[i][0]**2)+' '+str(2*invarianceJy[i][0]*invarianceJy[i][1])+' '+str(invarianceJx[i][0]/sqrt(betax_ratio))+' '+str(invarianceJx[i][1]/sqrt(betax_ratio))+' '+str(invarianceJy[i][0]/sqrt(betay_ratio))+' '+str(invarianceJy[i][1]/sqrt(betay_ratio))+' '+str(invarianceJx[i][0]**2/betax_ratio)+' '+str(2*invarianceJx[i][0]*invarianceJx[i][1]/betax_ratio)+' '+str(invarianceJy[i][0]**2/betay_ratio)+' '+str(2*invarianceJy[i][0]*invarianceJy[i][1]/betay_ratio)+'\n')


    fkick.close()


    if acswitch=='1':
        files=[ListOfZeroDPPX+ListOfNonZeroDPPX,ListOfZeroDPPY+ListOfNonZeroDPPY]

        fkickac.write('@ RescalingFactor_for_X %le '+str(betaxf_ratio)+'\n')
        fkickac.write('@ RescalingFactor_for_Y %le '+str(betayf_ratio)+'\n')
        fkickac.write('*  DPP  QX  QXRMS  QY  QYRMS  sqrt2JX  sqrt2JXSTD  sqrt2JY  sqrt2JYSTD  2JX  2JXSTD  2JY  2JYSTD  sqrt2JXRES  sqrt2JXSTDRES  sqrt2JYRES  sqrt2JYSTDRES  2JXRES  2JXSTDRES  2JYRES  2JYSTDRES\n')
        fkickac.write('$  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n')

        [invarianceJx,invarianceJy,tune,tuneRMS,dpp]=getkickac(MADTwiss_ac,files,Q1,Q2,Q1f,Q2f,acphasex_ac2bpmac,acphasey_ac2bpmac,bd,lhcphase)

        for i in range(0,len(dpp)):
            fkickac.write(str(dpp[i])+' '+str(tune[0][i])+' '+str(tuneRMS[0][i])+' '+str(tune[1][i])+' '+str(tuneRMS[1][i])+' '+str(invarianceJx[i][0])+' '+str(invarianceJx[i][1])+' '+str(invarianceJy[i][0])+' '+str(invarianceJy[i][1])+' '+str(invarianceJx[i][0]**2)+' '+str(2*invarianceJx[i][0]*invarianceJx[i][1])+' '+str(invarianceJy[i][0]**2)+' '+str(2*invarianceJy[i][0]*invarianceJy[i][1])+' '+str(invarianceJx[i][0]/sqrt(betax_ratio))+' '+str(invarianceJx[i][1]/sqrt(betax_ratio))+' '+str(invarianceJy[i][0]/sqrt(betay_ratio))+' '+str(invarianceJy[i][1]/sqrt(betay_ratio))+' '+str(invarianceJx[i][0]**2/betax_ratio)+' '+str(2*invarianceJx[i][0]*invarianceJx[i][1]/betax_ratio)+' '+str(invarianceJy[i][0]**2/betay_ratio)+' '+str(2*invarianceJy[i][0]*invarianceJy[i][1]/betay_ratio)+'\n')

        fkickac.close()




    ####### -------------- end

if __name__=="__main__":
    sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')
    options,args=parse_args()
    main(outputpath=options.output,
         dict_file=options.dict,
         files_to_analyse=options.files,
         twiss_model_file=options.Twiss,
         accel=options.ACCEL,
         lhcphase=options.lhcphase,
         BPMU=options.BPMUNIT,
         COcut=float(options.COcut),
         NBcpl=int(options.NBcpl),
         TBTana=options.TBTana,
         higher_order=options.higher)
