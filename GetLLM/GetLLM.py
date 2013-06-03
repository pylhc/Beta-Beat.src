'''
Created on 11/09/09

@author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)

@version: 2.38b

Python script to obtain Linear Lattice functions and More -> GetLLM



Usage1 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./
Usage2 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./




Some rules for variable name:
- Dictionary is used to contain the output of function
- Variable containing 'm' is a value directly obtained from measurment data
- Variable containing 'mdl' is a value related to model

Change history:
V1.0, 11/Feb/2008 by Masa. Aiba
V1.1, 18/Feb/2008:
- Debugging, add model phase and tunes to output
- add function to obtain DY
- add chromatic parameter (phase for non zero DPP)
V1.2, 22/Feb/2008:
- test version for beta with all BPM
V1.3, 29/Feb/2008:
- beta from phases is improved, averaging beta1, 2 and 3
V1.31, 12/Mar/2008:
- debugged on alpha3
V1.4, 12/Mar/2008:
- modify output to fit latest TSF format and to meet requests from Rogelio
- fix buggs in r.m.s. beta-beat and added to the output of getbetax/y.out
V1.5 Rogelio, 13 March 2008:
- Update to option parser, include BPMdictionary to filter BPMs not in Model
V1.51, 13/Mar/2008:
- Modify output to fit latest TSF format again. Add STD to beta.
V1.6, 15/Jul/2008:
- Add the integer part of tunes - assuming that the phase advance is always
  less than 1.0.
V1.71 27/Jul/2008:
- Add GetCO. Filter in dispersion calculation to exclude bad bpms.
V1.8, 13/Aug/2008 Ref. note by A. Franchi, R. T. Garcia, G. Vanbavinckhove:
- Add GetCoupling
- "Computation of the Coupling Resonance Driving term f1001 and the coupling
  coefficient C from turn-by-turn single-BPM data", 28/May/2008
- The GetCoupling.py is initiated by Glenn V. and imported into GetLLM /
  finalized by Masa. Aiba
- Some bugs are fixed - you can run even without liny file and find the
  results for horizontal plane only.
V1.81, 1/Sep/2008:
- For an accelerator in which the beam goes opposite direction to the model
  as in LHCB2, the beam direction parameter (bd) is added to GetPhases.
- Bug in the phi13 for the last monitor and the last but one is fixed.
V1.9, 21/Oct/2008:
- Add the beta from spectrum height.
V1.91, 08/Dec/2008:
- Add option - SUSSIX or SVD for input file
V1.99, 13/Feb/09:
- Add DPX! The MEMORIAL version for Version 1.** since all the linear lattice
  parameters become available!
- Add the option for the Harmonic analysis.
- Refine coding, add several lines of comment
V2.0, 17/Feb/2009:
- Major version up to output "More", that is, other than the linear lattice
  parameters.
- Add off-momentum lattice (dbeta/beta)/(dp/p)
- Add SPS coupling with "Pseudo-double plane BPM"-it need at this moment a
  list containing
- pre-paired H-V monitor. This should be replaced by a clever algorithm to
  find good pairs automatically.
V2.01, 10/Mar/2009:
- Fix bug on SPS double plane BPM monitor, in which missing BPM could cause
  an error.
- Modify BetaFromAmplitude to output invariant J (Rogelio / finalised by MA)
V2.02, 10/Mar/2009:
- Fix bug in getcoupling, bad input for phasex and phasey
V2.10, 13/Mar/2009, Glenn Vanbavinckhove:
- Added function for finding sextupole lines (amp and phases) + chiterms amp
V2.11. 26/Mar/2009
- Fix bug in Getphase (model phase advance from the last to first monitor).
V2.12, 28/03/2009
- Following bugs fixed : avoid negative square, change option -h to -l
- Add -r option as in correct.py
- Change the way to import BPM pair file for SPS. (import -> execfile)
V2.13, 06/Apl/2009
- Fix bug in weight function to accept negative dpp
V2.14, 08/Apl/2009
- Fix bug in Normalized dispersion to treat COcut correctly.
V2.15:
- Enable coupling in RHIC as in SPS
V2.16, 28/May/2009:
- Add STDBET for the beta from amplitude.
- Add option for off momentum beta-beating to choose algorithm, that is, beta
  from phase or amp
- Add a routine to detect wrong data having two lines in linx/y file with the
  same BPM name.
- Add a routine to avoid zero division due to exactly n*pi phase advance in
  beta from phase (see the last part of GetPhases).
V2.21, 23/June/2009:
- Add STDBET Model for off momentum beta beat phase.
V2.25:
- Fixed bd flag (must be -1 for beam2)
V2.25:
- Adding VERSION variable to be always output and modified in subsequent
  versions, do not forget!!!
V2.26:
- Adding F2000 (two different methods linear and non-linear)
V2.27:
- Adding new method for IP calculation
V2.28:
- Changing the rejection of bad BPM for the coupling phase - averaging the
  phase over the sets of data first , then cut if the q1 and q2 are very
  different.
- 24/Feb/2010 the change is not yet checked.
- Hyphens in the @ field of tfs files is not allowed:  The previous label
  "RMS-beta-beat" has been  moved to "RMSbetabeat"
V2.29 5/Mar/2010
- Change the default value for COcut from 1000 to 4000 as it was too small
V2.30 13/Sept/2010:
- Updating for AC-Dipole, implementing chromatic coupling (not done for RHIC
  or SPS)
V2.31 8/November/2010:
- Update for AC-Dipole: gives free beta,phase,coupling(global factor)
V2.32 15/January/2010:
- Taking models
V2.33 7/February/2011:
- implementing coupling correction for AC-dipole.
V2.34 7/04/2011:
- Updating to deal with chromatic twiss files.
V2.35 9/06/2011:
- Functions to cancel the AC dipole effect for beta, phase, and total phase,
  based on equations, are added.
- A function to calculate beta* from the phase advance between Q1s is added.
- Phase shift by tune is compensated for the LHC experiment data.
V2.36 30/09/2011:
- Rescaling algorithm for BetaFromAmp and action (by Andy) is implemented.
- 2nd function to calculate IP parameters from Q1s is added.
- The compensation of the phase shift by tune for LHC exp data is modified.
- A function to calculate action of the AC dipole excitation is added.
V2.38 08/03/2012:
- added main() function
- using raise ValueError() instead of sys.exit() some places
13/09/2012:
- merged in patch 2.35-2.37
V2.38b 03/dec/2012, tbach:
- reformatted comments
- changed all parts of code, where the program exits inside an exception to
  display the exception, because the messages are not always helpful
- removed ";" from all over the code (still some left)
- 207 errors, 983 warning, 574 infos left from static code analysis...
 - x.xxx, vimaier  16th Apr 2013:
    deleted functions function and GetCoupling2
    Changed GetCoupling2b to GetCoupling2
    Set some TODOs
    Changed wolin* and acswitch (* := x|y|x2|y2)
    Defined variables Q1,Q2,Q1f,Q2f,MUX,MUY,muxf,muyf with standard values
      Saves a lot of try/excepts while writing files
    Reformatted a lot of code
    Introduced tfs_file for all output files in main() --> Formatted output files

'''




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
import metaclass
import os
import string
import math

from numpy import sin,cos,tan
import numpy as np

import utils.tfs_file

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
            check_if_bpm_in_model = model.indx[bpm[1].upper()]  # @UnusedVariable
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
        return phase1ave%norm

def PhaseStd(phase0,norm):  #-- phases must be in [0,1) or [0,2*pi), norm = 1 or 2*pi

    phase0   =np.array(phase0)%norm
    phase1   =(phase0+0.5*norm)%norm-0.5*norm
    phase0ave=np.mean(phase0)
    phase1ave=np.mean(phase1)
    phase0std=math.sqrt(np.mean((phase0-phase0ave)**2))
    phase1std=math.sqrt(np.mean((phase1-phase1ave)**2))
    return min(phase0std,phase1std)

def GetPhasesTotal(MADTwiss,ListOfFiles,Q,plane,bd,oa,op):

    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms, MADTwiss)
    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and oa=="LHCB1": s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    if op=="1" and oa=="LHCB2": s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    bn1=string.upper(commonbpms[0][1])
    phaseT={}
    print "Reference BPM:", bn1, "Plane:", plane
    for i in range(0,len(commonbpms)):
        #bn2=string.upper(commonbpms[i+1][1]) ?
        bn2=string.upper(commonbpms[i][1])
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
        phi12=np.array(phi12)
        # for the beam circulating reversely to the model
        if bd==-1: phi12=1.0-phi12

        #phstd12=math.sqrt(np.average(phi12*phi12)-(np.average(phi12))**2.0+2.2e-15)
        #phi12=np.average(phi12)
        phstd12=PhaseStd(phi12,1.0)
        phi12  =PhaseMean(phi12,1.0)
        phaseT[bn2]=[phi12,phstd12,phmdl12,bn1]

    return [phaseT,commonbpms]


def GetPhases(MADTwiss,ListOfFiles,Q,plane,outputpath,beam_direction,accel,lhcphase):
    #print "Hi Get", len(ListOfFiles)

    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms, MADTwiss)
    length_commonbpms = len(commonbpms)
    #print len(commonbpms)
    #sys.exit()

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if lhcphase=="1" and accel=="LHCB1": 
        s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    if lhcphase=="1" and accel=="LHCB2": 
        s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    mu=0.0
    tunem=[]
    phase={} # Dictionary for the output containing [average phase, rms error]
    for i in range(0,length_commonbpms): # To find the integer part of tune as well, the loop is up to the last monitor
        bn1=string.upper(commonbpms[i%length_commonbpms][1])
        bn2=string.upper(commonbpms[(i+1)%length_commonbpms][1])
        bn3=string.upper(commonbpms[(i+2)%length_commonbpms][1])
        
        if bn1 == bn2 :
            print >> sys.stderr, "There seem two lines with the same BPM name "+bn1+" in linx/y file."
            print >> sys.stderr, "Please check your input data....leaving GetLLM."
            sys.exit(1)
            
        if plane == 'H':
            phmdl12 = MADTwiss.MUX[MADTwiss.indx[bn2]] - MADTwiss.MUX[MADTwiss.indx[bn1]]
            phmdl13 = MADTwiss.MUX[MADTwiss.indx[bn3]] - MADTwiss.MUX[MADTwiss.indx[bn1]]
        elif plane == 'V':
            phmdl12 = MADTwiss.MUY[MADTwiss.indx[bn2]] - MADTwiss.MUY[MADTwiss.indx[bn1]]
            phmdl13 = MADTwiss.MUY[MADTwiss.indx[bn3]] - MADTwiss.MUY[MADTwiss.indx[bn1]]
            
        if i == length_commonbpms-2:
            if plane == 'H':
                madtune = MADTwiss.Q1 % 1.0
            elif plane == 'V':
                madtune = MADTwiss.Q2 % 1.0
            
            if madtune > 0.5:
                madtune -= 1.0
            
            phmdl13 = phmdl13 % 1.0
            phmdl13 = phiLastAndLastButOne(phmdl13,madtune)
        elif i == length_commonbpms-1:
            if plane == 'H':
                madtune = MADTwiss.Q1 % 1.0
            elif plane == 'V':
                madtune = MADTwiss.Q2 % 1.0
            
            if madtune>0.5:
                madtune -= 1.0
            
            phmdl12 = phmdl12 % 1.0
            phmdl13 = phmdl13 % 1.0
            phmdl12 = phiLastAndLastButOne(phmdl12,madtune)
            phmdl13 = phiLastAndLastButOne(phmdl13,madtune)




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
                if MADTwiss.S[MADTwiss.indx[bn1]]<=s_lastbpm and MADTwiss.S[MADTwiss.indx[bn2]] >s_lastbpm: 
                    phm12 += beam_direction*Q
                if MADTwiss.S[MADTwiss.indx[bn1]]<=s_lastbpm and MADTwiss.S[MADTwiss.indx[bn3]] >s_lastbpm: 
                    phm13 += beam_direction*Q
                if MADTwiss.S[MADTwiss.indx[bn1]] >s_lastbpm and MADTwiss.S[MADTwiss.indx[bn2]]<=s_lastbpm: 
                    phm12 += -beam_direction*Q
                if MADTwiss.S[MADTwiss.indx[bn1]] >s_lastbpm and MADTwiss.S[MADTwiss.indx[bn3]]<=s_lastbpm: 
                    phm13 += -beam_direction*Q
            except: pass
            if phm12<0: phm12+=1
            if phm13<0: phm13+=1
            phi12.append(phm12)
            phi13.append(phm13)

        phi12=np.array(phi12)
        phi13=np.array(phi13)
        if beam_direction==-1: # for the beam circulating reversely to the model
            phi12=1.0-phi12
            phi13=1.0-phi13

        #if any(phi12)>0.9 and i !=len(commonbpms): # Very small phase advance could result in larger than 0.9 due to measurement error
        #       print 'Warning: there seems too large phase advance! '+bn1+' to '+bn2+' = '+str(phi12)+'in plane '+plane+', recommended to check.'
        phstd12=PhaseStd(phi12,1.0)
        phstd13=PhaseStd(phi13,1.0)
        phi12=PhaseMean(phi12,1.0)
        phi13=PhaseMean(phi13,1.0)
        #phstd12=math.sqrt(np.average(phi12*phi12)-(np.average(phi12))**2.0+2.2e-15)
        #phstd13=math.sqrt(np.average(phi13*phi13)-(np.average(phi13))**2.0+2.2e-15)
        #phi12=np.average(phi12)
        #phi13=np.average(phi13)
        tunemi=np.array(tunemi)
        if i<length_commonbpms-1 :
            tunem.append(np.average(tunemi))

        # Note that the phase advance between the last monitor and the first monitor should be find by taking into account the fractional part of tune.
        if i==length_commonbpms-2:
            tunem=np.array(tunem)
            tune=np.average(tunem)
            phi13=phiLastAndLastButOne(phi13,tune)
        elif i==length_commonbpms-1:
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
            bn1=string.upper(commonbpms[i][1])
            bn2=string.upper(commonbpms[i+1][1])
            bn3=string.upper(commonbpms[0][1])
        elif i==len(commonbpms)-1: # The last monitor
            bn1=string.upper(commonbpms[i][1])
            bn2=string.upper(commonbpms[0][1])
            bn3=string.upper(commonbpms[1][1])
        else : # Others
            bn1=string.upper(commonbpms[i][1])
            bn2=string.upper(commonbpms[i+1][1])
            bn3=string.upper(commonbpms[i+2][1])
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
        bn1=string.upper(commonbpms[i][1])
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
        bn1=string.upper(commonbpms[i][1])
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

    #print Amp2/Kick2


    delbeta=[]
    for i in range(0,len(commonbpms)):
        bn1=string.upper(commonbpms[i][1])
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

#-------------------------

def GetCO(MADTwiss, ListOfFiles):

    commonbpms=intersect(ListOfFiles)
    commonbpms=modelIntersect(commonbpms, MADTwiss)
    co={} # Disctionary for output
    for i in range(0,len(commonbpms)):
        bn1=string.upper(commonbpms[i][1])
        coi=0.0
        coi2=0.0
        for j in ListOfFiles:
            coi=coi + j.CO[j.indx[bn1]]
            coi2=coi2 + j.CO[j.indx[bn1]]**2
        coi=coi/len(ListOfFiles)
        corms=math.sqrt(coi2/len(ListOfFiles)-coi**2+2.2e-16)
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
    mydp=np.array(mydp)
    wf=np.array(abs(mydp))/sum(abs(mydp))*len(mydp) #Weight for the average depending on DPP


    # Find the global factor
    nd=[]
    ndmdl=[]
    badco=0
    for i in range(0,len(commonbpmsALL)):
        bn1=string.upper(commonbpmsALL[i][1])
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
        bn1=string.upper(commonbpmsALL[i][1])
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
        bn1=string.upper(commonbpms[i][1])
        if i==len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2=string.upper(commonbpms[0][1])
            phmdl12=2.*np.pi*(MADTwiss.Q1+MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
        else:
            bn2=string.upper(commonbpms[i+1][1])
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
        bn1=string.upper(commonbpms[i][1])
        if i==len(commonbpms)-1: # The first BPM is BPM2 for the last BPM
            bn2=string.upper(commonbpms[0][1])
            phmdl12=2.*np.pi*(MADTwiss.Q2+MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
        else:
            bn2=string.upper(commonbpms[i+1][1])
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
    commonbpmsALL=intersect(ALL)



    mydp=[]
    for j in ListOfNonZeroDPP:
        mydp.append(float(j.DPP))
    mydp=np.array(mydp)
    wf=np.array(abs(mydp))/sum(abs(mydp))*len(mydp) #Weitghs for the average


    dco={} # Dictionary for the output containing [(averaged)Disp, rms error]
    bpms=[]
    for i in range(0,len(commonbpmsALL)):
        bn1=string.upper(commonbpmsALL[i][1])
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


#-----------

def GetCoupling1(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2):

    # not applicable to db=-1 for the time being...

    tp=2.0*np.pi

    # find operation point
    try:
        #TODO: There is no global outputpath. Will always crash. See Github issue #9 (vimaier)
        fdi=open(outputpath+'Drive.inp','r')  # Drive.inp file is normally in the outputpath directory in GUI operation
        for line in fdi:
            if "TUNE X" in line:
                fracxinp = line.split("=")
                fracx = fracxinp[1]
            if "TUNE Y" in line:
                fracyinp = line.split("=")
                fracy = fracyinp[1]
        fdi.close()
    except:
        fracx = Q1 # Otherwise, the fractional parts are assumed to be below 0.5
        fracy = Q2

    if fracx < 0.0 :
        fracx = 1.0 - Q1
    else:
        fracx = Q1
    if fracy < 0.0 :
        fracx = 1.0 - Q2
    else:
        fracy=Q2

    if fracx > fracy:
        sign_QxmQy = 1.0
    else:
        sign_QxmQy = -1.0

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
        bn1=string.upper(dbpms[i][1])

        fij=[]
        q1j=[]
        q2j=[]
        badbpm=0
        for j in range(0,len(ListOfZeroDPPX)):
            jx=ListOfZeroDPPX[j]
            jy=ListOfZeroDPPY[j]
            C01ij=jx.AMP01[jx.indx[bn1]]
            C10ij=jy.AMP10[jy.indx[bn1]]
            fij.append(0.5*math.atan(math.sqrt(C01ij*C10ij)))

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
        q1j=np.array(q1j)
        q2j=np.array(q2j)
        q1=np.average(q1j)
        q2=np.average(q2j)

        if abs(q1-q2)<0.25:  # Very rough cut !!!!!!!!!!!!!!!!!!!
            qi=(q1+q2)/2.0
        elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            qi=q1 # Note that q1 and q2 are confined 0. to 1.
        else:
            badbpm=1
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1, "total so far", countBadPhase



        if badbpm==0:
            fij=np.array(fij)
            fi=np.average(fij)
            fistd=math.sqrt(np.average(fij*fij)-(np.average(fij))**2.0+2.2e-16)
            #qij=np.array(qij)
            #qi=np.average(qij)
            #qistd=math.sqrt(np.average(qij*qij)-(np.average(qij))**2.0+2.2e-16)
            qistd=math.sqrt(np.average(q1j*q1j)-(np.average(q1j))**2.0+2.2e-16) # Not very exact...
            fi=fi*complex(cos(tp*qi),sin(tp*qi))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])
            # Trailing "0,0" in following lists because of compatibility. 
            # See issue on github pylhc/Beta-Beat.src#3
            # --vimaier
            fwqw[bn1]=[[fi,fistd,0,0],[qi,qistd,0,0]]


    dbpms=dbpmt


    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)):
        jx=ListOfZeroDPPX[j]
        jy=ListOfZeroDPPY[j]
        bn1=string.upper(dbpms[i][1])
        CG=CG+math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(jx.MUX[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]])


    CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
    QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]


    return [fwqw,dbpms]

#-----------

def ComplexSecondaryLine(delta, cw, cw1, pw, pw1):
    tp=2.0*np.pi
    a1=complex(1.0,-tan(tp*delta))
    a2=cw*complex(cos(tp*pw),sin(tp*pw))
    a3=-1.0/cos(tp*delta)*complex(0.0,1.0)
    a4=cw1*complex(cos(tp*pw1),sin(tp*pw1))
    SL=a1*a2+a3*a4
    sizeSL=math.sqrt(SL.real**2+SL.imag**2)
    phiSL=(np.arctan2(SL.imag , SL.real)/tp) %1.0
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
    tp=2.0*np.pi
    C=cos(delta*tp)
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
    phase=(np.arctan2(sig.imag,sig.real)/tp) %1.0

    # computing error secondary line (h-)
    esig=(sig1*complex(1,-(1/C))-sig2*complex(0,1)*(SC))*edelta

    eamp=abs(esig)
    ephase=(np.arctan2(esig.imag,esig.real)/tp) %1.0

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

    tp=2.0*np.pi
    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)-1):
        bn1=string.upper(dbpms[i][1])
        bn2=string.upper(dbpms[i+1][1])

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
            [SA0p1ij,phi0p1ij]=ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]],
                    jx.PHASE01[jx.indx[bn1]], jx.PHASE01[jx.indx[bn2]])
            [SA0m1ij,phi0m1ij]=ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]],
                    -jx.PHASE01[jx.indx[bn1]], -jx.PHASE01[jx.indx[bn2]])
            [TBp10ij,phip10ij]=ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]],
                    jy.PHASE10[jy.indx[bn1]], jy.PHASE10[jy.indx[bn2]])
            [TBm10ij,phim10ij]=ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]],
                    -jy.PHASE10[jy.indx[bn1]], -jy.PHASE10[jy.indx[bn2]])


            #print SA0p1ij,phi0p1ij,SA0m1ij,phi0m1ij,TBp10ij,phip10ij,TBm10ij,phim10ij
            f1001ij.append(0.5*math.sqrt(TBp10ij*SA0p1ij/2.0/2.0))
            f1010ij.append(0.5*math.sqrt(TBm10ij*SA0m1ij/2.0/2.0))

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

        q1jd=np.array(q1jd)
        q2jd=np.array(q2jd)
        q1d=PhaseMean(q1jd,1.0)
        q2d=PhaseMean(q2jd,1.0)

        q1js=np.array(q1js)
        q2js=np.array(q2js)
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

            f1001ij=np.array(f1001ij)
            f1001i=np.average(f1001ij)
            f1001istd=math.sqrt(np.average(f1001ij*f1001ij)-(np.average(f1001ij))**2.0+2.2e-16)
            f1010ij=np.array(f1010ij)
            f1010i=np.average(f1010ij)
            f1010istd=math.sqrt(np.average(f1010ij*f1010ij)-(np.average(f1010ij))**2.0+2.2e-16)

            q1001i=PhaseMean(np.array([q1d,q2d]),1.0)
            q1010i=PhaseMean(np.array([q1s,q2s]),1.0)
            q1001istd=PhaseStd(np.append(q1jd,q2jd),1.0)
            q1010istd=PhaseStd(np.append(q1js,q2js),1.0)

            f1001i=f1001i*complex(cos(tp*q1001i),sin(tp*q1001i))
            f1010i=f1010i*complex(cos(tp*q1010i),sin(tp*q1010i))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])

            if bd==1:
                fwqw[bn1]=[[f1001i,f1001istd,f1010i,f1010istd],[q1001i,q1001istd,q1010i,q1010istd]]
            elif bd==-1:
                fwqw[bn1]=[[f1010i,f1010istd,f1001i,f1001istd],[q1010i,q1010istd,q1001i,q1001istd]]


    dbpms=dbpmt

    # possible correction ??
    #bn0=string.upper(dbpms[0][1])
    #up1=fwqw[bn0][0][0]
    #up2=fwqw[bn0][0][2]
    #for i in range(1,len(dbpms)):
        #bn0=string.upper(dbpms[i-1][1])
        #bn1=string.upper(dbpms[i][1])
        #df1001=math.sqrt(fwqw[bn0][0][0].real**2+fwqw[bn0][0][0].imag**2)/math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        #df1010=math.sqrt(fwqw[bn0][0][2].real**2+fwqw[bn0][0][2].imag**2)/math.sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2)
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
        bn1=string.upper(dbpms[i][1])
        CG=CG+math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
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

        # bpmhp will be used for storing in this section. Not used further. Replaced by 
        # 'tentative solution' with bpmpair. See some lines below.
        # --vimaier
#     bpmhp=[]
#     for i in range(0,len(bpmh)):
#         smin=1.0e10
#         jsave=0
#         for j in range (0,len(bpmv),10):
#             sdiff=abs(bpmh[i][0]-bpmv[j][0])
#             if sdiff<smin:
#                 smin=sdiff
#                 jsave=j

#         jlower=jsave-9
#         jupper=jsave+9
#         if jupper > len(bpmv):
#             jupper=len(bpmv)
#         for j in range (jlower,jupper):
#             sdiff=abs(bpmh[i][0]-bpmv[j][0])
#             if sdiff<smin:
#                 smin=sdiff
#                 jsave=j
# 
#         bpmhp.append([bpmh[i][0],bpmh[i][1],bpmv[jsave][1],0])


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
    # dbpms is replaced through 'tentative solution' (vimaier)
#     dbpms=bpmhp

    # tentative solution
    dbpms=bpmpair() # model BPM name
    countofmissingBPMs=0
    for i in range(0,len(dbpms)):
        wname=string.upper(dbpms[i][1]) # horizontal BPM basis of the pairing (model name)
        pname=string.upper(dbpms[i][2]) # vertical pair of the horizontal as in SPSBPMpairs (model name)
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
                # dphix is used only in commented out code beneath (vimaier)
#                 dphix=MADTwiss.MUX[MADTwiss.indx[string.upper(pname)]]-MADTwiss.MUX[MADTwiss.indx[string.upper(wname)]]
                dphiy=MADTwiss.MUY[MADTwiss.indx[string.upper(pname)]]-MADTwiss.MUY[MADTwiss.indx[string.upper(wname)]]
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
                    #dphix=MADTwiss.MUX[MADTwiss.indx[string.upper(pname)]]-MADTwiss.MUX[MADTwiss.indx[string.upper(wname)]]
                    #dphiy=MADTwiss.MUY[MADTwiss.indx[string.upper(pname)]]-MADTwiss.MUY[MADTwiss.indx[string.upper(wname)]]
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
        PseudoListX.append(metaclass.twiss(filex))
        PseudoListY.append(metaclass.twiss(filey))
        
        # Delete temp files again. (vimaier)
        os.remove(filex)
        os.remove(filey)


    return [PseudoListX,PseudoListY]

#----------------------- for finding the lines of the sextupoles (@ Glenn Vanbavinckhove)
def f2h(amp,ampphase,termj,factor,term,M2M):    # converts from f-term to h-term

    # conversion to include f2h
    tp=2.0*np.pi
    H=(amp/(termj*factor))*(1-np.e**complex(0,term*tp))

    Ampi=H.imag
    Ampr=H.real
    Amp=abs(H)/M2M
    phase=(math.atan2(Ampi,Ampr))/tp


    fh=[Ampi,Ampr,Amp,phase]

    return fh


def Getsextupole(MADTwiss,amp20list,phase,tune,j,k):
    '''
    function written to calculate resonance driving terms
    '''

    # constructing complex amplitude and phase using two BPM method

    bpms=intersect(amp20list)
    bpms=modelIntersect(bpms,MADTwiss)

    # Since beta,rmsbb(return_value[0:2]) is not used, slice return value([2:4])(vimaier)
    [bpms,invariantJx] = ( BetaFromAmplitude(MADTwiss,amp20list,'H') )[2:4]
    sqrt2jx=invariantJx[0]

    Q=tune+float(str(MADTwiss.Q1).split(".")[0])

    afactor=(1-cos(2*(j-k)*np.pi*Q))#(2*sin(np.pi*(j-k)*Q))
    #print (2*sin(np.pi*(j-k)*Q)),(1-cos(6*np.pi*Q))
    #sys.exit()
    pfactor=(np.pi*(j-k)*Q)

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
            # Since eampi,ephasei(return_value[2:4]) is not used, slice return value([0:1])(vimaier)
            ampi,phasei = ( ComplexSecondaryLineExtended(delta,edelta,amp_201,amp_202,phase_201,phase_202) )[0:2]


            if ampi!=0.0:

                amp_i_list.append(ampi)
                phase_i_list.append(phasei)

                if (j==3 and k==0):
                    factor=math.sqrt(2)### factor
                    fterm=ampi/(factor*2*j*sqrt2jx**2)
                    pterm=(phasei-phase[bpm.upper()][0]+0.25)%1

                    hterm=fterm/afactor

                    hpterm=(pterm-pfactor)%1


                elif (j==2 and k==1):
                    factor=math.sqrt(2)### factor
                    fterm=ampi/(factor*2*j*sqrt2jx**2)
                    pterm=(phasei-phase[bpm][0]+0.25)%1

                    hterm=fterm/afactor

                    hpterm=(pterm-pfactor)%1

                flist.append(fterm)
                fplist.append(pterm)
                hlist.append(hterm)
                hplist.append(hpterm)


        if len(amp_i_list)!=0.0:
            al=np.mean(amp_i_list)
            alstd=np.std(amp_i_list)

            pl=np.mean(phase_i_list)
            plstd=np.mean(phasei)

            fl=np.mean(flist)
            fstd=np.std(flist)

            fpl=np.mean(fplist)
            fpstd=np.std(fplist)

            hl=np.mean(hlist)
            hstd=np.std(hlist)

            hpl=np.mean(hplist)
            hpstd=np.std(hplist)


            htot[bpm]=[bpm,s,al,alstd,pl,plstd,fl,fstd,fpl,fpstd,hl,hstd,hpl,hpstd]


    return htot,afactor,pfactor


def Getoctopole(MADTwiss,plane,twiss_files,phaseI,Q,fname,fM,NAMES):
    '''
    for finding secondary lines of the octuple (@ Glenn Vanbavinckhove)
    '''

        # intersects BPMs
    dbpms=intersect(twiss_files[0])
    dbpms=modelIntersect(dbpms,MADTwiss)



    # value definition
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
    for j in range(0,len(twiss_files[0])):
        singleFilex=[twiss_files[0][j]]
        singleFiley=[twiss_files[1][j]]

        # Since beta,rmsbb,bpms(return_value[0:3]) is not used, slice return value([3])(vimaier)
        invariantJx = ( BetaFromAmplitude(MADTwiss,singleFilex,'H') )[3] 

        # Since beta,rmsbb,bpms(return_value[0:3]) is not used, slice return value([3])(vimaier)
        invariantJy = ( BetaFromAmplitude(MADTwiss,singleFiley,'V') )[3]

        invarianceJx.append(invariantJx)
        invarianceJy.append(invariantJy)


    # for the model
    for i in range(0,len(dbpms)):

        bpm=string.upper(dbpms[i][1])

        bpmC=MADTwiss.NAME[MADTwiss.indx[bpm]]


        for j in range(0,len(NAMES)):
            try:
                name=NAMES[j]

                if name==bpmC:

                    amp=abs(fM[j])
                    ampr=fM[i].real
                    ampi=fM[j].imag
                    phase=np.arctan2(ampi,ampr)%1

                    hMODELT.append(amp)
                    hMODELTr.append(ampr)
                    hMODELTi.append(ampi)
                    h_phase_MODELT.append(phase)



            except:
                print 'name '+str(NAMES[j])+' is not found in dictionary'
            hMODEL=[hMODELT,hMODELTi,hMODELTr,h_phase_MODELT]

    #calculation of f,q,h,qh
    for i in range(0,len(dbpms)-1):

        bn1=string.upper(dbpms[i][1])
        bn2=string.upper(dbpms[i+1][1])

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

        for j in range(0,len(twiss_files[0])):

            single_twiss = twiss_files[0][j]

            # for f4000
            if fname=='f4000':

                [A,phi]=ComplexSecondaryLine(dell, single_twiss.AMP_30[single_twiss.indx[bn1]], single_twiss.AMP_30[single_twiss.indx[bn2]], single_twiss.PHASE_30[single_twiss.indx[bn1]], single_twiss.PHASE_30[single_twiss.indx[bn2]])

                factor=float(8*invarianceJx[j][0]**1.5)   # 1 to fit with model
                term=float(4*Q[0])
                termj=4
                M2M=0.5

            #------ converting
            h=f2h(A,phi,termj,factor,term,M2M)

            #----- adding the terms
            AS.append(A)
            phaseS.append(phi)
            hSi.append(h[0])
            hSr.append(h[1])
            hS.append(h[2])
            h_phaseS.append(h[3])

        # array and taking average for all the input files for one BPM
        AS=np.array(AS)
        A_SRMS=math.sqrt(np.average(AS*AS)-(np.average(AS))**2+2.2e-16)

        phaseS=np.array(phaseS)
        try:
            phase_RMSS=math.sqrt(np.average(phaseS*phaseS)-(np.average(phaseS))**2+2.2e-16)
        except:
            phase_RMSS=0

        hS=np.array(hS)
        hSi=np.array(hSi)
        hSr=np.array(hSr)
        try:
            h_RMSS=math.sqrt(np.average(hS*hS)-(np.average(hS))**2+2.2e-16)
        except:
            h_RMSS=0

        h_phaseS=np.array(h_phaseS)
        try:
            phase_rms=np.average(h_phaseS*h_phaseS)-(np.average(h_phaseS))**2+2.2e-16
        except:
            phase_rms=0
        h_phase_RMSS=math.sqrt(phase_rms)

        # real output
        AT.append(np.average(AS))
        A_RMST.append(A_SRMS)

        phaseT.append(np.average(phaseS))
        phase_RMST.append(phase_RMSS)

        hT.append(np.average(hS))
        hTi.append(np.average(hSi))
        hTr.append(np.average(hSr))
        h_RMST.append(h_RMSS)

        h_phaseT.append(np.average(h_phaseS))
        h_phase_RMST.append(h_phase_RMSS)

        A=[AT,A_RMST,phaseT,phase_RMST]
        h=[hT,hTi,hTr,h_RMST,h_phaseT,h_phase_RMST]




    return [A,h,hMODEL,dbpms]


def computeChiTerms(amp,phase_20,phase,terms,J,plane,ima,rea):
    ''' for finding the chi terms '''

    #computes the chiterms for different inputs
    twoPi=2*np.pi

    delta1=((phase[1]-phase[0]-0.25)*twoPi)
    delta2=((phase[2]-phase[1]-0.25)*twoPi)

    inp=0.13 # ????
    #term1=((amp[0]*np.e**complex(0,phase_20[0]*twoPi)))/cos(delta1)
    #term2=((amp[1]*np.e**complex(0,phase_20[1]*twoPi)))*(tan(delta1)+tan(delta2))
    #term3=((amp[2]*np.e**complex(0,phase_20[2]*twoPi)))/cos(delta2)
    term1=((amp[0]*np.e**complex(0,(phase_20[0]+inp)*twoPi)))/cos(delta1)
    term2=((amp[1]*np.e**complex(0,(phase_20[1]+inp)*twoPi)))*(tan(delta1)+tan(delta2))
    term3=((amp[2]*np.e**complex(0,(phase_20[2]+inp)*twoPi)))/cos(delta2)
    chiTOT=(term1+term2+term3)

    chiAMP=abs(chiTOT)

    chiAMPi=chiTOT.imag
    chiAMPr=chiTOT.real
    #print chiTOT.imag,chiTOT.real
    chiPHASE=(((np.arctan2(chiTOT.imag,chiTOT.real)))/twoPi)%1
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
    files = filesF[0]

    dbpms = intersect(files)
    dbpms = modelIntersect(dbpms, MADTwiss)


    # initiliasing variables
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

        bn1=string.upper(dbpms[i][1])

        BPMS.append(bn1)

    #### invariance
    for j in range(0,len(ListOfZeroDPPX)):
        # Since betax,rmsbbx,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJX = ( BetaFromAmplitude(MADTwiss,ListOfZeroDPPX,'H') )[3]
        # Since betay,rmsbby,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJY= ( BetaFromAmplitude(MADTwiss,ListOfZeroDPPY,'V') )[3]
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

            phase=np.arctan2(MODEL[i].imag,MODEL[i].real)%1


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

        bn1=string.upper(dbpms[i][1])
        bn2=string.upper(dbpms[i+1][1])
        bn3=string.upper(dbpms[i+2][1])

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



        XI=np.array(XI)
        XIi=np.array(XIi)
        XIr=np.array(XIr)
        try:
            XIrms=math.sqrt(np.average(XI*XI)-np.average(XI)**2+2.2e-16)
        except:
            XIrms=0
        XI_phase=np.array(XI_phase)
        try:
            XI_phaseRMS=math.sqrt(np.average(XI_phase*XI_phase)-np.average(XI_phase)**2+2.2e-16)
        except:
            XI_phaseRMS=0


        XIT.append(np.average(XI))
        XITi.append(np.average(XIi))
        XITr.append(np.average(XIr))
        XIrmsT.append(XIrms)
        XI_phase_T.append(np.average(XI_phase))
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
    XIT=[]
    XIrmsT=[]
    XI_phase_T=[]
    XI_phaseRMS_T=[]

    BPMS=[]
    invarianceJx=[]
    invarianceJy=[]

    for i in range(0,len(dbpms)): # ask rogelio


        BPMS.append(bn1)

    #### invariance
    for j in range(0,len(files)):
        # Since betax,rmsbbx,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJX = ( BetaFromAmplitude(MADTwiss,ListOfZeroDPPX,'H') )[3]
        # Since betay,rmsbby,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJY = ( BetaFromAmplitude(MADTwiss,ListOfZeroDPPY,'V') )[3]
        invarianceJx.append(invariantJX[0])
        invarianceJy.append(invariantJY[0])


    for i in range(0,len(dbpms)):

        XI=[]
        XIrms=[]
        XI_phase=[]
        XI_phaseRMS=[]

        bn=string.upper(dbpms[i][1])
        bny=string.upper(dbpmsy[i][1])

        for j in range(0,len(files)):

            jx=files[j]

            jy=filesy[j]

            amp10x=jx.AMP01[jx.indx[bn]]
            amp10y=jy.AMP10[jy.indx[bny]]
            phase10x=jx.PHASE01[jx.indx[bn]]
            phasex=jx.MUX[jx.indx[bn]]

            XI1010=0.25*math.sqrt(amp10x*amp10y)
            phase1010=phase10x+phasex

            XI.append(XI1010)
            XI_phase.append(phase1010)

        XI=np.array(XI)
        XIrms=math.sqrt(np.average(XI*XI)-np.average(XI)**2+2.2e-16)
        XI_phase=np.array(XI_phase)
        XI_phaseRMS=math.sqrt(np.average(XI_phase*XI_phase)-np.average(XI_phase)**2+2.2e-16)


        XIT.append(np.average(XI))
        XIrmsT.append(XIrms)
        XI_phase_T.append(np.average(XI_phase))
        XI_phaseRMS_T.append(XI_phaseRMS)

        XItot=[XIT,XIrmsT,XI_phase_T,XI_phaseRMS_T]

    return [dbpms,XItot]

#---- construct offmomentum
def ConstructOffMomentumModel(MADTwiss,dpp, dictionary):

    j=MADTwiss
    bpms=intersect([MADTwiss])

    Qx=j.Q1+dpp*j.DQ1
    Qy=j.Q2+dpp*j.DQ2

    ftemp_name = "./TempTwiss.dat"
    ftemp=open(ftemp_name,"w")
    ftemp.write("@ Q1 %le "+str(Qx)+"\n")
    ftemp.write("@ Q2 %le "+str(Qy)+"\n")
    ftemp.write("@ DPP %le "+str(dpp)+"\n")
    ftemp.write("* NAME S BETX BETY ALFX ALFY MUX MUY\n")
    ftemp.write("$ %s %le %le  %le  %le  %le  %le %le\n")


    for i in range(0,len(bpms)):
        bn=string.upper(bpms[i][1])
        bns=bpms[i][0]

        # dbeta and dalpha will be extract via metaclass. As it is for the time being.
        ax=j.WX[j.indx[bn]]*cos(2.0*np.pi*j.PHIX[j.indx[bn]])
        bx=j.WX[j.indx[bn]]*sin(2.0*np.pi*j.PHIX[j.indx[bn]])
        bx1=bx+j.ALFX[j.indx[bn]]*ax
        NBETX=j.BETX[j.indx[bn]]*(1.0+ax*dpp)
        NALFX=j.ALFX[j.indx[bn]]+bx1*dpp
        NMUX=j.MUX[j.indx[bn]]+j.DMUX[j.indx[bn]]*dpp

        ay=j.WY[j.indx[bn]]*cos(2.0*np.pi*j.PHIY[j.indx[bn]])
        by=j.WY[j.indx[bn]]*sin(2.0*np.pi*j.PHIY[j.indx[bn]])
        by1=by+j.ALFY[j.indx[bn]]*ay
        NBETY=j.BETY[j.indx[bn]]*(1.0+ay*dpp)
        NALFY=j.ALFY[j.indx[bn]]+by1*dpp
        NMUY=j.MUY[j.indx[bn]]+j.DMUY[j.indx[bn]]*dpp

        ftemp.write('"'+bn+'" '+str(bns)+" "+str(NBETX)+" "+str(NBETY)+" "+str(NALFX)+" "+str(NALFY)+" "+str(NMUX)+" "+str(NMUY)+"\n")

    ftemp.close()
    dpptwiss=metaclass.twiss(ftemp_name,dictionary)
    
    # Delete temp file again(vimaier)
    os.remove(ftemp_name)
    

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

        # Since beta,rmsbb,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJx = ( BetaFromAmplitude(MADTwiss,[x],'H') )[3]
        # Since beta,rmsbb,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJy = ( BetaFromAmplitude(MADTwiss,[y],'V') )[3]

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

    bpml="null"
    bpmh="null"
    for ind in range(len(indxes)):
        #print ind
        name=model.NAME[ind]
        if "BPMSW.1L"+IP in name:
            bpml=name
            try:
                test = measured[0][bpml][0]  # @UnusedVariable
            except:
                bpml="null"
        if "BPMSW.1R"+IP in name:
            bpmh=name
            try:
                test = measured[0][bpmh][0]  # @UnusedVariable
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
        betastar=(2*math.sqrt(betaxl)*math.sqrt(betaxr)*sin(deltaphimodel*2*np.pi))/(betayl+betayr-2*math.sqrt(betaxl)*math.sqrt(betaxr)*cos(2*np.pi*deltaphimodel))*L
        location=((betaxl-betaxr)/(betaxl+betaxr-2*math.sqrt(betaxl)*math.sqrt(betaxr)*cos(2*np.pi*deltaphimodel)))*L

        deltaphi=(math.atan((L-location)/betastar)+math.atan((L+location)/betastar))/(2*np.pi)

        betahor=[IP,betastar,location,deltaphi,betaipx,deltaphimodel,0]


        print "horizontal betastar for ",IP," is ",str(betastar)," at location ",str(location), " of IP center with phase advance ",str(deltaphi)

        #vertical
        deltaphimodel=abs(model.MUY[model.indx[BPMright]]-model.MUY[model.indx[BPMleft]])


        betastar=(2*math.sqrt(betayl)*math.sqrt(betayr)*sin(deltaphimodel*2*np.pi))/(betayl+betayr-2*math.sqrt(betayl)*math.sqrt(betayr)*cos(2*np.pi*deltaphimodel))*L
        location=((betayl-betayr)/(betayl+betayr-2*math.sqrt(betayl)*math.sqrt(betayr)*cos(2*np.pi*deltaphimodel)))*L

        deltaphi=(math.atan((L-location)/betastar)+math.atan((L+location)/betastar))/(2*np.pi)

        betaver=[IP,betastar,location,deltaphi,betaipy,deltaphimodel,0]

        print "vertical betastar for ",IP," is ",str(betastar)," at location ",str(location), " of IP center with phase advance ",str(deltaphi)


    return [betahor,betaver]

def GetIP2(MADTwiss,Files,Q,plane,bd,oa,op):

    #-- Common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],string.upper(b[1])) for b in bpm]

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
                            dpsi=2*np.pi*bd*(Files[i].MUX[Files[i].indx[bpmr]]-Files[i].MUX[Files[i].indx[bpml]])
                        else:
                            dpsi=2*np.pi*(Q+bd*(Files[i].MUX[Files[i].indx[bpmr]]-Files[i].MUX[Files[i].indx[bpml]]))
                        #-- To compensate the phase shift by tune
                        if op=='1':
                            if (bd==1 and ip=='2') or (bd==-1 and ip=='8'): dpsi+=2*np.pi*Q
                    if plane=='V':
                        al=Files[i].AMPY[Files[i].indx[bpml]]
                        ar=Files[i].AMPY[Files[i].indx[bpmr]]
                        if list(zip(*bpm)[1]).index(bpmr)>list(zip(*bpm)[1]).index(bpml):
                            dpsi=2*np.pi*bd*(Files[i].MUY[Files[i].indx[bpmr]]-Files[i].MUY[Files[i].indx[bpml]])
                        else:
                            dpsi=2*np.pi*(Q+bd*(Files[i].MUY[Files[i].indx[bpmr]]-Files[i].MUY[Files[i].indx[bpml]]))
                        #-- To compensate the phase shift by tune
                        if op=='1':
                            if (bd==1 and ip=='2') or (bd==-1 and ip=='8'): dpsi+=2*np.pi*Q

                    #-- bet, alf, and math.sqrt(2J) from amp and phase advance
                    bet =L*(al**2+ar**2+2*al*ar*cos(dpsi))/(2*al*ar*sin(dpsi))
                    alf =(al**2-ar**2)/(2*al*ar*sin(dpsi))
                    bets=bet/(1+alf**2)
                    ds  =alf*bets
                    rt2J=math.sqrt(al*ar*sin(dpsi)/(2*L))
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

def getCandGammaQmin(fqwq,bpms,tunex,tuney,twiss):
    # Cut the fractional part of Q1 and Q2
    QQ1 = float( int(twiss.Q1) )
    QQ2 = float( int(twiss.Q2) )

    tunex=float(tunex)+QQ1
    tuney=float(tuney)+QQ2

    tunefactor=(cos(2*np.pi*tunex)-cos(2*np.pi*tuney))/(np.pi*(sin(2*np.pi*tunex)+sin(2*np.pi*tuney)))

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
            gamma=math.sqrt(1/(1/(1+4*(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2))))
            #print detC
            ffactor= 2*gamma*tunefactor*math.sqrt(abs(detC)) # cannot take abs
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

    Qmin=np.array(Qmin)

    Qminerr=math.sqrt(np.average(Qmin*Qmin)-(np.average(Qmin))**2+2.2e-16)
    Qminav=np.average(Qmin)




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
    factor_top_diff=math.sqrt(sin(np.pi*(tunedrivenx-tunefreey))*sin(np.pi*(tunefreex-tunedriveny)))
    factor_bottom_diff=sin(np.pi*(tunefreex-tunefreey))

    factor_diff=abs((factor_top_diff/factor_bottom_diff))

    print "Factor for coupling diff ",factor_diff

    # sum f1010
    factor_top_sum=math.sqrt(sin(np.pi*(tunedrivenx+tunefreey))*sin(np.pi*(tunefreex+tunedriveny)))
    factor_bottom_sum=sin(np.pi*(tunefreex+tunefreey))

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
    '''
    :Parameters:
        'phase': dict
            (bpm_name:string) --> (phase_list:[phi12,phstd12,phi13,phstd13,phmdl12,phmdl13,bn2])
            phi13, phstd13, phmdl12 and phmdl13 are note used.
        
    '''

    #print "Calculating free phase using model"

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

def getfreephaseTotal(phase,bpms,plane,MADTwiss,MADTwiss_ac):
    '''
    :Parameters:
        'phase': dict
            (bpm_name:string) --> (phase_list:[phi12,phstd12,phmdl12,bn1])
            phmdl12 and bn1 are note used.
        
    '''
    #print "Calculating free total phase using model"

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

    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))
    psid_ac2bpmac1=np.arctan((1+r)/(1-r)*tan(2*np.pi*psi_ac2bpmac1-np.pi*Q))%np.pi-np.pi+np.pi*Qd
    psid_ac2bpmac2=np.arctan((1+r)/(1-r)*tan(2*np.pi*psi_ac2bpmac2+np.pi*Q))%np.pi-np.pi*Qd

    return {bpmac1:psid_ac2bpmac1,bpmac2:psid_ac2bpmac2}


def GetFreePhaseTotal_Eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,op):

    #-- Select common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],string.upper(b[1])) for b in bpm]

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
    if plane=='H': psimdl=np.array([(MADTwiss.MUX[MADTwiss.indx[b[1]]]-MADTwiss.MUX[MADTwiss.indx[bpm[0][1]]])%1 for b in bpm])
    if plane=='V': psimdl=np.array([(MADTwiss.MUY[MADTwiss.indx[b[1]]]-MADTwiss.MUY[MADTwiss.indx[bpm[0][1]]])%1 for b in bpm])

    #-- Global parameters of the driven motion
    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    psiall=np.zeros((len(bpm),len(Files)))
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
        for k in range(len(bpm)): psiall[k][i]=psi[k]/(2*np.pi)  #-- phase range back to [0,1)

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
    bpm=[(b[0],string.upper(b[1])) for b in bpm]

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
    if plane=='H': psimdl=np.array([MADTwiss.MUX[MADTwiss.indx[b[1]]] for b in bpm])
    if plane=='V': psimdl=np.array([MADTwiss.MUY[MADTwiss.indx[b[1]]] for b in bpm])
    psi12mdl=(np.append(psimdl[1:],psimdl[0] +Q)-psimdl)%1
    psi13mdl=(np.append(psimdl[2:],psimdl[:2]+Q)-psimdl)%1

    #-- Global parameters of the driven motion
    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))

    #-- Loop for files, psid, Psi, Psid are w.r.t the AC dipole
    psi12all=np.zeros((len(bpm),len(Files)))
    psi13all=np.zeros((len(bpm),len(Files)))
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
        psi12=(np.append(psi[1:],psi[0] +2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
        #psi12=(np.append(psi[1:],psi[0])-psi)/(2*np.pi)  #-- phase range back to [0,1)
        psi13=(np.append(psi[2:],psi[:2]+2*np.pi*Q)-psi)/(2*np.pi)  #-- phase range back to [0,1)
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
    bpm = modelIntersect(intersect(Files),MADTwiss_ac)
    bpm = [(b[0],string.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    if op=="1" and bd==1: 
        s_lastbpm=MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L2.B1']]
    if op=="1" and bd==-1: 
        s_lastbpm=MADTwiss_ac.S[MADTwiss_ac.indx['BPMSW.1L8.B2']]

    #-- Determine the BPM closest to the AC dipole and its position
    for b in psid_ac2bpmac.keys():
        if '5L4' in b: 
            bpmac1=b
        if '6L4' in b: 
            bpmac2=b
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
    if plane=='H': betmdl=np.array([MADTwiss_ac.BETX[MADTwiss_ac.indx[b[1]]] for b in bpm])
    if plane=='V': betmdl=np.array([MADTwiss_ac.BETY[MADTwiss_ac.indx[b[1]]] for b in bpm])

    #-- Global parameters of the driven motion
    r=sin(np.pi*(Qd-Q))/sin(np.pi*(Qd+Q))

    #-- Loop for files
    betall=np.zeros((len(bpm),len(Files)))
    Adall=np.zeros((len(bpm),len(Files)))
    for i in range(len(Files)):
        if plane=='H':
            amp =np.array([2*Files[i].AMPX[Files[i].indx[b[1]]] for b in bpm])
            psid=bd*2*np.pi*np.array([Files[i].MUX[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        if plane=='V':
            amp =np.array([2*Files[i].AMPY[Files[i].indx[b[1]]] for b in bpm])
            psid=bd*2*np.pi*np.array([Files[i].MUY[Files[i].indx[b[1]]] for b in bpm])  #-- bd flips B2 phase to B1 direction
        for k in range(len(bpm)):
            try:
                if bpm[k][0]>s_lastbpm: psid[k]+=2*np.pi*Qd  #-- To fix the phase shift by Q
            except: pass
        Ad  =amp/map(math.sqrt,betmdl)
        psid=psid-(psid[k_bpmac]-psid_ac2bpmac[bpmac])
        Psid=psid+np.pi*Qd
        Psid[k_bpmac:]=Psid[k_bpmac:]-2*np.pi*Qd
        bet =(amp/np.mean(Ad))**2*(1+r**2+2*r*np.cos(2*Psid))/(1-r**2)
        for k in range(len(bpm)):
            betall[k][i]=bet[k]
            Adall[k][i]=Ad[k]

    #-- Output
    result={}
    bb=[]
    Adave=[]
    for k in range(len(bpm)):
        betave=np.mean(betall[k])
        betstd=math.sqrt(np.mean((betall[k]-betave)**2))
        bb.append((betave-betmdl[k])/betmdl[k])
        Adave.append(np.mean(Adall[k]))
        result[bpm[k][1]]=[betave,betstd,bpm[k][0]]
    bb=math.sqrt(np.mean(np.array(bb)**2))
    Ad=[np.mean(Adave),math.sqrt(np.mean((Adave-np.mean(Adave))**2))]

    return [result,bb,bpm,Ad]


def GetFreeCoupling_Eq(MADTwiss,FilesX,FilesY,Qh,Qv,Qx,Qy,psih_ac2bpmac,psiv_ac2bpmac,bd):

    #-- Details of this algorithms is in http://www.agsrhichome.bnl.gov/AP/ap_notes/ap_note_410.pdf

    #-- Check linx/liny files, may be redundant
    if len(FilesX)!=len(FilesY): return [{},[]]

    #-- Select common BPMs
    bpm=modelIntersect(intersect(FilesX+FilesY),MADTwiss)
    bpm=[(b[0],string.upper(b[1])) for b in bpm]

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
        f1001xArgAve=PhaseMean(f1001xArg[k],2*np.pi)
        f1001yArgAve=PhaseMean(f1001yArg[k],2*np.pi)
        f1010xArgAve=PhaseMean(f1010xArg[k],2*np.pi)
        f1010yArgAve=PhaseMean(f1010yArg[k],2*np.pi)
        if min(abs(f1001xArgAve-f1001yArgAve),2*np.pi-abs(f1001xArgAve-f1001yArgAve))>np.pi/2: badbpm=1
        if min(abs(f1010xArgAve-f1010yArgAve),2*np.pi-abs(f1010xArgAve-f1010yArgAve))>np.pi/2: badbpm=1

        #-- Output
        if badbpm==0:
            f1001AbsAve=np.mean(f1001Abs[k])
            f1010AbsAve=np.mean(f1010Abs[k])
            f1001ArgAve=PhaseMean(np.append(f1001xArg[k],f1001yArg[k]),2*np.pi)
            f1010ArgAve=PhaseMean(np.append(f1010xArg[k],f1010yArg[k]),2*np.pi)
            f1001Ave   =f1001AbsAve*np.exp(1j*f1001ArgAve)
            f1010Ave   =f1010AbsAve*np.exp(1j*f1010ArgAve)
            f1001AbsStd=math.sqrt(np.mean((f1001Abs[k]-f1001AbsAve)**2))
            f1010AbsStd=math.sqrt(np.mean((f1010Abs[k]-f1010AbsAve)**2))
            f1001ArgStd=PhaseStd(np.append(f1001xArg[k],f1001yArg[k]),2*np.pi)
            f1010ArgStd=PhaseStd(np.append(f1010xArg[k],f1010yArg[k]),2*np.pi)
            fwqw[bpm[k][1]]=[[f1001Ave          ,f1001AbsStd       ,f1010Ave          ,f1010AbsStd       ],
                             [f1001ArgAve/(2*np.pi),f1001ArgStd/(2*np.pi),f1010ArgAve/(2*np.pi),f1010ArgStd/(2*np.pi)]]  #-- Phases renormalized to [0,1)
            goodbpm.append(bpm[k])

    #-- Global parameters not implemented yet
    fwqw['Global']=['"null"','"null"']

    return [fwqw,goodbpm]

def GetFreeIP2_Eq(MADTwiss,Files,Qd,Q,psid_ac2bpmac,plane,bd,oa,op):

    #-- Common BPMs
    bpm=modelIntersect(intersect(Files),MADTwiss)
    bpm=[(b[0],string.upper(b[1])) for b in bpm]

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
        # Since beta,rmsbb,bpms(return_value[:3]) are not used, slice the return value([3]) (vimaier)
        invariantJx = ( GetFreeBetaFromAmp_Eq(MADTwiss_ac,[x],Qh,Qx,psih_ac2bpmac,'H',bd,op) )[3]
        # Since beta,rmsbb,bpms(return_value[:3]) are not used, slice the return value([3]) (vimaier)
        invariantJy = ( GetFreeBetaFromAmp_Eq(MADTwiss_ac,[y],Qv,Qy,psiv_ac2bpmac,'V',bd,op) )[3]
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
        fout.write('@ MAD_FILE %s "'+mad_file+'"'+'\n')
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

    outputpath = _fix_output(outputpath)

    if dict_file == "0":
        BPMdictionary = {}
    else:
        execfile(dict_file)
        BPMdictionary = dictionary   # temporaryly since presently name is not BPMdictionary

    listOfInputFiles = files_to_analyse.split(",")

    # Beam direction
    beam_direction = 1
    if accel == "LHCB2":
        beam_direction = -1 # THIS IS CORRECT, be careful with tune sign in SUSSIX and eigenmode order in SVD
    elif accel == "LHCB4":
        beam_direction = 1  # IS THIS CORRECT? I (rogelio) try for Simon...
        accel = "LHCB2" #This is because elements later are named B2 anyway, not B4

    #-- finding base model
    try:
        MADTwiss = metaclass.twiss(twiss_model_file,BPMdictionary) # MODEL from MAD
        print "Base model found!"
    except:
        print >> sys.stderr, "twiss file loading failed for:",twiss_model_file
        print >> sys.stderr, traceback.format_exc()
        sys.exit(1)

    #-- finding the ac dipole model
    with_ac_calc = False
    try:
        MADTwiss_ac = metaclass.twiss(twiss_model_file.replace(".dat","_ac.dat"))
        with_ac_calc = True
        print "Driven Twiss file found. AC dipole effects calculated with the effective model (get***_free2.out)"
    except:
        MADTwiss_ac = MADTwiss
        print "WARN: AC dipole effects not calculated. Driven twiss file does not exsist !"

    #-- Test if the AC dipole (MKQA) is in the model of LHC
    if with_ac_calc:
        if 'LHC' in accel:
            if 'MKQA.6L4.'+accel[3:] in MADTwiss.NAME:
                print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
            else:
                try:
                    MADTwissElem = metaclass.twiss(twiss_model_file.replace(".dat","_elements.dat"))
                    print "AC dipole found in the model. AC dipole effects calculated with analytic equations (get***_free.out)"
                except:
                    print 'WARN: AC dipoles not in the model. AC dipole effects not calculated with analytic equations !'
        else: print 'WARN: AC dipole effects calculated with analytic equations only for LHC for now'

    if BPMU == 'um': 
        COcut = COcut
    elif BPMU == 'mm': 
        COcut = COcut/1.0e3
    elif BPMU == 'cm': 
        COcut = COcut/1.0e4
    elif BPMU == 'm': 
        COcut = COcut/1.0e6




    if TBTana == "SUSSIX":
        Suffix_x = '_linx'
        Suffix_y = '_liny'
    elif TBTana == 'SVD':
        Suffix_x = '_svdx'
        Suffix_y = '_svdy'
    elif TBTana == 'HA':
        Suffix_x = '_hax'
        Suffix_y = '_hay'

    # Static variable of TfsFile to save the outputfile    
    utils.tfs_file.TfsFile.s_output_path = outputpath
    
    files_dict = {}
    files_dict['getphasex.out'] = utils.tfs_file.TfsFile('getphasex.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getphasey.out'] = utils.tfs_file.TfsFile('getphasey.out').add_getllm_header(VERSION, twiss_model_file)
    
    files_dict['getphasetotx.out'] = utils.tfs_file.TfsFile('getphasetotx.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getphasetoty.out'] = utils.tfs_file.TfsFile('getphasetoty.out').add_getllm_header(VERSION, twiss_model_file)

    if with_ac_calc:
        files_dict['getphasex_free.out'] = utils.tfs_file.TfsFile('getphasex_free.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getphasey_free.out'] = utils.tfs_file.TfsFile('getphasey_free.out').add_getllm_header(VERSION, twiss_model_file)
        
        files_dict['getphasex_free2.out'] = utils.tfs_file.TfsFile('getphasex_free2.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getphasey_free2.out'] = utils.tfs_file.TfsFile('getphasey_free2.out').add_getllm_header(VERSION, twiss_model_file)
       
        files_dict['getphasetotx_free.out'] = utils.tfs_file.TfsFile('getphasetotx_free.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getphasetoty_free.out'] = utils.tfs_file.TfsFile('getphasetoty_free.out').add_getllm_header(VERSION, twiss_model_file)
        
        files_dict['getphasetotx_free2.out'] = utils.tfs_file.TfsFile('getphasetotx_free2.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getphasetoty_free2.out'] = utils.tfs_file.TfsFile('getphasetoty_free2.out').add_getllm_header(VERSION, twiss_model_file)
        
    files_dict['getbetax.out'] = utils.tfs_file.TfsFile('getbetax.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getbetay.out'] = utils.tfs_file.TfsFile('getbetay.out').add_getllm_header(VERSION, twiss_model_file)

    if with_ac_calc:
        files_dict['getbetax_free.out'] = utils.tfs_file.TfsFile('getbetax_free.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getbetay_free.out'] = utils.tfs_file.TfsFile('getbetay_free.out').add_getllm_header(VERSION, twiss_model_file)
        
        files_dict['getbetax_free2.out'] = utils.tfs_file.TfsFile('getbetax_free2.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getbetay_free2.out'] = utils.tfs_file.TfsFile('getbetay_free2.out').add_getllm_header(VERSION, twiss_model_file)

    files_dict['getampbetax.out'] = utils.tfs_file.TfsFile('getampbetax.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getampbetay.out'] = utils.tfs_file.TfsFile('getampbetay.out').add_getllm_header(VERSION, twiss_model_file)

    if with_ac_calc:
        files_dict['getampbetax_free.out'] = utils.tfs_file.TfsFile('getampbetax_free.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getampbetay_free.out'] = utils.tfs_file.TfsFile('getampbetay_free.out').add_getllm_header(VERSION, twiss_model_file)

        files_dict['getampbetax_free2.out'] = utils.tfs_file.TfsFile('getampbetax_free2.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getampbetay_free2.out'] = utils.tfs_file.TfsFile('getampbetay_free2.out').add_getllm_header(VERSION, twiss_model_file)

    files_dict['getCOx.out'] = utils.tfs_file.TfsFile('getCOx.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getCOy.out'] = utils.tfs_file.TfsFile('getCOy.out').add_getllm_header(VERSION, twiss_model_file)

    files_dict['getNDx.out'] = utils.tfs_file.TfsFile('getNDx.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getDx.out'] = utils.tfs_file.TfsFile('getDx.out').add_getllm_header(VERSION, twiss_model_file)
    files_dict['getDy.out'] = utils.tfs_file.TfsFile('getDy.out').add_getllm_header(VERSION, twiss_model_file)

    files_dict['getcouple.out'] = utils.tfs_file.TfsFile('getcouple.out').add_getllm_header(VERSION, twiss_model_file)
    if with_ac_calc:
        files_dict['getcouple_free.out'] = utils.tfs_file.TfsFile('getcouple_free.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getcouple_free2.out'] = utils.tfs_file.TfsFile('getcouple_free2.out').add_getllm_header(VERSION, twiss_model_file)

  

    files_dict['getcoupleterms.out'] = utils.tfs_file.TfsFile('getcoupleterms.out').add_getllm_header(VERSION, twiss_model_file)

    if "LHC" in accel:
        files_dict['getIP.out'] = utils.tfs_file.TfsFile('getIP.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getIPx.out'] = utils.tfs_file.TfsFile('getIPx.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getIPy.out'] = utils.tfs_file.TfsFile('getIPy.out').add_getllm_header(VERSION, twiss_model_file)
        files_dict['getIPfromphase.out'] = utils.tfs_file.TfsFile('getIPfromphase.out').add_getllm_header(VERSION, twiss_model_file)
        if with_ac_calc:
            files_dict['getIPx_free.out'] = utils.tfs_file.TfsFile('getIPx_free.out').add_getllm_header(VERSION, twiss_model_file)
            files_dict['getIPy_free.out'] = utils.tfs_file.TfsFile('getIPy_free.out').add_getllm_header(VERSION, twiss_model_file)
            files_dict['getIPx_free2.out'] = utils.tfs_file.TfsFile('getIPx_free2.out').add_getllm_header(VERSION, twiss_model_file)
            files_dict['getIPy_free2.out'] = utils.tfs_file.TfsFile('getIPy_free2.out').add_getllm_header(VERSION, twiss_model_file)
            files_dict['getIPfromphase_free.out'] = utils.tfs_file.TfsFile('getIPfromphase_free.out').add_getllm_header(VERSION, twiss_model_file)
            files_dict['getIPfromphase_free2.out'] = utils.tfs_file.TfsFile('getIPfromphase_free2.out').add_getllm_header(VERSION, twiss_model_file)

    FileOfNonZeroDPPX = []
    FileOfNonZeroDPPY = []
    
    ListOfZeroDPPX = []
    ListOfNonZeroDPPX = []
    ListOfZeroDPPY = []
    ListOfNonZeroDPPY = []

    for file_in in listOfInputFiles:
        if file_in.endswith(".gz"):
            file_x = file_in.replace(".gz", Suffix_x+".gz")
        else:
            file_x = file_in+Suffix_x
        twiss_file_x = metaclass.twiss(file_x)
        try:
            dppi = twiss_file_x.DPP
        except:
            dppi = 0.0

        if type(dppi) != float:
            print 'Warning: DPP may not be given as a number in ',file_x ,'...trying to forcibly cast it as a number'
            try:
                dppi = float(dppi)
                print 'dppi= ',dppi
            except:
                print >> sys.stderr, 'but failing. DPP in ',file_x , ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)

        if dppi == 0.0:
            ListOfZeroDPPX.append(twiss_file_x)
            files_dict['getphasex.out'].add_filename_to_getllm_header(file_x)
            files_dict['getphasetotx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getbetax.out'].add_filename_to_getllm_header(file_x)
            files_dict['getampbetax.out'].add_filename_to_getllm_header(file_x)
            files_dict['getCOx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getNDx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getDx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getcouple.out'].add_filename_to_getllm_header(file_in)
            if "LHC" in accel:
                files_dict['getIPx.out'].add_filename_to_getllm_header(file_in)
                files_dict['getIPy.out'].add_filename_to_getllm_header(file_in)
                files_dict['getIPfromphase.out'].add_filename_to_getllm_header(file_in)
                if with_ac_calc:
                    files_dict['getIPx_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPy_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPx_free2.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPy_free2.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPfromphase_free.out'].add_filename_to_getllm_header(file_in)
                    files_dict['getIPfromphase_free2.out'].add_filename_to_getllm_header(file_in)
            if with_ac_calc:
                files_dict['getphasex_free.out'].add_filename_to_getllm_header(file_x)
                files_dict['getphasex_free2.out'].add_filename_to_getllm_header(file_x)
                files_dict['getphasetotx_free.out'].add_filename_to_getllm_header(file_x)
                files_dict['getphasetotx_free2.out'].add_filename_to_getllm_header(file_x)
                files_dict['getbetax_free.out'].add_filename_to_getllm_header(file_x)
                files_dict['getbetax_free2.out'].add_filename_to_getllm_header(file_x)
                files_dict['getampbetax_free.out'].add_filename_to_getllm_header(file_x)
                files_dict['getampbetax_free2.out'].add_filename_to_getllm_header(file_x)
                files_dict['getcouple_free.out'].add_filename_to_getllm_header(file_in)
                files_dict['getcouple_free2.out'].add_filename_to_getllm_header(file_in)
        else:
            ListOfNonZeroDPPX.append(twiss_file_x)
            FileOfNonZeroDPPX.append(file_x)
            files_dict['getNDx.out'].add_filename_to_getllm_header(file_x)
            files_dict['getDx.out'].add_filename_to_getllm_header(file_x)

        try:
            if file_in.endswith(".gz"):
                file_y = file_in.replace(".gz", Suffix_y+".gz")
            else:
                file_y = file_in+Suffix_y
            twiss_file_y = metaclass.twiss(file_y)
            try:
                dppi = twiss_file_y.DPP
            except:
                dppi = 0.0

            if type(dppi)!=float:
                print 'Warning: DPP may not be given as a number in ',file_y ,'...trying to forcibly cast it as a number'
                try:
                    dppi=float(dppi)
                    print 'dppi= ',dppi
                except:
                    print 'but failing. DPP in ',file_y , ' is something wrong. String?'
                    raise ValueError('leaving GetLLM')

            if dppi == 0.0:
                ListOfZeroDPPY.append(twiss_file_y)
                files_dict['getphasey.out'].add_filename_to_getllm_header(file_y)
                files_dict['getphasetoty.out'].add_filename_to_getllm_header(file_y)
                files_dict['getbetay.out'].add_filename_to_getllm_header(file_y)
                files_dict['getampbetay.out'].add_filename_to_getllm_header(file_y)
                files_dict['getCOy.out'].add_filename_to_getllm_header(file_y)
                files_dict['getDy.out'].add_filename_to_getllm_header(file_y)
                if with_ac_calc:
                    files_dict['getphasey_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasey_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasetoty_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getphasetoty_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getbetay_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getbetay_free2.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getampbetay_free.out'].add_filename_to_getllm_header(file_y)
                    files_dict['getampbetay_free2.out'].add_filename_to_getllm_header(file_y)
            else:
                ListOfNonZeroDPPY.append(twiss_file_y)
                FileOfNonZeroDPPY.append(file_y)
                files_dict['getDy.out'].add_filename_to_getllm_header(file_y)
        except:
            print 'Warning: There seems no '+str(file_y)+' file in the specified directory.'


    with_liny_for_zero_dppy = True #FLAG meaning there is _liny file for zero DPPY!
    with_liny_for_nonzero_dppy = True #FLAG meaning there is _liny file for non-zero DPPY!
    with_linx_for_zero_dppx = True
    with_linx_for_nonzero_dppx = True


    if len(ListOfZeroDPPY) == 0 :
        with_liny_for_zero_dppy = False  
    if len(ListOfNonZeroDPPY) == 0 :
        with_liny_for_nonzero_dppy = False  
    if len(ListOfZeroDPPX) == 0 :
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'
        with_linx_for_zero_dppx = False
    if len(ListOfNonZeroDPPX) == 0 :
        with_linx_for_nonzero_dppx = False

    if (len(ListOfNonZeroDPPX) != 0) and (len(ListOfZeroDPPX) == 0):
        ListOfZeroDPPX = ListOfNonZeroDPPX
        ListOfZeroDPPY = ListOfNonZeroDPPY
        with_linx_for_zero_dppx = True
        with_liny_for_zero_dppy = True
        with_liny_for_nonzero_dppy = False
        with_linx_for_nonzero_dppx = False
        print "Previous warning suppressed, running in chromatic mode"
        files_dict['getphasex.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getbetax.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getampbetax.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getCOx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getNDx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getDx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getcouple.out'].add_filename_to_getllm_header("chrommode")
        if with_ac_calc:
            files_dict['getcouple_free.out'].add_filename_to_getllm_header("chrommode")
            files_dict['getcouple_free2.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getphasey.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getbetay.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getampbetay.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getCOx.out'].add_filename_to_getllm_header("chrommode")
        files_dict['getDy.out'].add_filename_to_getllm_header("chrommode")

    Q1 = 0.0
    Q2 = 0.0
    MUX = 0.0
    MUY = 0.0
    
    Q1f = 0.0
    Q2f = 0.0
    muxf = 0.0
    muyf = 0.0
    
    if with_ac_calc:
        # Get fractional part: frac(62.23) = 0.23; 62.23 % 1 ==> 0.23 (vimaier)
        Q1f = abs(MADTwiss.Q1) % 1 #-- Free Q1 (tempolarlly, overwritten later)
        Q2f = abs(MADTwiss.Q2) % 1 #-- Free Q2 (tempolarlly, overwritten later)
        Q1 = abs(MADTwiss_ac.Q1) % 1 #-- Drive Q1 (tempolarlly, overwritten later)
        Q2 = abs(MADTwiss_ac.Q2) % 1 #-- Drive Q2 (tempolarlly, overwritten later)
        d1 =Q1-Q1f #-- Used later to calculate free Q1
        d2 =Q2-Q2f #-- Used later to calculate free Q2
    else:
        Q1f=ListOfZeroDPPX[0].Q1
        Q2f=ListOfZeroDPPY[0].Q2

    # Construct pseudo-double plane BPMs
    if (accel=="SPS" or "RHIC" in accel) and with_linx_for_zero_dppx and with_liny_for_zero_dppy:
        execfile(twiss_model_file.replace("twiss.dat","BPMpair.py"))
        [PseudoListX,PseudoListY] = PseudoDoublePlaneMonitors(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, BPMdictionary)


    #-------- Check monitor compatibility between data and model
    all_twiss_files = ListOfNonZeroDPPX+ListOfZeroDPPX+ListOfNonZeroDPPY+ListOfZeroDPPY
    for twiss_file in all_twiss_files:
        for bpm_name in twiss_file.NAME:
            #TODO: maybe easier with the usage of intersect?(vimaier)
            # Check if all BPMs are in the model(vimaier)
            try:
                MADTwiss.NAME[MADTwiss.indx[bpm_name]]
            except:
                try:
                    MADTwiss.NAME[MADTwiss.indx[string.upper(bpm_name)]]
                except:
                    print 'Monitor '+bpm_name+' cannot be found in the model!'
                    #exit()


    #-------- START Phase
    print 'Calculating phase'
    
    #---- Calling GetPhases first to save tunes
    #TODO:redundant code here? Better: if x: doX(); if y: doY();? Check it (vimaier)
    if with_linx_for_zero_dppx and with_liny_for_zero_dppy:
        #-- Calculate temp value of tune
        Q1temp = []
        Q2temp = []
        for twiss_f in ListOfZeroDPPX: 
            Q1temp.append(np.mean(twiss_f.TUNEX))
        for twiss_f in ListOfZeroDPPY: 
            Q2temp.append(np.mean(twiss_f.TUNEY))
        Q1temp = np.mean(Q1temp)
        Q2temp = np.mean(Q2temp)
        if len(ListOfZeroDPPX[0].NAME)==0:
            print "No BPMs in linx file"
            sys.exit(1)
        if  len(ListOfZeroDPPY[0].NAME)==0:
            print "No BPMs in liny file"
            sys.exit(1)
        [phasex,Q1,MUX,bpmsx] = GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1temp,'H',outputpath,beam_direction,accel,lhcphase)
        [phasey,Q2,MUY,bpmsy] = GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2temp,'V',outputpath,beam_direction,accel,lhcphase)
        #TODO: what is KK??(vimaier)
        print "KK end"
    elif with_linx_for_zero_dppx:
        #-- Calculate temp value of tune
        Q1temp = []
        for i in ListOfZeroDPPX: 
            Q1temp.append(np.mean(i.TUNEX))
        Q1temp = np.mean(Q1temp)

        [phasex,Q1,MUX,bpmsx] = GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1temp,'H',outputpath,beam_direction,accel,lhcphase)
        print 'liny missing and output x only ...'
    elif with_liny_for_zero_dppy:
        #-- Calculate temp value of tune
        Q2temp = []
        for i in ListOfZeroDPPY: 
            Q2temp.append(np.mean(i.TUNEY))
        Q2temp=np.mean(Q2temp)

        [phasey,Q2,MUY,bpmsy] = GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2temp,'V',outputpath,beam_direction,accel,lhcphase)
        print 'linx missing and output y only ...'

    #---- Re-run GetPhase to fix the phase shift by Q for exp data of LHC
    if lhcphase=="1":
        if with_linx_for_zero_dppx and with_liny_for_zero_dppy:
            [phasex,Q1,MUX,bpmsx] = GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',outputpath,beam_direction,accel,lhcphase)
            [phasey,Q2,MUY,bpmsy] = GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',outputpath,beam_direction,accel,lhcphase)
        elif with_linx_for_zero_dppx:
            [phasex,Q1,MUX,bpmsx] = GetPhases(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',outputpath,beam_direction,accel,lhcphase)
        elif with_liny_for_zero_dppy:
            [phasey,Q2,MUY,bpmsy] = GetPhases(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',outputpath,beam_direction,accel,lhcphase)


    #---- ac to free phase from eq and the model
    if with_ac_calc:
        if with_linx_for_zero_dppx:
            Q1f = Q1-d1  #-- Free H-tune
            try:    
                acphasex_ac2bpmac = GetACPhase_AC2BPMAC(MADTwissElem,Q1,Q1f,'H',accel)
            except: 
                acphasex_ac2bpmac = GetACPhase_AC2BPMAC(MADTwiss,Q1,Q1f,'H',accel)
                
            [phasexf,muxf,bpmsxf] = GetFreePhase_Eq(MADTwiss,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',beam_direction,lhcphase)
            [phasexf2,muxf2,bpmsxf2] = getfreephase(phasex,Q1,Q1f,bpmsx,MADTwiss_ac,MADTwiss,"H")
        if with_liny_for_zero_dppy:
            Q2f=Q2-d2  #-- Free V-tune
            try:    
                acphasey_ac2bpmac=GetACPhase_AC2BPMAC(MADTwissElem,Q2,Q2f,'V',accel)
            except: 
                acphasey_ac2bpmac=GetACPhase_AC2BPMAC(MADTwiss,Q2,Q2f,'V',accel)
            [phaseyf,muyf,bpmsyf]=GetFreePhase_Eq(MADTwiss,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',beam_direction,lhcphase)
            [phaseyf2,muyf2,bpmsyf2]=getfreephase(phasey,Q2,Q2f,bpmsy,MADTwiss_ac,MADTwiss,"V")

    #---- H plane result
    if with_linx_for_zero_dppx:

        phasexlist=[]
        phasex['DPP']=0.0
        phasexlist.append(phasex)
        tfs_file = files_dict['getphasex.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEX","STDPHX","PHXMDL","MUXMDL"])
        tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(len(bpmsx)):
            bn1=string.upper(bpmsx[i][1])
            bns1=bpmsx[i][0]
            phmdl=phasex[bn1][4]
            if i==len(bpmsx)-1:
                bn2=string.upper(bpmsx[0][1])
                bns2=bpmsx[0][0]
            else:
                bn2=string.upper(bpmsx[i+1][1])
                bns2=bpmsx[i+1][0]
            list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPX),phasex[bn1][0],phasex[bn1][1],phmdl,MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free phase
        if with_ac_calc:

            #-- from eq
            try:
                tfs_file = files_dict['getphasex_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEX","STDPHX","PHXMDL","MUXMDL"])
                tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(len(bpmsxf)):
                    bn1=string.upper(bpmsxf[i][1])
                    bns1=bpmsxf[i][0]
                    phmdlf=phasexf[bn1][4]
                    if i==len(bpmsxf)-1:
                        bn2=string.upper(bpmsxf[0][1])
                        bns2=bpmsxf[0][0]
                    else:
                        bn2=string.upper(bpmsxf[i+1][1])
                        bns2=bpmsxf[i+1][0]
                    list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPX),phasexf[bn1][0],phasexf[bn1][1],phmdlf,MADTwiss.MUX[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            tfs_file = files_dict['getphasex_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEX","STDPHX","PHXMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsxf2)):
                bn1=string.upper(bpmsxf2[i][1])
                bns1=bpmsxf2[i][0]
                phmdlf2=phasexf2[bn1][2]
                bn2=phasexf2[bn1][3]
                bns2=phasexf2[bn1][4]
                list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPX),phasexf2[bn1][0],phasexf2[bn1][1],phmdlf2,MADTwiss.MUX[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)

    #---- V plane result
    if with_liny_for_zero_dppy:

        phaseylist=[]
        phasey['DPP']=0.0
        phaseylist.append(phasey)
        tfs_file = files_dict['getphasey.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEY","STDPHY","PHYMDL","MUYMDL"])
        tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(len(bpmsy)):
            bn1=string.upper(bpmsy[i][1])
            bns1=bpmsy[i][0]
            phmdl=phasey[bn1][4]
            # TODO easier with modulo(vimaier)
            if i==len(bpmsy)-1:
                bn2=string.upper(bpmsy[0][1])
                bns2=bpmsy[0][0]
            else:
                bn2=string.upper(bpmsy[i+1][1])
                bns2=bpmsy[i+1][0]
            list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPY),phasey[bn1][0],phasey[bn1][1],phmdl,MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free phase
        if with_ac_calc:

            #-- from eq
            try:
                tfs_file = files_dict['getphasey_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEY","STDPHY","PHYMDL","MUYMDL"])
                tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(len(bpmsyf)):
                    bn1=string.upper(bpmsyf[i][1])
                    bns1=bpmsyf[i][0]
                    phmdlf=phaseyf[bn1][4]
                    if i==len(bpmsyf)-1:
                        bn2=string.upper(bpmsyf[0][1])
                        bns2=bpmsyf[0][0]
                    else:
                        bn2=string.upper(bpmsyf[i+1][1])
                        bns2=bpmsyf[i+1][0]
                    list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPY),phaseyf[bn1][0],phaseyf[bn1][1],phmdlf,MADTwiss.MUY[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            tfs_file = files_dict['getphasey_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEY","STDPHY","PHYMDL","MUYMDL"])
            tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsyf2)):
                bn1=string.upper(bpmsyf2[i][1])
                bns1=bpmsyf2[i][0]
                phmdlf2=phaseyf2[bn1][2]
                bn2=phaseyf2[bn1][3]
                bns2=phaseyf2[bn1][4]
                list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPY),phaseyf2[bn1][0],phaseyf2[bn1][1],phmdlf2,MADTwiss.MUY[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)


    #-------- START Total Phase
    print 'Calculating total phase'

    #---- H plane result
    if with_linx_for_zero_dppx:

        [phasexT,bpmsxT]=GetPhasesTotal(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',beam_direction,accel,lhcphase)
        tfs_file = files_dict['getphasetotx.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEX","STDPHX","PHXMDL","MUXMDL"])
        tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpmsxT)):
            bn1=string.upper(bpmsxT[i][1])
            bns1=bpmsxT[i][0]
            phmdl=phasexT[bn1][2]
            bn2=string.upper(bpmsxT[0][1])
            bns2=bpmsxT[0][0]
            list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPX),phasexT[bn1][0],phasexT[bn1][1],phmdl,MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free total phase
        if with_ac_calc:

            #-- from eq
            try:
                [phasexTf,bpmsxTf]=GetFreePhaseTotal_Eq(MADTwiss,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',beam_direction,lhcphase)
                tfs_file = files_dict['getphasetotx_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEX","STDPHX","PHXMDL","MUXMDL"])
                tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(0,len(bpmsxTf)):
                    bn1=string.upper(bpmsxTf[i][1])
                    bns1=bpmsxTf[i][0]
                    phmdlf=phasexTf[bn1][2]
                    bn2=string.upper(bpmsxTf[0][1])
                    bns2=bpmsxTf[0][0]
                    list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPX),phasexTf[bn1][0],phasexTf[bn1][1],phmdlf,MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                    tfs_file.add_table_row(list_row_entries)
            except: pass

            #-- from the model
            [phasexTf2,bpmsxTf2]=getfreephaseTotal(phasexT,bpmsxT,"H",MADTwiss,MADTwiss_ac)
            tfs_file = files_dict['getphasetotx_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEX","STDPHX","PHXMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsxTf2)):
                bn1=string.upper(bpmsxTf2[i][1])
                bns1=bpmsxTf2[i][0]
                phmdlf2=phasexTf2[bn1][2]
                bn2=string.upper(bpmsxTf2[0][1])
                bns2=bpmsxTf2[0][0]
                list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPX),phasexTf2[bn1][0],phasexTf2[bn1][1],phmdlf2,MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)


    #---- V plane result
    if with_liny_for_zero_dppy:

        [phaseyT,bpmsyT]=GetPhasesTotal(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',beam_direction,accel,lhcphase)
        tfs_file = files_dict['getphasetoty.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("MUX", "%le", str(MUX))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("MUY", "%le", str(MUY))
        tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEY","STDPHY","PHYMDL","MUYMDL"])
        tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpmsyT)):
            bn1=string.upper(bpmsyT[i][1])
            bns1=bpmsyT[i][0]
            phmdl=phaseyT[bn1][2]
            bn2=string.upper(bpmsyT[0][1])
            bns2=bpmsyT[0][0]
            list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPY),phaseyT[bn1][0],phaseyT[bn1][1],phmdl,MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free total phase
        if with_ac_calc:

            #-- from eq
            try:
                [phaseyTf,bpmsyTf]=GetFreePhaseTotal_Eq(MADTwiss,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',beam_direction,lhcphase)
                tfs_file = files_dict['getphasetoty_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("MUX", "%le", str(muxf))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("MUY", "%le", str(muyf))
                tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEY","STDPHY","PHYMDL","MUYMDL"])
                tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(0,len(bpmsyTf)):
                    bn1=string.upper(bpmsyTf[i][1])
                    bns1=bpmsyTf[i][0]
                    phmdlf=phaseyTf[bn1][2]
                    bn2=string.upper(bpmsyTf[0][1])
                    bns2=bpmsyTf[0][0]
                    list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPY),phaseyTf[bn1][0],phaseyTf[bn1][1],phmdlf,MADTwiss.MUY[MADTwiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            [phaseyTf2,bpmsyTf2]=getfreephaseTotal(phaseyT,bpmsyT,"V",MADTwiss,MADTwiss_ac)
            tfs_file = files_dict['getphasetoty_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("MUX", "%le", str(muxf2))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("MUY", "%le", str(muyf2))
            tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASEY","STDPHY","PHYMDL","MUYMDL"])
            tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsyTf2)):
                bn1=string.upper(bpmsyTf2[i][1])
                bns1=bpmsyTf2[i][0]
                phmdlf2=phaseyTf2[bn1][2]
                bn2=string.upper(bpmsyTf2[0][1])
                bns2=bpmsyTf2[0][0]
                list_row_entries = ['"'+bn1+'"','"'+bn2+'"',bns1,bns2,len(ListOfZeroDPPY),phaseyTf2[bn1][0],phaseyTf2[bn1][1],phmdlf2,MADTwiss.MUY[MADTwiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)


    #-------- START Beta
    print 'Calculating beta'

    #---- H plane
    if with_linx_for_zero_dppx:

        [betax,rmsbbx,alfax,bpms]=BetaFromPhase(MADTwiss_ac,ListOfZeroDPPX,phasex,'H')
        betax['DPP']=0
        betaxPhaseCopy=betax  #-- For Andy's BetaFromAmp re-scaling
        tfs_file = files_dict['getbetax.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbx))
        tfs_file.add_column_names(["NAME","S","COUNT","BETX","ERRBETX","STDBETX","ALFX","ERRALFX","STDALFX","BETXMDL","ALFXMDL","MUXMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpms)):
            bn1=string.upper(bpms[i][1])
            bns1=bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),betax[bn1][0],betax[bn1][1],betax[bn1][2],alfax[bn1][0],alfax[bn1][1],alfax[bn1][2],MADTwiss_ac.BETX[MADTwiss_ac.indx[bn1]],MADTwiss_ac.ALFX[MADTwiss_ac.indx[bn1]],MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free beta
        if with_ac_calc:

            #-- from eq
            try:
                [betaxf,rmsbbxf,alfaxf,bpmsf]=BetaFromPhase(MADTwiss,ListOfZeroDPPX,phasexf,'H')
                betaxPhaseCopyf=betaxf  #-- For Andy's BetaFromAmp re-scaling
                tfs_file = files_dict['getbetax_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf))
                tfs_file.add_column_names(["NAME","S","COUNT","BETX","ERRBETX","STDBETX","ALFX","ERRALFX","STDALFX","BETXMDL","ALFXMDL","MUXMDL"])
                tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(0,len(bpmsf)):
                    bn1=string.upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),betaxf[bn1][0],betaxf[bn1][1],betaxf[bn1][2],alfaxf[bn1][0],alfaxf[bn1][1],alfaxf[bn1][2],MADTwiss.BETX[MADTwiss.indx[bn1]],MADTwiss.ALFX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            [betaxf2,rmsbbxf2,alfaxf2,bpmsf2]=getFreeBeta(MADTwiss_ac,MADTwiss,betax,rmsbbx,alfax,bpms,'H')
            tfs_file = files_dict['getbetax_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf2))
            tfs_file.add_column_names(["NAME","S","COUNT","BETX","ERRBETX","STDBETX","ALFX","ERRALFX","STDALFX","BETXMDL","ALFXMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsf2)):
                bn1=string.upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),betaxf2[bn1][0],betaxf2[bn1][1],betaxf2[bn1][2],alfaxf2[bn1][0],alfaxf2[bn1][1],alfaxf2[bn1][2],MADTwiss.BETX[MADTwiss.indx[bn1]],MADTwiss.ALFX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)

    #---- V plane
    if with_liny_for_zero_dppy:

        [betay,rmsbby,alfay,bpms]=BetaFromPhase(MADTwiss_ac,ListOfZeroDPPY,phasey,'V')
        betay['DPP']=0
        betayPhaseCopy=betay  #-- For Andy's BetaFromAmp re-scaling
        tfs_file = files_dict['getbetay.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbby))
        tfs_file.add_column_names(["NAME","S","COUNT","BETY","ERRBETY","STDBETY","ALFY","ERRALFY","STDALFY","BETYMDL","ALFYMDL","MUYMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpms)):
            bn1=string.upper(bpms[i][1])
            bns1=bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),betay[bn1][0],betay[bn1][1],betay[bn1][2],alfay[bn1][0],alfay[bn1][1],alfay[bn1][2],MADTwiss_ac.BETY[MADTwiss_ac.indx[bn1]],MADTwiss_ac.ALFY[MADTwiss_ac.indx[bn1]],MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free beta
        if with_ac_calc:

            #-- from eq
            try:
                [betayf,rmsbbyf,alfayf,bpmsf]=BetaFromPhase(MADTwiss,ListOfZeroDPPY,phaseyf,'V')
                betayPhaseCopyf=betayf  #-- For Andy's BetaFromAmp re-scaling
                tfs_file = files_dict['getbetay_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf))
                tfs_file.add_column_names(["NAME","S","COUNT","BETY","ERRBETY","STDBETY","ALFY","ERRALFY","STDALFY","BETYMDL","ALFYMDL","MUYMDL"])
                tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(0,len(bpmsf)):
                    bn1=string.upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),betayf[bn1][0],betayf[bn1][1],betayf[bn1][2],alfayf[bn1][0],alfayf[bn1][1],alfayf[bn1][2],MADTwiss.BETY[MADTwiss.indx[bn1]],MADTwiss.ALFY[MADTwiss.indx[bn1]],MADTwiss.MUY[MADTwiss.indx[bn1]] ]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            [betayf2,rmsbbyf2,alfayf2,bpmsf2]=getFreeBeta(MADTwiss_ac,MADTwiss,betay,rmsbby,alfay,bpms,'V')
            tfs_file = files_dict['getbetay_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf2))
            tfs_file.add_column_names(["NAME","S","COUNT","BETY","ERRBETY","STDBETY","ALFY","ERRALFY","STDALFY","BETYMDL","ALFYMDL","MUYMDL"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsf2)):
                bn1=string.upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),betayf2[bn1][0],betayf2[bn1][1],betayf2[bn1][2],alfayf2[bn1][0],alfayf2[bn1][1],alfayf2[bn1][2],MADTwiss.BETY[MADTwiss.indx[bn1]],MADTwiss.ALFY[MADTwiss.indx[bn1]],MADTwiss.MUY[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)

    #------- Start beta from amplitude
    print 'Calculating beta from amplitude'
    betaxalist=[]
    betayalist=[]

    #---- H plane
    if with_linx_for_zero_dppx:

        [betax,rmsbbx,bpms,invJx]=BetaFromAmplitude(MADTwiss_ac,ListOfZeroDPPX,'H')
        betax['DPP']=0
        beta2_save=betax
        betaxalist.append(betax)

        #-- Rescaling
        betax_ratio=0
        skipped_bpmx=[]
        arcbpms=filterbpm(bpms)
        for bpm in arcbpms:
            name=string.upper(bpm[1]) # second entry is the name
            #Skip BPM with strange data
            if abs(betaxPhaseCopy[name][0]/betax[name][0])>100: skipped_bpmx.append(name)
            elif (betax[name][0]<0 or betaxPhaseCopy[name][0]<0): skipped_bpmx.append(name)
            else: betax_ratio=betax_ratio+(betaxPhaseCopy[name][0]/betax[name][0])
        try: betax_ratio=betax_ratio/(len(arcbpms)-len(skipped_bpmx))
        except: betax_ratio=1
        betax_rescale={}
        for bpm in map(string.upper,zip(*bpms)[1]): 
            betax_rescale[bpm]=[betax_ratio*betax[bpm][0],betax_ratio*betax[bpm][1],betax[bpm][2]]

        tfs_file = files_dict['getampbetax.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbx))
        tfs_file.add_descriptor("RescalingFactor", "%le", str(betax_ratio))
        tfs_file.add_column_names(["NAME","S","COUNT","BETX","BETXSTD","BETXMDL","MUXMDL","BETXRES","BETXSTDRES"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpms)):
            bn1=string.upper(bpms[i][1])
            bns1=bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),betax[bn1][0],betax[bn1][1],MADTwiss_ac.BETX[MADTwiss_ac.indx[bn1]],MADTwiss_ac.MUX[MADTwiss_ac.indx[bn1]],betax_rescale[bn1][0],betax_rescale[bn1][1] ]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free amp beta
        if with_ac_calc:

            #-- from eq
            try:
                
                # Since invJxf(return_value[3]) is not used, slice the return value([:3]) (vimaier)
                [betaxf,rmsbbxf,bpmsf] = ( GetFreeBetaFromAmp_Eq(MADTwiss_ac,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',beam_direction,lhcphase) )[:3]

                #-- Rescaling
                betaxf_ratio=0
                skipped_bpmxf=[]
                arcbpms=filterbpm(bpmsf)
                for bpm in arcbpms:
                    name=string.upper(bpm[1]) # second entry is the name
                    #Skip BPM with strange data
                    if abs(betaxPhaseCopyf[name][0]/betaxf[name][0])>10: skipped_bpmxf.append(name)
                    elif abs(betaxPhaseCopyf[name][0]/betaxf[name][0])<0.1: skipped_bpmxf.append(name)
                    elif (betaxf[name][0]<0 or betaxPhaseCopyf[name][0]<0): skipped_bpmxf.append(name)
                    else: betaxf_ratio=betaxf_ratio+(betaxPhaseCopyf[name][0]/betaxf[name][0])
                try: betaxf_ratio=betaxf_ratio/(len(arcbpms)-len(skipped_bpmxf))
                except: betaxf_ratio=1

                tfs_file = files_dict['getampbetax_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf))
                tfs_file.add_descriptor("RescalingFactor", "%le", str(betaxf_ratio))
                tfs_file.add_column_names(["NAME","S","COUNT","BETX","BETXSTD","BETXMDL","MUXMDL","BETXRES","BETXSTDRES"])
                tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(0,len(bpmsf)):
                    bn1=string.upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),betaxf[bn1][0],betaxf[bn1][1],MADTwiss.BETX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]],betaxf_ratio*betaxf[bn1][0],betaxf_ratio*betaxf[bn1][1] ]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            # Since invJxf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betaxf2,rmsbbxf2,bpmsf2] = (getFreeAmpBeta(betax,rmsbbx,bpms,invJx,MADTwiss_ac,MADTwiss,'H'))[:3]
            betaxf2_rescale = getFreeAmpBeta(betax_rescale,rmsbbx,bpms,invJx,MADTwiss_ac,MADTwiss,'H')[0]
            tfs_file = files_dict['getampbetax_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf2))
            tfs_file.add_descriptor("RescalingFactor", "%le", str(betax_ratio))
            tfs_file.add_column_names(["NAME","S","COUNT","BETX","BETXSTD","BETXMDL","MUXMDL","BETXRES","BETXSTDRES"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsf2)):
                bn1=string.upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),betaxf2[bn1][0],betaxf2[bn1][1],MADTwiss.BETX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]],betaxf2_rescale[bn1][0],betaxf2_rescale[bn1][1] ]
                tfs_file.add_table_row(list_row_entries)


    #---- V plane
    if with_liny_for_zero_dppy:

        [betay,rmsbby,bpms,invJy]=BetaFromAmplitude(MADTwiss_ac,ListOfZeroDPPY,'V')
        betay['DPP']=0
        betayalist.append(betay)

        #-- Rescaling
        betay_ratio=0
        skipped_bpmy=[]
        arcbpms=filterbpm(bpms)
        for bpm in arcbpms:
            name=string.upper(bpm[1]) # second entry is the name
            #Skip BPM with strange data
            if abs(betayPhaseCopy[name][0]/betay[name][0])>100: skipped_bpmy.append(name)
            elif (betay[name][0]<0 or betayPhaseCopy[name][0]<0): skipped_bpmy.append(name)
            else: betay_ratio=betay_ratio+(betayPhaseCopy[name][0]/betay[name][0])
        try: betay_ratio=betay_ratio/(len(arcbpms)-len(skipped_bpmy))
        except: betay_ratio=1
        betay_rescale={}
        for bpm in map(string.upper,zip(*bpms)[1]): betay_rescale[bpm]=[betay_ratio*betay[bpm][0],betay_ratio*betay[bpm][1],betay[bpm][2]]
        
        tfs_file = files_dict['getampbetay.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbby))
        tfs_file.add_descriptor("RescalingFactor", "%le", str(betay_ratio))
        tfs_file.add_column_names(["NAME","S","COUNT","BETY","BETYSTD","BETYMDL","MUYMDL","BETYRES","BETYSTDRES"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpms)):
            bn1=string.upper(bpms[i][1])
            bns1=bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),betay[bn1][0],betay[bn1][1],MADTwiss_ac.BETY[MADTwiss_ac.indx[bn1]],MADTwiss_ac.MUY[MADTwiss_ac.indx[bn1]],betay_rescale[bn1][0],betay_rescale[bn1][1] ]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free amp beta
        if with_ac_calc:

            #-- from eq
            try:
                # Since invJyf(return_value[3]) is not used, slice the return value([:3]) (vimaier)
                [betayf,rmsbbyf,bpmsf] = ( GetFreeBetaFromAmp_Eq(MADTwiss_ac,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',beam_direction,accel) )[:3]

                #-- Rescaling
                betayf_ratio=0
                skipped_bpmyf=[]
                arcbpms=filterbpm(bpmsf)
                for bpm in arcbpms:
                    name=string.upper(bpm[1]) # second entry is the name
                    #Skip BPM with strange data
                    if abs(betayPhaseCopyf[name][0]/betayf[name][0])>10: skipped_bpmyf.append(name)
                    elif (betayf[name][0]<0 or betayPhaseCopyf[name][0]<0): skipped_bpmyf.append(name)
                    elif abs(betayPhaseCopyf[name][0]/betayf[name][0])<0.1: skipped_bpmyf.append(name)
                    else: betayf_ratio=betayf_ratio+(betayPhaseCopyf[name][0]/betayf[name][0])
                try: betayf_ratio=betayf_ratio/(len(arcbpms)-len(skipped_bpmyf))
                except: betayf_ratio=1

                tfs_file = files_dict['getampbetay_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(Q1f))
                tfs_file.add_descriptor("Q2", "%le", str(Q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf))
                tfs_file.add_descriptor("RescalingFactor", "%le", str(betayf_ratio))
                tfs_file.add_column_names(["NAME","S","COUNT","BETY","BETYSTD","BETYMDL","MUYMDL","BETYRES","BETYSTDRES"])
                tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(0,len(bpmsf)):
                    bn1=string.upper(bpmsf[i][1])
                    bns1=bpmsf[i][0]
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),betayf[bn1][0],betayf[bn1][1],MADTwiss.BETY[MADTwiss.indx[bn1]],MADTwiss.MUY[MADTwiss.indx[bn1]],(betayf_ratio*betayf[bn1][0]),(betayf_ratio*betayf[bn1][1]) ]
                    tfs_file.add_table_row(list_row_entries)
            except: 
                pass

            #-- from the model
            # Since invJyf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)            
            [betayf2,rmsbbyf2,bpmsf2] = ( getFreeAmpBeta(betay,rmsbby,bpms,invJy,MADTwiss_ac,MADTwiss,'V') )[:3]
            betayf2_rescale=getFreeAmpBeta(betay_rescale,rmsbby,bpms,invJy,MADTwiss_ac,MADTwiss,'V')[0]
            tfs_file = files_dict['getampbetay_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(Q1f))
            tfs_file.add_descriptor("Q2", "%le", str(Q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf2))
            tfs_file.add_descriptor("RescalingFactor", "%le", str(betay_ratio))
            tfs_file.add_column_names(["NAME","S","COUNT","BETY","BETYSTD","BETYMDL","MUYMDL","BETYRES","BETYSTDRES"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpmsf2)):
                bn1=string.upper(bpmsf2[i][1])
                bns1=bpmsf2[i][0]
                list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),betayf2[bn1][0],betayf2[bn1][1],MADTwiss.BETY[MADTwiss.indx[bn1]],MADTwiss.MUY[MADTwiss.indx[bn1]],betayf2_rescale[bn1][0],betayf2_rescale[bn1][1] ]
                tfs_file.add_table_row(list_row_entries)

    #-------- START IP
    print 'Calculating IP'

    if "LHC" in accel:
        tfs_file = files_dict['getIP.out']
        tfs_file.add_column_names(["NAME","BETASTARH","BETASTARHMDL","H","PHIH","PHIXH","PHIHMDL","BETASTARV","BETASTARVMDL","V","PHIV","PHIYV","PHIVMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
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
            list_row_entries = ['"IP'+ip+'"',betahor[1],betahor[4],betahor[2],betahor[3],betahor[6],betahor[5],betaver[1],betaver[4],betaver[2],betaver[3],betaver[6],betaver[5] ]
            tfs_file.add_table_row(list_row_entries)


        #-- Parameters at IP1, IP2, IP5, and IP8
        IPx=GetIP2(MADTwiss_ac,ListOfZeroDPPX,Q1,'H',beam_direction,accel,lhcphase)
        IPy=GetIP2(MADTwiss_ac,ListOfZeroDPPY,Q2,'V',beam_direction,accel,lhcphase)
        tfs_file_x = files_dict['getIPx.out']
        tfs_file_x.add_column_names(["NAME","BETX","BETXSTD","BETXMDL","ALFX","ALFXSTD","ALFXMDL","BETX*","BETX*STD","BETX*MDL","SX*","SX*STD","SX*MDL","rt(2JX)","rt(2JX)STD" ])
        tfs_file_x.add_column_datatypes(["%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le"])

        tfs_file_y = files_dict['getIPy.out']
        tfs_file_y.add_column_names(["NAME","BETY","BETYSTD","BETYMDL","ALFY","ALFYSTD","ALFYMDL","BETY*","BETY*STD","BETY*MDL","SY*","SY*STD","SY*MDL","rt(2JY)","rt(2JY)STD" ])
        tfs_file_y.add_column_datatypes(["%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le"])

        for bpm_name in ('IP1','IP5','IP8','IP2'):
            try:
                list_row_entries = ['"'+bpm_name+'"']
                for k in IPx[bpm_name]:
                    list_row_entries.append(k) 
                tfs_file_x.add_table_row(list_row_entries)
            except: pass
            try:
                list_row_entries = ['"'+bpm_name+'"']
                for k in IPy[bpm_name]: 
                    list_row_entries.append(k) 
                tfs_file_y.add_table_row(list_row_entries)
            except: pass

        #-- ac to free parameters at IP1, IP2, IP5, and IP8
        if with_ac_calc:

            #-- From Eq
            IPxf=GetFreeIP2_Eq(MADTwiss,ListOfZeroDPPX,Q1,Q1f,acphasex_ac2bpmac,'H',beam_direction,accel,lhcphase)
            IPyf=GetFreeIP2_Eq(MADTwiss,ListOfZeroDPPY,Q2,Q2f,acphasey_ac2bpmac,'V',beam_direction,accel,lhcphase)
            
            tfs_file_x = files_dict['getIPx_free.out']
            tfs_file_x.add_column_names(["NAME","BETX","BETXSTD","BETXMDL","ALFX","ALFXSTD","ALFXMDL","BETX*","BETX*STD","BETX*MDL","SX*","SX*STD","SX*MDL","rt(2JXD)","rt(2JXD)STD" ])
            tfs_file_x.add_column_datatypes(["%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le"])
            
            tfs_file_y = files_dict['getIPy_free.out']
            tfs_file_y.add_column_names(["NAME","BETY","BETYSTD","BETYMDL","ALFY","ALFYSTD","ALFYMDL","BETY*","BETY*STD","BETY*MDL","SY*","SY*STD","SY*MDL","rt(2JYD)","rt(2JYD)STD" ])
            tfs_file_y.add_column_datatypes(["%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le"])
            
            for bpm_name in ('IP1','IP5','IP8','IP2'):
                try:
                    list_row_entries = ['"'+bpm_name+'"']
                    for k in IPxf[bpm_name]: 
                        list_row_entries.append(k) 
                    tfs_file_x.add_table_row(list_row_entries)
                except: pass
                try:
                    list_row_entries = ['"'+bpm_name+'"']
                    for k in IPyf[bpm_name]:
                        list_row_entries.append(k) 
                    tfs_file_y.add_table_row(list_row_entries)
                except: pass

            #-- From model
            IPxf2=GetFreeIP2(MADTwiss,MADTwiss_ac,IPx,'H',accel)
            IPyf2=GetFreeIP2(MADTwiss,MADTwiss_ac,IPy,'V',accel)
            
            tfs_file_x = files_dict['getIPx_free2.out']
            tfs_file_x.add_column_names(["NAME","BETX","BETXSTD","BETXMDL","ALFX","ALFXSTD","ALFXMDL","BETX*","BETX*STD","BETX*MDL","SX*","SX*STD","SX*MDL","rt(2JXD)","rt(2JXD)STD" ])
            tfs_file_x.add_column_datatypes(["%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le"])
            
            tfs_file_y = files_dict['getIPy_free2.out']
            tfs_file_y.add_column_names(["NAME","BETY","BETYSTD","BETYMDL","ALFY","ALFYSTD","ALFYMDL","BETY*","BETY*STD","BETY*MDL","SY*","SY*STD","SY*MDL","rt(2JYD)","rt(2JYD)STD" ])
            tfs_file_y.add_column_datatypes(["%s","%le","%le",      "%le",   "%le",  "%le"     ,"%le", "%le",   "%le",      "%le",  "%le",  "%le",  "%le",  "%le",     "%le"])
            for bpm_name in ('IP1','IP5','IP8','IP2'):
                try:
                    list_row_entries = ['"'+bpm_name+'"']
                    for k in IPxf2[bpm_name]: 
                        list_row_entries.append(k) 
                    tfs_file_x.add_table_row(list_row_entries)
                except: 
                    pass
                try:
                    list_row_entries = ['"'+bpm_name+'"']
                    for k in IPyf2[bpm_name]: 
                        list_row_entries.append(k) 
                    tfs_file_y.add_table_row(list_row_entries)
                except: 
                    pass
            
        #-- IP beta* and phase from phase only
        try:    IPfromphase=GetIPFromPhase(MADTwiss_ac,phasex,phasey,accel)
        except: print 'No output from IP from phase. H or V file missing?'
        tfs_file = files_dict['getIPfromphase.out']
        tfs_file.add_column_names(["NAME","2L","BETX*","BETX*STD","BETX*MDL","BETY*","BETY*STD","BETY*MDL","PHX","PHXSTD","PHXMDL","PHY","PHYSTD","PHYMDL" ])
        tfs_file.add_column_datatypes(["%s","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le",  "%le"])
        for bpm_name in ('IP1','IP5','IP8','IP2'):
            list_row_entries = ['"'+bpm_name+'"']
            try:
                for k in IPfromphase[bpm_name]: 
                    list_row_entries.append(k)
                tfs_file.add_table_row(list_row_entries)
            except: 
                pass

        #-- ac to free beta*
        if with_ac_calc:

            #-- from eqs
            try:    
                IPfromphasef=GetIPFromPhase(MADTwiss,phasexf,phaseyf,accel)
            except: 
                pass
            tfs_file = files_dict['getIPfromphase_free.out']
            tfs_file.add_column_names(["NAME","2L","BETX*","BETX*STD","BETX*MDL","BETY*","BETY*STD","BETY*MDL","PHX","PHXSTD","PHXMDL","PHY","PHYSTD","PHYMDL" ])
            tfs_file.add_column_datatypes(["%s","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le",  "%le"])
            for bpm_name in ('IP1','IP5','IP8','IP2'):
                list_row_entries = ['"'+bpm_name+'"']
                try:
                    for k in IPfromphasef[bpm_name]: 
                        list_row_entries.append(k)
                    tfs_file.add_table_row(list_row_entries)
                except: 
                    pass

            #-- from the model
            try:    
                IPfromphasef2=GetIPFromPhase(MADTwiss,phasexf2,phaseyf2,accel)
            except: 
                pass
            tfs_file = files_dict['getIPfromphase_free2.out']
            tfs_file.add_column_names(["NAME","2L","BETX*","BETX*STD","BETX*MDL","BETY*","BETY*STD","BETY*MDL","PHX","PHXSTD","PHXMDL","PHY","PHYSTD","PHYMDL" ])
            tfs_file.add_column_datatypes(["%s","%le","%le",   "%le",   "%le",  "%le"     ,"%le",     "%le",   "%le",  "%le",  "%le",  "%le",  "%le",  "%le"])
            for bpm_name in ('IP1','IP5','IP8','IP2'):
                list_row_entries = ['"'+bpm_name+'"']
                try:
                    for k in IPfromphasef2[bpm_name]: 
                        list_row_entries.append(k)
                    tfs_file.add_table_row(list_row_entries)
                except: 
                    pass

    #-------- START Orbit
    ListOfCOX=[]
    if with_linx_for_zero_dppx:

        [cox,bpms]=GetCO(MADTwiss, ListOfZeroDPPX)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOx.out']
        tfs_file.add_descriptor("TABLE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("TYPE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("SEQUENCE", "%05s", '"'+accel+'"')
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME","S","COUNT","X","STDX","XMDL","MUXMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpms)):
            bn1=string.upper(bpms[i][1])
            bns1=bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),cox[bn1][0],cox[bn1][1],MADTwiss.X[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)

        ListOfCOX.append(cox)




    ListOfCOY=[]
    if with_liny_for_zero_dppy:
        [coy,bpms]=GetCO(MADTwiss, ListOfZeroDPPY)
        # The output file can be directly used for orbit correction with MADX
        tfs_file = files_dict['getCOy.out']
        tfs_file.add_descriptor("TABLE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("TYPE", "%05s", '"ORBIT"')
        tfs_file.add_descriptor("SEQUENCE", "%05s", '"'+accel+'"')
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME","S","COUNT","Y","STDY","YMDL","MUYMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le"])
        for i in range(0,len(bpms)):
            bn1=string.upper(bpms[i][1])
            bns1=bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPY),coy[bn1][0],coy[bn1][1],MADTwiss.Y[MADTwiss.indx[bn1]],MADTwiss.MUY[MADTwiss.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)


        ListOfCOY.append(coy)


    #-------- Orbit for non-zero DPP
    if with_linx_for_nonzero_dppx:

        k = 0
        for j in ListOfNonZeroDPPX:
            list_with_single_twiss=[]
            list_with_single_twiss.append(j)
            filename = 'getCOx_dpp_'+str(k+1)+'.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPX[k])
            tfs_file.add_descriptor("DPP", "%le", str(float(j.DPP)))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            [codpp,bpms]=GetCO(MADTwiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME","S","COUNT","X","STDX","XMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpms)):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]
                list_row_entries = ['"'+bn1+'"',bns1, len(ListOfZeroDPPX), codpp[bn1][0], codpp[bn1][1], MADTwiss.X[MADTwiss.indx[bn1]], MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)
            ListOfCOX.append(codpp)
            k += 1

    if with_liny_for_nonzero_dppy:
        k=0
        for j in ListOfNonZeroDPPY:
            list_with_single_twiss=[]
            list_with_single_twiss.append(j)
            filename = 'getCOy_dpp_'+str(k+1)+'.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPY[k])
            tfs_file.add_descriptor("DPP", "%le", str(float(j.DPP)))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            [codpp,bpms]=GetCO(MADTwiss, list_with_single_twiss)
            tfs_file.add_column_names(["NAME","S","COUNT","Y","STDY","YMDL","MUYMDL"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpms)):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]
                #TODO: why ListOfZeroDPPY.. above used ListOfNonZeroDPPY(vimaier)
                list_row_entries = ['"'+bn1+'"',bns1, len(ListOfZeroDPPY), codpp[bn1][0], codpp[bn1][1], MADTwiss.Y[MADTwiss.indx[bn1]], MADTwiss.MUY[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)
            ListOfCOY.append(codpp)
            k+=1



    #-------- START Dispersion

    if with_linx_for_zero_dppx and with_linx_for_nonzero_dppx:


        [nda,Dx,DPX,bpms] = NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, beta2_save, COcut)
        tfs_file = files_dict['getNDx.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME","S","COUNT","NDX","STDNDX","DX","DPX","NDXMDL","DXMDL","DPXMDL","MUXMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            ndmdl = MADTwiss.DX[MADTwiss.indx[bn1]] / math.sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfNonZeroDPPX),nda[bn1][0],nda[bn1][1],Dx[bn1][0],DPX[bn1],ndmdl,MADTwiss.DX[MADTwiss.indx[bn1]],MADTwiss.DPX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)




        [dxo,bpms] = DispersionfromOrbit(ListOfZeroDPPX,ListOfNonZeroDPPX,ListOfCOX,COcut,BPMU)

        DPX = GetDPX(MADTwiss,dxo,bpms)
        tfs_file = files_dict['getDx.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME","S","COUNT","DX","STDDX","DPX","DXMDL","DPXMDL","MUXMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
        for i in range(len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfNonZeroDPPX),dxo[bn1][0],dxo[bn1][1],DPX[bn1],MADTwiss.DX[MADTwiss.indx[bn1]],MADTwiss.DPX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)


    if with_liny_for_zero_dppy and with_liny_for_nonzero_dppy:
        [dyo,bpms] = DispersionfromOrbit(ListOfZeroDPPY,ListOfNonZeroDPPY,ListOfCOY,COcut,BPMU)
        DPY = GetDPY(MADTwiss,dyo,bpms)
        tfs_file = files_dict['getDy.out']
        tfs_file.add_descriptor("Q1", "%le", str(Q1))
        tfs_file.add_descriptor("Q2", "%le", str(Q2))
        tfs_file.add_column_names(["NAME","S","COUNT","DY","STDDY","DPY","DYMDL","DPYMDL","MUYMDL"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])

        for i in range(len(bpms)):
            bn1 = string.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfNonZeroDPPY),dyo[bn1][0],dyo[bn1][1],DPY[bn1],MADTwiss.DY[MADTwiss.indx[bn1]],MADTwiss.DPY[MADTwiss.indx[bn1]],MADTwiss.MUY[MADTwiss.indx[bn1]] ]
            tfs_file.add_table_row(list_row_entries)


    #-------- START coupling.
    print "Calculating coupling"

    if with_linx_for_zero_dppx and with_liny_for_zero_dppy:

        #-- Coupling in the model
        try:    MADTwiss.Cmatrix()
        except: pass

        #-- Main part
        if NBcpl == 1:
            
            # Avoids crashing the programm(vimaier)
            fwqwf = None
            fwqwf2 = None
            
            [fwqw,bpms] = GetCoupling1(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2)
            tfs_file = files_dict['getcouple.out']
            tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
            tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
            tfs_file.add_column_names(["NAME","S","COUNT","F1001W","FWSTD","Q1001W","QWSTD","MDLF1001R","MDLF1001I"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(len(bpms)):
                bn1 = string.upper(bpms[i][1])
                bns1 = bpms[i][0]
                try:    
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),(math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)),fwqw[bn1][0][1],fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag,MADTwiss.f1001[MADTwiss.indx(bn1)].real,MADTwiss.f1001[MADTwiss.indx(bn1)].imag,MADTwiss_ac.f1010[MADTwiss_ac.indx(bn1)].real,MADTwiss_ac.f1010[MADTwiss_ac.indx(bn1)].imag ]
                #-- Output zero if the model does not have couping parameters
                except: 
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),(math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)),fwqw[bn1][0][1],fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag, 0.0, 0.0]
                tfs_file.add_table_row(list_row_entries)

        elif NBcpl == 2:

            if accel=="SPS" or "RHIC" in accel:
                [phasexp,Q1,MUX,bpmsx]=GetPhases(MADTwiss,PseudoListX,'H',outputpath,beam_direction,accel,lhcphase)
                [phaseyp,Q2,MUY,bpmsy]=GetPhases(MADTwiss,PseudoListY,'V',outputpath,beam_direction,accel,lhcphase)
                [fwqw,bpms]=GetCoupling2(MADTwiss,PseudoListX,PseudoListY,Q1,Q2,phasexp,phaseyp,beam_direction,accel)
            else:
                [fwqw,bpms]=GetCoupling2(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2,phasexlist[0],phaseylist[0],beam_direction,accel)
            tfs_file = files_dict['getcouple.out']
            tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
            tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
            tfs_file.add_column_names(["NAME","S","COUNT","F1001W","FWSTD1","F1001R","F1001I","F1010W","FWSTD2","F1010R","F1010I","MDLF1001R","MDLF1001I","MDLF1010R","MDLF1010I"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(len(bpms)):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]
                try:    
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),(math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)),fwqw[bn1][0][1],fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag, math.sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2),fwqw[bn1][0][3],fwqw[bn1][0][2].real,fwqw[bn1][0][2].imag,MADTwiss_ac.f1001[MADTwiss_ac.indx[bn1]].real,MADTwiss_ac.f1001[MADTwiss_ac.indx[bn1]].imag,MADTwiss_ac.f1010[MADTwiss_ac.indx[bn1]].real,MADTwiss_ac.f1010[MADTwiss_ac.indx[bn1]].imag ]
                #-- Output zero if the model does not have couping parameters
                except: 
                    list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX), math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2),fwqw[bn1][0][1],fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag,math.sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2),fwqw[bn1][0][3],fwqw[bn1][0][2].real,fwqw[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0 ]
                tfs_file.add_table_row(list_row_entries)

            #-- ac to free coupling
            if with_ac_calc:

                #-- analytic eqs
                try:
                    [fwqwf,bpmsf]=GetFreeCoupling_Eq(MADTwiss,ListOfZeroDPPX,ListOfZeroDPPY,Q1,Q2,Q1f,Q2f,acphasex_ac2bpmac,acphasey_ac2bpmac,beam_direction)
                    tfs_file = files_dict['getcouple_free.out']
                    tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
                    tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
                    tfs_file.add_column_names(["NAME","S","COUNT","F1001W","FWSTD1","F1001R","F1001I","F1010W","FWSTD2","F1010R","F1010I","MDLF1001R","MDLF1001I","MDLF1010R","MDLF1010I"])
                    tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
                    for i in range(len(bpmsf)):
                        bn1=string.upper(bpmsf[i][1])
                        bns1=bpmsf[i][0]
                        try:    
                            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX), math.sqrt(fwqwf[bn1][0][0].real**2+fwqwf[bn1][0][0].imag**2),fwqwf[bn1][0][1],fwqwf[bn1][0][0].real,fwqwf[bn1][0][0].imag,math.sqrt(fwqwf[bn1][0][2].real**2+fwqwf[bn1][0][2].imag**2),fwqwf[bn1][0][3],fwqwf[bn1][0][2].real,fwqwf[bn1][0][2].imag,MADTwiss.f1001[MADTwiss.indx[bn1]].real,MADTwiss.f1001[MADTwiss.indx[bn1]].imag,MADTwiss.f1010[MADTwiss.indx[bn1]].real,MADTwiss.f1010[MADTwiss.indx[bn1]].imag ]
                        #-- Output zero if the model does not have couping parameters
                        except: 
                            list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX), math.sqrt(fwqwf[bn1][0][0].real**2+fwqwf[bn1][0][0].imag**2),fwqwf[bn1][0][1],fwqwf[bn1][0][0].real,fwqwf[bn1][0][0].imag,math.sqrt(fwqwf[bn1][0][2].real**2+fwqwf[bn1][0][2].imag**2),fwqwf[bn1][0][3],fwqwf[bn1][0][2].real,fwqwf[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0 ]
                        tfs_file.add_table_row(list_row_entries)
                except: 
                    traceback.print_exc()
                    pass

                #-- global factor
                [fwqwf2,bpmsf2]=getFreeCoupling(Q1f,Q2f,Q1,Q2,fwqw,MADTwiss,bpms)
                tfs_file = files_dict['getcouple_free2.out']
                tfs_file.add_descriptor("CG", "%le", str(fwqw['Global'][0]))
                tfs_file.add_descriptor("QG", "%le", str(fwqw['Global'][1]))
                tfs_file.add_column_names(["NAME","S","COUNT","F1001W","FWSTD1","F1001R","F1001I","F1010W","FWSTD2","F1010R","F1010I","MDLF1001R","MDLF1001I","MDLF1010R","MDLF1010I"])
                tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
                for i in range(len(bpmsf2)):
                    bn1=string.upper(bpmsf2[i][1])
                    bns1=bpmsf2[i][0]
                    try:    
                        list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),math.sqrt(fwqwf2[bn1][0][0].real**2+fwqwf2[bn1][0][0].imag**2),fwqwf2[bn1][0][1],fwqwf2[bn1][0][0].real,fwqwf2[bn1][0][0].imag,math.sqrt(fwqwf2[bn1][0][2].real**2+fwqwf2[bn1][0][2].imag**2),fwqwf2[bn1][0][3],fwqwf2[bn1][0][2].real,fwqwf2[bn1][0][2].imag, MADTwiss.f1001[MADTwiss.indx[bn1]].real,MADTwiss.f1001[MADTwiss.indx[bn1]].imag,MADTwiss.f1010[MADTwiss.indx[bn1]].real,MADTwiss.f1010[MADTwiss.indx[bn1]].imag ]
                    #-- Output zero if the model does not have couping parameters
                    except: 
                        list_row_entries = ['"'+bn1+'"',bns1,len(ListOfZeroDPPX),math.sqrt(fwqwf2[bn1][0][0].real**2+fwqwf2[bn1][0][0].imag**2),fwqwf2[bn1][0][1],fwqwf2[bn1][0][0].real,fwqwf2[bn1][0][0].imag,math.sqrt(fwqwf2[bn1][0][2].real**2+fwqwf2[bn1][0][2].imag**2),fwqwf2[bn1][0][3],fwqwf2[bn1][0][2].real,fwqwf2[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0 ]
                    tfs_file.add_table_row(list_row_entries)

        else:
            raise ValueError('Number of monitors for coupling analysis should be 1 or 2 (option -n)')

        #-- Convert to C-matrix:
        if with_ac_calc and (fwqwf is not None or fwqwf2 is not None):
            try:    
                [coupleterms,Qminav,Qminerr,bpms] = getCandGammaQmin(fwqwf,bpmsf,Q1f,Q2f,MADTwiss)
            except: 
                [coupleterms,Qminav,Qminerr,bpms] = getCandGammaQmin(fwqwf2,bpmsf2,Q1f,Q2f,MADTwiss)
        else: 
            [coupleterms,Qminav,Qminerr,bpms] = getCandGammaQmin(fwqw,bpms,Q1f,Q2f,MADTwiss)
        
        tfs_file = files_dict['getcoupleterms.out']
        tfs_file.add_descriptor("DQMIN", "%le", Qminav)
        tfs_file.add_descriptor("DQMINE", "%le", Qminerr)
        tfs_file.add_column_names(["NAME","S","DETC","DETCE","GAMMA","GAMMAE","C11","C12","C21","C22"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
        for bpm in bpms:
            bps = bpm[0]
            bpmm = bpm[1].upper()
            list_row_entries = [bpmm,bps,coupleterms[bpmm][0],coupleterms[bpmm][1],coupleterms[bpmm][2],coupleterms[bpmm][3],coupleterms[bpmm][4],coupleterms[bpmm][5],coupleterms[bpmm][6],coupleterms[bpmm][7] ]
            tfs_file.add_table_row(list_row_entries)

        #-- For chromatic coupling
        fwqw['DPP'] = 0

    #-------- Phase, Beta and coupling for non-zero DPP

    print "Phase and Beta for non-zero DPP"
    #TODO: what is a thingie??(vimaier)
    print "lenght of zerothingie "+ str(len(ListOfNonZeroDPPX))
    print "lenght of zerothingie "+ str(len(ListOfNonZeroDPPY))


    if with_linx_for_nonzero_dppx:
        plane='H'
        k=0
        for j in ListOfNonZeroDPPX:
            dpop=float(j.DPP)
            list_with_single_twiss=[]
            list_with_single_twiss.append(j)
            filename = 'getphasex_dpp_'+str(k+1)+'.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPX[k])
            tfs_file.add_descriptor("DPP", "%le", str(dpop))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
                
            DPPTwiss=ConstructOffMomentumModel(MADTwiss,dpop,BPMdictionary)
            [phasex,Q1DPP,MUX,bpms]=GetPhases(DPPTwiss,list_with_single_twiss,Q1,plane,outputpath,beam_direction,accel,lhcphase)
            phasex['DPP']=dpop
            phasexlist.append(phasex)

            tfs_file.add_descriptor("Q1DPP", "%le", str(Q1DPP))
            tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASE","STDPH","PHXMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
            
            length_bpms = len(bpms)
            for i in range(0,length_bpms):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]
               
                index = (i+1) % length_bpms
                bn2 = string.upper(bpms[index][1])
                bns2 = bpms[index][0]
                
                try:
                    phmdl=phasexlist[0][bn1][4]
                except:
                    phmdl=0.0
                #phmdl=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
                list_row_entries = ['"'+bn1+'" ','"'+bn2+'"',bns1,bns2,1,phasex[bn1][0],phasex[bn1][1],phmdl,MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)

            betax={}
            alfax={}
            rmsbbx=0.
            [betax,rmsbbx,alfax,bpms]=BetaFromPhase(MADTwiss,list_with_single_twiss,phasex,plane)
            betax['DPP']=dpop
            betaxa={}
            [betaxa,rmsbbx,bpms,invJx]=BetaFromAmplitude(MADTwiss,list_with_single_twiss,plane)
            betaxa['DPP']=dpop
            betaxalist.append(betaxa)
            
            filename = 'getbetax_dpp_'+str(k+1)+'.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPX[k])
            tfs_file.add_descriptor("DPP", "%le", str(dpop))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            tfs_file.add_column_names(["NAME","S","COUNT","BETX","ERRBETX","STDBETX","ALFX","ERRALFX","STDALFX","BETXMDL","ALFXMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpms)):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]
                list_row_entries = ['"'+bn1+'" ',bns1,len(ListOfZeroDPPX),betax[bn1][0],betax[bn1][1],betax[bn1][2],alfax[bn1][0],alfax[bn1][1],alfax[bn1][2],MADTwiss.BETX[MADTwiss.indx[bn1]],MADTwiss.ALFX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)
            k+=1


    if with_liny_for_nonzero_dppy:
        plane='V'
        k=0

        for j in ListOfNonZeroDPPY:
            dpop = float(j.DPP)
            list_with_single_twiss = []
            list_with_single_twiss.append(j)
            filename = 'getphasey_dpp_'+str(k+1)+'.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPY[k])
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            
            DPPTwiss = ConstructOffMomentumModel(MADTwiss,dpop,BPMdictionary)
            [phasey,Q2DPP,MUY,bpms] = GetPhases(DPPTwiss,list_with_single_twiss,Q2,plane,outputpath,beam_direction,accel,lhcphase)
            phasey['DPP'] = dpop
            phaseylist.append(phasey)
            
            tfs_file.add_descriptor("Q2DPP", "%le", str(Q2DPP))
            tfs_file.add_column_names(["NAME","NAME2","S","S1","COUNT","PHASE","STDPH","PHYMDL","MUYMDL"])
            tfs_file.add_column_datatypes(["%s","%s","%le","%le","%le","%le","%le","%le","%le"])
            
            length_bpms = len(bpms)
            for i in range(0,length_bpms):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]

                index = (i+1) % length_bpms
                bn2 = string.upper(bpms[index][1])
                bns2 = bpms[index][0]

                try:
                    phmdl=phaseylist[0][bn1][4]
                except:
                    phmdl=0.0
                #phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
                list_row_entries = ['"'+bn1+'" ','"'+bn2+'" ',bns1,bns2,1,phasey[bn1][0],phasey[bn1][1],phmdl,MADTwiss.MUY[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)

            betay={}
            alfay={}
            rmsbby=0.
            [betay,rmsbby,alfay,bpms]=BetaFromPhase(DPPTwiss,list_with_single_twiss,phasey,plane)
            betay['DPP']=dpop
            betaya={}
            [betaya,rmsbby,bpms,invJy]=BetaFromAmplitude(DPPTwiss,list_with_single_twiss,plane)
            betaya['DPP']=dpop
            betayalist.append(betaya)
            
            filename = 'getbetay_dpp_'+str(k+1)+'.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_filename_to_getllm_header(FileOfNonZeroDPPY[k])
            tfs_file.add_descriptor("DPP", "%le", str(dpop))
            tfs_file.add_descriptor("Q1", "%le", str(Q1))
            tfs_file.add_descriptor("Q2", "%le", str(Q2))
            #TODO: check if it should be Y instead of X in the column names since it is getbetaY_dpp...out (vimaier)
            tfs_file.add_column_names(["NAME","S","COUNT","BETX","ERRBETX","STDBETX","ALFX","ERRALFX","STDALFX","BETXMDL","ALFXMDL","MUXMDL"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])
            for i in range(0,len(bpms)):
                bn1=string.upper(bpms[i][1])
                bns1=bpms[i][0]
                list_row_entries = ['"'+bn1+'" ',bns1,len(ListOfZeroDPPY),betay[bn1][0],betay[bn1][1],betay[bn1][2],alfay[bn1][0],alfay[bn1][1],alfay[bn1][2],MADTwiss.BETX[MADTwiss.indx[bn1]],MADTwiss.ALFX[MADTwiss.indx[bn1]],MADTwiss.MUX[MADTwiss.indx[bn1]] ]
                tfs_file.add_table_row(list_row_entries)

            k += 1


    if with_liny_for_nonzero_dppy and with_linx_for_nonzero_dppx:

        if len(ListOfNonZeroDPPX)!=len(ListOfNonZeroDPPY):
            raise ValueError("list of dppx is not equal list of dppy")


        for j in range(len(ListOfNonZeroDPPX)):
            dpop=float(ListOfNonZeroDPPX[j].DPP)
            list_with_single_twiss_x = []
            list_with_single_twiss_y = []
            list_with_single_twiss_x.append(ListOfNonZeroDPPX[j])
            list_with_single_twiss_y.append(ListOfNonZeroDPPY[j])
            ### coupling
            try:
                MADTwiss.Cmatrix()
            except:
                0.0

            if accel=="SPS" or "RHIC" in accel:
                plane='H'
                [phasexp,Q1,MUX,bpmsx] = GetPhases(MADTwiss,PseudoListX,plane,outputpath,beam_direction,accel,lhcphase)
                plane = 'V'
                [phaseyp,Q2,MUY,bpmsy] = GetPhases(MADTwiss,PseudoListY,plane,outputpath,beam_direction,accel,lhcphase)
                [fwqw,bpms] = GetCoupling2(MADTwiss, PseudoListX, PseudoListY, Q1, Q2, phasexp, phaseyp, beam_direction, accel)
            elif NBcpl == 1:
                [fwqw,bpms] = GetCoupling1(MADTwiss, list_with_single_twiss_x, list_with_single_twiss_y, Q1, Q2)
            elif NBcpl == 2:
                print phasexlist[j+1]['DPP'],dpop
                [fwqw,bpms] = GetCoupling2(MADTwiss, list_with_single_twiss_x, list_with_single_twiss_y, Q1, Q2, phasexlist[j+1], phaseylist[j+1], beam_direction, accel)
                if with_ac_calc:
                    [fwqw,bpms] = getFreeCoupling(Q1f,Q2f,Q1,Q2,fwqw,MADTwiss,bpms)

            else:
                raise ValueError('Number of monitors for coupling analysis (option -n) should be 1 or 2.')

            fwqw['DPP'] = dpop

            ####
    #---------------------------------------- Start getsextupoles @ Glenn Vanbavinckhove

    if not higher_order:
        print "Not analysing higher order..."
        return 0


    # fsex1200 and 2100 is never used and ouputname is wrong(vimaier)
#     fsex1200 = open(outputpath+'getsex1200.out','w')
#     fsex1200.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')
#     fsex2100 = open(outputpath+'getsex1200.out','w')
#     fsex2100.write('@ MAD_FILE %s "'+twiss_model_file+'"'+'\n')


    if TBTana=="SUSSIX":
    #-> 1) f3000 line (-2,0)
    #-> 2) f1200 line  (2,0)
    #-> 3) f2100 line  (0,0)

    # global stuff

    # 1)

        htot, afactor, pfactor = Getsextupole(MADTwiss,ListOfZeroDPPX,phasexlist[0],Q1f,3,0)
        
        filename = 'getsex3000.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_descriptor("f2h_factor", "%le", str(afactor))
        tfs_file.add_descriptor("p_f2h_factor", "%le", str(pfactor))



        tfs_file.add_column_names(["NAME","S","AMP_20","AMP_20std","PHASE_20","PHASE_20std","f3000","f3000std","phase_f_3000","phase_f_3000std","h3000","h3000_std","phase_h_3000","phase_h_3000_std"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])

        for bpm in htot:
            li=htot[bpm]
            list_row_entries = [li[0],li[1],li[2],li[3],li[4],li[5],li[6],li[7],li[8],li[9],li[10],li[11],li[12],li[13] ]
            tfs_file.add_table_row(list_row_entries)




        # --------------------------------------- end getsextupoles
        #---------------------------------------- begin getchiterms @ Glenn Vanbavinckhove
        #-> 1) chi3000
        #-> 2) chi1010
        #-> 2) chi4000

        # 1) chi3000
        filename = 'getchi3000.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_column_names(["NAME","S","S1","S2","X3000","X3000i","X3000r","X3000RMS","X3000PHASE","X3000PHASERMS","X3000M","X3000Mi","X3000Mr","X3000MPHASE"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])

        files=[ListOfZeroDPPX,ListOfZeroDPPY]
        name='chi3000'
        plane='H'

        [dbpms,POS,XItot,XIMODEL]=getChiTerms(MADTwiss,files,plane,name,ListOfZeroDPPX,ListOfZeroDPPY)

        for i in range(0,len(dbpms)-2):
            bn=string.upper(dbpms[i][1])

            list_row_entries = ['"'+bn+'"',POS[0][i],POS[1][i],POS[2][i],XItot[0][i],XItot[1][i],XItot[2][i],XItot[3][i],XItot[4][i],XItot[5][i],XIMODEL[0][i],XIMODEL[1][i],XIMODEL[2][i],XIMODEL[3][i] ]
            tfs_file.add_table_row(list_row_entries)

    # 2) chi1010

        if  accel != 'SPS':
            filename = 'getchi1010.out'
            files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
            tfs_file = files_dict[filename]
            tfs_file.add_column_names(["NAME","S","X1010","X1010RMS","X1010PHASE","X1010PHASERMS","X1010M","X1010MPHASE"])
            tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le"])

            files = [ListOfZeroDPPX,ListOfZeroDPPY]
            name = 'chi1010'
            plane = 'H'

            [dbpms,XItot] = getchi1010(MADTwiss,files,plane,name,bn1,ListOfZeroDPPX,ListOfZeroDPPY)

            for i in range(len(dbpms)-2):
                bn = string.upper(dbpms[i][1])
                bns = dbpms[i][0]
                list_row_entries = ['"'+bn+'"',bns,XItot[0][i],XItot[1][i],XItot[2][i],XItot[3][i],'0','0' ]
                tfs_file.add_table_row(list_row_entries)

    # 1) chi4000


        
        filename = 'getchi4000.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_column_names(["NAME","S","S1","S2","X4000","X4000i","X4000r","X4000RMS","X4000PHASE","X4000PHASERMS","X4000M","X4000Mi","X4000Mr","X4000MPHASE"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])

#         files=[ListOfZeroDPPX,ListOfZeroDPPY]
#         name='chi4000'
#         plane='H'
# 
#         [dbpms,POS,XItot,XIMODEL]=getChiTerms(MADTwiss,files,plane,name,ListOfZeroDPPX,ListOfZeroDPPY)
# 
#         for i in range(0,len(dbpms)-2):
#                    bn=string.upper(dbpms[i][1])
#                 list_row_entries = ['"'+bn+'"',POS[0][i],POS[1][i],POS[2][i],XItot[0][i],XItot[1][i],XItot[2][i],XItot[3][i],XItot[4][i],XItot[5][i],XIMODEL[0][i],XIMODEL[1][i],XIMODEL[2][i],XIMODEL[3][i] ]
#                 tfs_file.add_table_row(list_row_entries)





    #---------------------------------------- end chiterms
    #-----------------------------------------begin octupole
    #->  1) f4000 (-3,0)
    #
        filename = 'getoct4000.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_column_names(["NAME","S","AMP_30","AMP_30RMS","PHASE_30","PHASE_30RMS","H4000","H4000I","H4000R","H4000RMS","H4000PHASE","H4000PHASERMS","H4000M","H4000MI","H4000MR","HMPHASE4000"])
        tfs_file.add_column_datatypes(["%s","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])

        files=[ListOfZeroDPPX,ListOfZeroDPPY]
        plane='H'
        name='f4000'

        #TODO: f2100M and NAMES not defined befor. What to do? Delete section? (vimaier)
#         [A,h,hMODEL,dbpms]=Getoctopole(MADTwiss,plane,files,phasexlist[0],Q,name,f2100M,NAMES)
# 
#         for i in range(0,len(dbpms)-1):
#     
#                 bn=string.upper(dbpms[i][1])
#                 bns=dbpms[i][0]
#                 list_row_entries = ['"'+bn+'"',bns,A[0][i],A[1][i],A[2][i],A[3][i],h[0][i],h[1][i],h[2][i],h[3][i],h[4][i],h[5][i],hMODEL[0][i],hMODEL[1][i],hMODEL[2][i],hMODEL[3][i] ]
#                 tfs_file.add_table_row(list_row_entries)

    #-----------------------------------------end octupole

    #----------------------------- begin get Q,JX,delta

    files=[ListOfZeroDPPX+ListOfNonZeroDPPX,ListOfZeroDPPY+ListOfNonZeroDPPY]

    filename = 'getkick.out'
    files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
    tfs_file = files_dict[filename]
    tfs_file.add_descriptor("RescalingFactor_for_X", "%le", str(betax_ratio))
    tfs_file.add_descriptor("RescalingFactor_for_Y", "%le", str(betay_ratio))
    tfs_file.add_column_names(["DPP","QX","QXRMS","QY","QYRMS","sqrt2JX","sqrt2JXSTD","sqrt2JY","sqrt2JYSTD","2JX","2JXSTD","2JY","2JYSTD","sqrt2JXRES","sqrt2JXSTDRES","sqrt2JYRES","sqrt2JYSTDRES","2JXRES","2JXSTDRES","2JYRES","2JYSTDRES"])
    tfs_file.add_column_datatypes(["%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])

    [invarianceJx,invarianceJy,tune,tuneRMS,dpp]=getkick(files,MADTwiss)

    for i in range(0,len(dpp)):
        list_row_entries = [dpp[i],tune[0][i],tuneRMS[0][i],tune[1][i],tuneRMS[1][i],invarianceJx[i][0],invarianceJx[i][1],invarianceJy[i][0],invarianceJy[i][1],(invarianceJx[i][0]**2),(2*invarianceJx[i][0]*invarianceJx[i][1]),(invarianceJy[i][0]**2),(2*invarianceJy[i][0]*invarianceJy[i][1]),(invarianceJx[i][0]/math.sqrt(betax_ratio)),(invarianceJx[i][1]/math.sqrt(betax_ratio)),(invarianceJy[i][0]/math.sqrt(betay_ratio)),(invarianceJy[i][1]/math.sqrt(betay_ratio)),(invarianceJx[i][0]**2/betax_ratio),(2*invarianceJx[i][0]*invarianceJx[i][1]/betax_ratio),(invarianceJy[i][0]**2/betay_ratio),(2*invarianceJy[i][0]*invarianceJy[i][1]/betay_ratio) ]
        tfs_file.add_table_row(list_row_entries)




    if with_ac_calc:
        files=[ListOfZeroDPPX+ListOfNonZeroDPPX,ListOfZeroDPPY+ListOfNonZeroDPPY]
        
        filename = 'getkickac.out'
        files_dict[filename] = utils.tfs_file.TfsFile(filename).add_getllm_header(VERSION, twiss_model_file)
        tfs_file = files_dict[filename]
        tfs_file.add_descriptor("RescalingFactor_for_X", "%le", str(betaxf_ratio))
        tfs_file.add_descriptor("RescalingFactor_for_Y", "%le", str(betayf_ratio))
        tfs_file.add_column_names(["DPP","QX","QXRMS","QY","QYRMS","sqrt2JX","sqrt2JXSTD","sqrt2JY","sqrt2JYSTD","2JX","2JXSTD","2JY","2JYSTD","sqrt2JXRES","sqrt2JXSTDRES","sqrt2JYRES","sqrt2JYSTDRES","2JXRES","2JXSTDRES","2JYRES","2JYSTDRES"])
        tfs_file.add_column_datatypes(["%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le","%le"])

        [invarianceJx,invarianceJy,tune,tuneRMS,dpp]=getkickac(MADTwiss_ac,files,Q1,Q2,Q1f,Q2f,acphasex_ac2bpmac,acphasey_ac2bpmac,beam_direction,lhcphase)

        for i in range(0,len(dpp)):
            list_row_entries = [dpp[i],tune[0][i],tuneRMS[0][i],tune[1][i],tuneRMS[1][i],invarianceJx[i][0],invarianceJx[i][1],invarianceJy[i][0],invarianceJy[i][1],(invarianceJx[i][0]**2),(2*invarianceJx[i][0]*invarianceJx[i][1]),(invarianceJy[i][0]**2),(2*invarianceJy[i][0]*invarianceJy[i][1]),(invarianceJx[i][0]/math.sqrt(betax_ratio)),(invarianceJx[i][1]/math.sqrt(betax_ratio)),(invarianceJy[i][0]/math.sqrt(betay_ratio)),(invarianceJy[i][1]/math.sqrt(betay_ratio)),(invarianceJx[i][0]**2/betax_ratio),(2*invarianceJx[i][0]*invarianceJx[i][1]/betax_ratio),(invarianceJy[i][0]**2/betay_ratio),(2*invarianceJy[i][0]*invarianceJy[i][1]/betay_ratio) ]
            tfs_file.add_table_row(list_row_entries)


    for tfsfile in files_dict.itervalues():
        tfsfile.write_to_file(formatted=True)
        
    ####### -------------- end

if __name__=="__main__":
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
