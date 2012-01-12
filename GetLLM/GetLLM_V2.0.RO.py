## Python script to obtain Linear Lattice function and More -> GetLLM
## Version 1.51
## Version-up history:V1.0, 11/Feb/2008 by Masa. Aiba
##                    V1.1, 18/Feb/2008 Debugging, add model phase and tunes to output
##                                      add function to obtain DY
##                                      add chromatic parameter (phase for non zero DPP)
##                    V1.2, 22/Feb/2008 test version for beta with all BPM
##                    V1.3, 29/Feb/2008 beta from phases is improved, averaging beta1, 2 and 3
##                    V1.31, 12/Mar/2008 debugged on alpha3
##                    V1.4, 12/Mar/2008 modify output to fit latest TSF format and to meet requests from Rogelio
##                                      fix buggs in r.m.s. beta-beat and added to the output of getbetax/y.out
##                    V1.5 Update to option parser, include BPMdictionary to filter BPMs not in Model
##                                  Rogelio, 13 March 2008
##                    V1.51, 13/Mar/2008 Modify output to fit latest TSF format again. Add STD to beta.
##                    V1.6, 15/Jul/2008 Add the integer part of tunes - assuming that the phase advance is always less than 1.0.
##                    V1.71 27/Jul/2008 Add GetCO. Filter in dispersion calculation to exclude bad bpms.
##                    V1.8, 13/Aug/2008 Add GetCoupling. - Ref. note by A. Franchi, R. T. Garcia, G. Vanbavinckhove
##                                      "Computation of the Coupling Resonance Driving term f1001 and the coupling coefficient C
##                                       from turn-by-turn single-BPM data", 28/May/2008
##                                      The GetCoupling.py is initiated by Glenn V.
##                                      and imported into GetLLM / finalized by Masa. Aiba
##                                      Some bugs are fixed - you can run even without liny file and
##                                      find the results for horizontal plane only.
##                    V1.81, 1/Sep/2008 For an accelerator in which the beam goes opposite direction to the model as in LHCB2,
##                                       the beam direction parameter (bd) is added to GetPhases.
##                                       Bug in the phi13 for the last monitor and the last but one is fixed.
##                    V1.9, 21/Oct/2008 Add the beta from spectrum height.
##                    V1.91, 08/Dec/2008 Add option - SUSSIX or SVD for input file

##                    V1.99, 13/Feb/09 Add DPX! The MEMORIAL version for Version 1.** since all the linear lattice parameters become available!
##                                     Add the option for the Harmonic analysis.
##                                     Refine coding, add several lines of comment

##                    V2.0, 17/Feb/2009 Major version up to output "More", that is, other than the linear lattice parameters.
##                                      Add off-momentum lattice (dbeta/beta)/(dp/p)
##                                      Add SPS coupling with "Pseudo-double plane BPM"-it need at this moment a list containing
##                                      pre-paired H-V monitor. This should be replaced by a clever algorithm to find good pairs automatically.

## Usage1 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./
## Usage2 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./


## Some rules for variable name: Dictionary is used to contain the output of function
##                               Valiable containing 'm' is a value directly obtained from measurment data
##                               Valiable containing 'mdl' is a value related to model




from metaclass import *
from Numeric import *
from math import *
import cmath
import sys, pickle,os
#import operator
from string import *

# tentative solution for SPS pseudo double plane BPM
from SPSBPMpair import *


#######################################################
#                 Functions                           #
#######################################################

#------------

def modelIntersect(expbpms, model):
	bpmsin=[]
	for bpm in expbpms:
		try:
			check=model.indx[bpm[1].upper()]
			bpmsin.append(bpm)
		except:
			print bpm, "Not in Model"
	if len(bpmsin)==0:
		print "Zero intersection of Exp and Model"
		print "Please, provide a good Dictionary"
		print "Now we better leave!"
		sys.exit()			
	return bpmsin


def intersect(ListOfFile): 
	'''Pure intersection of all bpm names in all files '''
	if len(ListOfFile)==0:
		print "Nothing to intersect!!!!"
		sys.exit()
	z=ListOfFile[0].NAME
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



def GetPhases(MADTwiss,ListOfFiles,plane,outputpath,bd):

	commonbpms=intersect(ListOfFiles)
	commonbpms=modelIntersect(commonbpms, MADTwiss)
	
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
			phmdl13=phiLastAndLastButOne(phmdl13,madtune)
		elif i==len(commonbpms)-1:
			if plane=='H':
				madtune=MADTwiss.Q1 % 1.0
			elif plane=='V':
				madtune=MADTwiss.Q2 % 1.0
			if madtune>0.5:
				madtune=madtune-1.0
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
			if phm12<0: phm12+=1
			if phm13<0: phm13+=1
			phi12.append(phm12)
			phi13.append(phm13)
		phi12=array(phi12)
		phi13=array(phi13)
		if bd==-1: # for the beam circulating reversely to the model
			phi12=1.0-phi12
			phi13=1.0-phi13
		if phi12>0.9 and i !=len(commonbpms): # Very small phase advance could result in larger than 0.9 due to measurement error
			print 'Warning: there seems too large phase advance! '+bn1+' to '+bn2+' = '+str(phi12)+'in plane '+plane+', recommended to check.'
		phstd12=sqrt(average(phi12*phi12)-(average(phi12))**2.0)
		phstd13=sqrt(average(phi13*phi13)-(average(phi13))**2.0)
		phi12=average(phi12)
		phi13=average(phi13)
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
		phase[bn1]=[phi12,phstd12,phi13,phstd13,phmdl12,phmdl13]

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
		# Find beta1 from phases assuming model transfer matrix
		# Matrix M: BPM1-> BPM2
		# Matrix N: BPM1-> BPM3
		M11=sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
		M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
		N11=sqrt(betmdl3/betmdl1)*(cos(phmdl13)+alpmdl1*sin(phmdl13))
		N12=sqrt(betmdl1*betmdl3)*sin(phmdl13)
		denom=M11/M12-N11/N12
		numer=1/tan(ph2pi12)-1/tan(ph2pi13)
		beti1=numer/denom
		beterr1=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
		beterr1=beterr1+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
		beterr1=beterr1/abs(denom)
		denom=M12/M11-N12/N11
		numer=-M12/M11/tan(ph2pi12)+N12/N11/tan(ph2pi13)
		alfi1=numer/denom
		alferr1=abs(M12/M11*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
		alferr1=alferr1+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
		alferr1=alferr1/abs(denom)

		# Find beta2 from phases assuming model transfer matrix
		# Matrix M: BPM1-> BPM2
		# Matrix N: BPM2-> BPM3
		M22=sqrt(betmdl1/betmdl2)*(cos(phmdl12)-alpmdl2*sin(phmdl12))
		M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
		N11=sqrt(betmdl3/betmdl2)*(cos(phmdl23)+alpmdl2*sin(phmdl23))
		N12=sqrt(betmdl2*betmdl3)*sin(phmdl23)
		denom=M22/M12+N11/N12
		numer=1/tan(ph2pi12)+1/tan(ph2pi23)
		beti2=numer/denom
		beterr2=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
		beterr2=beterr2+abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
		beterr2=beterr2/abs(denom)
		denom=M12/M22+N12/N11
		numer=M12/M22/tan(ph2pi12)-N12/N11/tan(ph2pi23)
		alfi2=numer/denom
		alferr2=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
		alferr2=alferr2+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
		alferr2=alferr2/abs(denom)
		

		# Find beta3 from phases assuming model transfer matrix
		# Matrix M: BPM2-> BPM3
		# Matrix N: BPM1-> BPM3
		M22=sqrt(betmdl2/betmdl3)*(cos(phmdl23)-alpmdl3*sin(phmdl23))
		M12=sqrt(betmdl2*betmdl3)*sin(phmdl23)
		N22=sqrt(betmdl1/betmdl3)*(cos(phmdl13)-alpmdl3*sin(phmdl13))
		N12=sqrt(betmdl1*betmdl3)*sin(phmdl13)
		denom=M22/M12-N22/N12
		numer=1/tan(ph2pi23)-1/tan(ph2pi13)
		beti3=numer/denom
		beterr3=abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
		beterr3=beterr3+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
		beterr3=beterr3/abs(denom)
		denom=M12/M22-N12/N22
		numer=M12/M22/tan(ph2pi23)-N12/N22/tan(ph2pi13)
		alfi3=numer/denom
		alferr3=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
		alferr3=alferr3+abs(N22/N12*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
		alferr3=alferr3/abs(denom)

		betii.append([beti1,beterr1,beti2,beterr2,beti3,beterr3])
		alfii.append([alfi1,alferr1,alfi2,alferr2,alfi3,alferr3])


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
		beterr=sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/sqrt(3.)
		betstd=sqrt((betii[ib1][0]**2.+betii[ib2][2]**2.+betii[ib3][4]**2.)/3.-((betii[ib1][0]+betii[ib2][2]+betii[ib3][4])/3.)**2.)
		alfi=(alfii[ib1][0]+alfii[ib2][2]+alfii[ib3][4])/3.
		alferr=sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/sqrt(3.)
		alfstd=sqrt((alfii[ib1][0]**2.+alfii[ib2][2]**2.+alfii[ib3][4]**2.)/3.-((alfii[ib1][0]+alfii[ib2][2]+alfii[ib3][4])/3.)**2.)

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
	for i in range(0,len(commonbpms)):
		bn1=upper(commonbpms[i][1])
		Ampi=0.0
                root2Ji=0.0
		for j in ListOfFiles:
			if plane=='H':
				Ampi+=j.AMPX[j.indx[bn1]]
                                root2Ji+=j.PK2PK[j.indx[bn1]]/2.     
			elif plane=='V':
				Ampi+=j.AMPY[j.indx[bn1]]
                                root2Ji+=j.PK2PK[j.indx[bn1]]/2.
		Ampi=Ampi/len(ListOfFiles)
                root2Ji=root2Ji/len(ListOfFiles)
		Amp.append(Ampi)
                
		if plane=='H':
                        tembetax=MADTwiss.BETX[MADTwiss.indx[bn1]]
			SumA+=Ampi**2/tembetax
                        root2J.append(root2Ji/sqrt(tembetax))
		if plane=='V':
                        tembetay=MADTwiss.BETY[MADTwiss.indx[bn1]]
			SumA+=Ampi**2/tembetay
                        root2J.append(root2Ji/sqrt(tembetay))
		
	Kick=SumA/len(commonbpms) # Assuming the average of beta is constant
        root2J=array(root2J)
        root2Jave=average(root2J)
        root2Jrms=sqrt(average(root2J*root2J)-root2Jave**2)
	
	delbeta=[]
	for i in range(0,len(commonbpms)):
		bn1=upper(commonbpms[i][1])
		beta[bn1]=Amp[i]**2/Kick
		if plane=='H':
			betmdl=MADTwiss.BETX[MADTwiss.indx[bn1]]
		elif plane=='V':
			betmdl=MADTwiss.BETY[MADTwiss.indx[bn1]]
		delbeta.append((beta[bn1]-betmdl)/betmdl)

	delbeta=array(delbeta)
	rmsbb=sqrt(average(delbeta*delbeta))
	return [beta,rmsbb,commonbpms, root2Jave, root2Jrms]


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
		corms=sqrt(coi2/len(ListOfFiles)-coi**2)
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
		mydp.append(j.DPP)
		gf.append(0.0) # just to initialize the array, may not be a clever way...
	mydp=array(mydp)
	wf=array(mydp)/sum(mydp)*len(mydp) #Weight for the average depending on DPP


	# Find the global factor
	nd=[]
	ndmdl=[]
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
			coi=0

	ndmdl=array(ndmdl)
	avemdl=average(ndmdl)

	gf=array(gf)
	gf=gf/avemdl/len(commonbpmsALL)


	# Find normalized dispersion and Dx construction
	nd=array(nd)
	Dx={}
	dummy=0.0 # dummy for DXSTD
	bpms=[]
	for i in range(0,len(commonbpmsALL)):
		ndi=[]
		bn1=upper(commonbpmsALL[i][1])
		bns1=commonbpmsALL[i][0]
		try:
			coac[bn1]
			for j in range(0,nzdpp): # the range(0,nzdpp) instead of ListOfZeroDPPX is used because the index j is used in the loop
				ndi.append(nd[i][j]/gf[j])
			ndi=array(ndi)
			ndstd=sqrt(average(ndi*ndi)-(average(ndi))**2.0)
			ndas=average(wf*ndi)
			nda[bn1]=[ndas,ndstd]
			Dx[bn1]=[nda[bn1][0]*sqrt(betax[bn1]),dummy]
			bpms.append([bns1,bn1])
		except:
			0
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



#-----------

def DispersionfromOrbit(ListOfZeroDPP,ListOfNonZeroDPP,ListOfCO,COcut):


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
		mydp.append(j.DPP)
	mydp=array(mydp)
	wf=array(mydp)/sum(mydp)*len(mydp) #Weitghs for the average

	dco={} # Dictionary for the output containing [(averaged)Disp, rms error]
	bpms=[]
	for i in range(0,len(commonbpmsALL)):
		bn1=upper(commonbpmsALL[i][1])
		bns1=commonbpmsALL[i][0]

		try:
			coi=coac[bn1]
			dcoi=[]
			for j in ListOfNonZeroDPP:
				dcoi.append((j.CO[j.indx[bn1]]-coi[0])/j.DPP)
			dcoi=array(dcoi)
			dcostd=sqrt(average(dcoi*dcoi)-(average(dcoi))**2.0)
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
	for i in range(0,len(dbpms)):
		bn1=upper(dbpms[i][1])

		fij=[]
		qij=[]
		badbpm=0
		for j in range(0,len(ListOfZeroDPPX)):
			jx=ListOfZeroDPPX[j]
			jy=ListOfZeroDPPY[j]
			C01ij=jx.AMP01[jx.indx[bn1]]
			C10ij=jy.AMP10[jy.indx[bn1]]
			fij.append(0.5*atan(sqrt(C01ij*C10ij)))

			q1=(jx.MUX[jx.indx[bn1]]-jy.PHASE10[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
			q2=(jx.PHASE01[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]]-0.25)%1.0
			q1=(0.5-q1)%1.0 # This sign change in the real part is to comply with MAD output
			q2=(0.5-q2)%1.0

			if abs(q1-q2)<0.25: 
				qij.append((q1+q2)/2.0)
			elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
				qij.append(q1) # Note that q1 and q2 are confined 0. to 1.
			else:
				badbpm=1

		
		if badbpm==0:
			fij=array(fij)
			fi=average(fij)
			fistd=sqrt(average(fij*fij)-(average(fij))**2.0)
			qij=array(qij)
			qi=average(qij)
			qistd=sqrt(average(qij*qij)-(average(qij))**2.0)
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

def GetCoupling2(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, phasex, phasey, bd):


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
	for i in range(0,len(dbpms)-1):
		bn1=upper(dbpms[i][1])
		bn2=upper(dbpms[i+1][1])

		delx= phasex[bn1][0] - 0.25
		dely= phasey[bn1][0] - 0.25
		
		f1001ij=[]
		q1001ij=[]
		f1010ij=[]
		q1010ij=[]
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
				q1=(phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
				q2=(-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0
			elif bd==-1:
				q1=(phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
				q2=-(-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0
			#print q1,q2
			q1=(0.5-q1)%1.0 # This sign change in the real part is to comply with MAD output
			q2=(0.5-q2)%1.0
				

			if abs(q1-q2)<0.25: 
				q1001ij.append((q1+q2)/2.0)
			elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
				q1001ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
			else:
				badbpm=1
				q1001ij.append(q1)

			if bd==1:
				q1=(phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
				q2=(phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0
			if bd==-1:
				q1=(phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
				q2=-(phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0
			#print q1,q2
			q1=(0.5-q1)%1.0 # This sign change in the real part is to comply with MAD output
			q2=(0.5-q2)%1.0

			if abs(q1-q2)<0.25: 
				q1010ij.append((q1+q2)/2.0)
			elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
				q1010ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
			else:
				badbpm=1
				q1010ij.append(q1)

		if badbpm==0:
			f1001ij=array(f1001ij)
			f1001i=average(f1001ij)
			f1001istd=sqrt(average(f1001ij*f1001ij)-(average(f1001ij))**2.0)
			f1010ij=array(f1010ij)
			f1010i=average(f1010ij)
			f1010istd=sqrt(average(f1010ij*f1010ij)-(average(f1010ij))**2.0)
			q1001ij=array(q1001ij)
			q1001i=average(q1001ij)
			q1001istd=sqrt(average(q1001ij*q1001ij)-(average(q1001ij))**2.0)
			q1010ij=array(q1010ij)
			q1010i=average(q1010ij)
			q1010istd=sqrt(average(q1010ij*q1010ij)-(average(q1010ij))**2.0)
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
		
	
	CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
	QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0	
	fwqw['Global']=[CG,QG]


	return [fwqw,dbpms]


#--------------
def GetOffMomentumLattice(MADTwiss, ListOfFiles, betalist):

	bpms=intersect(ListOfFiles)
        bpmsl=[]

	slope={}
	for i in range(0,len(bpms)):
		bn=upper(bpms[i][1])
		check=0
		slopei=0.0
		for j in range(1,len(betalist)):
			try:
				#slopei+=(betalist[j][bn][0]/betalist[0][bn][0]-1.0)/betalist[j]['DPP']
				slopei+=(betalist[j][bn]/betalist[0][bn]-1.0)/betalist[j]['DPP']
			except:
				check=1
		if check==0:
			slopei=slopei/(len(betalist)-1)
			slope[bn]=slopei
			bpmsl.append([bpms[i][0],bpms[i][1]])
			

	return [slope,bpmsl]


	

#---------------------------------


def PseudoDoublePlaneMonitors(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY):



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
		print i
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
	dbpms=bpmpair()

 	for i in range(0,len(dbpms)):
		wname=dbpms[i][1]
		pname=dbpms[i][2]
		ws=dbpms[i][0]
		for j in range(0,len(ListOfZeroDPPX)):
			jx=ListOfZeroDPPX[j]
			jy=ListOfZeroDPPY[j]
			if dbpms[i][3]==0:
				wampx=jx.AMPX[jx.indx[wname]]
				wampy=jy.AMPY[jy.indx[pname]]
				wamp01=jx.AMP01[jx.indx[wname]]
				wamp10=jy.AMP10[jy.indx[pname]]
				wtunex=jx.TUNEX[jx.indx[wname]]
				wtuney=jy.TUNEY[jy.indx[pname]]
				dphix=MADTwiss.MUX[MADTwiss.indx[upper(pname)]]-MADTwiss.MUX[MADTwiss.indx[upper(wname)]]
				dphiy=MADTwiss.MUY[MADTwiss.indx[upper(pname)]]-MADTwiss.MUY[MADTwiss.indx[upper(wname)]]
				wmux=jx.MUX[jx.indx[wname]]
				wmuy=(jy.MUY[jy.indx[pname]]-dphiy)%1.0
				if (wmuy > 0.5): wmuy=wmuy-1
				wphase01=jx.PHASE01[jx.indx[wname]]
				wphase10=(jy.PHASE10[jy.indx[pname]]-dphiy)%1.0
				if (wphase10 > 0.5): wphase10=wphase10-1
			elif dbpms[i][3]==1:
				wampx=jx.AMPX[jx.indx[pname]]
				wampy=jy.AMPY[jy.indx[wname]]
				wamp01=jx.AMP01[jx.indx[pname]]
				wamp10=jy.AMP10[jy.indx[wname]]
				wtunex=jx.TUNEX[jx.indx[pname]]
				wtuney=jy.TUNEY[jy.indx[wname]]
				dphix=MADTwiss.MUX[MADTwiss.indx[upper(pname)]]-MADTwiss.MUX[MADTwiss.indx[upper(wname)]]
				dphiy=MADTwiss.MUY[MADTwiss.indx[upper(pname)]]-MADTwiss.MUY[MADTwiss.indx[upper(wname)]]
				wmux=(jx.MUX[jx.indx[pname]]-dphix)%1.0
				if (wmux > 0.5): wmux=wmux-1
				wmuy=jy.MUY[jy.indx[wname]]
				wphase01=(jx.PHASE01[jx.indx[pname]]-dphix)%1.0
				wphase10=jy.PHASE10[jy.indx[wname]]
				if (wphase01 > 0.5): wphase01=wphase01-1
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

	



#######################################################
#                   Main part                         #
#######################################################


# Path to accelerator setting file
# This path shuold be changed to be compatible to your working system.
accpath="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/CoreFiles/"


#-- Find index of python command in the system call
#i=0
#for entry in sys.argv:
#	if '.py' in entry: 
#		indpy=i
#		break
#	i=i+1

#-- Reading sys.argv
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
                metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-d", "--dictionary",
                help="File with the BPM dictionary",
                metavar="DICT", default="0", dest="dict")
parser.add_option("-m", "--model",
                help="Twiss File",
                metavar="TwissFile", default="0", dest="Twiss")
parser.add_option("-f", "--files",
                help="Files from analysis, separated by comma",
                metavar="TwissFile", default="0", dest="files")
parser.add_option("-o", "--output",
                help="Output Path",
                metavar="OUT", default="./", dest="output")
parser.add_option("-c", "--cocut",
                help="Cut for closed orbit measurement [um]",
                metavar="COCUT", default=1000, dest="COcut")
parser.add_option("-n", "--nbcpl",
                help="Analysis option for couplng, 1 bpm or 2 bpms",
                metavar="NBCPL", default=2, dest="NBcpl")
parser.add_option("-t", "--tbtana",
                help="Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA",
                metavar="TBTANA", default="SUSSIX", dest="TBTana")



(options, args) = parser.parse_args()


listOfInputFiles=options.files.split(",")
file0=options.Twiss
outputpath=options.output
if options.dict=="0":
	BPMdictionary={}
else:
	execfile(options.dict)
	BPMdictionary=dictionary   # temporaryly since presently name is not BPMdictionary

#file0=sys.argv[indpy+2]
MADTwiss=twiss(file0, BPMdictionary) # MODEL from MAD


try:
	facc=accpath+str(options.ACCEL)+'/accelerator.dat'
	faccs=twiss(facc)
	BPMU=faccs.BPMUNIT
except:
	BPMU='um'

COcut= float(options.COcut)

if BPMU=='um' or BPMU=='um': COcut=COcut
elif BPMU=='mm' or BPMU=='mm': COcut=COcut/1.0e3
elif BPMU=='cm' or BPMU=='cm': COcut=COcut/1.0e4
elif BPMU=='m' or BPMU=='m': COcut=COcut/1.0e6

# For selecting the coupling measurement method
NBcpl= int(options.NBcpl)


# Beam direction
bd=1
if options.ACCEL=="LHCB2":
	bd=-1 # note that the x axis has the same direction to BPM data. Otherwise another treatment should be done.


if options.TBTana=="SUSSIX":
	Suffix1='_linx'
	Suffix2='_liny'
elif options.TBTana=='SVD':
	Suffix1='_svdx'
	Suffix2='_svdy'
elif options.TBTana=='HA':
	Suffix1='_hax'
	Suffix2='_hay'
	

fphasex=open(outputpath+'getphasex.out','w')
fphasey=open(outputpath+'getphasey.out','w')

fphasex.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fphasey.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fphasex.write('@ FILES %s "')
fphasey.write('@ FILES %s "')

fbetax=open(outputpath+'getbetax.out','w')
fbetay=open(outputpath+'getbetay.out','w')
fbetax.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fbetay.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fbetax.write('@ FILES %s "')
fbetay.write('@ FILES %s "')

fabetax=open(outputpath+'getampbetax.out','w')
fabetay=open(outputpath+'getampbetay.out','w')
fabetax.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fabetay.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fabetax.write('@ FILES %s "')
fabetay.write('@ FILES %s "')

fcox=open(outputpath+'getCOx.out','w')
fcoy=open(outputpath+'getCOy.out','w')
fcox.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fcoy.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fcox.write('@ FILES %s "')
fcoy.write('@ FILES %s "')

fNDx=open(outputpath+'getNDx.out','w')
fDx=open(outputpath+'getDx.out','w')
fDy=open(outputpath+'getDy.out','w')
fNDx.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fDx.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fDy.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fNDx.write('@ FILES %s "')
fDx.write('@ FILES %s "')
fDy.write('@ FILES %s "')

fcouple=open(outputpath+'getcouple.out','w')
fcouple.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fcouple.write('@ FILES %s "')


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
	file1=filein+Suffix1
	x1=twiss(file1)
	try:
		dppi=x1.DPP
	except:
		dppi=0.0

	if float(dppi)==0.0:
		ListOfZeroDPPX.append(twiss(file1))
		FileOfZeroDPPX.append(file1)
		fphasex.write(file1+' ')
		fbetax.write(file1+' ')
		fabetax.write(file1+' ')
		fcox.write(file1+' ')
		fNDx.write(file1+' ')
		fDx.write(file1+' ')
		fcouple.write(filein+' ')
	else:
		ListOfNonZeroDPPX.append(twiss(file1))
		FileOfNonZeroDPPX.append(file1)
		fNDx.write(file1+' ')
		fDx.write(file1+' ')

	try:
		file1=filein+Suffix2
		y1=twiss(file1)
		try:
			dppi=y1.DPP
		except:
			dppi=0.0
		if dppi==0.0:
			ListOfZeroDPPY.append(twiss(file1))
			FileOfZeroDPPY.append(file1)
			fphasey.write(file1+' ')
			fbetay.write(file1+' ')
			fabetay.write(file1+' ')
			fcoy.write(file1+' ')
			fDy.write(file1+' ')
		else:
			ListOfNonZeroDPPY.append(twiss(file1))
			FileOfNonZeroDPPY.append(file1)
			fDy.write(file1+' ')
	except:
		print 'Warning: There seems no '+str(file1)+' file in the specified directory.' 


fphasex.write('"'+'\n')
fphasey.write('"'+'\n')
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


#-------- START Phases
if wolinx!=1:
	plane='H'
	[phasex,Q1,MUX,bpmsx]=GetPhases(MADTwiss,ListOfZeroDPPX,plane,outputpath,bd)
	phasex_save=phasex
if woliny!=1:
	plane='V'
	[phasey,Q2,MUY,bpmsy]=GetPhases(MADTwiss,ListOfZeroDPPY,plane,outputpath,bd)
	phasey_save=phasey
	fphasey.write('@ Q1 %le '+str(Q1)+'\n')
	fphasey.write('@ MUX %le '+str(MUX)+'\n')
	fphasey.write('@ Q2 %le '+str(Q2)+'\n')
	fphasey.write('@ MUY %le '+str(MUY)+'\n')
	fphasey.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHYMDL MUYMDL\n')
	fphasey.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpmsy)):
		bn1=upper(bpmsy[i][1])
		bns1=bpmsy[i][0]
		phmdl=phasey[bn1][4]
		if i==len(bpmsy)-1:
			bn2=upper(bpmsy[0][1])
			bns2=bpmsy[0][0]
			
		else:
			bn2=upper(bpmsy[i][1])
			bns2=bpmsy[i][0]	
		fphasey.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phasey[bn1][0])+' '+str(phasey[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

fphasey.close()


if wolinx!=1:
	fphasex.write('@ Q1 %le '+str(Q1)+'\n')
	fphasex.write('@ MUX %le '+str(MUX)+'\n')

	try:
		fphasex.write('@ Q2 %le '+str(Q2)+'\n')
		fphasex.write('@ MUY %le '+str(MUY)+'\n')
	except:
		fphasex.write('@ Q2 %le '+'0.0'+'\n')
		fphasex.write('@ MUY %le '+'0.0'+'\n')
	fphasex.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHXMDL MUXMDL\n')
	fphasex.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpmsx)):
		bn1=upper(bpmsx[i][1])
		bns1=bpmsx[i][0]
		phmdl=phasex[bn1][4]
		if i==len(bpmsx)-1:
			bn2=upper(bpmsx[0][1])
			bns2=bpmsx[0][0]
		else:
			bn2=upper(bpmsx[i][1])
			bns2=bpmsx[i][0]	
		fphasex.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasex[bn1][0])+' '+str(phasex[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
	
fphasex.close()




#-------- START Beta

betaxlist=[]
betaylist=[]

if wolinx!=1:
	plane='H'
	betax={}
	alfax={}
	rmsbbx=0.
	[betax,rmsbbx,alfax,bpms]=BetaFromPhase(MADTwiss,ListOfZeroDPPX,phasex,plane)
	betax['DPP']=0
	#beta_save=betax
	betaxlist.append(betax)
	fbetax.write('@ Q1 %le '+str(Q1)+'\n')
	try:
		fbetax.write('@ Q2 %le '+str(Q2)+'\n')
	except:
		fbetax.write('@ Q2 %le '+'0.0'+'\n')
	fbetax.write('@ RMS-beta-beat %le '+str(rmsbbx)+'\n')
	fbetax.write('* NAME   POS    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
	fbetax.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fbetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1][0])+' '+str(betax[bn1][1])+' '+str(betax[bn1][2])+' '+str(alfax[bn1][0])+' '+str(alfax[bn1][1])+' '+str(alfax[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )

fbetax.close()

if woliny!=1:
	plane='V'
	betay={}
	alfay={}
	rmsbby=0.
	[betay,rmsbby,alfay,bpms]=BetaFromPhase(MADTwiss,ListOfZeroDPPY,phasey,plane)
	betay['DPP']=0
	betaylist.append(betay)
	fbetay.write('@ Q1 %le '+str(Q1)+'\n')
	fbetay.write('@ Q2 %le '+str(Q2)+'\n')
	fbetay.write('@ RMS-beta-beat %le '+str(rmsbby)+'\n')
	fbetay.write('* NAME   POS    COUNT  BETY   ERRBETY STDBETY ALFY   ERRALFY STDALFY BETYMDL ALFYMDL MUYMDL\n')
	fbetay.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fbetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(betay[bn1][2])+' '+str(alfay[bn1][0])+' '+str(alfay[bn1][1])+' '+str(alfay[bn1][2])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

fbetay.close()


#------- Start beta from amplitude


betaxalist=[]
betayalist=[]

if wolinx!=1:
	plane='H'
	betax={}
	rmsbbx=0.
	[betax,rmsbbx,bpms, r2jx, r2jxerr]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPX,plane)
	betax['DPP']=0
	beta2_save=betax
	betaxalist.append(betax)
	fabetax.write('@ Q1 %le '+str(Q1)+'\n')
        fabetax.write('@ r2Jx %le '+str(r2jx)+'\n')
        fabetax.write('@ r2Jxerr %le '+str(r2jxerr)+'\n')
	try:
		fabetax.write('@ Q2 %le '+str(Q2)+'\n')
	except:
		fabetax.write('@ Q2 %le '+'0.0'+'\n')
	fabetax.write('@ RMS-beta-beat %le '+str(rmsbbx)+'\n')
	fabetax.write('* NAME   POS    COUNT  BETX   BETXMDL MUXMDL\n')
	fabetax.write('$ %s     %le    %le    %le    %le     %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fabetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n')

fabetax.close()

if woliny!=1:
	plane='V'
	betay={}
	rmsbby=0.
	[betay,rmsbby,bpms, r2jy, rsjyerr]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPY,plane)
	betay['DPP']=0
	betayalist.append(betay)
	fabetay.write('@ Q1 %le '+str(Q1)+'\n')
	fabetay.write('@ Q2 %le '+str(Q2)+'\n')
	fabetay.write('@ RMS-beta-beat %le '+str(rmsbby)+'\n')
	fabetay.write('* NAME   POS    COUNT  BETY   BETYMDL MUYMDL\n')
	fabetay.write('$ %s     %le    %le    %le    %le     %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fabetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n')

fabetay.close()



#-------- START Orbit
ListOfCOX=[]
if wolinx!=1:
	
	[cox,bpms]=GetCO(MADTwiss, ListOfZeroDPPX)
	# The output file can be directly used for orbit correction with MADX
	fcox.write('@ TABLE %05s "ORBIT"\n')
	fcox.write('@ TYPE %05s "ORBIT"\n')
	fcox.write('@ SEQUENCE %05s "'+options.ACCEL+'"\n')
	fcox.write('@ Q1 %le '+str(Q1)+'\n')

	try:
		fcox.write('@ Q2 %le '+str(Q2)+'\n')
	except:
		fcox.write('@ Q2 %le '+'0.0'+'\n')
	fcox.write('* NAME   POS1   COUNT  X      STDX   XMDL   MUXMDL\n')
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
	fcoy.write('@ SEQUENCE %05s "'+options.ACCEL+'"\n')
	fcoy.write('@ Q1 %le '+str(Q1)+'\n')
	fcoy.write('@ Q2 %le '+str(Q2)+'\n')
	fcoy.write('* NAME   POS1   COUNT  Y      STDY   YMDL   MUYMDL\n')
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
		fcoDPP.write('@ MAD_FILE: %s "'+file0+'"'+'\n')
		fcoDPP.write('@ FILE %s "')
		fcoDPP.write(FileOfNonZeroDPPX[k]+' "'+'\n')
		fcoDPP.write('@ DPP %le '+str(j.DPP)+'\n')
		try:
			fcoDPP.write('@ Q1 %le '+str(Q1)+'\n')
		except:
			fcoDPP.write('@ Q1 %le '+'0.0'+'\n')
		try:
			fcoDPP.write('@ Q2 %le '+str(Q2)+'\n')
		except:
			fcoDPP.write('@ Q2 %le '+'0.0'+'\n')
		[codpp,bpms]=GetCO(MADTwiss, SingleFile)
		fcoDPP.write('* NAME   POS1   COUNT  X      STDX   XMDL   MUXMDL\n')
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
		fcoDPP.write('@ MAD_FILE: %s "'+file0+'"'+'\n')
		fcoDPP.write('@ FILE %s "')
		fcoDPP.write(FileOfNonZeroDPPY[k]+' "'+'\n')
		fcoDPP.write('@ DPP %le '+str(j.DPP)+'\n')
		try:
			fcoDPP.write('@ Q1 %le '+str(Q1)+'\n')
		except:
			fcoDPP.write('@ Q1 %le '+'0.0'+'\n')
		try:
			fcoDPP.write('@ Q2 %le '+str(Q2)+'\n')
		except:
			fcoDPP.write('@ Q2 %le '+'0.0'+'\n')
		[codpp,bpms]=GetCO(MADTwiss, SingleFile)
		fcoDPP.write('* NAME   POS1   COUNT  Y      STDY   YMDL   MUYMDL\n')
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
	fNDx.write('* NAME   POS    COUNT  NDX    STDNDX DX     DPX    NDXMDL DXMDL  DPXMDL MUXMDL\n')
	fNDx.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		ndmdl=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
		fNDx.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPX))+' '+str(nda[bn1][0])+' '+str(nda[bn1][1])+' '+str(Dx[bn1][0])+' '+str(DPX[bn1])+' '+str(ndmdl)+' '+str(MADTwiss.DX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.DPX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )




	[dxo,bpms]=DispersionfromOrbit(ListOfZeroDPPX,ListOfNonZeroDPPX,ListOfCOX,COcut)
	DPX=GetDPX(MADTwiss,dxo,bpms)
	fDx.write('@ Q1 %le '+str(Q1)+'\n')
	fDx.write('@ Q2 %le '+str(Q2)+'\n')
	fDx.write('* NAME   POS    COUNT  DX     STDDX  DPX    DXMDL  DPXMDL MUXMDL\n')
	fDx.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fDx.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPX))+' '+str(dxo[bn1][0])+' '+str(dxo[bn1][1])+' '+str(DPX[bn1])+' '+str(MADTwiss.DX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.DPX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )


fNDx.close()
fDx.close()



if woliny!=1 and woliny2!=1:
	[dyo,bpms]=DispersionfromOrbit(ListOfZeroDPPY,ListOfNonZeroDPPY,ListOfCOY,COcut)
	fDy.write('@ Q1 %le '+str(Q1)+'\n')
	fDy.write('@ Q2 %le '+str(Q2)+'\n')
	fDy.write('* NAME   POS    COUNT  DY     STDDY  MUYMDL\n')
	fDy.write('$ %s     %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fDy.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPY))+' '+str(dyo[bn1][0])+' '+str(dyo[bn1][1])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

fDy.close()

	

#-------- Phase and Beta for non-zero DPP
if wolinx2!=1:
	plane='H'
	k=0
	for j in ListOfNonZeroDPPX:
		dpop=j.DPP
		SingleFile=[]
		SingleFile.append(j)
		file1=outputpath+'getphasex_dpp_'+str(k+1)+'.out'
		fphDPP=open(file1,'w')
		fphDPP.write('@ MAD_FILE: %s "'+file0+'"'+'\n')
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
		[phasex,Q1DPP,MUX,bpms]=GetPhases(MADTwiss,SingleFile,plane,outputpath,bd)
		fphDPP.write('@ Q1DPP %le '+str(Q1DPP)+'\n')
		fphDPP.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHXMDL MUXMDL\n')
		fphDPP.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
		for i in range(0,len(bpms)):
			bn1=upper(bpms[i][1])
			bns1=bpms[i][0]
			if i==len(bpms)-1:
				bn2=upper(bpms[0][1])
				bns2=bpms[0][0]
			
			else:
				bn2=upper(bpms[i][1])
				bns2=bpms[i][0]
			try:
				phmdl=phasex_save[bn1][4]
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
		[betaxa,rmsbbx,bpms, r2Jx, r2Jxerr]=BetaFromAmplitude(MADTwiss,SingleFile,plane)
		betaxa['DPP']=dpop
		betaxalist.append(betaxa)
		file2=outputpath+'getbetax_dpp_'+str(k+1)+'.out'
		fbetaxDPP=open(file2,'w')
		fbetaxDPP.write('@ MAD_FILE: %s "'+file0+'"'+'\n')
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
		#fbetaxDPP.write('@ RMS-beta-beat %le '+str(rmsbbx)+'\n')
		fbetaxDPP.write('* NAME   POS    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
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
		dpop=j.DPP
		SingleFile=[]
		SingleFile.append(j)
		file1=outputpath+'getphasey_dpp_'+str(k+1)+'.out'
		fphDPP=open(file1,'w')
		fphDPP.write('@ MAD_FILE %s "'+file0+'"'+'\n')
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
		[phasey,Q2DPP,MUY,bpms]=GetPhases(MADTwiss,SingleFile,plane,outputpath,bd)
		fphDPP.write('@ Q2DPP %le '+str(Q2DPP)+'\n')
		fphDPP.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHYMDL MUYMDL\n')
		fphDPP.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
		for i in range(0,len(bpms)):
			bn1=upper(bpms[i][1])
			bns1=bpms[i][0]
			if i==len(bpms)-1:
				bn2=upper(bpms[0][1])
				bns2=bpms[0][0]
			
			else:
				bn2=upper(bpms[i][1])
				bns2=bpms[i][0]
			try:
				phmdl=phasey_save[bn1][4]
			except:
				phmdl=0.0
			#phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
			fphDPP.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(1)+' '+str(phasey[bn1][0])+' '+str(phasey[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
		fphDPP.close()

		betay={}
		alfay={}
		rmsbby=0.
		[betay,rmsbby,alfay,bpms]=BetaFromPhase(MADTwiss,SingleFile,phasey,plane)
		betay['DPP']=dpop
		betaylist.append(betay)
		betaya={}
		[betaya,rmsbby,bpms, r2Jy, r2Jyerr]=BetaFromAmplitude(MADTwiss,SingleFile,plane)
		betaya['DPP']=dpop
		betayalist.append(betaya)
		file2=outputpath+'getbetay_dpp_'+str(k+1)+'.out'
		fbetayDPP=open(file2,'w')
		fbetayDPP.write('@ MAD_FILE: %s "'+file0+'"'+'\n')
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
		#fbetayDPP.write('@ RMS-beta-beat %le '+str(rmsbbx)+'\n')
		fbetayDPP.write('* NAME   POS    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
		fbetayDPP.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
		for i in range(0,len(bpms)):
			bn1=upper(bpms[i][1])
			bns1=bpms[i][0]
			fbetayDPP.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(betay[bn1][2])+' '+str(alfay[bn1][0])+' '+str(alfay[bn1][1])+' '+str(alfay[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
		fbetayDPP.close()

		k+=1


#-------- Find db/b vs dP/P

if wolinx!=1 and wolinx2!=1:
	file1=outputpath+'getdpplatticex.out'
	fdppx=open(file1,'w')
	slopex={}
	bpms=[]
	ListOfFiles=ListOfZeroDPPX+ListOfNonZeroDPPX
	#[slopex,bpms]=GetOffMomentumLattice(MADTwiss, ListOfFiles, betaxlist)
	[slopex,bpms]=GetOffMomentumLattice(MADTwiss, ListOfFiles, betaxalist)
	fdppx.write('* NAME   POS    COUNT  SBETX\n')
	fdppx.write('$ %s     %le    %le    %le\n')
	for i in range(0,len(bpms)):
		bn=upper(bpms[i][1])
		bns=bpms[i][0]
		fdppx.write('"'+bn+'" '+str(bns)+' '+str(len(ListOfNonZeroDPPX))+' '+str(slopex[bn])+'\n')
	fdppx.close()

if woliny!=1 and woliny2!=1:
	file1=outputpath+'getdpplatticey.out'
	fdppy=open(file1,'w')
	slopey={}
	bpms=[]
	ListOfFiles=ListOfZeroDPPY+ListOfNonZeroDPPY
	#[slopey,bpms]=GetOffMomentumLattice(MADTwiss, ListOfFiles, betaylist)
	[slopey,bpms]=GetOffMomentumLattice(MADTwiss, ListOfFiles, betayalist)
	fdppy.write('* NAME   POS    COUNT  SBETY\n')
	fdppy.write('$ %s     %le    %le    %le\n')
	for i in range(0,len(bpms)):
		bn=upper(bpms[i][1])
		bns=bpms[i][0]
		fdppy.write('"'+bn+'" '+str(bns)+' '+str(len(ListOfNonZeroDPPY))+' '+str(slopey[bn])+'\n')
	fdppy.close()


#-------- START coupling. If there is no double plane BPMs, nothing to do and leave the program. Will be updated to have pseudo-double plane BPM

sys.exit()
# Construct pseudo-double plane BPMs
if options.ACCEL=="SPS" and wolinx!=1 and woliny!=1 :
          [PseudoListX,PseudoListY]=PseudoDoublePlaneMonitors(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY)


          


if wolinx!=1 and woliny!=1:
	try:
		MADTwiss.Cmatrix()
	except:
		0.0

	if options.ACCEL=="SPS":
		plane='H'
		[phasexp,Q1,MUX,bpmsx]=GetPhases(MADTwiss,PseudoListX,plane,outputpath,bd)
		plane='V'
		[phaseyp,Q2,MUY,bpmsy]=GetPhases(MADTwiss,PseudoListY,plane,outputpath,bd)
		[fwqw,bpms]=GetCoupling2(MADTwiss, PseudoListX, PseudoListY, Q1, Q2, phasexp, phaseyp, bd)
	elif NBcpl==1:
		[fwqw,bpms]=GetCoupling1(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2)
	elif NBcpl==2:
		[fwqw,bpms]=GetCoupling2(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, phasex, phasey, bd)
	else:
		print 'Number of monitors for coupling analysis (option -n) should be 1 or 2.'
		print 'Leaving the coupling analysis...'
		sys.exit()
	

	try:
		fcouple.write('@ CG %le '+str(fwqw['Global'][0])+'\n')
		fcouple.write('@ QG %le '+str(fwqw['Global'][1])+'\n')
		if NBcpl==1:
				fcouple.write('* NAME   POS    COUNT  F1001W FWSTD  Q1001W QWSTD MDLF1001R MDLF1001I\n')
				fcouple.write('$ %s     %le    %le    %le    %le    %le    %le   %le       %le\n')
		elif NBcpl==2:
				fcouple.write('* NAME   POS    COUNT  F1001W FWSTD1 F1001R F1001I F1010W FWSTD2 F1010R F1010I MDLF1001R MDLF1001I MDLF1010R MDLF1010I\n')
				fcouple.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le    %le    %le    %le       %le       %le       %le\n')


		for i in range(0,len(bpms)):
			bn1=upper(bpms[i][1])
			bns1=bpms[i][0]
			
			if NBcpl==1:
				try:
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(MADTwiss.f1001[MADTwiss.indx(bn1)].real)+' '+str(MADTwiss.f1001[MADTwiss.indx(bn1)].imag)+' '+str(MADTwiss.f1010[MADTwiss.indx(bn1)].real)+' '+str(MADTwiss.f1010[MADTwiss.indx(bn1)].imag)+'\n')
				except: # Output zero if the model does not have couping parameters
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' 0.0 0.0'+'\n')
				

			elif NBcpl==2:
				try:
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2))+' '+str(fwqw[bn1][0][3])+' '+str(fwqw[bn1][0][2].real)+' '+str(fwqw[bn1][0][2].imag)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].imag)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].imag)+'\n')
				except: # Output zero if the model does not have couping parameters
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2))+' '+str(fwqw[bn1][0][3])+' '+str(fwqw[bn1][0][2].real)+' '+str(fwqw[bn1][0][2].imag)+' 0.0 0.0 0.0 0.0'+'\n')
					
			
	except:
		0.0
		
