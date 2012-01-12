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


## Usage1 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./
## Usage2 >pythonafs ../GetLLM_V1.8.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./


## lin files as many as you like

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
		#print "len",len(z)
	#SORT by S
	result=[]
	x0=ListOfFile[0]
	for bpm in z:
		result.append((x0.S[x0.indx[bpm]], bpm))	
		
	result.sort()
	return result


#------------ Get phases

def phiLastAndLastButOne(phi,frac,tune):
	if frac > 0.0:
		phit=phi+tune
		if phit>1.0: phit=phit-1.0
	elif frac <= 0.0:
		phit=phi+(1.0+tune)
		if phit>1.0: phit=phit-1.0
	return phit



def GetPhases(MADTwiss,ListOfFiles,plane,outputpath,bd):

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
		fracx=0.1 # Otherwise, the fractional parts are assumed to be below 0.5
		fracy=0.1 # Tentatively set 0.1

	


	commonbpms=intersect(ListOfFiles)
	commonbpms=modelIntersect(commonbpms, MADTwiss)
	#zdpp=len(ListOfFiles)
	
	mu=0.0
        tunem=[]
	phase={} # Dictionary for the output containing [average phase, rms error]
	for i in range(1,len(commonbpms)+1): # To find the integer part of tune as well, the loop is up to the last monitor
		if i==len(commonbpms)-1:
			bn1=upper(commonbpms[i-1][1])
			bn2=upper(commonbpms[i][1])
			bn3=upper(commonbpms[0][1])
			bns1=commonbpms[i-1][0]
			bns2=commonbpms[i][0]
		elif i==len(commonbpms):
			bn1=upper(commonbpms[i-1][1])
			bn2=upper(commonbpms[0][1])
			bn3=upper(commonbpms[1][1])
			bns1=commonbpms[i-1][0]
			bns2=commonbpms[0][0]
		else :
			bn1=upper(commonbpms[i-1][1])
			bn2=upper(commonbpms[i][1])
			bn3=upper(commonbpms[i+1][1])
			bns1=commonbpms[i-1][0]
			bns2=commonbpms[i][0]	
		phi12=[]
		phi13=[]
		tunemi=[]
		for j in ListOfFiles:
			# Phase has units of 2pi
			if plane=='H':
				phm12=(j.MUX[j.indx[bn2]]-j.MUX[j.indx[bn1]])
				phm13=(j.MUX[j.indx[bn3]]-j.MUX[j.indx[bn1]])
				tunemi.append(j.TUNEX[j.indx[bn1]])
			elif plane=='V':
				phm12=(j.MUY[j.indx[bn2]]-j.MUY[j.indx[bn1]])
				phm13=(j.MUY[j.indx[bn3]]-j.MUY[j.indx[bn1]])
				tunemi.append(j.TUNEY[j.indx[bn1]])
			if phm12<0: phm12+=1
			if phm13<0: phm13+=1
			phi12.append(phm12)
			phi13.append(phm13)
		phi12=array(phi12)
		phi13=array(phi13)
		if bd==-1:
			phi12=1.0-phi12
			phi13=1.0-phi13
		if phi12>0.9 and i !=len(commonbpms):
			print 'Warning: there seems too large phase advance! '+bn1+' to '+bn2+' = '+str(phi12)
			print plane
		phstd12=sqrt(average(phi12*phi12)-(average(phi12))**2.0)
		phstd13=sqrt(average(phi13*phi13)-(average(phi13))**2.0)
		phi12=average(phi12)
		phi13=average(phi13)
		tunemi=array(tunemi)
		if i<len(commonbpms)-1 :
			tunem.append(average(tunemi))
		if i==len(commonbpms)-1:
			tunem=array(tunem)
			tune=average(tunem)
			if plane=='H':
				phi13=phiLastAndLastButOne(phi13,fracx,tune)
			else:
				phi13=phiLastAndLastButOne(phi13,fracy,tune)
		elif i==len(commonbpms):
			if plane=='H':
				phi12=phiLastAndLastButOne(phi12,fracx,tune)
			else:
				phi12=phiLastAndLastButOne(phi12,fracy,tune)
			if plane=='H':
				phi13=phiLastAndLastButOne(phi13,fracx,tune)
			else:
				phi13=phiLastAndLastButOne(phi13,fracy,tune)
		mu=mu+phi12
		phase[bn1]=[phi12,phstd12,phi13,phstd13]

	return [phase,tune,mu,commonbpms]


#-------- Beta from pahses

def BetaFromPhase(MADTwiss,ListOfZeroDPP,phase,plane):

	alfa={}
	beta={}
	commonbpms=intersect(ListOfZeroDPP)
	commonbpms=modelIntersect(commonbpms,MADTwiss)
	alfii=[]
        betii=[]
	delbeta=[]
	for i in range(0,len(commonbpms)-3):
		bn1=upper(commonbpms[i][1])
		bn2=upper(commonbpms[i+1][1])
		bn3=upper(commonbpms[i+2][1])
		ph2pi12=2.*pi*phase[bn1][0]
		ph2pi23=2.*pi*phase[bn2][0]
		ph2pi13=2.*pi*phase[bn1][2]
		# Find the model transfer matrices for beta1
		if plane=='H':
			phmdl12=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
			phmdl13=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn3]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
			phmdl23=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn3]]-MADTwiss.MUX[MADTwiss.indx[bn2]])
			betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
			betmdl2=MADTwiss.BETX[MADTwiss.indx[bn2]]
			betmdl3=MADTwiss.BETX[MADTwiss.indx[bn3]]
			alpmdl1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
			alpmdl2=MADTwiss.ALFX[MADTwiss.indx[bn2]]
			alpmdl3=MADTwiss.ALFX[MADTwiss.indx[bn3]]
		elif plane=='V':
			phmdl12=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
			phmdl13=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn3]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
			phmdl23=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn3]]-MADTwiss.MUY[MADTwiss.indx[bn2]])
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



	i=0
	bn1=upper(commonbpms[i][1])
	beti=betii[i][0]
	beterr=betii[i][1]
	betstd=0.0
	beta[bn1]=(beti,beterr,betstd)
	alfi=alfii[i][0]
	alferr=alfii[i][1]
	alfstd=0.0
	alfa[bn1]=(alfi,alferr,alfstd)
	if plane=='H':
		betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
	elif plane=='V':
		betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
	delbeta.append((beti-betmdl1)/betmdl1)


	i=1
	bn1=upper(commonbpms[i][1])
	beti=(betii[i][0]+betii[i-1][2])/2.
	beterr=sqrt(betii[i][1]**2.+betii[i-1][3]**2.)/sqrt(2.)
	betstd=sqrt((betii[i][0]**2.+betii[i-1][2]**2.)/2.-((betii[i][0]+betii[i-1][2])/2.)**2.)
	beta[bn1]=(beti,beterr,betstd)
	alfi=(alfii[i][0]+alfii[i-1][2])/2.
	alferr=sqrt(alfii[i][1]**2.+alfii[i-1][3]**2.)/sqrt(2.)
	alfstd=sqrt((alfii[i][0]**2.+alfii[i-1][2]**2.)/2.-((alfii[i][0]+alfii[i-1][2])/2.)**2.)
	alfa[bn1]=(alfi,alferr,alfstd)
	if plane=='H':
		betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
	elif plane=='V':
		betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
	delbeta.append((beti-betmdl1)/betmdl1)

	for i in range(2,len(commonbpms)-3):
		bn1=upper(commonbpms[i][1])
		beti=(betii[i][0]+betii[i-1][2]+betii[i-2][4])/3.
		# error?
		beterr=sqrt(betii[i][1]**2.+betii[i-1][3]**2.+betii[i-2][5]**2.)/sqrt(3.)
		betstd=sqrt((betii[i][0]**2.+betii[i-1][2]**2.+betii[i-2][4]**2.)/3.-((betii[i][0]+betii[i-1][2]+betii[i-2][4])/3.)**2.)
		beta[bn1]=(beti,beterr,betstd)
		alfi=(alfii[i][0]+alfii[i-1][2]+alfii[i-2][4])/3.
		alferr=sqrt(alfii[i][1]**2.+alfii[i-1][3]**2.+alfii[i-2][5]**2.)/sqrt(3.)
		alfstd=sqrt((alfii[i][0]**2.+alfii[i-1][2]**2.+alfii[i-2][4]**2.)/3.-((alfii[i][0]+alfii[i-1][2]+alfii[i-2][4])/3.)**2.)
		alfa[bn1]=(alfi,alferr,alfstd)
		if plane=='H':
			betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
		elif plane=='V':
			betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
		delbeta.append((beti-betmdl1)/betmdl1)
		

	delbeta=array(delbeta)
	rmsbb=sqrt(average(delbeta*delbeta))
	return [beta,rmsbb,alfa,commonbpms]

#--------------------------------

def BetaFromAmplitude(MADTwiss,ListOfZeroDPP,plane):

	beta={}

	commonbpms=intersect(ListOfZeroDPP)
	commonbpms=modelIntersect(commonbpms,MADTwiss)
	SumA=0.0
	Amp=[]
	for i in range(0,len(commonbpms)):
		bn1=upper(commonbpms[i][1])
		Ampi=0.0
		for j in ListOfZeroDPP:
			if plane=='H':
				Ampi+=j.AMPX[j.indx[bn1]]
			elif plane=='V':
				Ampi+=j.AMPY[j.indx[bn1]]
		Ampi=Ampi/len(ListOfZeroDPP)
		Amp.append(Ampi)
		if plane=='H':
			SumA+=Ampi**2/MADTwiss.BETX[MADTwiss.indx[bn1]]
		if plane=='V':
			SumA+=Ampi**2/MADTwiss.BETY[MADTwiss.indx[bn1]]
		
	Kick=SumA/len(commonbpms)

	
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
	return [beta,rmsbb,commonbpms]


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


#--------------
def NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX, COcut):

	nzdpp=len(ListOfNonZeroDPPX) # How many non zero dpp
	zdpp=len(ListOfZeroDPPX)  # How many zero dpp
	if zdpp==0 or nzdpp ==0:
		print 'Error: No data for dp/p=0 or for dp/p!=0.'
		dum0={}
		dum1=[]
		return [dum0, dum1]

		#sys.exit() #!?

	coac=ListOfCOX[0] # COX dictionary after cut bad BPMs
	coact={}
	for i in coac:
		if (coac[i][1] < COcut):
			coact[i]=coac[i]
	coac=coact


	#coact={}
	#for i in coac:
		#Dt=[]
		#for j in range(1,len(ListOfCOX)):
			#codpp=ListOfCOX[j]
			#Dj=(codpp[i][0]-coac[i][0])/ListOfNonZeroDPPX[j-1].DPP
			#Dt.append(Dj)
		#Dt=array(Dt)
		#Dt=average(Dt) # Tentative dispersion
		#coj=0.0
		#coj2=0.0
		#for j in range(1,len(ListOfCOX)):
			#codpp=ListOfCOX[j]
			#coj=coj+(codpp[i][0]-coac[i][0])-Dt*ListOfNonZeroDPPX[j-1].DPP
			#coj2=coj2 + ((codpp[i][0]-coac[i][0])-Dt*ListOfNonZeroDPPX[j-1].DPP)**2
		#coj=coj/nzdpp
		#corms=sqrt(coj2/nzdpp-coj**2)
		#if (corms < COcut):
			#coact[i]=coac[i]

	#coac=coact



	nda={} # Dictionary for the output containing [average Disp, rms error]

	ALL=ListOfZeroDPPX+ListOfNonZeroDPPX
	commonbpmsALL=intersect(ALL)
        commonbpmsALL=modelIntersect(commonbpmsALL, MADTwiss)
	
	mydp=[]
	gf=[]
	for j in ListOfNonZeroDPPX:
		mydp.append(j.DPP)
		gf.append(0.0)
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
			ampi=ampi/zdpp # Note that ampi is averaged with kick size weighting


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




	# Find normalized dispersion and its rms
	nd=array(nd)
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
			bpms.append([bns1,bn1]) 
		except:
			0
	return [nda,bpms]



#-----------

def DYfromOrbit(ListOfZeroDPPY,ListOfNonZeroDPPY,ListOfCOY,COcut):


	coac=ListOfCOY[0] # COX dictionary after cut bad BPMs
	coact={}
	for i in coac:
		if (coac[i][1] < COcut):
			coact[i]=coac[i]

	coac=coact

        ALL=ListOfZeroDPPY+ListOfNonZeroDPPY
	commonbpmsALL=intersect(ALL)

	nzdpp=len(ListOfNonZeroDPPY) # How many non zero dpp
	zdpp=len(ListOfZeroDPPY)  # How many zero dpp


	mydp=[]
	for j in ListOfNonZeroDPPY:
		mydp.append(j.DPP)
	mydp=array(mydp)
	wf=array(mydp)/sum(mydp)*len(mydp) #Weitghs for the average

	dyo={} # Dictionary for the output containing [average Disp, rms error]
	bpms=[]
	for i in range(0,len(commonbpmsALL)):
		bn1=upper(commonbpmsALL[i][1])
		bns1=commonbpmsALL[i][0]

		try:
			coi=coac[bn1]
			dyoi=[]
			for j in ListOfNonZeroDPPY:
				dyoi.append((j.CO[j.indx[bn1]]-coi[0])/j.DPP)
			dyoi=array(dyoi)
			dyostd=sqrt(average(dyoi*dyoi)-(average(dyoi))**2.0)
			dyos=average(wf*dyoi)
			dyo[bn1]=[dyos,dyostd]
			bpms.append([bns1,bn1]) 
		except:
			coi=0
	return [dyo,bpms]

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

	# possibel correction
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


#######################################################
#                   Main part                         #
#######################################################


# Path to accelerator setting file
# This path shuold be changed to be compatible to your wourking system.
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
                help="Turn-by-turn data analysis algorithm: SUSSIX SVD",
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


# Beam deirection
bd=1
if options.ACCEL=="LHCB2":
	bd=-1 # note that the x axis has the same direction to BPM data. Otherwise another treatment should be done.


if options.TBTana=="SUSSIX":
	Suffix1='_linx'
	Suffix2='_liny'
elif options.TBTana=='SVD':
	Suffix1='_svdx'
	Suffix2='_svdy'


#outputpath=sys.argv[indpy+1]

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

fDx=open(outputpath+'getDx.out','w')
fDy=open(outputpath+'getDy.out','w')
fDx.write('@ MAD_FILE %s "'+file0+'"'+'\n')
fDy.write('@ MAD_FILE %s "'+file0+'"'+'\n')
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

	if dppi==0.0:
		ListOfZeroDPPX.append(twiss(file1))
		FileOfZeroDPPX.append(file1)
		fphasex.write(file1+' ')
		fbetax.write(file1+' ')
		fabetax.write(file1+' ')
		fcox.write(file1+' ')
		fDx.write(file1+' ')
		fcouple.write(filein+' ')
	else:
		ListOfNonZeroDPPX.append(twiss(file1))
		FileOfNonZeroDPPX.append(file1)
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
fDx.write('"'+'\n')
fDy.write('"'+'\n')
fcouple.write('"'+'\n')




if len(ListOfZeroDPPY)==0 :
	woliny=1  #FLAG meaning there is no _liny file for zero DPPY!
if len(ListOfNonZeroDPPY)==0 :
	woliny2=1  #FLAG meaning there is no _liny file for non-zero DPPY!



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

plane='H'
[phasex,Q1,MUX,bpmsx]=GetPhases(MADTwiss,ListOfZeroDPPX,plane,outputpath,bd)

if woliny!=1:
	plane='V'
	[phasey,Q2,MUY,bpmsy]=GetPhases(MADTwiss,ListOfZeroDPPY,plane,outputpath,bd)
	fphasey.write('@ Q1 %le '+str(Q1)+'\n')
	fphasey.write('@ MUX %le '+str(MUX)+'\n')
	fphasey.write('@ Q2 %le '+str(Q2)+'\n')
	fphasey.write('@ MUY %le '+str(MUY)+'\n')
	fphasey.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHYMDL MUYMDL\n')
	fphasey.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
	#bpms=intersect(ListOfZeroDPPY)
	#bpms=modelIntersect(bpms, MADTwiss)
	for i in range(1,len(bpmsy)+1):
		if i==len(bpmsy):
			bn1=upper(bpmsy[i-1][1])
			bn2=upper(bpmsy[0][1])
			bns1=bpmsy[i-1][0]
			bns2=bpmsy[0][0]
			phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]+MADTwiss.Q2-MADTwiss.MUY[MADTwiss.indx[bn1]]
		else:
			bn1=upper(bpmsy[i-1][1])
			bn2=upper(bpmsy[i][1])
			bns1=bpmsy[i-1][0]
			bns2=bpmsy[i][0]	
			phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
		fphasey.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phasey[bn1][0])+' '+str(phasey[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

fphasey.close()


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
#bpms=intersect(ListOfZeroDPPX)
#bpms=modelIntersect(bpms, MADTwiss)
for i in range(1,len(bpmsx)+1):
	if i==len(bpmsx):
		bn1=upper(bpmsx[i-1][1])
		bn2=upper(bpmsx[0][1])
		bns1=bpmsx[i-1][0]
		bns2=bpmsx[0][0]
		phmdl=MADTwiss.MUX[MADTwiss.indx[bn2]]+MADTwiss.Q1-MADTwiss.MUX[MADTwiss.indx[bn1]]
	else:
		bn1=upper(bpmsx[i-1][1])
		bn2=upper(bpmsx[i][1])
		bns1=bpmsx[i-1][0]
		bns2=bpmsx[i][0]	
		phmdl=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
	fphasex.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasex[bn1][0])+' '+str(phasex[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
	
fphasex.close()



#-------- START Beta

plane='H'
betax={}
alfax={}
rmsbbx=0.
[betax,rmsbbx,alfax,bpms]=BetaFromPhase(MADTwiss,ListOfZeroDPPX,phasex,plane)
fbetax.write('@ Q1 %le '+str(Q1)+'\n')
try:
	fbetax.write('@ Q2 %le '+str(Q2)+'\n')
except:
	fbetax.write('@ Q2 %le '+'0.0'+'\n')
fbetax.write('@ RMS-beta-beat %le '+str(rmsbbx)+'\n')
fbetax.write('* NAME   POS    COUNT  BETX   ERRBETX STDBETX ALFX   ERRALFX STDALFX BETXMDL ALFXMDL MUXMDL\n')
fbetax.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
#bpms=intersect(ListOfZeroDPPX)
#bpms=modelIntersect(bpms, MADTwiss)
for i in range(1,len(bpms)-3):
	bn1=upper(bpms[i-1][1])
	bn2=upper(bpms[i][1])
	bns1=bpms[i-1][0]
	bns2=bpms[i][0]
	fbetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1][0])+' '+str(betax[bn1][1])+' '+str(betax[bn1][2])+' '+str(alfax[bn1][0])+' '+str(alfax[bn1][1])+' '+str(alfax[bn1][2])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )

fbetax.close()

if woliny!=1:
	plane='V'
	betay={}
	alfay={}
	rmsbby=0.
	[betay,rmsbby,alfay,bpms]=BetaFromPhase(MADTwiss,ListOfZeroDPPY,phasey,plane)
	fbetay.write('@ Q1 %le '+str(Q1)+'\n')
	fbetay.write('@ Q2 %le '+str(Q2)+'\n')
	fbetay.write('@ RMS-beta-beat %le '+str(rmsbby)+'\n')
	fbetay.write('* NAME   POS    COUNT  BETY   ERRBETY STDBETY ALFY   ERRALFY STDALFY BETYMDL ALFYMDL MUYMDL\n')
	fbetay.write('$ %s     %le    %le    %le    %le     %le     %le    %le     %le     %le     %le     %le\n')
	#bpms=intersect(ListOfZeroDPPY)
	#bpms=modelIntersect(bpms, MADTwiss)
	for i in range(1,len(bpms)-3):
		bn1=upper(bpms[i-1][1])
		bn2=upper(bpms[i][1])
		bns1=bpms[i-1][0]
		bns2=bpms[i][0]
		fbetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(betay[bn1][2])+' '+str(alfay[bn1][0])+' '+str(alfay[bn1][1])+' '+str(alfay[bn1][2])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.ALFY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

fbetay.close()


#------- Start beta from amplitude

plane='H'
betax={}
rmsbbx=0.
[betax,rmsbbx,bpms]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPX,plane)
fabetax.write('@ Q1 %le '+str(Q1)+'\n')
try:
	fabetax.write('@ Q2 %le '+str(Q2)+'\n')
except:
	fabetax.write('@ Q2 %le '+'0.0'+'\n')
fabetax.write('@ RMS-beta-beat %le '+str(rmsbbx)+'\n')
fabetax.write('* NAME   POS    COUNT  BETX   BETXMDL MUXMDL\n')
fabetax.write('$ %s     %le    %le    %le    %le     %le\n')
#bpms=intersect(ListOfZeroDPPX)
#bpms=modelIntersect(bpms, MADTwiss)
for i in range(1,len(bpms)-3):
	bn1=upper(bpms[i-1][1])
	bn2=upper(bpms[i][1])
	bns1=bpms[i-1][0]
	bns2=bpms[i][0]
	fabetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n')

fabetax.close()

if woliny!=1:
	plane='V'
	betay={}
	rmsbby=0.
	[betay,rmsbby,bpms]=BetaFromAmplitude(MADTwiss,ListOfZeroDPPY,plane)
	fabetay.write('@ Q1 %le '+str(Q1)+'\n')
	fabetay.write('@ Q2 %le '+str(Q2)+'\n')
	fabetay.write('@ RMS-beta-beat %le '+str(rmsbby)+'\n')
	fabetay.write('* NAME   POS    COUNT  BETY   BETYMDL MUYMDL\n')
	fabetay.write('$ %s     %le    %le    %le    %le     %le\n')
	#bpms=intersect(ListOfZeroDPPY)
	#bpms=modelIntersect(bpms, MADTwiss)
	for i in range(1,len(bpms)-3):
		bn1=upper(bpms[i-1][1])
		bn2=upper(bpms[i][1])
		bns1=bpms[i-1][0]
		bns2=bpms[i][0]
		fabetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n')

fabetay.close()



#-------- START Orbit

[cox,bpms]=GetCO(MADTwiss, ListOfZeroDPPX)

fcox.write('@ Q1 %le '+str(Q1)+'\n')

try:
	fcox.write('@ Q2 %le '+str(Q2)+'\n')
except:
	fcox.write('@ Q2 %le '+'0.0'+'\n')
fcox.write('* NAME   POS1   COUNT  COX    STDCOX COXMDL MUXMDL\n')
fcox.write('$ %s     %le    %le    %le    %le    %le    %le\n')
for i in range(0,len(bpms)):
	bn1=upper(bpms[i][1])
	bns1=bpms[i][0]
	fcox.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(cox[bn1][0])+' '+str(cox[bn1][1])+' '+str(MADTwiss.X[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )

fcox.close()

ListOfCOX=[]
ListOfCOX.append(cox)


if woliny!=1:
	[coy,bpms]=GetCO(MADTwiss, ListOfZeroDPPY)
	fcoy.write('@ Q1 %le '+str(Q1)+'\n')
	fcoy.write('@ Q2 %le '+str(Q2)+'\n')
	fcoy.write('* NAME   POS1   COUNT  COY    STDCOY COYMDL MUYMDL\n')
	fcoy.write('$ %s     %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fcoy.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(coy[bn1][0])+' '+str(coy[bn1][1])+' '+str(MADTwiss.Y[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )


	ListOfCOY=[]
	ListOfCOY.append(coy)
fcoy.close()



#-------- Orbit for non-zero DPP

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
	fcoDPP.write('@ Q1 %le '+str(Q1)+'\n')
	try:
		fcoDPP.write('@ Q2 %le '+str(Q2)+'\n')
	except:
		fcoDPP.write('@ Q2 %le '+'0.0'+'\n')
	[codpp,bpms]=GetCO(MADTwiss, SingleFile)
	fcoDPP.write('* NAME   POS1   COUNT  COX    STDCOX COYMDL MUYMDL\n')
	fcoDPP.write('$ %s     %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpms)-1):
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
		fcoDPP.write('@ Q1 %le '+str(Q1)+'\n')
		try:
			fcoDPP.write('@ Q2 %le '+str(Q2)+'\n')
		except:
			fcoDPP.write('@ Q2 %le '+'0.0'+'\n')
		[codpp,bpms]=GetCO(MADTwiss, SingleFile)
		fcoDPP.write('* NAME   POS1   COUNT  COX    STDCOX COYMDL MUYMDL\n')
		fcoDPP.write('$ %s     %le    %le    %le    %le    %le    %le\n')
		for i in range(0,len(bpms)-1):
			bn1=upper(bpms[i][1])
			bns1=bpms[i][0]
			fcoDPP.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(codpp[bn1][0])+' '+str(codpp[bn1][1])+' '+str(MADTwiss.Y[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
		fcoDPP.close()
		ListOfCOY.append(codpp)
		k+=1



#-------- START Dispersion

[nda,bpms]=NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, ListOfCOX,COcut)
fDx.write('@ Q1 %le '+str(Q1)+'\n')
try:
	fDx.write('@ Q2 %le '+str(Q2)+'\n')
except:
	fDx.write('@ Q2 %le '+'0.0'+'\n')
fDx.write('* NAME   POS    COUNT  NDX    STDNDX DX     NDXMDL DXMDL  MUXMDL\n')
fDx.write('$ %s     %le    %le    %le    %le    %le    %le    %le    %le\n')
ALL=ListOfZeroDPPX+ListOfNonZeroDPPX
for i in range(0,len(bpms)):
	bn1=upper(bpms[i][1])
	bns1=bpms[i][0]
	ndmdl=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
	try:
		dxi=nda[bn1][0]*sqrt(betax[bn1][0])
	except:
		dxi=0.0
	fDx.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPX))+' '+str(nda[bn1][0])+' '+str(nda[bn1][1])+' '+str(dxi)+' '+str(ndmdl)+' '+str(MADTwiss.DX[MADTwiss.indx[bn1]])+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )

fDx.close()



if woliny!=1 and woliny2!=1:
	[dyo,bpms]=DYfromOrbit(ListOfZeroDPPY,ListOfNonZeroDPPY,ListOfCOY,COcut)
	fDy.write('@ Q1 %le '+str(Q1)+'\n')
	fDy.write('@ Q2 %le '+str(Q2)+'\n')
	fDy.write('* NAME   POS    COUNT  DY     STDDY  MUYMDL\n')
	fDy.write('$ %s     %le    %le    %le    %le    %le\n')
	ALL=ListOfZeroDPPY+ListOfNonZeroDPPY
	for i in range(0,len(bpms)):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		fDy.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPY))+' '+str(dyo[bn1][0])+' '+str(dyo[bn1][1])+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )

fDy.close()


#-------- START coupling

if woliny!=1:
	try:
		MADTwiss.Cmatrix()
	except:
		0.0
	if NBcpl==1:
		[fwqw,bpms]=GetCoupling1(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2)
	elif NBcpl==2:
		[fwqw,bpms]=GetCoupling2(MADTwiss, ListOfZeroDPPX, ListOfZeroDPPY, Q1, Q2, phasex, phasey, bd)
	else:
		print 'Number of monitors for coupling analysis (option -n) should be 1 or 2.'
		print 'Leaving the coupling analysis...'
	

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
				except:
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' 0.0 0.0'+'\n')
				

			elif NBcpl==2:
				try:
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2))+' '+str(fwqw[bn1][0][3])+' '+str(fwqw[bn1][0][2].real)+' '+str(fwqw[bn1][0][2].imag)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1001[MADTwiss.indx[bn1]].imag)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].real)+' '+str(MADTwiss.f1010[MADTwiss.indx[bn1]].imag)+'\n')
				except:
					fcouple.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2))+' '+str(fwqw[bn1][0][1])+' '+str(fwqw[bn1][0][0].real)+' '+str(fwqw[bn1][0][0].imag)+' '+str(sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2))+' '+str(fwqw[bn1][0][3])+' '+str(fwqw[bn1][0][2].real)+' '+str(fwqw[bn1][0][2].imag)+' 0.0 0.0 0.0 0.0'+'\n')
					
			
	except:
		0.0
		
	

#-------- Phase for non-zero DPP

k=0
for j in ListOfNonZeroDPPX:
	SingleFile=[]
	SingleFile.append(j)
	file1=outputpath+'getphasex_dpp_'+str(k+1)+'.out'
	fphDPP=open(file1,'w')
	fphDPP.write('@ MAD_FILE: %s "'+file0+'"'+'\n')
	fphDPP.write('@ FILE %s "')
	fphDPP.write(FileOfNonZeroDPPX[k]+' ')
	fphDPP.write('"'+'\n')
	fphDPP.write('@ Q1 %le '+str(Q1)+'\n')
	try:
		fphDPP.write('@ Q2 %le '+str(Q2)+'\n')
	except:
		fphDPP.write('@ Q2 %le '+'0.0'+'\n')
	plane='H'
	[phase,Q1DPP,MUX,bpms]=GetPhases(MADTwiss,SingleFile,plane,outputpath,bd)
	fphDPP.write('@ Q1DPP %le '+str(Q1DPP)+'\n')
	fphDPP.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHXMDL MUXMDL\n')
	fphDPP.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
	for i in range(0,len(bpms)-1):
		bn1=upper(bpms[i][1])
		bns1=bpms[i][0]
		phmdl=MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]]
		fphDPP.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(1)+' '+str(phase[bn1][0])+' '+str(phase[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUX[MADTwiss.indx[bn1]])+'\n' )
	fphDPP.close()
	k+=1


if woliny2!=1:
	k=0
	for j in ListOfNonZeroDPPY:
		SingleFile=[]
		SingleFile.append(j)
		file1=outputpath+'getphasey_dpp_'+str(k+1)+'.out'
		fphDPP=open(file1,'w')
		fphDPP.write('@ MAD_FILE %s "'+file0+'"'+'\n')
		fphDPP.write('@ FILE %s "')
		fphDPP.write(FileOfNonZeroDPPY[k]+' ')
		fphDPP.write('"'+'\n')
		fphDPP.write('@ Q1 %le '+str(Q1)+'\n')
		try:
			fphDPP.write('@ Q2 %le '+str(Q2)+'\n')
		except:
			fphDPP.write('@ Q2 %le '+'0.0'+'\n')
		plane='V'
		[phase,Q2DPP,MUY,bpms]=GetPhases(MADTwiss,SingleFile,plane,outputpath,bd)
		fphDPP.write('@ Q2DPP %le '+str(Q2DPP)+'\n')
		fphDPP.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH  PHYMDL MUYMDL\n')
		fphDPP.write('$ %s     %s     %le    %le    %le    %le    %le    %le    %le\n')
		for i in range(0,len(bpms)-1):
			bn1=upper(bpms[i][1])
			bns1=bpms[i][0]
			phmdl=MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]]
			fphDPP.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(1)+' '+str(phase[bn1][0])+' '+str(phase[bn1][1])+' '+str(phmdl)+' '+str(MADTwiss.MUY[MADTwiss.indx[bn1]])+'\n' )
		fphDPP.close()
		k+=1

