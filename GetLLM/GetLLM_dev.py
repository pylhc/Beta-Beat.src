## Python script to obtain Linear Lattice function and More -> GetLMM
## Version 1.1 (still under development)
## Version-up history:V1.0, 11/Feb/2008 by Masa. Aiba

## Usage python2.5 GetPhaseDisp25_V1.0.py <ouput path> <lin file1 w/o _linx,y> <lin file2 w/o _linx,y> ... # lin files as many as you like

## Some rules for variable name: Dictionary is used to contain the output of function
##                               Valiable containing 'm' is a value directly obtained from measurment data
##                               Valiable containing 'mdl' is a value related to model



from metaclass25 import *
from numpy import *
import sys, pickle
#import operator
from string import *


#######################################################
#                 Functions                           #
#######################################################

#------------
def intersect(ListOfFile): 
	'''Pure intersection of all bpm names in all files '''
	if len(ListOfFile)==0:
		print "Nothing to intersect!!!!"
		exit()
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

def GetPhases(ListOfZeroDPP,plane):


	commonbpms=intersect(ListOfZeroDPP)
	zdpp=len(ListOfZeroDPP)

	
	phase={} # Dictionary for the output containing [average phase, rms error]
	for i in range(1,len(commonbpms)):
		bn1=upper(commonbpms[i-1][1])
		bn2=upper(commonbpms[i][1])
		bns1=commonbpms[i-1][0]
		bns2=commonbpms[i][0]
		phi=[]
		for j in ListOfZeroDPP:
			# Phase has unit of 2pi
			if plane=='H':
				phm=(j.MUX[j.indx[bn2]]-j.MUX[j.indx[bn1]])
			elif plane=='V':
				phm=(j.MUY[j.indx[bn2]]-j.MUY[j.indx[bn1]])
			if phm<0: phm+=1
			phi.append(phm)
		phi=array(phi)
		phstd=sqrt(average(phi*phi)-(average(phi))**2.0)
		phi=average(phi)
		phase[bn1]=[phi,phstd]
		
	return phase


#-------- Beta from pahses

def BetaFromPhase(MADTwiss,ListOfZeroDPP,phase,plane):

	beta={}
	commonbpms=intersect(ListOfZeroDPP)

	delbeta=[]
	for i in range(0,len(commonbpms)-2):
		bn1=upper(commonbpms[i][1])
		bn2=upper(commonbpms[i+1][1])
		bn3=upper(commonbpms[i+2][1])
		ph2pi12=2.*pi*phase[bn1][0]
		ph2pi13=2.*pi*(phase[bn1][0]+phase[bn2][0])
		# Find the model transfer matrices
		if plane=='H':
			phmdl12=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
			phmdl13=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn3]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
			betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
			betmdl2=MADTwiss.BETX[MADTwiss.indx[bn2]]
			betmdl3=MADTwiss.BETX[MADTwiss.indx[bn3]]
			alpmdl1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
		elif plane=='V':
			phmdl12=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
			phmdl13=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn3]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
			betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
			betmdl2=MADTwiss.BETY[MADTwiss.indx[bn2]]
			betmdl3=MADTwiss.BETY[MADTwiss.indx[bn3]]
			alpmdl1=MADTwiss.ALFY[MADTwiss.indx[bn1]]
		M11=sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
		M12=sqrt(betmdl1*betmdl2)*sin(phmdl12)
		N11=sqrt(betmdl3/betmdl1)*(cos(phmdl13)+alpmdl1*sin(phmdl13))
		N12=sqrt(betmdl1*betmdl3)*sin(phmdl13)

		# Find beta from phases assuming model transfer matrix  
		denom=M11/M12-N11/N12
		numer=1/tan(ph2pi12)-1/tan(ph2pi13)
		beti=numer/denom
		##### Factor 1/2 for phase13(=phase12+phase23) error is a bit optimistic side.
		betstd=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2)
		betstd=betstd+(abs(phase[bn1][1])+abs(phase[bn2][1]))/2.*abs((1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2)
		betstd=betstd/abs(denom)
		beta[bn1]=(beti,betstd)
		delbeta.append(beti-betmdl1)
	delbeta=array(delbeta)
	rmsbb=sqrt(average(delbeta*delbeta))
	#print beta
	return [beta,rmsbb]



#--------------
def NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX):

	nda={} # Dictionary for the output containing [average Disp, rms error]
        ALL=ListOfZeroDPPX+ListOfNonZeroDPPX
	bpmsMODEL=MADTwiss.NAME
	commonbpmsALL=intersect(ALL)

	nzdpp=len(ListOfNonZeroDPPX) # How many non zero dpp
	zdpp=len(ListOfZeroDPPX)  # How many zero dpp
	if zdpp==0 or nzdpp ==0:
		print 'Error: No data for dp/p=0 or for dp/p!=0.'
		exit() #!?


	mydp=[]
	gf=[]
	for j in range(0,nzdpp):
		mydp.append(ListOfNonZeroDPPX[j].DPP)
		gf.append(0.0)
	mydp=array(mydp)
	wf=array(mydp)/sum(mydp)*len(mydp) #Weitghs for the average


	co0=[]
	amp=[]
	nd=[]
	ndmdl=[]
	for i in range(0,len(commonbpmsALL)):
		bn1=upper(commonbpmsALL[i][1])
		bns1=commonbpmsALL[i][0]
		ndmdli=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
		ndmdl.append(ndmdli)

		coi=0.0
		ampi=0.0
		for j in ListOfZeroDPPX:
			bn1=upper(commonbpmsALL[i][1])
			coi+=j.CO[j.indx[bn1]]
			ampi+=j.AMPX[j.indx[bn1]]
                coi=coi/zdpp
		ami=ampi/zdpp

		ndi=[]
		for j in range(0,nzdpp):
			bn1=upper(commonbpmsALL[i][1])
			bns1=commonbpmsALL[i][0]
			orbit=ListOfNonZeroDPPX[j].CO[ListOfNonZeroDPPX[j].indx[bn1]]-coi
			ndm=orbit/ampi
			gf[j]+=ndm
			ndi.append(ndm)
		nd.append(ndi)

	ndmdl=array(ndmdl)
	avex0=average(ndmdl)

	gf=array(gf)
	gf=gf/avex0/len(commonbpmsALL)


	nd=array(nd)
	for i in range(0,len(commonbpmsALL)):
		ndi=[]
		bn1=upper(commonbpmsALL[i][1])
		bns1=commonbpmsALL[i][0]
		for j in range(0,nzdpp):
			ndi.append(nd[i][j]/gf[j])
		ndi=array(ndi)
		ndstd=sqrt(average(ndi*ndi)-(average(ndi))**2.0)
		ndas=average(wf*ndi)
		nda[bn1]=[ndas,ndstd]
	return nda









#######################################################
#                   Main part                         #
#######################################################

#-- Find index of python command in the system call
i=0
for entry in sys.argv:
	if '.py' in entry: 
		indpy=i
		break
	i=i+1

#-- Reading sys.argv

outputpath=sys.argv[indpy+1]

fphasex=open(outputpath+'getphasex.out','w')
fphasey=open(outputpath+'getphasey.out','w')
fphasex.write('@ FILES %s ')
fphasey.write('@ FILES %s ')

fbetax=open(outputpath+'getbetax.out','w')
fbetay=open(outputpath+'getbetay.out','w')
fbetax.write('@ FILES %s ')
fbetay.write('@ FILES %s ')

fDx=open(outputpath+'getDx.out','w')
#fDy=open(outputpath+'getbetay.out','w')
fDx.write('@ FILES %s ')
#fDy.write('@ FILES %s ')


file0=sys.argv[indpy+2]
MADTwiss=twiss(file0) # MODEL from MAD
fphasex.write('@ MAD_FILE %s '+file0+' '+'\n')
fphasey.write('@ MAD_FILE %s '+file0+' '+'\n')
fbetax.write('@ MAD_FILE %s '+file0+' '+'\n')
fbetay.write('@ MAD_FILE %s '+file0+' '+'\n')
fDx.write('@ MAD_FILE %s '+file0+' '+'\n')


ListOfZeroDPPX=[]
ListOfNonZeroDPPX=[]
ListOfZeroDPPY=[]
ListOfNonZeroDPPY=[]
woliny=0  # Let's assume there is liny for the moment
for i in range(indpy+3,len(sys.argv)):
	file1=sys.argv[i]+'_linx'
	x1=twiss(file1)
	try:
		dppi=x1.DPP
	except:
		dppi=0.0

	if dppi==0.0:
		ListOfZeroDPPX.append(twiss(file1))
		fphasex.write(file1+' ')
		fbetax.write(file1+' ')
		fDx.write(file1+' ')
	else:
		ListOfNonZeroDPPX.append(twiss(file1))
		fDx.write(file1+' ')

	try:
		file1=sys.argv[i]+'_liny'
		y1=twiss(file1)
		fphasex.write(file1+' ')
		fbetay.write(file1+' ')
		try:
			dppi=y1.DPP
		except:
			dppi=0.0
		if dppi==0.0:
			ListOfZeroDPPY.append(twiss(file1))
			fphasey.write(file1+' ')
			fbetay.write(file1+' ')
		else:
			ListOfNonZeroDPPY.append(twiss(file1))
	except:
		print 'Warning: There seems no '+str(file1)+' file in the specified directory.' 


fphasex.write('\n')
fphasey.write('\n')
fbetax.write('\n')
fbetay.write('\n')
fDx.write('\n')

if len(ListOfZeroDPPY)==0 :
	woliny=1  #FLAG meaning there is no _liny file!



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
				exit()





#-------- START Phases

plane='H'
phasex=GetPhases(ListOfZeroDPPX,plane)
fphasex.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH\n')
fphasex.write('$ %s     %s     %le    %le    %le    %le    %le\n')
bpms=intersect(ListOfZeroDPPX)
for i in range(1,len(bpms)):
	bn1=upper(bpms[i-1][1])
	bn2=upper(bpms[i][1])
	bns1=bpms[i-1][0]
	bns2=bpms[i][0]
	fphasex.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPX))+' '+str(phasex[bn1][0])+' '+str(phasex[bn1][1])+'\n' )


if woliny!=1:
	plane='V'
	phasey=GetPhases(ListOfZeroDPPY,plane)
	fphasey.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH\n')
	fphasey.write('$ %s     %s     %le    %le    %le    %le    %le\n')
	bpms=intersect(ListOfZeroDPPY)
	for i in range(1,len(bpms)):
		bn1=upper(bpms[i-1][1])
		bn2=upper(bpms[i][1])
		bns1=bpms[i-1][0]
		bns2=bpms[i][0]
		fphasey.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(ListOfZeroDPPY))+' '+str(phasey[bn1][0])+' '+str(phasey[bn1][1])+'\n' )


#-------- START Beta

plane='H'
betax={}
rmsbbx=0.
[betax,rmsbbx]=BetaFromPhase(MADTwiss,ListOfZeroDPPX,phasex,plane)
fbetax.write('* NAME   POS    COUNT  BETX   STDBT  BETXMDL\n')
fbetax.write('$ %s     %le    %le    %le    %le    %le\n')
bpms=intersect(ListOfZeroDPPX)
for i in range(1,len(bpms)-2):
	bn1=upper(bpms[i-1][1])
	bn2=upper(bpms[i][1])
	bns1=bpms[i-1][0]
	bns2=bpms[i][0]
	fbetax.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPX))+' '+str(betax[bn1][0])+' '+str(betax[bn1][1])+' '+str(MADTwiss.BETX[MADTwiss.indx[bn1]])+'\n' )


if woliny!=1:
	plane='V'
	betay={}
	rmsbby=0.
	[betay,rmsbby]=BetaFromPhase(MADTwiss,ListOfZeroDPPY,phasey,plane)
	fbetay.write('* NAME   POS    COUNT  BETX   STDBT  BETXMDL\n')
	fbetay.write('$ %s     %le    %le    %le    %le    %le\n')
	bpms=intersect(ListOfZeroDPPY)
	for i in range(1,len(bpms)-2):
		bn1=upper(bpms[i-1][1])
		bn2=upper(bpms[i][1])
		bns1=bpms[i-1][0]
		bns2=bpms[i][0]
		fbetay.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfZeroDPPY))+' '+str(betay[bn1][0])+' '+str(betay[bn1][1])+' '+str(MADTwiss.BETY[MADTwiss.indx[bn1]])+'\n' )


#-------- START Dispersion

nda=NormDispX(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX)
fDx.write('* NAME   POS    COUNT  NDX    STDNDX NDXMDL\n')
fDx.write('$ %s     %le    %le    %le    %le    %le\n')
ALL=ListOfZeroDPPX+ListOfNonZeroDPPX
bpms=intersect(ALL)
for i in range(0,len(bpms)):
	bn1=upper(bpms[i][1])
	bns1=bpms[i][0]
	ndm=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
	fDx.write('"'+bn1+'" '+str(bns1)+' '+str(len(ListOfNonZeroDPPX))+' '+str(nda[bn1][0])+' '+str(nda[bn1][1])+' '+str(ndm)+'\n' )


