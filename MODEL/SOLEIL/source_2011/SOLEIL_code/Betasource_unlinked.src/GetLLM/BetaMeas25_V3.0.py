## Python script to obtain beta-function via phase advances
## Version 3.0 (Translated from Fortran version)
## 31/Jan/2008 by Masa. Aiba

## Usage python2.5 GetPhaseDisp25_V1.0.py <ouput path> <lin file1 w/o _lin> <lin file2 w/o _lin> ...


from metaclass25 import *
from numpy import *
import sys, pickle
#import operator
from string import *

#------------
def intersect(xx): 
	'''Pure intersection of all bpm names in all files '''
	if len(xx)==0:
		print "Nothing to intersect!!!!"
		exit()
	z=xx[0].NAME
	for b in xx:
		z=filter(lambda x: x in z   , b.NAME)	
		#print "len",len(z)
	#SORT by S
	result=[]
	x0=xx[0]
	for bpm in z:
		result.append((x0.S[x0.indx[bpm]], bpm))	
		
	result.sort()
	return result

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
			phasem12=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn2]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
			phasem13=2.*pi*(MADTwiss.MUX[MADTwiss.indx[bn3]]-MADTwiss.MUX[MADTwiss.indx[bn1]])
			betam1=MADTwiss.BETX[MADTwiss.indx[bn1]]
			betam2=MADTwiss.BETX[MADTwiss.indx[bn2]]
			betam3=MADTwiss.BETX[MADTwiss.indx[bn3]]
			alpham1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
		elif plane=='V':
			phasem12=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn2]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
			phasem13=2.*pi*(MADTwiss.MUY[MADTwiss.indx[bn3]]-MADTwiss.MUY[MADTwiss.indx[bn1]])
			betam1=MADTwiss.BETY[MADTwiss.indx[bn1]]
			betam2=MADTwiss.BETY[MADTwiss.indx[bn2]]
			betam3=MADTwiss.BETY[MADTwiss.indx[bn3]]
			alpham1=MADTwiss.ALFY[MADTwiss.indx[bn1]]
		M11=sqrt(betam2/betam1)*(cos(phasem12)+alpham1*sin(phasem12))
		M12=sqrt(betam1*betam2)*sin(phasem12)
		N11=sqrt(betam3/betam1)*(cos(phasem13)+alpham1*sin(phasem13))
		N12=sqrt(betam1*betam3)*sin(phasem13)

		# Find beta from phases assuming model transfer matrix  
		denom=M11/M12-N11/N12
		numer=1/tan(ph2pi12)-1/tan(ph2pi13)
		beti=numer/denom
		##### Factor 1/2 for phase13(=phase12+phase23) error is a bit optimistic side.
		betstd=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2)
		betstd=betstd+(abs(phase[bn1][1])+abs(phase[bn2][1]))/2.*abs((1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2)
		betstd=abs(betstd/denom)
		beta[bn1]=(beti,betstd)
		delbeta.append(beti-betam1)
	delbeta=array(delbeta)
	rmsbb=sqrt(average(delbeta*delbeta))
	#print beta
	return [beta,rmsbb]


	



#-------- START 
x=[]
y=[]


#-- Find index of python command in the system call
i=0
for entry in sys.argv:
	if '.py' in entry: 
		indpy=i
		break
	i=i+1

#Reading sys.argv

outputpath=sys.argv[indpy+1]

fx=open(outputpath+'getphasex_Dx.out','w')
fy=open(outputpath+'getphasey.out','w')
fx.write('@ FILES %s ');fy.write('@ FILES %s ')

file0=sys.argv[indpy+2]
MADTwiss=twiss(file0) # MODEL from MAD
fx.write('@ MAD_FILE %s '+file0+' '+'\n')
fx.write('@ FILES %s ')


ListOfZeroDPPX=[]
ListOfNonZeroDPPX=[]
ListOfZeroDPPY=[]
ListOfNonZeroDPPY=[]
woliny=0  # Let's assume there is liny for the moment
for i in range(indpy+3,len(sys.argv)):
	file1=sys.argv[i]+'_linx'
	x1=twiss(file1)
	fx.write(file1+' ')
	try:
		dppi=x1.DPP
	except:
		dppi=0.0
	if dppi==0.0:
		ListOfZeroDPPX.append(twiss(file1))
	else:
		ListOfNonZeroDPPX.append(twiss(file1))

	try:
		file1=sys.argv[i]+'_liny'
		y1=twiss(file1)
		fy.write(file1+' ')
		try:
			dppi=y1.DPP
		except:
			dppi=0.0
		if dppi==0.0:
			ListOfZeroDPPY.append(twiss(file1))
		else:
			ListOfNonZeroDPPY.append(twiss(file1))
	except:
		woliny=1  #FLAG meaning there is no _liny file!




fx.write('\n')
fy.write('\n')

if len(ListOfZeroDPPY)==0:
	print 'Warning: There seems to be no LINY file in the specified directory.' 


fx.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH\n')
fy.write('* NAME   NAME2  POS1   POS2   COUNT  PHASE  STDPH\n')
fx.write('$ %s     %s     %le    %le    %le    %le    %le\n')
fy.write('$ %s     %s     %le    %le    %le    %le    %le\n')



#-------- START Phases

commonbpmsx=intersect(ListOfZeroDPPX)
#if woliny==0:
commonbpmsy=intersect(ListOfZeroDPPY)

phasex={}
for i in range(1,len(commonbpmsx)):
	bn1=upper(commonbpmsx[i-1][1])
	bn2=upper(commonbpmsx[i][1])
	bns1=commonbpmsx[i-1][0]
	bns2=commonbpmsx[i][0]
	a=[]
	for c in ListOfZeroDPPX:
		adv=(c.MUX[c.indx[bn2]]-c.MUX[c.indx[bn1]])
		if adv<0: adv+=1
		a.append(adv)
	a=array(a)
	ave=average(a)
	ave2=average(a**2)
	phstd=sqrt(abs(ave2-ave**2))
	phasex[bn1]=[ave,phstd]
		
	fx.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(a))+' '+str(ave)+' '+str(phstd)+'\n' )


fx.close()


#
# Now Vertical
#

if len(ListOfZeroDPPY)==0:
	exit()

phasey={}
for i in range(1,len(commonbpmsy)):
	bn1=upper(commonbpmsy[i-1][1])
	bn2=upper(commonbpmsy[i][1])
	bns1=commonbpmsy[i-1][0]
	bns2=commonbpmsy[i][0]
	a=[]
	for c in ListOfZeroDPPY:
		adv=(c.MUY[c.indx[bn2]]-c.MUY[c.indx[bn1]])
		if adv<0: adv+=1
		a.append(adv)
	a=array(a)
	ave=average(a)
	ave2=average(a**2)
	phstd=sqrt(abs(ave2-ave**2))
	phasey[bn1]=[ave,phstd]

	fy.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(a))+' '+str(ave)+' '+str(phstd)+'\n' )


fy.close()


#-------- START Beta

plane='H'
betax={}
rmsbbx=0.
[betax,rmsbbx]=BetaFromPhase(MADTwiss,ListOfZeroDPPX,phasex,plane)

plane='V'
betay={}
rmsbby=0.
[betay,rmsbby]=BetaFromPhase(MADTwiss,ListOfZeroDPPY,phasey,plane)


print betax
print betay
