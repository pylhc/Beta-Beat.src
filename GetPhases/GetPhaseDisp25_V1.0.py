## Python script to obtain phase advances between BPMs and normalized dispersion
## Version 1.0 (still under development)
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



#--------------
def ComputeNormDisp(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, nda):


        x1=ListOfZeroDPPX+ListOfNonZeroDPPX
	bpmsx0=MADTwiss.NAME
	bpmsx1=intersect(x1)

	nzdpp=len(ListOfNonZeroDPPX) # How many non zero dpp
	zdpp=len(ListOfZeroDPPX)  # How many zero dpp
	if zdpp==0 or nzdpp ==0:
		print 'Error: No data for dp/p=0 or for dp/p!=0.'
		exit() #!?


	mydp=[]
	gf=[]
	for i in range(0,nzdpp):
		mydp.append(ListOfNonZeroDPPX[i].DPP)
		gf.append(0.0)
	mydp=array(mydp)
	wf=array(mydp)/sum(mydp) #Weitghs for the average


	co0=[]
	amp=[]
	dd=[]
	a=[]
	for i in range(0,len(bpmsx1)):
		bn1=upper(bpmsx1[i][1])
		bns1=bpmsx1[i][0]
		nd=MADTwiss.DX[MADTwiss.indx[bn1]]/sqrt(MADTwiss.BETX[MADTwiss.indx[bn1]])
		a.append(nd)

		coi=0.0
		ampi=0.0
		for j in range(0,zdpp):
			bn1=upper(bpmsx1[i][1])
			coi+=x1[j].CO[ListOfZeroDPPX[j].indx[bn1]]
			ampi+=x1[j].AMPX[ListOfZeroDPPX[j].indx[bn1]]
                coi=coi/zdpp
		ami=ampi/zdpp

		d=[]
		for j in range(0,nzdpp):
			bn1=upper(bpmsx1[i][1])
			bns1=bpmsx1[i][0]
			orbit=ListOfNonZeroDPPX[j].CO[ListOfNonZeroDPPX[j].indx[bn1]]-coi
			nd=orbit/ampi
			gf[j]+=nd
			d.append(nd)
		dd.append(d)

	a=array(a)
	avex0=average(a)

	gf=array(gf)
	gf=gf/avex0/len(bpmsx1)


	dd=array(dd)
	for i in range(0,len(bpmsx1)):
		ndai=0.0
		bn1=bpmsx1[i][1]
		bns1=bpmsx1[i][0]
		for j in range(0,nzdpp):
			ndai+=dd[i][j]*wf[j]/gf[j]
		nda[bn1]=ndai
	return nda




#-------- START Reading sys.argv
x=[];y=[]


#-- Find index of python command in the system call
i=0
for entry in sys.argv:
	if '.py' in entry: 
		indpy=i
		break
	i=i+1

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
for i in range(indpy+3,len(sys.argv)):
	file1=sys.argv[i]+'_linx'
	x1=twiss(file1)
	fx.write(file1+' ')
	try:
		dppi=x1.DPP
	except:
		dppi=0.0
	if dppi==0 or dppi==0.0:
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
		if dppi==0 or dppi==0.0:
			ListOfZeroDPPY.append(twiss(file1))
		else:
			ListOfNonZeroDPPY.append(twiss(file1))
	except:
		woliny=1




fx.write('\n')
fy.write('\n')

if len(ListOfZeroDPPY)==0:
	print 'Warning: There seems no LINY file in your system call.' 


fx.write('* NAME   NAME2  POS1   POS2  COUNT  PHASE  RMS    NDX\n')
fy.write('* NAME   NAME2  POS1   POS2  COUNT  PHASE  RMS\n')
fx.write('$ %s     %s     %le    %le   %le    %le    %le    %le\n')
fy.write('$ %s     %s     %le    %le   %le    %le    %le\n')



#-------- START Dispersion


nda={}
ComputeNormDisp(MADTwiss, ListOfZeroDPPX, ListOfNonZeroDPPX, nda)



#-------- START Phases

commonbpmsx=intersect(ListOfZeroDPPX)
#commonbpmsy=intersect(ListOfZeroDPPY)
commonbpmsy=[]


for i in range(1,len(commonbpmsx)):
	bn1=commonbpmsx[i-1][1]
	bn2=commonbpmsx[i][1]
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
	try:
		ndai=nda[bn1]
	except:
		ndai=-100.
	fx.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(a))+' '+str(ave)+' '+str( sqrt(abs(ave2-ave**2)))+' '+str(ndai)+'\n' )


fx.close()


#
# Now Vertical
#

if len(ListOfZeroDPPY)==0:
	exit()


for i in range(1,len(commonbpmsy)):
	bn1=commonbpmsy[i-1][1]
	bn2=commonbpmsy[i][1]
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
	fy.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(a))+' '+str(ave)+' '+str( sqrt(abs(ave2-ave**2)))+'\n' )


fy.close()
