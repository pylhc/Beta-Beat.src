from metaclass25 import *
from numpy import *
from string import *
import sys, pickle
#import operator

#------------
def intersect(xx): 
	'''Pure intersection of all bpm names in all files '''
	z=xx[0].NAME
	#print xx
	for b in xx:
		z=filter(lambda x: x in z   , b.NAME)
	#SORT by S
	result=[]
	x0=xx[0]
	for bpm in z:
		result.append((x0.S[x0.indx[bpm]], upper(bpm)))	
		
	result.sort()
	#print result
	return result


#-------- START
x1=[]

#---- Find index of python command in the system call
i=0
for entry in sys.argv:
	if '.py' in entry: 
		indpy=i
		break
	i=i+1

outputpath=sys.argv[indpy+1]

fx=open(outputpath+'test.out','w')

file0=sys.argv[indpy+2]
x0=twiss(file0) # MODEL from MAD
fx.write('@ MAD_FILE %s '+file0+' '+'\n')
fx.write('@ FILES %s ')



for i in range(indpy+3,len(sys.argv)):
	file1=sys.argv[i]
	x1.append(twiss(file1))
	fx.write(file1+' ')
	
fx.write('\n')
fx.write('* NAME   S      Norm_DX\n')
fx.write('$ %s     %le    %le\n')

bpmsx0=x0.NAME
bpmsx1=intersect(x1)

mydp=[]
for i in range(0,len(x1)):
	try:
		mydp.append(x1[i].DPP)
	except:
		mydp.append(0.0)

mydp=array(mydp)
nzdpp=len(filter(lambda x: x!=0, mydp)) # How many non zero dpp 
wf=array(mydp)/sum(mydp) #Weitghs for the average
zdpp=len(x1)-nzdpp  # How many zero dpp

if zdpp==0 or nzdpp ==0:
	print 'Error: No data for dp/p=0 or for dp/p!=0.'
	excit() #!?


a=[]
for i in range(0,len(bpmsx1)):
	bn1=upper(bpmsx1[i][1])
	bns1=bpmsx1[i][0]
        nd=x0.DX[x0.indx[bn1]]/sqrt(x0.BETX[x0.indx[bn1]])
	a.append(nd)
a=array(a)
avex0=average(a)
#print avex0/

co0=[]
amp=[]
for i in range(0,len(bpmsx1)):
	coi=0.0
	ampi=0.0
	for j in range(0,len(x1)):
		if mydp[j]==0.0:
			bn1=upper(bpmsx1[i][1])
			coi+=x1[j].CO[x1[j].indx[bn1]]
			ampi+=x1[j].AMPX[x1[j].indx[bn1]]
	co0.append(coi)
	amp.append(ampi)
co0=array(co0)/zdpp
amp=array(amp)/zdpp


dd=[]
gf=[]
for j in range(0,len(x1)):
	if mydp[j]!=0.0:
		d=[]
		for i in range(0,len(bpmsx1)):
			bn1=upper(bpmsx1[i][1])
			bns1=bpmsx1[i][0]
			orbit=x1[j].CO[x1[j].indx[bn1]]-co0[i]
			nd=orbit/amp[i]
			d.append(nd)
		d=array(d)
		gf.append(avex0/average(d))
		dd.append(d)
	else :
		d=[]
		for i in range(0,len(bpmsx1)):
			d.append(0.0)
		d=array(d)
		gf.append(0.0)
		dd.append(d)

dd=array(dd)


for i in range(0,len(d)):
	nda=0.0
	bn1=bpmsx1[i][1]
	bns1=bpmsx1[i][0]
	for j in range(0,len(dd)):
		if mydp[j]!=0.0 or mydp[j]!=0 :
			nda+=gf[j]*dd[j][i]*wf[j]
	fx.write('"'+bn1+'" '+str(bns1)+' '+str(nda)+'\n' )


fx.close()
