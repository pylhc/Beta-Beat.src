from metaclass import *
from Numeric import *
import sys
import pickle
#import operator



####################################
def intersect(xx): 
####################################
	'''Pure intersection of all bpm names in all files '''
	z=xx[0].NAME
	for b in xx:
		z=filter(lambda x: x in z   , b.NAME)	
		print "len",len(z)
	#SORT by S
	result=[]
	x0=xx[0]
	for bpm in z:
		result.append((x0.S[x0.indx[bpm]], bpm))	
		
	result.sort()
	return result



#####################################
def intersection(xx, level=3):
####################################
	# DO NOT USE !!!!!!!
	'''Intersect the BPMs that appear at least <level> times 
	   Of course if not enough file to meet level, level=# of files  
	'''
	bpmCount={}
	bpmPos={}
	result=[]
	if len(xx)<level:
		level=len(xx)
	for b in xx:
		for i in range(b.NAME):
			nam=b.NAME[i]
			if nam not in bpmCount:
				bpmCount[nam]=0
			 	bpmPos[nam]=b.S[i]
			bpmCount[nam]=bpmCount[nam]+1
			if bpmCount[nam]==level:
				result.append((nam,bpmPos[nam]))		
	
	return sorted(result, key=operator.itemgetter(1))			


#################### START #######################

x=[]
y=[]


#Find index of python command in the system call
i=0
for entry in sys.argv:
	if '.py' in entry: 
		indpy=i
		break
	i=i+1

outputpath=sys.argv[indpy+1]

fx=open(outputpath+'getphasex.out','w')
fy=open(outputpath+'getphasey.out','w')
fx.write('@ FILES %s ')
fy.write('@ FILES %s ')

tunesx=[]
tunesy=[]

for files in sys.argv[indpy+2:]:
    #if 'sdds' in files:
        x.append(twiss(files+'_linx'))
        y.append(twiss(files+'_liny'))
	fx.write(files+' ')
        fy.write(files+' ')
	tunesx.append(x[-1].Q1)
	tunesy.append(y[-1].Q2)
fx.write('\n')
fy.write('\n')
fx.write('@ Q1 %le '+str(average(tunesx))+'\n')
fy.write('@ Q2 %le '+str(average(tunesy))+'\n')





fx.write('* NAME   NAME2  POS1   POS2  COUNT  PHASE  RMS\n')
fy.write('* NAME   NAME2  POS1   POS2  COUNT  PHASE  RMS\n')
fx.write('$ %s     %s     %le    %le   %le    %le    %le\n')
fy.write('$ %s     %s     %le    %le   %le    %le    %le\n')



commonbpmsx=intersect(x)
commonbpmsy=intersect(y)
#print commonbpmsx, commonbpmsy


for i in range(1,len(commonbpmsx)):
	bn1=commonbpmsx[i-1][1]
	bn2=commonbpmsx[i][1]
	bns1=commonbpmsx[i-1][0]
	bns2=commonbpmsx[i][0]
	a=[]
	for c in x:
            adv=c.MUX[c.indx[bn2]]-c.MUX[c.indx[bn1]]
            if adv<0: adv+=1
            a.append(adv)
	a=array(a)
	ave=average(a)
	ave2=average(a**2)
	fx.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(a))+' '+str(ave)+' '+str( sqrt(abs(ave2-ave**2)))+'\n' )




fx.close()


#
# Now Vertical
#

for i in range(1,len(commonbpmsy)):
	bn1=commonbpmsy[i-1][1]
	bn2=commonbpmsy[i][1]
	bns1=commonbpmsy[i-1][0]
	bns2=commonbpmsy[i][0]
	a=[]
	for c in y:
		adv=c.MUY[c.indx[bn2]]-c.MUY[c.indx[bn1]]
		if adv<0: adv+=1
		a.append(adv)
	a=array(a)
	ave=average(a)
	ave2=average(a**2)
	fy.write('"'+bn1+'" '+'"'+bn2+'" '+str(bns1)+' '+str(bns2)+' '+str(len(a))+' '+str(ave)+' '+str( sqrt(abs(ave2-ave**2)))+'\n' )


fy.close()

