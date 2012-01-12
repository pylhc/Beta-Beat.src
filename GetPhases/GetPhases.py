from metaclass import *
from Numeric import *
import sys
import pickle
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


for files in sys.argv[indpy+2:]:
    #if 'sdds' in files:
        x.append(twiss(files+'_linx'))
        y.append(twiss(files+'_liny'))
        fx.write(files+' ')
        fy.write(files+' ')
fx.write('\n')
fy.write('\n')



fx.write('* NAME   NAME2  POS1   POS2  COUNT  PHASE  RMS\n')
fy.write('* NAME   NAME2  POS1   POS2  COUNT  PHASE  RMS\n')
fx.write('$ %s     %s     %le    %le   %le    %le    %le\n')
fy.write('$ %s     %s     %le    %le   %le    %le    %le\n')



a={}
loc={}

for c in x:
    for i in range(1,len(c.NAME)):
        if c.SLABEL[i-1]==1 and c.SLABEL[i]==1:
            twobpms=c.NAME[i-1]+' '+c.NAME[i]
            twolocations=str(c.S[i-1])+' '+str(c.S[i])
            if twobpms not in a:
                a[twobpms]=[]
                loc[twobpms]=twolocations
            adv=c.MUX[i]-c.MUX[i-1]
            if adv<0: adv+=1
            a[twobpms].append(adv)


for (k,v) in a.iteritems():
    v=array(v)
    ave=average(v)
    ave2=average(v**2)
    k1,k2=k.split()[0], k.split()[1]
    #print "\""+k1+"\"","\""+k2+"\"", loc[k], len(v), ave, sqrt(ave2-ave**2)
    fx.write('"'+k1+'" '+'"'+k2+'" '+loc[k]+' '+str(len(v))+' '+str(ave)+' '+str( sqrt(abs(ave2-ave**2)))+'\n' )

fx.close()


#
# Now Vertical
#

a={}
loc={}

for c in y:
    for i in range(1,len(c.NAME)):
        if c.SLABEL[i-1]==1 and c.SLABEL[i]==1:
            twobpms=c.NAME[i-1]+' '+c.NAME[i]
            twolocations=str(c.S[i-1])+' '+str(c.S[i])
            if twobpms not in a:
                a[twobpms]=[]
                loc[twobpms]=twolocations
            adv=c.MUY[i]-c.MUY[i-1]
            if adv<0: adv+=1
            a[twobpms].append(adv)


for (k,v) in a.iteritems():
    v=array(v)
    ave=average(v)
    ave2=average(v**2)
    k1,k2=k.split()[0], k.split()[1]
    #print "\""+k1+"\"","\""+k2+"\"", loc[k], len(v), ave, sqrt(ave2-ave**2)
    fy.write('"'+k1+'" '+'"'+k2+'" '+loc[k]+' '+str(len(v))+' '+str(ave)+' '+str( sqrt(abs(ave2-ave**2)))+'\n' )


fy.close()




#pickle.dump(a,open('dict.baseline2','w'))
#pickle.dump(a,open('dict.bi8-tq4_bo7-tq5_bo11-tq4++-0.005','w'))
#pickle.dump(loc,open('dict.baseline2.loc','w'))
