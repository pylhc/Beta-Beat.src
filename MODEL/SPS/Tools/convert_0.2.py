try:
	from Numeric import *
except:
	from numpy import *
from string import split, replace
import sys
import os.path

f22=sys.argv[1]

filename=f22
f = open(filename,'r');x=f.readlines();f.close()
header=x[8]

f2=open(sys.argv[2],'r');x2=f2.readlines();f2.close()

scoord=[]
for kk in range(len(x2)):
    s=split(x2[kk]);scoord.append((s[0],float(s[1])))
scoord = dict(scoord)

def nchg(name):
    if '.H' in name: return replace(name,'.H','')
    if '.V' in name: return replace(name,'.V','')

indx=[];indy=[];name=[];data=[]
for i in range(len(x)):
    s = split(x[i])
    if len(s) >= 2:
        if s[0]=='#DATA-SIZE':
            totBpms = int(s[1]); nturns = int(s[2])
        if s[0]=='#MONITOR':
            name.append(s[1])
	
	if s[0]=='#ACQ-TIME':
	    timestamp=s[1]

	if s[0]=='#ACQ-DATE':
	    datestamp=s[1]

        if s[0] == '#TURN-POS': # 2nd array element is turn number
            data.append(s[2:])
#data=array(data)

if os.path.exists(f22+'.sdds'):
	f = open(f22+'.new','w')
else:
	f = open(f22+'.sdds.new','w')
f.write(header);indx=0;indy=0;indNone=0
f.write('#NTURNS '+str(nturns)+'\n')
f.write('#ACQ-TIME '+timestamp+'\n')
f.write('#ACQ-DATE '+datestamp+'\n')
for ff in range(len(name)):
    if '.H' in name[ff]:
        plane='0';indx=indx+1
    elif '.V' in name[ff]:
        plane='1';indy=indy+1
    else: indNone=indNone+1
    
    f.write(plane+' '+name[ff]+'     '+str(scoord[nchg(name[ff])])+' ')
    for gg in range(nturns):
        f.write(data[gg][ff]+'  ')
    f.write('\n')
    
f.close()
print '# of Turns', nturns
print 'Hor Bpms:',indx, 'Ver Bpms:', indy
if indNone != 0:
    print indNone, "BPMs found with unknown plane"
