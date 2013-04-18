
from YASPmetaclassLHC import *
import sys, gzip, os
from metaclass import *



dire='./'
output1="multiturn1.dat"
output2="multiturn2.dat"
fout1=open(output1, "w")
fout2=open(output2, "w")

files = os.listdir(dire)
files = filter(lambda x: 'Dataset' in x    , files)
turns = len(files)

print turns, "files in dir:", dire

turnFiles={}

ListFiles=[]
for el in files:
    turn=int(el.split('-')[1])
    #print el, turn
    turnFiles[turn]=YASPtwiss(el)
    if turn<90:
        ListFiles.append(turnFiles[turn])
    


lhcb1=twiss('/afs/cern.ch/user/r/rtomas/lintrack/Beta-Beat.src/MODEL/LHCB1/twiss.dat')
lhcb2=twiss('/afs/cern.ch/user/r/rtomas/lintrack/Beta-Beat.src/MODEL/LHCB2/twiss.dat')



fout1.write('#title\n')
fout2.write('#title\n')

turns=90

def intersectYASP(ListOfFile): 
	'''Pure intersection of all bpm names in all files '''
	if len(ListOfFile)==0:
		print "Nothing to intersect!!!!"
		sys.exit()
	h=ListOfFile[0].HNAME
	for b in ListOfFile:
		h=filter(lambda x: x in h   , b.HNAME)
        v=ListOfFile[0].VNAME
	for b in ListOfFile:
		v=filter(lambda x: x in v   , b.VNAME)
        
	
	return h, v

    

h , v = intersectYASP(ListFiles)

for el in h:
    if '.B1' in el:
        lhc=lhcb1
        f=fout1
    elif '.B2' in el:
        lhc=lhcb2
        f=fout2
    else:
        print "BPM ", el, "Not B1 Not B2 !!!!!!!!!"
        sys.exit()
    s=lhc.S[lhc.indx[el]]
    #print el, s
    f.write("0 "+el+" "+str(s)+" ")
    for i in range(1,turns+1):
        tF=turnFiles[i]
        f.write(str(tF.HCO[tF.Hindx[el]]*1000)+" ") #The *1000 is to convert to mum
    f.write("\n")



for el in v:
    if '.B1' in el:
        lhc=lhcb1
        f=fout1
    elif '.B2' in el:
        lhc=lhcb2
        f=fout2
    else:
        print "BPM ", el, "Not B1 Not B2 !!!!!!!!!"
        sys.exit()
    s=lhc.S[lhc.indx[el]]
    #print el, s
    f.write("1 "+el+" "+str(s)+" ")
    for i in range(1,turns+1):
        tF=turnFiles[i]
        f.write(str(tF.VCO[tF.Vindx[el]]*1000)+" ") #The *1000 is to convert to mum
    f.write("\n")


fout1.close()

fout2.close()


        
    









