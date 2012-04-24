
from metaclass import twiss
from string import *

twisss=twiss('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/RHICY//250GeV_1d0.opt/twiss.dat')

names=twisss.NAME

fil=open("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/RHICY//mydictionary.py","w")

bpmpair=open("BPMpair.py","w")

print >> fil,"dictionary={"

print >> bpmpair,"def bpmpair():"
bpmpair.write("    variables=[")

for name in range(len(names)):

    print >>fil,'    "'+names[name]+'":["'+names[name].replace(" ","").replace("_","-")+'","'+names[name].replace(" ","").replace("_","-")+'"],'



print >> fil,'"c":"c"}'

for name in range(len(names)):
    temp=name
    lo=lower(names[temp])
    if "h" in lo:
        print >>bpmpair,"["+str(twisss.S[twisss.indx[names[temp]]])+",'"+lo+"','"+lower(names[temp+1])+"', 0 ],"
        #print "before ", temp
        #temp=temp+1
        #print "after ", temp
    elif "v" in lo:
        x=0
    else:
        print >>bpmpair,"["+str(twisss.S[twisss.indx[names[temp]]])+",'"+lo+"','"+lo+"', 0 ],"
        

print >>bpmpair,"return variables"


fil.close
