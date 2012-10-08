from allmightyclass import *


data="/afs/cern.ch/user/g/gvanbavi/scratch0/NON_LINEAR/LHC/F2000/output/getphasex.out"

sbss=sbs(data)

sbss.getsegmentbpms('BPM.8L2.B1','BPMYB.5L2.B1')

print sbss.segmentbpms

sbss.filterbpms(sbss.segmentbpms)

print sbss.fbpms

