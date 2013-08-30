import sys, os
from optparse import OptionParser
from metaclass import *


parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC",
                metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-s", "--start",
                help="Start point s=0 of the Accel, default are predetermined",
                metavar="start", default="0", dest="start")
parser.add_option("-q", "--quads",
                help="Quadrupoles to apply DK=1e-5 separated by commas",
                metavar="quads", default="", dest="quads")
parser.add_option("-v", "--value",
                help="Value to change quadrupoles",
                metavar="value", default="1e-5", dest="change")

(options, args) = parser.parse_args()
accel=options.ACCEL
start=options.start
value=float(options.change)
if start=="0":
    if accel=="LCHB1":
        start="MSIA.EXIT.B1"
    elif accel=="LHCB2":
        start="MKI.A5R8.B2"


dic={"Q1R1":["MQXA.1R1"],
     "Q1L1":["MQXA.1L1"],
     "Q1R2":["MQXA.1R2"],
     "Q1L2":["MQXA.1L2"],
     "Q1R5":["MQXA.1R5"],
     "Q1L5":["MQXA.1L5"],
     "Q1R8":["MQXA.1R8"],
     "Q1L8":["MQXA.1L8"],
     "Q2R1":["MQXB.A2R1","MQXB.B2R1"],
     "Q2L1":["MQXB.A2L1","MQXB.B2L1"],
     "Q2R2":["MQXB.A2R2","MQXB.B2R2"],
     "Q2L2":["MQXB.A2L2","MQXB.B2L2"],
     "Q2R5":["MQXB.A2R5","MQXB.B2R5"],
     "Q2L5":["MQXB.A2L5","MQXB.B2L5"],
     "Q2R8":["MQXB.A2R8","MQXB.B2R8"],
     "Q2L8":["MQXB.A2L8","MQXB.B2L8"],
     "Q3R8":["MQXA.3R8"],
     "Q3L8":["MQXA.3L8"]
     }
try:
    quads=dic[options.quads]
except:
    print "Quad -in -q option not defined in dictionary"
    sys.exit()

print "Quads:", quads, len(quads)
fs=open('SelectQuad.madx', "w")
fc=open('ChangeQuad.madx', "w")
for i in range(len(quads)):

    print >>fs, "SELECT, FLAG=ERROR, RANGE=",quads[i],";"
    print >>fc, quads[i]+"->K1=",quads[i]+"->K1+ (",value,");"
fc.close()
fs.close()

cpath="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/"

path='./'
filename=path+'/var4mad.sh'
file4nad=open(filename,'w')

file4nad.write('sed -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%ACCEL/\''+str(options.ACCEL)+'\'/g\' \\\n')
file4nad.write('<'+cpath+'/job.orbits.mask > '+path+'/t.madx \n')
file4nad.close()

os.system("chmod +x "+path+'/var4mad.sh')
os.system(path+'/var4mad.sh')
os.system("madx < "+path+'/t.madx')


ref=twiss("orbit.ref","r")
dk=twiss("orbit.+dk","r")

f=open("diffOrbit."+options.quads+"."+options.ACCEL,"w")
print >>f, "@ Qxdiff %le ",dk.Q1-ref.Q1
print >>f, "@ Qydiff %le ",dk.Q2-ref.Q2
print >>f, "* NAME S X Y"
print >>f, "$ %s %le %le %le"
for i in range(len(ref.S)):
    print >>f,ref.NAME[i],ref.S[i],dk.X[i]-ref.X[i],dk.Y[i]-ref.Y[i]
f.close()
