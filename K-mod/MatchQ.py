import sys, os
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')
from optparse import OptionParser
from metaclass25 import *


parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 ",
                metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-s", "--start",
                help="Start point s=0 of the Accel, default are predetermined",
                metavar="start", default="0", dest="start")
parser.add_option("-I", "--IR",
                help="IR number (1,2,5,8) for the beta* measurement, default 1, ",
                metavar="IR", default="1", dest="IR")
parser.add_option("-p", "--plane",
                help="horizontal or vertical (Q1 or Q2)",
                metavar="Q12", default="Q1", dest="Q12")
parser.add_option("-d", "--delta",
                help="DQ from Q1R and from Q1L separated by commas as measured",
                metavar="deltas", default="0,0", dest="deltas")
parser.add_option("-k", "--deltaK",
                help="DK of Q1s used in the K-mod measurement",
                metavar="deltak", default="1e-5", dest="deltak")
parser.add_option("-F", "--FAST",
                help="FAST==1 no madx matching, only formula",
                metavar="FAST", default="0", dest="FAST")
parser.add_option("-m", "--modifiers",
                help="modifiers MADX file to define beta*",
                metavar="", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/modifiers.madx", dest="modifiers")


(options, args) = parser.parse_args()

accel=options.ACCEL
start=options.start
B12=options.ACCEL[-2:]
IR=options.IR
Q12=options.Q12
deltar=options.deltas.split(',')[0]
deltal=options.deltas.split(',')[1]
deltak=options.deltak

if Q12=="Q1":
    XY="x"
elif Q12=="Q2":
    XY="y"


if options.FAST=="1":
    deltar=float(deltar)
    deltal=float(deltal)
    deltak=float(deltak)
    pi=3.1415926
    
    q1r=abs(deltar/(deltak*6.37)*4*pi);
    q1l=abs(deltal/(deltak*6.37)*4*pi);
    betastar=2*26.15**2/(q1r+q1l);
    betastar26_5=2*26.5**2/(q1r+q1l);
    print "beta"+XY+"*=", betastar26_5, betastar
    print "beta"+XY+" @ q1r=",q1r
    print "beta"+XY+" @ q1l=",q1l
    sys.exit()
    

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
     "Q2R5":["MQXB.A2R1","MQXB.B2R5"],
     "Q2L5":["MQXB.A2L1","MQXB.B2L5"],
     "Q2R8":["MQXB.A2R2","MQXB.B2R8"],
     "Q2L8":["MQXB.A2L2","MQXB.B2L8"]
     }


cpath="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/K-mod/"

path='./'
filename=path+'/var4mad.sh'
file4nad=open(filename,'w')

file4nad.write('sed -e \'s/%START/\''+str(start)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%ACCEL/\''+str(options.ACCEL)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%IR/\''+str(IR)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%B12/\''+str(B12)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%Q12/\''+str(Q12)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%deltar/\''+str(deltar)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%deltal/\''+str(deltal)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%deltak/\''+str(deltak)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%XY/\''+str(XY)+'\'/g\' \\\n')
file4nad.write('    -e \'s/%MODIFIERS/'+options.modifiers.replace('/','\/')+'/g\' \\\n')
file4nad.write('<'+cpath+'/job.Qmatch.mask > '+path+'/t.madx \n')
file4nad.close()

os.system("chmod +x "+path+'/var4mad.sh')
os.system(path+'/var4mad.sh')

os.system("madx < "+path+'/t.madx')

sys.exit()

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
