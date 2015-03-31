import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD/")

#from metaclass import *
from metaclass25 import *


t=twiss('./twiss.C.dat')
t.Cmatrix()

f=open("C.madx","w")
print >>f, "cminusreal=", t.F1001R[0], ";"
print >>f, "cminusimag=", t.F1001I[0], ";"
print >>f, "cplusreal=", t.F1010R[0], ";"
print >>f, "cplusimag=", t.F1010I[0], ";"
f.close()



f=open("f.dat","w")
for i in range(len(t.F1001R)):
    print >>f ,t.NAME[i], t.S[i], t.F1001R[i], t.F1001I[i], t.F1010R[i],  t.F1010I[i]


f.close()
