
import sys
sys.path.append('/afs/cern.ch/user/r/rtomas/lintrack/Scientific-bad/')

from mycsv import *


#p="./"
#time=635


p=sys.argv[1]
time=float(sys.argv[2])


LSDA=csv(p+'LSDA.csv')
LSDB=csv(p+'LSDB.csv')
LSFA=csv(p+'LSFA.csv')
LSFB=csv(p+'LSFB.csv')
LSFC=csv(p+'LSFC.csv')

print "!Sexts at time ", time, ";"
print "KSDA = ",LSDA.interfunc(time), ";"
print "KSDB = ",LSDB.interfunc(time), ";"
print "KSFA = ",LSFA.interfunc(time), ";"
print "KSFB = ",LSFB.interfunc(time), ";"
print "KSFC = ",LSFC.interfunc(time), ";"
