import sys
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Python_Classes4MAD/')

from optparse import OptionParser
from metaclass import twiss
import matplotlib.pyplot as plt


parser = OptionParser()

parser.add_option("-f", "--file",
		 help="Name of measured ORM data file",
		 metavar="FILE", default="measuredORM.dat",dest="file")

parser.add_option("-p", "--path",
		 help="Path to experimental data files",
		 metavar="PATH", default="./",dest="path")

parser.add_option("-n", "--nbpms",
		 help="Number of BPMs",
		 metavar="NBPMS", default=16,dest="nbpms")


(options, args) = parser.parse_args()


datafilename = options.path+options.file
print 'ORM file = ', datafilename
ORMmeas=twiss(datafilename)

nbpms = options.nbpms
ncorrs = len(ORMmeas.s) / nbpms

print(' ncorrs = ', ncorrs)

plt.figure(1)


for i in range(ncorrs):
    plt.subplot(211)  
    startidx = i*nbpms
    stopidx = (i+1)*nbpms - 1
    plt.plot(ORMmeas.s[startidx:stopidx], ORMmeas.MX[startidx:stopidx])
    
    plt.subplot(212)     
    plt.plot(ORMmeas.s[startidx:stopidx], ORMmeas.MY[startidx:stopidx])
    
    

plt.show()
