from metaclass import *
from optparse import OptionParser
import sys


parser = OptionParser()
parser.add_option("-l", "--label",
                help="Label: IP1 IP2 IP3 IP4...",
                metavar="LABEL", default="",dest="label")

parser.add_option("-p", "--path", 
                help="Path to model files and output",
                metavar="PATH", default="./", dest="path")

parser.add_option("-e", "--exp", 
                help="path to experimental files, only used if fast!=1",
                metavar="EXP", default="./", dest="exp")

parser.add_option("-s", "--start", 
                help="Start BPM",
                metavar="START", default="./", dest="start")

#parser.add_option("-w", "--w", # Path to Chromaticity functions
               # help="Path to  chromaticity functions, by default this is skiped",
                #metavar="wpath", default="0", dest="wpath")



(options, args) = parser.parse_args()



Model=twiss(options.path+"/twiss_"+options.label+".dat")





Model.Cmatrix()

names=Model.NAME

#
# Write initvals.dat with f terms
#
f=open(options.path+"/initvals.dat","w")
f1001=Model.f1001[Model.indx[names[0]]]
f1010=Model.f1010[Model.indx[names[0]]]
print >>f, "f1001r=",f1001.real,";"
print >>f, "f1001i=",f1001.imag,";"
print >>f, "f1010r=",f1010.real,";"
print >>f, "f1010i=",f1010.imag,";"

f.close()










