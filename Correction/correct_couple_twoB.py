import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy

import pickle
from metaclass import twiss
try:
	from Numeric import *
	from LinearAlgebra import *
except:
	from numpy import *

#from numpy.oldnumeric.linear_algebra import generalized_inverse
from os import system
from metaclass import twiss
import random,re,sys
#from AllLists_couple import *

from optparse import OptionParser
from GenMatrix_couple_twoB import *



########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel",
		 help="What accelerator: LHCB1 LHCB2",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-b", "--accel2",
		 help="What accelerator for another beam: LHCB1 LHCB2",
		 metavar="ACCEL2", default="LHCB2",dest="ACCEL2")
parser.add_option("-t", "--tech",
		 help="Which algorithm: SVD MICADO",
		 metavar="TECH", default="SVD",dest="TECH")
parser.add_option("-p", "--path",
		 help="Path to experimental files",
		 metavar="PATH", default="./",dest="path")
parser.add_option("-q", "--path2",
		 help="Path to experimental files for the other beam",
		 metavar="PATH2", default="./",dest="path2")
parser.add_option("-c", "--cut",
                  help="Singular value cut for the generalized inverse",
                  metavar="CUT", default=0.1 , dest="cut")
parser.add_option("-e", "--errorcut",
                  help="Maximum error allowed for the phase measurement",
                  metavar="ERRORCUT", default="0.013,0.2" , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase",
                  metavar="MODELCUT", default="0.1,0.4" , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to BetaBeat repository (default is the afs repository)",
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
parser.add_option("-o", "--optics",
                  help="optics",
                  metavar="OPT", default="nominal.opt" , dest="OPT")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction (default is 1e-6)",
                  metavar="MinStr", default=0.000001 , dest="MinStr")
parser.add_option("-v", "--Variables",
                  help="variables split with ,",
                  metavar="var", default="MQTb1" , dest="var")
parser.add_option("-j", "--JustOneBeam",
                  help="0 only 1 beam corrrection, 1 for 2 beams correction (default is 0)",
                  metavar="JUSTONEBEAM", default='0' , dest="JustOneBeam")


parser.add_option("-w", "--weight",
                help="Weighting factor (f1001r, f1001i, f1010r, f1010i, dispersiony)",
                metavar="WGT", default="1,1,1,1,1", dest="WGT")


(options, args) = parser.parse_args()



print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath = options.rpath
print "Path to Repository:", betapath
accelpath=betapath+'/MODEL/'+options.ACCEL+'/'+options.OPT+'/'
accelpath2=betapath+'/MODEL/'+options.ACCEL2+'/'+options.OPT+'/'
print "Path to Accelerator model", accelpath
print "Selected algorithm:", options.TECH

MinStr = float(options.MinStr)
print "Minimum corrector strength", MinStr


j=int(options.JustOneBeam)
##### implementing cuts for

modelcutC=float(options.modelcut.split(",")[0])
errorcutC=float(options.errorcut.split(",")[0])
modelcutD=float(options.modelcut.split(",")[1])
errorcutD=float(options.errorcut.split(",")[1])
cut= float(options.cut)
print "Model, error and SVD cuts:", modelcutC, errorcutC, cut
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/Coupling/FullResponse_couple','r'))


if j==1:
    FullResponse2=pickle.load(open(accelpath2+'/Coupling/FullResponse_couple','r'))
print "Loading ended"


couple=twiss(options.path+'/getcouple.out')
if j==1:
    couple2=twiss(options.path2+'/getcouple.out')

wei=[]
for i in range(0,5):
    wei.append(int(options.WGT.split(",")[i]))
print "weight=",wei

if wei[4]=="1":
	dispy=twiss(options.path+'/getDy.out')
        if j==1:
            dispy2=twiss(options.path2+'/getDy.out')
else:
	dispy=[]
        dispy2=[]


execfile(accelpath+'/Coupling/AllLists_couple.py')

if j==1:
    execfile(accelpath2+'/Coupling/AllLists_couple.py')
# extra depdency to be able to handle to different magnets group
listvar=options.var.split(",")
varslist=[]
for var in listvar:
    exec('variable='+var+'()')
    varslist=varslist+variable

variables=varslist


MADTwiss=FullResponse['0']
MADTwiss.Cmatrix()
mode='C'
couplelist=MakeListCD(couple, MADTwiss, modelcutC, errorcutC, mode)
#print couplelist
mode='D'
dispylist=MakeListCD(dispy, MADTwiss, modelcutD, errorcutD, mode)


if j==1:
    MADTwiss2=FullResponse2['0']
    MADTwiss2.Cmatrix()
    mode='C'
    couplelist2=MakeListCD(couple2, MADTwiss2, modelcutC, errorcutC, mode)
    mode='D'
    dispylist2=MakeListCD(dispy2, MADTwiss2, modelcutD, errorcutD, mode)

print "Input ready"


if j==0:
    couple_inp=couple_input(varslist, couplelist, dispy, wei)
    sensitivity_matrix=couple_inp.computeSensitivityMatrix(FullResponse)
    [deltas, varslist ] = correctcouple(couple, dispy, couple_inp, cut=cut, app=0, path=options.path)
if j==1:
    couple_inp2=couple_input2(varslist, couplelist, couplelist2, dispy, dispy2, wei)
    sensitivity_matrix=couple_inp2.computeSensitivityMatrix2(FullResponse,FullResponse2)
    [deltas, varslist ] = correctcouple2(couple, couple2, dispy, dispy2, couple_inp2, cut=cut, app=0, path=options.path)


print deltas


if "LHC" in options.ACCEL:   #.knob should always exist to be sent to LSA!
    system("cp "+options.path+"/changeparameters.tfs "+options.path+"/changeparameters.knob")

    # madx table
    b=twiss(options.path+"/changeparameters.tfs")
    mad=open(options.path+"/changeparameters.madx", "w")
    names=b.NAME
    delta=b.DELTA

    for i in range(len(names)):

        if cmp(delta[i],0)==1:
            mad.write(names[i]+" = "+names[i]+" + "+str(delta[i])+";\n");
        else:
            mad.write(names[i]+" = "+names[i]+" "+str(delta[i])+";\n");


    mad.write("return;");

    mad.close()

