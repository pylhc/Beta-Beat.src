

#!/usr/bin/env pythonafs


# Just to make sure that the path to the libraires is defined 
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy

import pickle
from Numeric import *
#from numpy.oldnumeric.linear_algebra import generalized_inverse
from os import system
from metaclass import twiss
import random,re,sys
#from AllLists_couple import *
from LinearAlgebra import *
from optparse import OptionParser
from GenMatrix_couple import *
from BCORR import *

########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel", 
		 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
#parser.add_option("-t", "--tech", 
		 #help="Which algorithm: SVD MICADO",
		 #metavar="TECH", default="SVD",dest="TECH")
#parser.add_option("-n", "--ncorr", 
		 #help="Number of Correctors for MICADO",
		 #metavar="NCORR", default=5,dest="ncorr")
parser.add_option("-p", "--path",
		 help="Path to experimental files",
		 metavar="PATH", default="./",dest="path")
parser.add_option("-c", "--cut",
                  help="Singular value cut for the generalized inverse",
                  metavar="CUT", default=0.1 , dest="cut")
parser.add_option("-e", "--errorcut",
                  help="Maximum error allowed for the phase measurement",
                  metavar="ERRORCUT", default=0.013 , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase",
                  metavar="MODELCUT", default=0.02 , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to BetaBeat repository (default is the afs repository)",
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction",
                  metavar="MinStr", default=0.000001 , dest="MinStr")


(options, args) = parser.parse_args()



print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath = options.rpath
print "Path to Repository:", betapath
accelpath=betapath+'MODEL/'+options.ACCEL
print "Path to Accelerator model", accelpath

#if options.TECH=="MICADO":
    #ncorr=int(options.ncorr)
    #print "Number of Correctors Used:", options.ncorr
#else:
    #MinStr = options.MinStr
    #print "Minimum corrector strength", MinStr
# Do not need MICADO for coupling correction because the number of knobs is small
MinStr = options.MinStr
print "Minimum corrector strength", MinStr
#modelcut=float(options.modelcut)
#errorcut=float(options.errorcut)
cut= float(options.cut)
#print "Model, error and SVD cuts:", modelcut, errorcut, cut
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/Coupling/FullResponse_couple.Numeric','r'))
print "Loading ended"


couple=twiss(options.path+'/getcouple.out')

execfile(accelpath+'/Coupling/AllLists_couple.py')
exec('varslist=squadvars'+options.ACCEL+'()')

#if options.ACCEL=="LHCB1":
#	varslist=squadvarsb1()
#        print varslist
#        #intqx=64
#        #intqy=59              # corrected 06/May/08 by M.A.
#if options.ACCEL=="LHCB2":    # added 
#	varslist=squadvarsb2() # added 
#        #intqx=64              # added 
#        #intqy=59              # added 06/May/08 by M.A.

variables=varslist

couplelist=MakeList(couple, FullResponse['0'])

wei=[10,10,1,1] # Weights of f1001 and f1010

couple_inp=couple_input(varslist, couplelist, wei)

sensitivity_matrix=couple_inp.computeSensitivityMatrix(FullResponse)


[deltas, varslist ] = correctcouple(couple, couple_inp, cut=cut, app=0, path=options.path)



if "LHC" in options.ACCEL:   #.knob should always exist to be sent to LSA!
    system("cp "+options.path+"/changeparameters.tfs "+options.path+"/changeparameters.knob")

