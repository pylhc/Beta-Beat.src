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
from AllLists import *
from LinearAlgebra import *
from optparse import OptionParser
from GenMatrixRHIC import *


def mymatch(cases,x):
	res=0
	for case in cases:
		res=(res or re.match('.*'+case, x))
	return res

def RemoveMatch(li,cases):
	return filter(lambda x: not (mymatch(cases,x))   , li)







########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel", 
		 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-p", "--path",
		 help="",
		 metavar="PATH", default="./",dest="path")

(options, args) = parser.parse_args()

print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/"
print "Path to Repository:", betapath
accelpath=betapath+'MODEL/'+options.ACCEL
print "Path to Accelerator model", accelpath
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/FullResponse.Numeric','r'))
print "Loading ended"


x=twiss(options.path+'x.out')
y=twiss(options.path+'y.out')

if options.ACCEL=="LHCB1":
	varslist=quadvarsb1()
if options.ACCEL=="RHIC":
	allquads=pickle.load(open(accelpath+'/qCirc.pic','r'))
	varslist=pickle.load(open(accelpath+'/qCirc.pic','r'))
	#print varslist
	print "Initial number of quads:", len(varslist)
	varslist=RemoveMatch(varslist,['w456m','k4m','k56m','kd','kf'])
	print "Number of quads after filtering:", len(varslist)
	print varslist
variables=varslist
phasexlist=MakePairs(x, FullResponse['0'])
phasexlist=RemoveMatch(phasexlist, ['rbpm.bo6-bh3.1'])
phaseylist=MakePairs(y, FullResponse['0'])
phaseylist=RemoveMatch(phaseylist, ['rbpm.bo6-bh3.1'])
betaxlist=[]
betaylist=betaxlist
displist=[]
wei=[1,1,1,1,1,1] # Weights of phasex phasey betax betay disp and tunes

beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)

sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)

correctbeatEXP(x,y, beat_inp, cut=0.04, app=0,allvars=allquads)


#--- evaluate beta-beat for output
bbt=betabeatEXP(x,y, FullResponse['0'])
print "Initial beatings:"
print "rmsx  rmsy  peakx  peaky  rmsphix  rmsphiy  peakphix  peakphiy   rmsdx  peakdx"
print bbt


#------ Apply svd correction
#correctbeat(x, beat_inp, cut=0.04, app=0)
    
#----- compute twiss after correction
#system('madx < job.random.madx > scum  ')
#z=twiss('twiss.dat')
#bb1=betabeat(x,y);bb=betabeat(z,y)
    

    
    

    
