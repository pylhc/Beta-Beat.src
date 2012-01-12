

#!/usr/bin/env pythonafs

def mad2dev(mad_name):
    #name_map = {'MDH.11605':'logical.MDHD11832' , 'MDH.61804':'logical.MDHA6803'} # THIS LINE!!!! ASK ILYA!!
    mad_name=mad_name.replace('"','')
    #if mad_name in name_map:
    #        return name_map[mad_name]
    return 'logical.' + mad_name

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
from GenMatrix import *
from BCORR import *

########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel", 
		 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
		 metavar="ACCEL", default="LHCB1",dest="ACCEL")
parser.add_option("-t", "--tech", 
		 help="Which algorithm: SVD MICADO",
		 metavar="TECH", default="SVD",dest="TECH")
parser.add_option("-n", "--ncorr", 
		 help="Number of Correctors for MICADO",
		 metavar="NCORR", default=5,dest="ncorr")
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
                  help="Minimum strength of correctors in SVD correction (defaul is 1e-6)",
                  metavar="MinStr", default=0.000001 , dest="MinStr")
parser.add_option("-j", "--JustOneBeam",
                  help="0 Just quads from one beam are used, 1 all quads (defaul is 0)",
                  metavar="JUSTONEBEAM", default=0 , dest="JustOneBeam")


(options, args) = parser.parse_args()



print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath = options.rpath
print "Path to Repository:", betapath
accelpath=betapath+'MODEL/'+options.ACCEL
print "Path to Accelerator model", accelpath
print "Selected algorithm:", options.TECH

JustOneBeam=int(options.JustOneBeam)
print "Going to use quads from one beam only (0 yes, 1 no):", JustOneBeam


if options.TECH=="MICADO":
    ncorr=int(options.ncorr)
    print "Number of Correctors Used:", options.ncorr
else:
    MinStr = float(options.MinStr)
    print "Minimum corrector strength", MinStr
modelcut=float(options.modelcut)
errorcut=float(options.errorcut)
cut= float(options.cut)
print "Model, error and SVD cuts:", modelcut, errorcut, cut
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/FullResponse.Numeric','r'))
print "Loading ended"




x=twiss(options.path+'/getphasex.out')
y=twiss(options.path+'/getphasey.out')
try:
	dx=twiss(options.path+'/getDx.out')
except:
	print "WARNING: No good dispersion or inexistent file getDx"
	print "WARNING: Correction will not take into account NDx"
	dx=[]
#dy=twiss(options.path+'/getDy.out')    # Should we correct Dy at some point???
#bx=twiss(options.path+'/getbetax.out')  # let's ignore betas for the momment...
#by=twiss(options.path+'/getbetay.out')


intqx=0
intqy=0

if options.ACCEL=="LHCB1":
	varslist=quadvarsb1()
        if JustOneBeam==0:
            varslist=filter(lambda x: 'B1' in x.upper(),varslist )
        intqx=64
        intqy=59              # corrected 06/May/08 by M.A.
if options.ACCEL=="LHCB2":    # added 
	varslist=quadvarsb2() # added
        if JustOneBeam==0:
            varslist=filter(lambda x: 'B2' in x.upper(), varslist)
        intqx=64              # added 
        intqy=59              # added 06/May/08 by M.A.
if options.ACCEL=="RHIC":
	varslist=pickle.load(open(accelpath+'/qCirc.pic','r'))
if options.ACCEL=="SPS":
	varslistBAD=["b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20", "b21", "b22", "b23", "b24", "b25", "b26", "b27", "b28", "b29", "b30", "b31", "b32", "b33", "b34", "b35", "b36", "b37", "b38", "b39", "b40", "b41", "b42", "b43", "b44", "b45", "b46", "b47", "b48", "b49", "b50", "b51", "b52", "b53", "b54", "b55", "b56", "b57", "b58", "b59", "b60", "b61", "b62", "b63", "b64", "b65", "b66", "b67", "b68","b69", "b70", "b71", "b72", "b73", "b74", "b75", "b76", "b77", "b78", "b79", "b80", "b81", "b82", "b83", "b84", "b85", "b86", "b87", "b88", "b89", "b90", "b91", "b92", "b93", "b94", "b95", "b96", "b97", "b98", "b99", "b100", "b101", "b102", "b103", "b104", "b105", "b106", "b107", "b108", "b109"]
	#These kickers: MDH11605 MDHB61804
        #do not work to LSA so they are removed by removing b7 b9 b98 b100
        varslist=["b2", "b3", "b4", "b5", "b6", "b8", "b10", "b11", "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20", "b21", "b22", "b23", "b24", "b25", "b26", "b27", "b28", "b29", "b30", "b31", "b32", "b33", "b34", "b35", "b36", "b37", "b38", "b39", "b40", "b41", "b42", "b43", "b44", "b45", "b46", "b47", "b48", "b49", "b50", "b51", "b52", "b53", "b54", "b55", "b56", "b57", "b58", "b59", "b60", "b61", "b62", "b63", "b64", "b65", "b66", "b67", "b68","b69", "b70", "b71", "b72", "b73", "b74", "b75", "b76", "b77", "b78", "b79", "b80", "b81", "b82", "b83", "b84", "b85", "b86", "b87", "b88", "b89", "b90", "b91", "b92", "b93", "b94", "b95", "b96", "b97", "b99", "b101", "b102", "b103", "b104", "b105", "b106", "b107", "b108", "b109"]
        intqx=26
        intqy=26

# remember to add integer part of the tunes to exp data!!!
x.Q1=x.Q1+intqx
y.Q2=y.Q2+intqy

variables=varslist
phasexlist=MakePairs(x, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
phaseylist=MakePairs(y, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
betaxlist=[]
betaylist=betaxlist
displist=MakeList(dx, FullResponse['0'])

wei=[1,1,1,1,1,10] # Weights of phasex phasey betax betay disp and tunes

beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)

sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)

#f=open("svd","w")
#print >>f, singular_value_decomposition(sensitivity_matrix)[1]
#f.close()

if options.TECH=="SVD":
    [deltas, varslist ] = correctbeatEXP(x,y,dx, beat_inp, cut=cut, app=0, path=options.path)
    #print varslist
    if "LHC" in options.ACCEL :
        iteration=0
        # Let's remove too low useless correctors
        while (len(filter(lambda x: abs(x)< MinStr, deltas))>0):
            iteration=1+iteration
            il=len(varslist)
            varslist_t=[]
            for i in range(0,il):
              #print MinStr, deltas[i]
              if (abs(deltas[i]) > MinStr):
                varslist_t.append(varslist[i])
            varslist=varslist_t
            if len(varslist)==0:
                print "You want to correct with too high cut on the corrector strength"
                sys.exit()
            beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)
            sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
            [deltas, varslist ] = correctbeatEXP(x,y,dx, beat_inp, cut=cut, app=0, path=options.path)
            print "Initial correctors:", il, ". Current: ",len(varslist), ". Removed for being lower than:", MinStr, "Iteration:", iteration

    
if options.TECH=="MICADO":
    bNCorrNumeric(x,y,dx,beat_inp, cut=cut,ncorr=ncorr,app=0,path=options.path)
    
if options.ACCEL=="SPS":
	b=twiss(options.path+"/changeparameters.tfs")
	execfile(accelpath+'/Bumps.py')    # LOADS corrs
	execfile(accelpath+'/BumpsYASP.py') # LOADS corrsYASP
        #Output for YASP...
	f=open(options.path+"/changeparameters.yasp", "w")
        #Output for Knob...
        g=open(options.path+"/changeparameters.knob", "w")
        f.write("#PLANE H\n")
	f.write("#UNIT RAD\n")
	
        g.write("* NAME  DELTA \n")
        g.write("$ %s    %le   \n")
	plane = 'H'
	beam = '1'
	for corr in corrsYASP:
		print >>f, "#SETTING", corr,  corrsYASP[corr]
	for corr in corrs:
                print >>g, "K"+corr, corrs[corr]
	f.close()
        g.close()
if "LHC" in options.ACCEL:   #.knob should always exist to be sent to LSA!
    system("cp "+options.path+"/changeparameters.tfs "+options.path+"/changeparameters.knob")


#--- evaluate beta-beat for output
#bbt=betabeatEXP(x,y, FullResponse['0'])
#print "Initial beatings:"
#print "rmsx  rmsy  peakx  peaky  rmsphix  rmsphiy  peakphix  peakphiy   rmsdx  peakdx"
#print bbt


#------ Apply svd correction
#correctbeat(x, beat_inp, cut=0.04, app=0)
    
#----- compute twiss after correction
#system('madx < job.random.madx > scum  ')
#z=twiss('twiss.dat')
#bb1=betabeat(x,y);bb=betabeat(z,y)
    

    
    
    
