

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
from GenMatrix_coupleDy import *
#from BCORR import *

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
                  help="Maximum error allowed for the coupling maesurement and Dy measurement",
                  metavar="ERRORCUT", default="0.1,0.1" , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase",
                  metavar="MODELCUT", default="0.1,0.1" , dest="modelcut")
parser.add_option("-r", "--rpath",
                  help="Path to BetaBeat repository (default is the afs repository)",
                  metavar="RPATH", default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/" , dest="rpath")
parser.add_option("-s", "--MinStr",
                  help="Minimum strength of correctors in SVD correction",
                  metavar="MinStr", default=0.000001 , dest="MinStr")
parser.add_option("-d", "--Dy",
                  help="Turn-on=1 or Turn-off=0 vertical dispersion correction",
                  metavar="Dy", default=0, dest="Dy")
parser.add_option("-o", "--opt",
                  help="To specify the optics",
                  metavar="OPT", default="nominal.opt", dest="opt")
parser.add_option("-v", "--Variables",
                  help="variables split with ,",
                  metavar="var", default="MQTb1" , dest="var")



(options, args) = parser.parse_args()



print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath = options.rpath
print "Path to Repository:", betapath
accelpath=betapath
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

modelcuts=options.modelcut.split(",")
modelcutC=float(modelcuts[0])
modelcutD=float(modelcuts[1])
errorcuts=options.errorcut.split(",")
errorcutC=float(errorcuts[0])
errorcutD=float(errorcuts[1])

#modelcut=float(options.modelcut)
#errorcut=float(options.errorcut)
cut= float(options.cut)
#print "Model, error and SVD cuts:", modelcut, errorcut, cut
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/Coupling/FullResponse_couple.Numeric','r'))

print "Loading ended"



couple=twiss(options.path+'/getcouple.out')
if int(options.Dy)==1:
	dispy=twiss(options.path+'/getDy.out')
else:
	dispy=[]

execfile(accelpath+'/Coupling/AllLists_couple.py')
listvar=options.var.split(",")
varslist=[]
for var in listvar:
    
    exec('variable='+var+'()')
    varslist=varslist+variable

print varslist


#if options.ACCEL=="LHCB1":
#	varslist=squadvarsb1()
#if options.ACCEL=="LHCB2":    
#	varslist=squadvarsb2()
#if options.ACCEL=="SPS":
#	#varslist=SPSbumps()
#	varslist=["v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11", "v12", "v13", "v14", "v15", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31", "v32", "v33", "v34", "v35", "v36", "v37", "v38", "v39", "v40", "v41", "v42", "v43", "v44", "v45", "v46", "v47", "v48", "v49", "v50", "v51", "v52", "v53", "v54", "v55", "v56", "v57", "v58", "v59", "v60", "v61","v62", "v63", "v64", "v65", "v66", "v67", "v68", "v69", "v70", "v71", "v72", "v73", "v74", "v75", "v76", "v77", "v78", "v79", "v80", "v81", "v82", "v83", "v84", "v85", "v86", "v87", "v88", "v89", "v90", "v91", "v92", "v93", "v94", "v95", "v96", "v97", "v98", "v99", "v100", "v101","v102", "v103", "v104", "v105", "v106", "v107"]


variables=varslist



MADTwiss=FullResponse['0']
MADTwiss.Cmatrix()
mode='C'
couplelist=MakeList(couple, MADTwiss, modelcutC, errorcutC, mode)
#print couplelist
mode='D'
dispylist=MakeList(dispy, MADTwiss, modelcutD, errorcutD, mode)


if int(options.Dy)==0:
	wei=[1,1,0,0,0] # Weights of f1001 and f1010
elif int(options.Dy)==1:
	wei=[1,1,0,0,1] # Weights of f1001, f1010 and Dy
else:
	print "Specify the vertical dispersion correction option 1(turn-on) or 0(turn-off)"




print "the weight option is "+str(wei)

optDy=options.Dy
print "entering couple input"
couple_inp=couple_input(varslist, couplelist, dispylist, wei)
print "computing the sensitivity matrix"
sensitivity_matrix=couple_inp.computeSensitivityMatrix(FullResponse)
#print "sensitivity matrix", sensitivity_matrix

print "computing correct coupling"
[deltas, varslist ] = correctcouple(couple, dispy, couple_inp, cut=cut, app=0, path=options.path)

#print deltas



#print deltas, varslist
print "handling data"
if options.ACCEL=="SPS":
	v=twiss(options.path+"/changeparameters_couple.tfs")
	print '\nFor SPS, vertical bump file is imported:'
	print accelpath+'/Coupling/VBumps.py'
	execfile(accelpath+'/Coupling/VBumps.py')    # LOADS corrs
	execfile(accelpath+'/Coupling/VBumpsYASP.py') # LOADS corrsYASP
        #Output for YASP...
	f=open(options.path+"/changeparameters_couple.yasp", "w")
        #Output for Knob...
	#print options.path
        g=open(options.path+"/changeparameters_couple.knob", "w")
	h=open(options.path+"/changeparameters_couple.madx","w")
        f.write("#PLANE V\n")
	f.write("#UNIT RAD\n")
	
        g.write("* NAME  DELTA \n")
        g.write("$ %s    %le   \n")
	
	plane = 'V'
	beam = '1'
	for vcorr in vcorrsYASP:
		print >>f, "#SETTING", vcorr,  vcorrsYASP[vcorr]
		
	for vcorr in vcorrs:
                print >>g, "K"+vcorr, vcorrs[vcorr]
		print >>h, vcorr,"->KICK:=",vcorrs[vcorr],";"



	f.close()
	g.close()
	h.write('return;')
	h.close()


if "LHC" in options.ACCEL:   #.knob should always exist to be sent to LSA!
    system("cp "+options.path+"/changeparameters_couple.tfs "+options.path+"/changeparameters_couple.knob")

  # madx table
    b=twiss(options.path+"/changeparameters_couple.tfs")
    mad=open(options.path+"/changeparameters_couple.madx", "w")
    names=b.NAME
    delta=b.DELTA

    for i in range(len(names)):

        if cmp(delta[i],0)==1:
            mad.write(names[i]+" = "+names[i]+" + "+str(delta[i])+";\n");
        else:
            mad.write(names[i]+" = "+names[i]+" "+str(delta[i])+";\n");


    mad.write("return;");

    mad.close()

print "Correcting couple Dy finished"
