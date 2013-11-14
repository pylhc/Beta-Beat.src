

#!/usr/bin/env pythonafs


# Just to make sure that the path to the libraires is defined
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy

import pickle
from metaclass import twiss
#try:
#	from Numeric import *
#	from LinearAlgebra import *
#except:
from numpy import *

#from numpy.oldnumeric.linear_algebra import generalized_inverse
from os import system
#from metaclass import twiss
import random,re,sys
#from AllLists_couple import *

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
                  help="weight on corrections (f1001.re, f1001.im, f1010.re, f1010.im, Dy)",
                  metavar="Dy", default="1,1,0,0,0", dest="Dy")
parser.add_option("-o", "--opt",
                  help="To specify the optics",
                  metavar="OPT", default="nominal.opt", dest="opt")
parser.add_option("-v", "--Variables",
                  help="variables split with ,",
                  metavar="var", default="MQSb1" , dest="var")



(options, args) = parser.parse_args()



print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath = options.rpath
print "Path to Repository:", betapath
accel=options.ACCEL

if "LHC" in accel:
	main="LHCB"
	accelpath=betapath+'/MODEL/'+main+'/fullresponse/'+options.ACCEL+'/'
else:
	main="SPS"
	accelpath=betapath+'/MODEL/'+main+'/fullresponse/'+options.ACCEL+'/'

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
FullResponse=pickle.load(open('./MODEL/LHCB1/1.5m.opt/Coupling/FullResponse_couple','r'))

print "Loading ended"



couple=twiss(options.path+'/getcouple.out')
weights=options.Dy.split(',')

if weights[4]=="1":
	dispy=twiss(options.path+'/getDy.out')
else:
	dispy=[]

execfile(accelpath+'/AllLists_couple.py')
print accelpath+'/AllLists_couple.py'
listvar=options.var.split(",")
print listvar
varslist=[]
for var in listvar:

    exec('variable='+var+'()')
    varslist=varslist+variable

variables=varslist



MADTwiss=FullResponse['0']
MADTwiss.Cmatrix()
mode='C'
couplelist=make_list(couple, MADTwiss, modelcutC, errorcutC, mode)
#print couplelist
mode='D'
dispylist=make_list(dispy, MADTwiss, modelcutD, errorcutD, mode)


#if int(options.Dy)==0:
	#wei=[1,1,0,0,0] # Weights of f1001 and f1010
#elif int(options.Dy)==1:
	#wei=[1,1,0,0,1] # Weights of f1001, f1010 and Dy
#else:
	#print "Specify the vertical dispersion correction option 1(turn-on) or 0(turn-off)"
wei=[int(weights[0]),int(weights[1]),int(weights[2]),int(weights[3]),int(weights[4])]




print "the weight option is "+str(wei)

optDy=options.Dy
print "entering couple input",len(couplelist)
couple_inp=CoupleInput(varslist, couplelist, dispylist, wei)
print "computing the sensitivity matrix"
sensitivity_matrix=couple_inp.computeSensitivityMatrix(FullResponse)
#print "sensitivity matrix", sensitivity_matrix

print "computing correct coupling "
[deltas, varslist ] = correctcouple(couple, dispy, couple_inp, cut=cut, app=0, path=options.path)

print deltas



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

    ##### for bumps
    if "bumps" in listvar:
	    print "passing trough bumps loop"
	    v=twiss(options.path+"/changeparameters_couple.tfs")
	    #print v.NAME
	    system('rm '+options.path+"/changeparameters_couple.tfs")
	    execfile(accelpath+'/Coupling/Bumps.py')
	    execfile(accelpath+'/Coupling/mydictionary.py')
	    filefile=open(options.path+"/changeparameters_couple.tfs","w")
	    filefile.write("* NAME  DELTA\n")
	    filefile.write("$ %s     %le\n")
	    for vcorr in corrs:
		          #print vcorr
			  #filefile.write(dictionary[vcorr]+" "+ str(corrs[vcorr])+"\n")
			  filefile.write(vcorr+" "+ str(corrs[vcorr])+"\n")

	    filefile.close()
    #####

    v=twiss(options.path+"/changeparameters_couple.tfs")
    mad=open(options.path+"/changeparameters_couple.madx", "w")
    names=v.NAME
    delta=v.DELTA

    for i in range(len(names)):

	    if "bumps" in listvar:

		    if cmp(delta[i],0)==1:
			    mad.write(names[i]+"->KICK:="+str(delta[i])+";\n");
		    else:
			    mad.write(names[i]+"->KICK:="+str(delta[i])+";\n");


	    else:

		    if cmp(delta[i],0)==1:
			    mad.write(names[i]+" = "+names[i]+" + "+str(delta[i])+";\n");
		    else:
			    mad.write(names[i]+" = "+names[i]+" "+str(delta[i])+";\n");


    mad.write("return;");

    mad.close()



print "Correcting couple Dy finished with weight "+ str(wei)
