"""
See docstring of iterative_correction.py (vimaier)
"""


##
## correct.py with two beam correction

import sys
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src///Numeric-23_p2.3/lib/python2.3/site-packages/Numeric/')

import pickle
try:
	from Numeric import *
	from LinearAlgebra import *
except:
	from numpy import *
from os import system
from metaclass import twiss
import random,re,sys
from optparse import OptionParser
from GenMatrix_twoB import *
from BCORR import *



########### START ###############







def  correctpy(options, args):

    print "Selected accelerator:", options.ACCEL
    print "Path to measurements:", options.path
    betapath = options.rpath
    print "Path to Repository:", betapath
    accelpath=options.OPT.split(",")[0]+"/"
    #accelpath=options.OPT+"/"
    #accelpath2=options.OPT+"/"
    print "Path to Accelerator model", accelpath
    print "Selected algorithm:", options.TECH


    if options.TECH=="MICADO":
        ncorr=int(options.ncorr)
        print "Number of Correctors Used:", options.ncorr
    else:
        if options.ACCEL=="SPS":
            MinStr =0.0
        else:
            MinStr = float(options.MinStr)
        print "Minimum corrector strength", MinStr


    j=int(options.JustOneBeam)
##### implementing cuts for
    
    modelcut=float(options.modelcut.split(",")[0])
    errorcut=float(options.errorcut.split(",")[0])
    modelcutdx=float(options.modelcut.split(",")[1])
    errorcutdx=float(options.errorcut.split(",")[1])
    cut= float(options.cut)
    print "Model, error and SVD cuts:", modelcut, errorcut, cut
    print "Starting loading Full Response optics"
    FullResponse=pickle.load(open(accelpath+'/FullResponse','r'))
    if j==1:
    	accelpath2=options.OPT.split(",")[1]+'/'
        FullResponse2=pickle.load(open(accelpath2+'/FullResponse','r'))
    print "Loading ended"


    x=twiss(options.path+'/getphasex.out.copy')
    y=twiss(options.path+'/getphasey.out.copy')
    #xbet=twiss(options.path+'/getphasex.out.copy')
    #ybet=twiss(options.path+'/getphasey.out.copy')
    
    if j==1:
        x2=twiss(options.path2+'/getphasex.out.copy')
        y2=twiss(options.path2+'/getphasey.out.copy')

    try:
        if options.WGT.split(",")[4]=='1':
            dx=twiss(options.path+'/getNDx.out.copy') # changed by Glenn Vanbavinckhove (26/02/09)
        else:
            dx=[]
        if j==1 and options.WGT.split(",")[4]=='1':
            dx2=twiss(options.path2+'/getNDx.out.copy')
        else:
            dx2=[]
    except:
        print "WARNING: No good dispersion or inexistent file getDx"
        print "WARNING: Correction will not take into account NDx"
        dx=[]
        dx2=[]
        

    execfile(options.rpath+'/MODEL/LHCB/fullresponse/'+options.ACCEL+'/AllLists.py')
    if j==1:
        execfile(options.rpath+'/MODEL/LHCB/fullresponse/'+options.ACCEL2+'/AllLists.py')
        # extra depdency to be able to handle to different magnets group
    listvar=options.var.split(",")
    varslist=[]
    for var in listvar:
    
        #exec('variable='+var+'()')
        #varslist=varslist+variable
        varslist=varslist+eval(var+'()')
    

    intqx=int(FullResponse['0'].Q1)
    intqy=int(FullResponse['0'].Q2)
    if j==1:
        intqx2=int(FullResponse2['0'].Q1)
        intqy2=int(FullResponse2['0'].Q2)

    print "Integer part of tunes: ", intqx, intqy

    # remember to add integer part of the tunes to exp data!!!
    if x.Q1 > 0.0:
        x.Q1=x.Q1+intqx
    else:
        x.Q1=x.Q1+intqx+1.0
    if y.Q2 > 0.0:
        y.Q2=y.Q2+intqy
    else:
        y.Q2=y.Q2+intqy+1.0

    if j==1:
        if x2.Q1 > 0.0:
            x2.Q1=x2.Q1+intqx2
        else:
            x2.Q1=x2.Q1+intqx2+1.0
        if y2.Q2 > 0.0:
            y2.Q2=y2.Q2+intqy2
        else:
            y2.Q2=y2.Q2+intqy2+1.0


    if j==0:
        print "Experiment tunes: ", x.Q1, y.Q2
    if j==1:
        print "Experiment tunes: ", x.Q1, y.Q2, x2.Q1, y2.Q2


    variables=varslist
    phasexlist=MakePairs(x, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
    phaseylist=MakePairs(y, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
#   betaxlist=MakeBetaList(x, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
 #  betaylist=MakeBetaList(y, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
    betaxlist=[]
    betaylist=betaxlist
    displist=MakeList(dx, FullResponse['0'])
    if j==1:
        phasexlist2=MakePairs(x2, FullResponse2['0'], modelcut=modelcut, errorcut=errorcut)
        phaseylist2=MakePairs(y2, FullResponse2['0'], modelcut=modelcut, errorcut=errorcut)
        betaxlist2=[]
        betaylist2=betaxlist2
        displist2=MakeList(dx2, FullResponse2['0'])

    print "Input ready"

    wei=[]
    for i in range(0,6):
        wei.append(float(options.WGT.split(",")[i]))#wei=[1,1,1,1,1,10] # Weights of phasex phasey betax betay disp and tunes
        #wei.append(int(options.WGT.split(",")[i]))#wei=[1,1,1,1,1,10] # Weights of phasex phasey betax betay disp and tunes

    print "weight=",wei

    if j==0:
        beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)
        sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
    if j==1:
        beat_inp=beat_input2(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, phasexlist2, phaseylist2, betaxlist2, betaylist2, displist2, wei)
        sensitivity_matrix=beat_inp.computeSensitivityMatrix2(FullResponse,FullResponse2)



    if options.TECH=="SVD":
        if j==0:
            [deltas, varslist ] = correctbeatEXP(x,y,dx, beat_inp, cut=cut, app=0, path=options.path)
        if j==1:
            [deltas, varslist ] = correctbeatEXP2(x,x2,y,y2,dx,dx2, beat_inp, cut=cut, app=0, path=options.path)
        if 1:#All accelerators
            iteration=0 # Let's remove too low useless correctors
            while (len(filter(lambda x: abs(x)< MinStr, deltas))>0):
                iteration=1+iteration
                il=len(varslist)
                varslist_t=[]
                for i in range(0,il):
                    if (abs(deltas[i]) > MinStr):
                        varslist_t.append(varslist[i])
                varslist=varslist_t
                if len(varslist)==0:
                    print "You want to correct with too high cut on the corrector strength"
                    sys.exit()
                if j==0:
                    beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)
                    sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
                    [deltas, varslist ] = correctbeatEXP(x,y,dx, beat_inp, cut=cut, app=0, path=options.path)
                if j==1:
                    beat_inp=beat_input2(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, phasexlist2, phaseylist2, betaxlist2, betaylist2, displist2, wei)
                    sensitivity_matrix=beat_inp.computeSensitivityMatrix2(FullResponse,FullResponse2)
                    [deltas, varslist ] = correctbeatEXP2(x,x2,y,y2,dx,dx2, beat_inp, cut=cut, app=0, path=options.path)
                    print "Initial correctors:", il, ". Current: ",len(varslist), ". Removed for being lower than:", MinStr, "Iteration:", iteration
        print deltas
    
    
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
        for corr in corrsYASP:
            print >>f, "#SETTING", corr,  corrsYASP[corr]
        for corr in corrs:
            print >>g, "K"+corr, corrs[corr]
        f.close()
        g.close()
        
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

    return 0    





    

if __name__ == '__main__':


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
    parser.add_option("-n", "--ncorr", 
		 help="Number of Correctors for MICADO",
		 metavar="NCORR", default=5,dest="ncorr")
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
                  metavar="MODELCUT", default="0.02,0.2" , dest="modelcut")
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
                help="Weighting factor (phasex, phasey, betax, betay, dispersion, tunes)",
                metavar="WGT", default="1,1,0,0,0,10", dest="WGT")


    (options, args) = parser.parse_args()
    correctpy(options, args)
    
