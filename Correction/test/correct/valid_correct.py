##
## version November updated for LHC by Glenn Vanbavinckhove + removing mad2dev (handled in GUI now)
## + also removed option j (one beam or both beam, was not used for anything)
##################################################################################################
##
##
##
##

import sys, os
import pickle
from optparse import OptionParser
import json



import __init__ # @UnusedImport init will include paths
from Python_Classes4MAD.GenMatrix import *
from Python_Classes4MAD.BCORR import *
import Python_Classes4MAD.metaclass
import Utilities.iotools

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
                  help="Maximum error allowed for the phase and dispersion measurements, separated by commas; e.g. -e 0.013,0.2",
                  metavar="ERRORCUT", default="0.013,0.2" , dest="errorcut")
parser.add_option("-m", "--modelcut",
                  help="Maximum difference allowed between model and measured phase and dispersion, separated by commas; e.g. -e 0.02,0.2",
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
                  help="0 Just quads from one beam are used, 1 all quads (default is 0)",
                  metavar="JUSTONEBEAM", default=0 , dest="JustOneBeam")

parser.add_option("-w", "--weight",
                help="Weighting factor (phasex, phasey, betax, betay, dispersion, tunes)",
                metavar="WGT", default="1,1,0,0,1,10", dest="WGT")


(options, args) = parser.parse_args()





def  MakeBetaList(x, m, modelcut=40, errorcut=20):   # Errors are in meters (
    t=[]
    cou=0
    keys=x.__dict__.keys()
    if "BETY" in keys:
        bmdl="BETYMDL"
        STD=x.STDBETY
        BET=x.BETY
    else:
        bmdl="BETXMDL"
        STD=x.STDBETX
        BET=x.BETX
    print "Number of x BPMs",len(x.NAME)
    for i in range(len(x.NAME)):
        bm=x.__dict__[bmdl][i]
        if (STD[i] < errorcut and abs(BET[i]-bm) < modelcut):
            try:
                m.indx[x.NAME[i].upper()]
            except:
                print "Not in Response:", x.NAME[i].upper()
                cou=cou+1
            else:
                t.append(x.NAME[i])
        else:
            cou=cou+1
    if cou > 0:
        print "Warning in MakeBetaList: ", cou, " BPM  removed from data for not beeing in the model or having too large error deviations: ", bmdl, modelcut, "STDPH",errorcut, "LEN", len(t)
    return t



print "Selected accelerator:", options.ACCEL
options.path = os.path.normpath(options.path)
# in future versions, one should always use os.path.join, then we do not need to take care of the ending slash (tbach)
print "Path to measurements:", options.path
options.rpath = os.path.normpath(options.rpath) + os.sep
betapath = options.rpath
print "Path to Repository:", betapath
options.OPT = os.path.normpath(options.OPT) + os.sep
accelpath=options.OPT
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

##### implementing cuts for

modelcut=float(options.modelcut.split(",")[0])
errorcut=float(options.errorcut.split(",")[0])
modelcutdx=float(options.modelcut.split(",")[1])
errorcutdx=float(options.errorcut.split(",")[1])
cut= float(options.cut)
print "Model, error and SVD cuts:", modelcut, errorcut, cut
print "Starting loading Full Response optics"
FullResponse=pickle.load(open(accelpath+'/FullResponse','rb'))
print "Loading ended"



try:
    x= Python_Classes4MAD.metaclass.twiss(options.path+'/getphasex_free.out')
    y= Python_Classes4MAD.metaclass.twiss(options.path+'/getphasey_free.out')
    print "Loading free"
except:
    x= Python_Classes4MAD.metaclass.twiss(options.path+'/getphasex.out')
    y= Python_Classes4MAD.metaclass.twiss(options.path+'/getphasey.out')

try:
    xbet= Python_Classes4MAD.metaclass.twiss(options.path+'/getbetax_free.out')
    ybet= Python_Classes4MAD.metaclass.twiss(options.path+'/getbetay_free.out')
except:
    xbet= Python_Classes4MAD.metaclass.twiss(options.path+'/getbetax.out')
    ybet= Python_Classes4MAD.metaclass.twiss(options.path+'/getbetay.out')

try:
    dx= Python_Classes4MAD.metaclass.twiss(options.path+'/getNDx.out') # changed by Glenn Vanbavinckhove (26/02/09)
except:
    print "WARNING: No good dispersion or inexistent file getDx"
    print "WARNING: Correction will not take into account NDx"
    dx=[]


accel=options.ACCEL

if "LHC" in accel:
    main="LHCB"
    accelpath=betapath+'/MODEL/'+main+'/fullresponse/'+options.ACCEL+'/'
else:
    main="SPS"
    accelpath=betapath+'/MODEL/'+main+'/fullresponse/'+options.ACCEL+'/'

knobsdict=json.load(file(accelpath + '/AllLists.json','r'))
# extra depdency to be able to handle to different magnets group
listvar=options.var.split(",")
varslist=[]
for var in listvar:
    variable=knobsdict[var]
    varslist=varslist+variable



intqx=int(FullResponse['0'].Q1)
intqy=int(FullResponse['0'].Q2)

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


print "Experiment tunes: ", x.Q1, y.Q2

variables=varslist
phasexlist=MakePairs(x, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
phaseylist=MakePairs(y, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
betaxlist=MakeBetaList(xbet, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
betaylist=MakeBetaList(ybet, FullResponse['0'], modelcut=modelcut, errorcut=errorcut)
displist=MakeList(dx, FullResponse['0'])
print "Input ready"

wei=[]
for i in range(0,6):
    wei.append(int(options.WGT.split(",")[i]))
#wei=[1,1,1,1,1,10] # Weights of phasex phasey betax betay disp and tunes
print "weight=",wei
beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)

sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)



if options.TECH=="SVD":
    [deltas, varslist ] = correctbeatEXP(x,y,dx, beat_inp, cut=cut, app=0, path=options.path, xbet=xbet, ybet=ybet)
    if 1:                           #All accelerators
        iteration=0
        # Let's remove too low useless correctors
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
            beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)
            sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
            [deltas, varslist ] = correctbeatEXP(x,y,dx, beat_inp, cut=cut, app=0, path=options.path, xbet=xbet, ybet=ybet)
            print "Initial correctors:", il, ". Current: ",len(varslist), ". Removed for being lower than:", MinStr, "Iteration:", iteration
    print deltas


if options.TECH=="MICADO":
    bNCorrNumeric(x,y,dx,beat_inp, cut=cut,ncorr=ncorr,app=0,path=options.path, beta_x=xbet, beta_y=ybet)

if options.ACCEL=="SPS":
    b= Python_Classes4MAD.metaclass.twiss(options.path+"/changeparameters.tfs")
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
    src = os.path.join(os.path.join(options.path, "changeparameters.tfs"))
    dst = os.path.join(os.path.join(options.path, "changeparameters.knob"))
    Utilities.iotools.copy_item(src, dst)

    # madx table
    b= Python_Classes4MAD.metaclass.twiss(options.path+"/changeparameters.tfs")
    mad=open(options.path+"/changeparameters.madx", "w")
    names=b.NAME
    delta=b.DELTA

    for i in range(len(names)):
        if cmp(delta[i],0)==1:
            mad.write(names[i]+" = "+names[i]+" + "+str(delta[i])+";\n")
        else:
            mad.write(names[i]+" = "+names[i]+" "+str(delta[i])+";\n")

    mad.write("return;")
    mad.close()
