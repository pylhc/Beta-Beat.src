
#--- beta beat for store with numpy

import pickle
import json

import os,shutil

import __init__ # @UnusedImport init will include paths
from optparse import OptionParser
import Python_Classes4MAD.metaclass
import Python_Classes4MAD.GenMatrix_chromcouple

########### START ###############


parser = OptionParser()
parser.add_option("-a", "--accel",
                 help="What accelerator: LHCB1 LHCB2 SPS RHIC",
                 metavar="ACCEL", default="LHCB1",dest="ACCEL")
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
                  help="weight on corrections (f1001.re, f1001.im, f1010.re, f1010.im)",
                  metavar="Dy", default="1,1,0,0,0", dest="Dy")
parser.add_option("-o", "--opt",
                  help="Path to fullresponse_chromcouple",
                  metavar="OPT", default="./", dest="opt")
parser.add_option("-v", "--Variables",
                  help="variables split with ,",
                  metavar="var", default="kss" , dest="var")



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
FullResponse=pickle.load(open(options.opt+'/FullResponse_chromcouple','rb'))

print "Loading ended"


couple = Python_Classes4MAD.metaclass.twiss(os.path.join(options.path,'chromcoupling.out'))



weights=options.Dy.split(',')


path_all_lists_json_file = os.path.join(accelpath, "AllLists_chromcouple.json")
knobsdict = json.load(file(path_all_lists_json_file, 'r'))
print "Loaded json file: " + path_all_lists_json_file
listvar=options.var.split(",")
print listvar
varslist=[]
for var in listvar:
    varslist = varslist + knobsdict[var]

variables=varslist



MADTwiss=FullResponse['0']
MADTwiss.Cmatrix()
mode='C'
couplelist = Python_Classes4MAD.GenMatrix_chromcouple.make_list(couple, MADTwiss, modelcutC, errorcutC, mode)
if len(couplelist)==0:
    raise ValueError("No valid BPM measurements, maybe your model-/errorcuts are too strict?")
#print couplelist
#mode='D'
#dispylist=MakeList(dispy, MADTwiss, modelcutD, errorcutD, mode)


#if int(options.Dy)==0:
    #wei=[1,1,0,0,0] # Weights of f1001 and f1010
#elif int(options.Dy)==1:
    #wei=[1,1,0,0,1] # Weights of f1001, f1010 and Dy
#else:
    #print "Specify the vertical dispersion correction option 1(turn-on) or 0(turn-off)"
wei=[int(weights[0]),int(weights[1]),int(weights[2]),int(weights[3]),int(weights[4])]




print "the weight option is "+str(wei)


print "entering couple input",len(couplelist)
chromcouple_inp = Python_Classes4MAD.GenMatrix_chromcouple.ChromCoupleInput(varslist, couplelist, wei)
print "computing the sensitivity matrix"
sensitivity_matrix=chromcouple_inp.computeSensitivityMatrix(FullResponse)
#print "sensitivity matrix", sensitivity_matrix

print "computing correct coupling "
[deltas, varslist ] = Python_Classes4MAD.GenMatrix_chromcouple.correctcouple(couple, chromcouple_inp, cut=cut, app=0, path=options.path)

print deltas



#print deltas, varslist
print "handling data"



if "LHC" in options.ACCEL:   #.knob should always exist to be sent to LSA!
    #TODO: check if gui sends changeparameters_chromcouple.knob or changeparameters_couple.knob(vimaier)
    shutil.copy(options.path+"/changeparameters_chromcouple.tfs", options.path+"/changeparameters_chromcouple.knob")


    #####

    v = Python_Classes4MAD.metaclass.twiss(options.path+"/changeparameters_chromcouple.tfs")
    mad=open(options.path+"/changeparameters_chromcouple.madx", "w")
    names=v.NAME
    delta=v.DELTA

    for i in range(len(names)):



        if cmp(delta[i],0)==1:
            mad.write(names[i]+" = "+names[i]+" + "+str(delta[i])+";\n");
        else:
            mad.write(names[i]+" = "+names[i]+" "+str(delta[i])+";\n");


    mad.write("return;");

    mad.close()



print "Correcting chrom couple finished "
