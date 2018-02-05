#!/usr/bin/env pythonafs


# Just to make sure that the path to the libraires is defined
import sys
import pickle
import os
import optparse
import json

import numpy
import __init__ # @UnusedImport init will include paths
import Utilities.iotools
import Python_Classes4MAD.GenMatrix_coupleDy as GenMatrix_coupleDy
import Python_Classes4MAD.metaclass as metaclass

########### START ###############
# /usr/bin/python /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src//Correction/correct_coupleDy.py
# -a LHCB1
# -p /afs/cern.ch/work/t/tbach/public/betabeatGui/temp/2012-09-28/LHCB1/Results/22-50-04_NORMALANALYSIS_SUSSIX_1
# -c 0.01
# -e 0.01,0.1
# -m 0.02,0.1
# -r /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/
# -s 0.00003
# -d 1,1,0,0,0
# -o /afs/cern.ch/work/t/tbach/public/betabeatGui/temp/2012-09-28/models/LHCB1/test2fr/
# -v coupling_knobs

# "accel", "path", "cut", "errorcut", "modelcut", "MinStr", "Dy", "opt", "Variables", "rpath"

parser = optparse.OptionParser()
parser.add_option("-a", "--accelerator",
                  help="What accelerator: LHCB1 LHCB2 SPS RHIC",
                  metavar="ACCEL", default="LHCB1", dest="ACCEL")
parser.add_option("-p", "--path",
                  help="Path to experimental files",
                  metavar="PATH", default="./", dest="path")
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
                  help="variables split with , - if coupling_knobs is given, the modelcut will be calculated automatically",
                  metavar="var", default="MQSb1" , dest="var")

(options, args) = parser.parse_args()

# internal options
printDebug = False  # If True, internal debug information will be printed (tbach)

modelcuts = options.modelcut.split(",")
modelcutC = float(modelcuts[0])
modelcutD = float(modelcuts[1])
errorcuts = options.errorcut.split(",")
errorcutC = float(errorcuts[0])
errorcutD = float(errorcuts[1])
weights = options.Dy.split(',')

print "Selected accelerator:", options.ACCEL
print "Path to measurements:", options.path
betapath = options.rpath
print "Path to Repository:", betapath
accelerator = options.ACCEL

if "LHC" in accelerator:
    main = "LHCB"
    accelpath = betapath + '/MODEL/' + main + '/fullresponse/' + accelerator + '/'
else:
    main = "SPS"
    accelpath = betapath + '/MODEL/' + main + '/fullresponse/' + accelerator + '/'

print "Path to Accelerator model", accelpath

MinStr = options.MinStr
print "Minimum corrector strength", MinStr

cut = float(options.cut)
print "Starting loading Full Response optics"
FullResponse = pickle.load(open(options.opt + '/FullResponse_couple', 'rb'))

print "Loading ended"

if os.path.isfile(options.path + "/getcouple_free.out"):
    couple = metaclass.twiss(options.path + '/getcouple_free.out')
    print "Will use free coupling"
elif os.path.isfile(options.path + "/getcouple.out"):
    couple = metaclass.twiss(options.path + '/getcouple.out')
    print "WARN: Free coupling not found!"
else:
    print "path", options.path, " does not contain a getcouple_free.out or getcouple.out file"
    sys.exit(1)


if "coupling_knobs" in options.var:
    if 0 != couple.F1001W.size:
        print "coupling_knobs mode. Trying to do automatic correcting"
        # we want to have a value which indicates the worst 5 percent, so we sort and get the value from the 95% index (tbach)
        sortedF1001W = numpy.sort(couple.F1001W)
        fivePercentIndex = int(numpy.floor(sortedF1001W.size * 0.95)) # no requirement to be accurate (tbach)

        fivePercentValue = sortedF1001W[fivePercentIndex];

        if printDebug: print "fivePercentIndex:", fivePercentIndex

        if numpy.allclose(fivePercentValue, 0):
            print "calculated value for modelcut is 0, exit"
            sys.exit(1)

        print "modelcutC will be changed from:", modelcutC, "to:", fivePercentValue
        modelcutC = fivePercentValue


dispy = []
if weights[4] == "1":
    dispy = metaclass.twiss(options.path + '/getDy.out')

path_all_lists_json_file = os.path.join(accelpath, "AllLists_couple.json")
knobsdict=json.load(file(path_all_lists_json_file, 'r'))
print "Loaded json file: " + path_all_lists_json_file
listvar = options.var.split(",")
varslist = []
for var in listvar:
    varslist = varslist + knobsdict[var]

MADTwiss = FullResponse['0']
MADTwiss.Cmatrix()
mode = 'C'
couplelist = GenMatrix_coupleDy.make_list(couple, MADTwiss, modelcutC, errorcutC, mode)
mode = 'D'
dispylist = GenMatrix_coupleDy.make_list(dispy, MADTwiss, modelcutD, errorcutD, mode)

wei = [int(weights[0]), int(weights[1]), int(weights[2]), int(weights[3]), int(weights[4])]

print "the weight option is " + str(wei)

optDy = options.Dy
print "entering couple input", len(couplelist)
couple_inp = GenMatrix_coupleDy.CoupleInput(varslist, couplelist, dispylist, wei)
print "computing the sensitivity matrix"
sensitivity_matrix = couple_inp.computeSensitivityMatrix(FullResponse)

print "computing correct coupling "
[deltas, varslist ] = GenMatrix_coupleDy.correctcouple(couple, dispy, couple_inp, cut=cut, app=0, path=options.path)

print deltas

print "handling data"
if options.ACCEL == "SPS":
    v = metaclass.twiss(options.path + "/changeparameters_couple.tfs")
    print '\nFor SPS, vertical bump file is imported:'
    print accelpath + '/Coupling/VBumps.py'
    execfile(accelpath + '/Coupling/VBumps.py')    # LOADS corrs
    execfile(accelpath + '/Coupling/VBumpsYASP.py') # LOADS corrsYASP
        #Output for YASP...
    f = open(options.path + "/changeparameters_couple.yasp", "w")
        #Output for Knob...
    #print options.path
    g = open(options.path + "/changeparameters_couple.knob", "w")
    h = open(options.path + "/changeparameters_couple.madx", "w")
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
    src = os.path.join(options.path, "changeparameters_couple.tfs")
    dst = os.path.join(options.path, "changeparameters_couple.knob")
    Utilities.iotools.copy_item(src, dst)

    ##### for bumps
    if "bumps" in listvar:
        print "passing trough bumps loop"
        v = metaclass.twiss(options.path + "/changeparameters_couple.tfs")
        Utilities.iotools.delete_item(os.path.join(options.path, "changeparameters_couple.tfs"))
        execfile(options.opt + '//Bumps.py')
        execfile(options.opt + '//mydictionary.py')
        filefile = open(options.path + "/changeparameters_couple.tfs", "w")
        filefile.write("* NAME  DELTA\n")
        filefile.write("$ %s     %le\n")
        for vcorr in corrs:
            filefile.write(vcorr+" "+ str(corrs[vcorr])+"\n")

        filefile.close()
    #####

    v = metaclass.twiss(options.path + "/changeparameters_couple.tfs")
    mad = open(options.path + "/changeparameters_couple.madx", "w")
    names = v.NAME
    delta = v.DELTA

    for i in range(len(names)):
        if "bumps" in listvar:
            if cmp(delta[i], 0) == 1:
                mad.write(names[i] + "->KICK:=" + str(delta[i]) + ";\n")
            else:
                mad.write(names[i] + "->KICK:=" + str(delta[i]) + ";\n")

        else:
            if cmp(delta[i], 0) == 1:
                mad.write(names[i] + " = " + names[i] + " + " + str(delta[i]) + ";\n")
            else:
                mad.write(names[i] + " = " + names[i] + " " + str(delta[i]) + ";\n")

    mad.write("return;")
    mad.close()

print "Correcting couple Dy finished with weight " + str(wei)
