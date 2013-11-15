"""
.. module: MODEL.LHCB.fullresponse.generateFullResponse_parallel

Created on ??

Creates the fullresponse and stores it in the following 'pickled' files:
 - FullResponse
 - FullResponse_couple
 - FullResponse_chromcouple
 The files are saved in option -p. They are used in the correction scripts.

Options::

  -a ACCEL, --accel=ACCEL
                        Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL
  -p PATH, --path=PATH  path to save. Have to contain 'job.iterate.madx' and 'modifiers.madx'
  -c CORE, --core=CORE  core files
  -k K, --deltak=K      delta k to be applied to quads for sensitivity matrix


.. moduleauthor:: Unknown
"""


import os
import math
import time
import cPickle
import optparse
import multiprocessing
import json

import numpy

import __init__ # @UnusedImport init will include paths
import Python_Classes4MAD.metaclass as metaclass
import Python_Classes4MAD.madxrunner as madxrunner

def shell_command(cmd):
#    print 'process id:', os.getpid()
    ret = os.system(cmd)
    if ret is not 0:
        raise ValueError("COMMAND: %s finished with exit value %i" % (cmd, ret))

devNull = open(os.devnull, "w")

def callMadx(pathToInputFile, attemptsLeft=5):
    result = madxrunner.runForInputFile(pathToInputFile, stdout=devNull)
    if result is not 0: # then madx failed for whatever reasons, lets try it again (tbach)
        print "madx failed. result:", result, "pathToInputFile:", pathToInputFile, "attempts left:", attemptsLeft
        if attemptsLeft is 0:
            raise Exception("madx finally failed, can not continue")
        print "lets wait 0.5s and try again..."
        time.sleep(0.5)
        return callMadx(pathToInputFile, attemptsLeft - 1)
    return result
    
def dump(pathToDump, content):
    dumpFile = open(pathToDump, 'wb')
    cPickle.Pickler(dumpFile, -1).dump(content)
    dumpFile.close()


def parallel_command(period,numberOfCases):
    global path
    iterfile=open(path+'/iter.madx', 'r')
    lines=iterfile.readlines()
    iterfile.close()
    casesperprocess=int(math.ceil(numberOfCases*1.0/numberofCPUs))
    linesperprocess=casesperprocess*period
    
    iterFilePaths = []
    for i in range(len(lines)):      #split the iter.madx using in final number of processes 
        if (i % linesperprocess == 0):
            proid = i / linesperprocess + 1
            iterFilePath = path + "/iter." + str(proid) + ".madx"
            iterFile = open(iterFilePath, 'w')
            iterFilePaths.append(iterFilePath)
        iterFile.write(lines[i])
        if i == len(lines) - 1 or (i % linesperprocess == linesperprocess - 1):
            iterFile.close()
            
    
    # Prepare copies of the job.iterate.madx and all the shell commands
    madxFilePaths = []
    for i in range(1,proid+1):
        cmd='sed \'s/iter.madx/iter.'+str(i)+'.madx/g\' '+path+'/job.iterate.madx > '+path+'/job.iterate.'+str(i)+'.madx'
        shell_command(cmd)
        madxFilePaths.append(path+'/job.iterate.'+str(i)+'.madx')
    
#    print "send jobs to madx in parallel, number of jobs:", len(madxFilePaths)
    pool.map(callMadx, madxFilePaths)
    
    # clean up again (tbach)
    for madxFilePathsItem in madxFilePaths:
        os.remove(madxFilePathsItem)
    for iterFilePathsItem in iterFilePaths:
        os.remove(iterFilePathsItem)


def loadtwiss_beta(varandpath):
    (var, path) = varandpath
#    print "Reading twiss." + var
    x = 0
    try:
        x = metaclass.twiss(path + "/twiss." + var)
        os.remove(path + "/twiss." + var)
    except IOError as e:
        print e
        return []
    return var, x

def loadtwiss_coup(varandpath):
    var, path = varandpath
#    print "Reading twiss." + var
    x = 0
    try:
        x = metaclass.twiss(path + "/twiss." + var)
        x.Cmatrix()
        os.remove(path + "/twiss." + var)
    except IOError as e:
        print e
        return []
    return var, x

def loadtwiss_chrom_coup(varandpathanddpp):
    var, path, dpp = varandpathanddpp
#    print  "Reading twiss.dp+." + var
    xp = 0
    xm = 0
    try:
        xp = metaclass.twiss(path + "/twiss.dp+." + var)
        xp.Cmatrix()
#    print  "Reading twiss.dp-." + var
        xm = metaclass.twiss(path + "/twiss.dp-." + var)
        xm.Cmatrix()
    except IOError as e:
        print e
        return []
    # Initializing and Calculating chromatic coupling for every BPM
    xp.Cf1001r = []
    xp.Cf1001i = []
    xp.Cf1010r = []
    xp.Cf1010i = []
    for j in range(len(xp.NAME)):

        vvv = (xp.F1001R[j] - xm.F1001R[j]) / (2 * dpp)
        xp.Cf1001r.append(vvv)

        vvv = (xp.F1001I[j] - xm.F1001I[j]) / (2 * dpp)
        xp.Cf1001i.append(vvv)

        vvv = (xp.F1001R[j] - xm.F1001R[j]) / (2 * dpp)
        xp.Cf1010r.append(vvv)

        vvv = (xp.F1010I[j] - xm.F1010I[j]) / (2 * dpp)
        xp.Cf1010i.append(vvv)

    os.remove(path + "/twiss.dp+." + var)
    os.remove(path + "/twiss.dp-." + var)
    return var, xp



if __name__ == '__main__':
    timeStartGlobal = time.time()
    numberofCPUs = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=numberofCPUs)
    
    
    ##### optionparser
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                      help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                      default="LHCB1", dest="accel")
    parser.add_option("-p", "--path",
                      help="path to save",
                      default="./", dest="path")
    parser.add_option("-c", "--core",
                      help="core files",
                      default="/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB/fullresponse/", dest="core")
    parser.add_option("-k", "--deltak",
                      help="delta k to be applied to quads for sensitivity matrix",
                      default="0.00002", dest="k")
    
    
    (options, args) = parser.parse_args()
    
    # paths
    corepath = options.core
    accel = options.accel
    path = options.path
    
    variables = 0
    #
    # Chromatic coupling
    #
    
    FullResponse = {}   #Initialize FullResponse
    knobsdict = json.load(file(corepath + "/" + accel + '/AllLists_chromcouple.json', 'r'))
    variables = knobsdict["kss"]
    delta1 = numpy.zeros(len(variables)) * 1.0   #Zero^th of the variables
    incr = numpy.ones(len(variables)) * 0.05    #increment of variables
    dpp = 0.0001
    FullResponse['incr'] = incr           #Store this info for future use
    FullResponse['delta1'] = delta1
    
    
    ######## loop over normal variables
    f = open(path + '/iter.madx', 'w')
    for i in range(0, len(delta1)) : #Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        print >> f, var, "=", var, "+(", delta[i], ");"
        print >> f, "twiss, deltap= " + str(dpp) + ",file=\"" + path + "/twiss.dp+." + var + "\";"
        print >> f, "twiss, deltap=-" + str(dpp) + ",file=\"" + path + "/twiss.dp-." + var + "\";"
        print >> f, var, "=", var, "-(", delta[i], ");"
    
    
    print >> f, "twiss, deltap= " + str(dpp) + ",file=\"" + path + "/twiss.dp+.0\";"
    print >> f, "twiss, deltap=-" + str(dpp) + ",file=\"" + path + "/twiss.dp-.0\";"
    f.close()
    print "Running MADX"
    #shell_command('/afs/cern.ch/group/si/slap/bin/madx < '+path+'/job.iterate.madx')
    parallel_command(period=4, numberOfCases=len(delta1) + 1) # period=4 since there are 4 lines in iter.madx per case, numberofcases has +1 since there is the 0 case
    
    
    
    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, path, dpp])
    a = pool.map(loadtwiss_chrom_coup, newvarsforloop)
    for key, value in a:
        FullResponse[key] = value
    
    dump(path + '/FullResponse_chromcouple', FullResponse)
    
    
    
    #
    # Coupling
    #
    
    FullResponse = {}   #Initialize FullResponse
    knobsdict = json.load(file(corepath + "/" + accel + '/AllLists_couple.json', 'r'))
    variables = knobsdict["Qs"]
    delta1 = numpy.zeros(len(variables)) * 1.0   #Zero^th of the variables
    incr = numpy.ones(len(variables)) * 0.0001    #increment of variables
    
    
    FullResponse['incr'] = incr           #Store this info for future use
    FullResponse['delta1'] = delta1       #"     "     "
    
    ######## loop over normal variables
    f = open(path + '/iter.madx', 'w')
    for i in range(0, len(delta1)) : #Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        print >> f, var, "=", var, "+(", delta[i], ");"
        print >> f, "twiss, file=\"" + path + "/twiss." + var + "\";"
        print >> f, var, "=", var, "-(", delta[i], ");"
    
    print >> f, "twiss, file=\"" + path + "/twiss.0\";"
    f.close()
    #Sending the mad jobs in parallel
    parallel_command(period=3, numberOfCases=len(delta1) + 1)
    #Loading the twiss files into fullresp in parallel
    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, path])
    a = pool.map(loadtwiss_coup, newvarsforloop)
    for key, value in a:
        FullResponse[key] = value
        
    dump(path + '/FullResponse_couple', FullResponse)
    
    
    
    #
    #
    # Beta
    #
    #
    FullResponse = {}   #Initialize FullResponse
    knobsdict = json.load(file(corepath + "/" + accel + '/AllLists.json', 'r'))
    variables = knobsdict["Q"]
    delta1 = numpy.zeros(len(variables)) * 1.0   #Zero^th of the variables
    #incr=ones(len(variables))*0.00005    #increment of variables    #### when squeeze low twiss fails because of to big delta
    incr = numpy.ones(len(variables)) * float(options.k)
    
    
    FullResponse['incr'] = incr           #Store this info for future use
    FullResponse['delta1'] = delta1       #"     "     "
    
    ######## loop over normal variables
    f = open(path + '/iter.madx', 'w')
    for i in range(0, len(delta1)) : #Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        print >> f, var, "=", var, "+(", delta[i], ");"
        print >> f, "twiss, file=\"" + path + "/twiss." + var + "\";"
        print >> f, var, "=", var, "-(", delta[i], ");"
    
    print >> f, "twiss,file=\"" + path + "/twiss.0\";"
    f.close()
    
    parallel_command(period=3, numberOfCases=len(delta1) + 1)
    #
    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, path])
    
    a = pool.map(loadtwiss_beta, newvarsforloop)
    for key, value in a:
        FullResponse[key] = value
    
    dump(path + '/FullResponse', FullResponse)
    
    timeGlobal = time.time() - timeStartGlobal
    print "Duration:", timeGlobal, "s"
