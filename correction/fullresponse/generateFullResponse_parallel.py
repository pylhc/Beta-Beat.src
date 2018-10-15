r"""
.. module: correction.fullresponse.generateFullResponse_parallel

Created on ??

This Python script produces the responses of the beta, phase and horizontal dispersion on the quadrupole strengths or
on other variables as horizontal orbit bumps at sextupoles.
The fullresponse.couple does the corresponding calculation for the coupling and the vertical dispersion.

The response matrices are stored it in the following 'pickled' files:
 - FullResponse
 - FullResponse_couple
 - FullResponse_chromcouple

 The files are saved in option -p(output_path).

These response matrices are used to calculate the corrections by the correction scripts.


Parallel execution:
*multiprocessing* from Python stdlib is used to spwan subprocesses and run MADX in parallel. Further the loading of the
MADX outputfiles will be executed in parallel.

Possible improvements:
 - The order of the three _gen*_for_*() functions are irrelevant and they do not depend on each other. Thus the execution
   could be parallelized. Therefore the created iter.madx and twiss files need an unique prefix. E.g.: 'beta_iter.madx'.
   All three functions would need own process pools and the handling of the input data have to be considered.
   For test_fullresponse_parallel.py the comparing of \*iter.madx should then be avoided since the 'valid' script does not produce them.
 - The three functions in the main function do basically the same. A good pattern to handle this is the template method:
   http://en.wikipedia.org/wiki/Template_method_pattern
 - Removing of the exec statements. This is very ugly and could be replaced with a function *get_variables_for_accel_and_fullresponse_type(accel, type)*.
   The function could import the corresponding module by checking function parameters in if/else statements and deliver the
   appropriate list.

Usage example::

    python generateFullResponse_parallel.py --accel=LHCB1 --path=/afs/cern.ch/work/v/vimaier/public/betabeatGui/temp/2013-10-22/models/LHCB1/a_mod --core=/afs/cern.ch/work/v/vimaier/public/Beta-Beat.src/MODEL/LHCB1/fullresponse --deltak=0.00002

Options::

  -a ACCEL, --accel=ACCEL
                        Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL
  -p PATH, --path=PATH  path to save. Have to contain 'job.iterate.madx' and 'modifiers.madx'
  -c CORE, --core=CORE  core files
  -k K, --deltak=K      delta k to be applied to quads for sensitivity matrix


.. moduleauthor:: Unknown
"""


import math
import time
import cPickle
import optparse
import multiprocessing
import json

import numpy
import sys
import os
from os.path import abspath, join, dirname, pardir
new_path = abspath(join(dirname(abspath(__file__)), pardir, pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

import Python_Classes4MAD.metaclass as metaclass
import utils.iotools
import correction.correction_helpers

import madx_wrapper


UNWANTED_CATEGORIES = ["LQ", "MQX", "MQXT", "Q", "QIP15", "QIP2", "getListsByIR"]


#===================================================================================================
# _parse_args()-function
#===================================================================================================
def _parse_args():
    ''' Parses the arguments, checks for valid input and returns tupel '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                      help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                      default="LHCB1", dest="accel")
    parser.add_option("-p", "--path",
                      help="path to save",
                      default="./", dest="path")
    parser.add_option("-c", "--core",
                      help="core files",
                      default="../correction/fullresponse/", dest="core")
    parser.add_option("-k", "--deltak",
                      help="delta k to be applied to quads for sensitivity matrix",
                      default="0.00002", dest="k")
    parser.add_option("-l", "--deltakl",
                      help="delta kl to be applied to quads for sensitivity matrix",
                      default="0.0", dest="kl")
    parser.add_option("-f", "--f", action="store_false", dest="fullresponse")
    parser.add_option("-t","--t", action="store_true", dest="fullresponse")

    options, _ = parser.parse_args()

    return options


class _InputData(object):
    """ Static class to access user input parameter and num_of_cpus and process pool"""
    output_path = ""
    delta_k = 0.0
    delta_kl = 0.0
    core_path_with_accel = ""

    number_of_cpus = 0
    process_pool = None
    fullresponse = True

    @staticmethod
    def static_init(accel, output_path, path_to_core_files_without_accel, delta_k, delta_kl, fullresponse):
        if accel not in ("LHCB1", "LHCB2", "SPS", "PSBOOSTER", "RHIC", "SOLEIL","PS"):
            raise ValueError("Unknown accelerator: " + accel)
        if not utils.iotools.dirs_exist(output_path):
            raise ValueError("Output path does not exists. It has to contain job.iterator.madx and modifiers.madx.")
        if not correction.correction_helpers.can_str_be_parsed_to_number(delta_k):
            raise ValueError("Delta k is not a number: " + delta_k)
        if not correction.correction_helpers.can_str_be_parsed_to_number(delta_kl):
            raise ValueError("Delta kl is not a number: " + delta_kl)
        if not utils.iotools.dirs_exist(os.path.join(path_to_core_files_without_accel, accel)):
            raise ValueError("Core path does not exist: " + _InputData.core_path_with_accel)

        _InputData.output_path = output_path
        _InputData.delta_k = float(delta_k)
        _InputData.delta_kl = float(delta_kl)
        _InputData.core_path_with_accel = os.path.join(path_to_core_files_without_accel, accel)
        _InputData.number_of_cpus = multiprocessing.cpu_count()
        _InputData.process_pool = multiprocessing.Pool(processes=_InputData.number_of_cpus)
        _InputData.fullresponse = fullresponse

    def __init__(self):
        raise NotImplementedError("static class _InputData cannot be instantiated")


#=======================================================================================================================
# main()-function
#=======================================================================================================================
def main(accel, output_path, path_to_core_files_without_accel, delta_k, delta_kl, fullresponse):

    _InputData.static_init(accel, output_path, path_to_core_files_without_accel, delta_k, delta_kl, fullresponse)
    if fullresponse == True:
        _generate_fullresponse_for_coupling(fullresponse)
        _generate_fullresponse_for_chromatic_coupling()
        _generate_fullresponse_for_beta()
    if fullresponse == False:
        _generate_fullresponse_for_coupling(fullresponse)


def _generate_fullresponse_for_chromatic_coupling():
    print "_generate_fullresponse_for_chromatic_coupling"
    path_all_lists_json_file = os.path.join(_InputData.core_path_with_accel, "AllLists_chromcouple.json")
    knobsdict = json.load(file(path_all_lists_json_file, 'r'))
    print "Loaded json file: " + path_all_lists_json_file
    variables = knobsdict["kss"]
    delta1 = numpy.zeros(len(variables)) * 1.0   # Zero^th of the variables
    incr = numpy.ones(len(variables)) * 0.05    # increment of variables
    incr_dict = {}
    dpp = 0.0001

    FullResponse = {}   # Initialize FullResponse
    FullResponse['incr'] = incr           # Store this info for future use
    FullResponse['incr_dict'] = incr_dict # "     "     "
    FullResponse['delta1'] = delta1

    ######## loop over normal variables
    f = open(_join_with_output("iter.madx"), "w")
    for i in range(0, len(delta1)):  # Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        incr_dict[var] = incr[i]
        print >> f, var, "=", var, "+(", delta[i], ");"
        print >> f, "twiss, deltap= " + str(dpp) + ",file=\"" + _join_with_output("twiss.dp+.") + var + "\";"
        print >> f, "twiss, deltap=-" + str(dpp) + ",file=\"" + _join_with_output("twiss.dp-.") + var + "\";"
        print >> f, var, "=", var, "-(", delta[i], ");"

    print >> f, "twiss, deltap= " + str(dpp) + ",file=\"" + _join_with_output("twiss.dp+.0") + "\";"
    print >> f, "twiss, deltap=-" + str(dpp) + ",file=\"" + _join_with_output("twiss.dp-.0") + "\";"
    f.close()
    print "Running MADX"
    _parallel_command(period=4, number_of_cases=len(delta1) + 1)  # period=4 since there are 4 lines in iter.madx per case, number_of_cases has +1 since there is the 0 case
    print "Finished MADX"

    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, _InputData.output_path, dpp])
    a = _InputData.process_pool.map(_loadtwiss_chrom_coup, newvarsforloop)
    for key, value in a:
        FullResponse[key] = value

    _dump(_join_with_output('FullResponse_chromcouple'), FullResponse)

    utils.iotools.copy_item(_join_with_output("iter.madx"), _join_with_output("chromcouple_iter.madx"))


def _generate_fullresponse_for_coupling(fullresponse):

    print "_generate_fullresponse_for_coupling"
    path_all_lists_json_file = os.path.join(_InputData.core_path_with_accel, "AllLists_couple.json")
    knobsdict = json.load(file(path_all_lists_json_file, 'r'))
    print "Loaded json file: " + path_all_lists_json_file
    if fullresponse == True:
        print knobsdict
        variables = knobsdict["Qs"]
        
    if fullresponse == False:
        print "Small coupling"
        variables = knobsdict["coupling_knobs"]
    delta1 = numpy.zeros(len(variables)) * 1.0   # Zero^th of the variables
    incr = numpy.ones(len(variables)) * 0.0001    # increment of variables
    incr_dict = {}

    FullResponse = {}   # Initialize FullResponse
    FullResponse['incr'] = incr           # Store this info for future use
    FullResponse['incr_dict'] = incr_dict # "     "     "
    FullResponse['delta1'] = delta1       # "     "     "

    ######## loop over normal variables
    f = open(_join_with_output('iter.madx'), 'w')
    for i in range(0, len(delta1)):  # Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        incr_dict[var] = incr[i]
        print >> f, var, "=", var, "+(", delta[i], ");"
        print >> f, "twiss, file=\"" + _join_with_output("twiss." + var) + "\";"
        print >> f, var, "=", var, "-(", delta[i], ");"

    print >> f, "twiss, file=\"" + _join_with_output("twiss.0") + "\";"
    f.close()
    #Sending the mad jobs in parallel
    print "Running MADX"
    _parallel_command(period=3, number_of_cases=len(delta1) + 1)
    print "Finished MADX"

    #Loading the twiss files into fullresp in parallel
    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, _InputData.output_path])
    a = _InputData.process_pool.map(_loadtwiss_coup, newvarsforloop)
    for key, value in a:
        FullResponse[key] = value

    _dump(_join_with_output("FullResponse_couple"), FullResponse)
    utils.iotools.copy_item(_join_with_output("iter.madx"), _join_with_output("couple_iter.madx"))


def _generate_fullresponse_for_beta():
    print "_generate_fullresponse_for_beta"
    path_all_lists_json_file = os.path.join(_InputData.core_path_with_accel, "AllLists.json")
    knobsdict = json.load(file(path_all_lists_json_file, 'r'))
    print "Loaded json file: " + path_all_lists_json_file
    variables = []
    for category in knobsdict.keys():
        if category not in UNWANTED_CATEGORIES:
            variables += knobsdict[category]
    variables = list(set(variables))
    delta1 = numpy.zeros(len(variables)) * 1.0   # Zero^th of the variables
    #incr=ones(len(variables))*0.00005    #increment of variables    #### when squeeze low twiss fails because of to big delta
    incr = numpy.ones(len(variables)) * _InputData.delta_k
    incr_dict = {}
    if _InputData.delta_kl > 0.0:
        variables2 = knobsdict["LQ"]
        delta2 = numpy.zeros(len(variables2)) * 1.0   # Zero^th of the variables
        incr2 = numpy.ones(len(variables2)) * _InputData.delta_kl
        variables = variables + variables2
        delta1 = numpy.concatenate((delta1,delta2))
        incr = numpy.concatenate((incr,incr2))
    FullResponse = {}   # Initialize FullResponse
    FullResponse['incr'] = incr           # Store this info for future use
    FullResponse['incr_dict'] = incr_dict # "     "     "
    FullResponse['delta1'] = delta1       # "     "     "
    if _InputData.core_path_with_accel[-5:] == "LHCB1":
        path_B1tfs = os.path.join(_InputData.core_path_with_accel, "varsKLvsK_B1.tfs")
        twiss_kl = metaclass.twiss(path_B1tfs)
    if _InputData.core_path_with_accel[-5:] == "LHCB2":
        path_B2tfs = os.path.join(_InputData.core_path_with_accel, "varsKLvsK_B2.tfs")
        twiss_kl = metaclass.twiss(path_B2tfs)
    
    ######## loop over normal variables
    f = open(_join_with_output("iter.madx"), "w")
    for i in range(0, len(delta1)):  # Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        incr_dict[var] = incr[i]
        if var[0] == "l":
            print >> f, var, "=", var, "+(", delta[i], ");", var[1:], "=", var[1:], "+(", var, "/", str(twiss_kl.LENGTH[twiss_kl.indx[var]]), ");"
            print >> f, "twiss, file=\"" + _join_with_output("twiss." + var) + "\";"
            print >> f, var[1:], "=", var[1:], "-(", var, "/", str(twiss_kl.LENGTH[twiss_kl.indx[var]]), ");", var, "=", var, "-(", delta[i], ");"
        else:
            print >> f, var, "=", var, "+(", delta[i], ");"
            print >> f, "twiss, file=\"" + _join_with_output("twiss." + var) + "\";"
            print >> f, var, "=", var, "-(", delta[i], ");"

    print >> f, "twiss,file=\"" + _join_with_output("twiss.0") + "\";"
    f.close()

    print "Running MADX"
    _parallel_command(period=3, number_of_cases=len(delta1) + 1)
    print "Finished MADX"

    #Loading the twiss files into fullresp in parallel
    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, _InputData.output_path])
    a = _InputData.process_pool.map(_loadtwiss_beta, newvarsforloop)
    for key, value in a:
        FullResponse[key] = value

    _dump(_join_with_output("FullResponse"), FullResponse)


#=======================================================================================================================
# helper functions
#=======================================================================================================================

def _join_with_output(*path_tokens):
    return os.path.join(_InputData.output_path, *path_tokens)


DEV_NULL = os.devnull


def _callMadx(pathToInputFile, attemptsLeft=5):
    try:
        madx_wrapper.resolve_and_run_file(pathToInputFile, log_file=DEV_NULL)
    except madx_wrapper.MadxError:  # then madx failed for whatever reasons, lets try it again (tbach)
        print "madx failed. pathToInputFile:", pathToInputFile, "attempts left:", attemptsLeft
        if attemptsLeft is 0:
            raise Exception("madx finally failed, can not continue")
        print "lets wait 0.5s and try again..."
        time.sleep(0.5)
        _callMadx(pathToInputFile, attemptsLeft - 1)


def _dump(pathToDump, content):
    dumpFile = open(pathToDump, 'wb')
    cPickle.Pickler(dumpFile, -1).dump(content)
    dumpFile.close()


def _parallel_command(period, number_of_cases):
    iterfile = open(_join_with_output("iter.madx"), 'r')
    lines = iterfile.readlines()
    iterfile.close()
    casesperprocess = int(math.ceil(number_of_cases*1.0 / _InputData.number_of_cpus))
    linesperprocess = casesperprocess * period

    iterFilePaths = []
    for i in range(len(lines)):   # split the iter.madx using in final number of processes
        if (i % linesperprocess == 0):
            proid = i / linesperprocess + 1
            iterFilePath = _join_with_output("iter." + str(proid) + ".madx")
            iterFile = open(iterFilePath, 'w')
            iterFilePaths.append(iterFilePath)
        iterFile.write(lines[i])
        if i == len(lines) - 1 or (i % linesperprocess == linesperprocess - 1):
            iterFile.close()

    # Prepare copies of the job.iterate.madx and all the shell commands
    madxFilePaths = []
    for i in range(1, proid + 1):
        cmd = 'sed \'s/iter.madx/iter.'+str(i)+'.madx/g\' '+_InputData.output_path+'/job.iterate.madx > '+_InputData.output_path+'/job.iterate.'+str(i)+'.madx'
        _shell_command(cmd)
        madxFilePaths.append(_InputData.output_path + '/job.iterate.' + str(i) + '.madx')

#    print "send jobs to madx in parallel, number of jobs:", len(madxFilePaths)
    _InputData.process_pool.map(_callMadx, madxFilePaths)

    # clean up again (tbach)
    for madxFilePathsItem in madxFilePaths:
        os.remove(madxFilePathsItem)
    for iterFilePathsItem in iterFilePaths:
        os.remove(iterFilePathsItem)


def _shell_command(cmd):
#    print 'process id:', os.getpid()
    ret = os.system(cmd)
    if ret is not 0:
        raise ValueError("COMMAND: %s finished with exit value %i" % (cmd, ret))


def _loadtwiss_beta(varandpath):
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


def _loadtwiss_coup(varandpath):
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


def _loadtwiss_chrom_coup(varandpathanddpp):
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


#=======================================================================================================================
# main invocation
#=======================================================================================================================
def _start():
    timeStartGlobal = time.time()
    options = _parse_args()

    main(
         accel=options.accel,
         output_path=options.path,
         path_to_core_files_without_accel=options.core,
         delta_k=options.k,
         delta_kl=options.kl,
         fullresponse=options.fullresponse
         )

    timeGlobal = time.time() - timeStartGlobal
    print "Duration:", timeGlobal, "s"

if __name__ == '__main__':
    _start()
