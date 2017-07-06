"""
This Python script produces the responses of the beta, phase and horizontal dispersion on the quadrupole strengths or
on other variables as horizontal orbit bumps at sextupoles.

The response matrix is stored in the following 'pickled' file:
 - FullResponse
The files are saved in option -p(output_path).

These response matrices are used to calculate the corrections by the correction scripts.
"""

import os
import sys
import math
import time
import cPickle
import optparse
import multiprocessing
import json
from collections import OrderedDict
import numpy
import pandas
sys.path.append("/afs/cern.ch/work/j/jcoellod/public/Beta-Beat.src/")
#import metaclass
from Utilities import tfs_pandas 
from madx import madx_wrapper

#===================================================================================================
# _parse_args()-function
#===================================================================================================
def _parse_args():
    ''' Parses the arguments, checks for valid input and returns tupel '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                      help="Which accelerator: LHCB1 LHCB2",
                      default="LHCB1", dest="accel")
    parser.add_option("-p", "--path",
                      help="path to save",
                      default="./", dest="path")
    parser.add_option("-c", "--core",
                      help="core files",
                      default="/afs/cern.ch/work/l/lmalina/online_model/", dest="core")
    parser.add_option("-k", "--deltak",
                      help="delta K1L to be applied to quads for sensitivity matrix",
                      default="0.00002", dest="k")

    options, _ = parser.parse_args()

    return options


class _InputData(object):
    """ Static class to access user input parameter and num_of_cpus and process pool"""
    output_path = ""
    delta_k = 0.0
    core_path_with_accel = ""
    core_path_without_accel = ""
    accel = ""

    number_of_cpus = 0
    process_pool = None

    @staticmethod
    def static_init(accel, output_path, path_to_core_files_without_accel, delta_k):
        if accel not in ("LHCB1", "LHCB2"):
            raise ValueError("Unknown accelerator: " + accel)
        _InputData.output_path = output_path
        _InputData.delta_k = float(delta_k)
        _InputData.accel = accel
        _InputData.core_path_with_accel = os.path.join(path_to_core_files_without_accel, accel)
        _InputData.core_path_without_accel = path_to_core_files_without_accel
        _InputData.number_of_cpus = multiprocessing.cpu_count()
        _InputData.process_pool = multiprocessing.Pool(processes=_InputData.number_of_cpus)

    def __init__(self):
        raise NotImplementedError("static class _InputData cannot be instantiated")


#=======================================================================================================================
# response()-function
#=======================================================================================================================
def response(accel, output_path, path_to_core_files_without_accel, delta_k):

    _InputData.static_init(accel, output_path, path_to_core_files_without_accel, delta_k)
    _generate_fullresponse_for_beta()

def _generate_fullresponse_for_beta():
    print "_generate_fullresponse_for_beta"
    path_all_lists_json_file = os.path.join(_InputData.core_path_without_accel, "AllLists.json") #"vqs.json" or AllLists.json
    knobsdict = json.load(file(path_all_lists_json_file, 'r'))
    print "Loaded json file: " + path_all_lists_json_file
    var_key = ""
    if _InputData.accel == "LHCB1":
        var_key = "VQ1"
    if _InputData.accel == "LHCB2":
        var_key = "VQ2"
    var_key="Q" #TODO
    
    variables = knobsdict[var_key] 
    delta1 = numpy.zeros(len(variables)) * 1.0   # Zero^th of the variables
    incr = numpy.ones(len(variables)) * _InputData.delta_k
    incr_dict = {}
    
    ######## loop over normal variables
    f = open(_join_with_output("iter.madx"), "w")
    for i in range(0, len(delta1)):  # Loop over variables
        delta = numpy.array(delta1)
        delta[i] = delta[i] + incr[i]
        var = variables[i]
        incr_dict[var] = incr[i]
        print >> f, var, "=", var, "+(", delta[i], ");"
        print >> f, "twiss, file=\"" + _join_with_output("twiss." + var) + "\";"
        print >> f, var, "=", var, "-(", delta[i], ");"
    incr_dict['0'] = 0.0
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
    FullResponse = {}
    for key, value in a:
        value['incr']=incr_dict[key]
        FullResponse[key] = value
    
    resp=pandas.Panel.from_dict(FullResponse)
    resp=resp.transpose(2,0,1)
    # After transpose e.g: resp[NDX, kqt3, bpm12l1.b1]
    # The magnet called "0" is no change (nominal model)
    resp['NDX'] = resp.xs('DX',axis=0) / numpy.sqrt(resp.xs('BETX',axis=0))
    resp['BBX'] = resp.xs('BETX',axis=0) / resp.loc['BETX','0',:]
    resp['BBY'] = resp.xs('BETY',axis=0) / resp.loc['BETY','0',:]
    resp = resp.subtract(resp.xs('0'),axis=1)
    # Remove beta-beating of nominal model with itself (bunch of zeros)
    resp.drop('0',axis=1,inplace=True)
    resp = resp.div(resp.loc['incr',:,:])
    full = {'MUX': resp.xs('MUX',axis=0),
            'MUY': resp.xs('MUY',axis=0),
            'BBX': resp.xs('BBX',axis=0),
            'BBY': resp.xs('BBY',axis=0),
            'NDX': resp.xs('NDX',axis=0),
            'Q' : resp.loc[['Q1','Q2'],:,resp.minor_axis[0]]
            }
    _dump(_join_with_output("FullResponse"), full)


#=======================================================================================================================
# helper functions
#=======================================================================================================================

def _join_with_output(*path_tokens):
    return os.path.join(_InputData.output_path, *path_tokens)

#def _add_type_and_set_index(x,ind):
#     x['TYPE']=ind
#     x.set_index('TYPE', inplace=True)
 #    return x

DEV_NULL = os.devnull


def _callMadx(pathToInputFile):
    return madx_wrapper.resolve_and_run_file(pathToInputFile, log_file=DEV_NULL)
    

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
    ret = os.system(cmd)
    if ret is not 0:
        raise ValueError("COMMAND: %s finished with exit value %i" % (cmd, ret))


def _loadtwiss_beta(varandpath):
    (var, path) = varandpath
    x = 0
    try:
        x = tfs_pandas.read_tfs(path + "/twiss." + var)
        x = x.set_index('NAME').drop_duplicates()
        x['Q1']=x.headers['Q1']
        x['Q2']=x.headers['Q2']
        os.remove(path + "/twiss." + var)
        print x.headers['Q2']
    except IOError as e:
        print e
        return []
    return var, x


#=======================================================================================================================
# main invocation
#=======================================================================================================================
def _start():
    timeStartGlobal = time.time()
    options = _parse_args()

    response(
         accel=options.accel,
         output_path=options.path,
         path_to_core_files_without_accel=options.core,
         delta_k=options.k
         )

    timeGlobal = time.time() - timeStartGlobal
    print "Duration:", timeGlobal, "s"

if __name__ == '__main__':
    _start()
