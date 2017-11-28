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
import argparse
import multiprocessing
import json
import numpy
import pandas
sys.path.append("/afs/cern.ch/work/j/jcoellod/public/Beta-Beat.src/")
from Utilities import tfs_pandas
from madx import madx_wrapper
from Utilities import logging_tools
from Utilities.contexts import timeit

UNWANTED_VARIABLE_CATEGORIES = ["LQ", "MQX", "MQXT", "Q", "QIP15", "QIP2", "getListsByIR"]

LOG = logging_tools.get_logger(__name__, level_console=0)

"""
================================= Parse Arguments ===========================================
"""


def _parse_args():
    """ Parses the arguments, checks for valid input and returns options """
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--accel",
                        help="Which accelerator: LHCB1 LHCB2",
                        default="LHCB1",
                        dest="accel",
                        type=str)
    parser.add_argument("-o", "--out",
                        help="Path to outputfolder",
                        default="./",
                        dest="output_path",
                        type=str)
    parser.add_argument("-f", "--filename",
                        help="Name of fullresponse file.",
                        default="FullResponse",
                        dest="output_file",
                        type=str)
    parser.add_argument("-v", "--variables",
                        help="Path to file containing variable lists",
                        default="./all_lists/LHCB1/AllLists.json",
                        dest="variable_filepath",
                        type=str)
    parser.add_argument("-k", "--deltak",
                        help="delta K1L to be applied to quads for sensitivity matrix",
                        default=0.00002,
                        dest="delta_k",
                        type=float)

    options, remainder = parser.parse_known_args()

    if len(remainder) > 0:
        LOG.warning("'{file:s}' does not know the arguments '{args:s}'".format(
            file=os.path.split(__file__)[1],
            args=', '.join(remainder),
        ))

    return options


"""
================================= FullResponse ===========================================
"""

def generate_fullresponse(opt):
    LOG.debug("Generating Fullresponse.")
    variables = _get_variables_from_file(opt.variable_filepath)
    process_pool, incr_dict = _generate_madx_jobs(variables, opt.output_path, opt.delta_k)
    fullresponse = _load_madx_results(variables, opt.output_path, process_pool, incr_dict)
    _save_fullresponse(os.path.join(opt.output_path, opt.output_file), fullresponse)


def _get_variables_from_file(variables_filepath):
    """ Load variables list from json file """
    LOG.debug("Loading variables from file {:s}".format(variables_filepath))

    with open(variables_filepath, 'r') as var_file:
        var_dict = json.load(var_file)

    variables = []
    for category in var_dict.keys():
        # TODO: Make it more general and without hardcoded unwanted categories
        if category not in UNWANTED_VARIABLE_CATEGORIES:
            variables += var_dict[category]
    return list(set(variables))


def _generate_madx_jobs(variables, outputpath, delta_k):
    """ Generates madx job-files and executes them in parallel """
    LOG.debug("Generating MADX jobs.")
    incr_dict = {'0': 0.0}

    # loop over normal variables
    with open(os.path.join(outputpath, "iter.madx"), "w") as f:
        for i in range(0, len(variables)):  # Loop over variables
            var = variables[i]
            incr_dict[var] = delta_k
            f.write("{var:s}={var:s}+({delta:f});\n".format(var=var, delta=delta_k))
            f.write("twiss, file='{:s}';\n".format(os.path.join(outputpath, "twiss." + var)))
            f.write("{var:s}={var:s}-({delta:f});\n".format(var=var, delta=delta_k))

        f.write("twiss, file='{:s}';\n".format(os.path.join(outputpath, "twiss.0")))

    process_pool = _parallel_command(outputpath, period=3, number_of_cases=len(variables) + 1)
    return process_pool, incr_dict


def _load_madx_results(variables, outputpath, process_pool, incr_dict):
    LOG.debug("Loading Madx Results.")
    varsforloop = variables + ['0']
    newvarsforloop = []
    for value in varsforloop:
        newvarsforloop.append([value, outputpath])
    a = process_pool.map(_loadtwiss_beta, newvarsforloop)
    FullResponse = {}
    for key, value in a:
        value['incr'] = incr_dict[key]
        FullResponse[key] = value

    resp = pandas.Panel.from_dict(FullResponse)
    resp = resp.transpose(2, 0, 1)
    # After transpose e.g: resp[NDX, kqt3, bpm12l1.b1]
    # The magnet called "0" is no change (nominal model)
    resp['NDX'] = resp.xs('DX', axis=0) / numpy.sqrt(resp.xs('BETX', axis=0))
    resp['BBX'] = resp.xs('BETX', axis=0) / resp.loc['BETX', '0', :]
    resp['BBY'] = resp.xs('BETY', axis=0) / resp.loc['BETY', '0', :]
    resp = resp.subtract(resp.xs('0'), axis=1)
    # Remove beta-beating of nominal model with itself (bunch of zeros)
    resp.drop('0', axis=1, inplace=True)
    resp = resp.div(resp.loc['incr', :, :])
    full = {'MUX': resp.xs('MUX', axis=0),
            'MUY': resp.xs('MUY', axis=0),
            'BBX': resp.xs('BBX', axis=0),
            'BBY': resp.xs('BBY', axis=0),
            'NDX': resp.xs('NDX', axis=0),
            'Q': resp.loc[['Q1', 'Q2'], :, resp.minor_axis[0]]
            }
    return full


def _save_fullresponse(outputfile, fullresponse):
    """ Dumping the FullResponse file """
    LOG.debug("Saving Fullresponse into file '{:s}'".format(outputfile))
    with open(outputfile, 'wb') as dump_file:
        cPickle.Pickler(dump_file, -1).dump(fullresponse)


"""
================================= Helpers ===========================================
"""

def _join_with_output(output_root, *path_tokens):
    return os.path.join(output_root, *path_tokens)


def _callMadx(pathToInputFile):
    return madx_wrapper.resolve_and_run_file(pathToInputFile, log_file=os.devnull)


def _parallel_command(outputpath, period, number_of_cases):
    with open(os.path.join(outputpath, "iter.madx"), 'r') as iterfile:
        lines = iterfile.readlines()

    number_of_cpus = multiprocessing.cpu_count()
    process_pool = multiprocessing.Pool(processes=number_of_cpus)

    casesperprocess = int(math.ceil(number_of_cases*1.0 / number_of_cpus))
    linesperprocess = casesperprocess * period

    iterFilePaths = []
    for i in range(len(lines)):   # split the iter.madx using in final number of processes
        if (i % linesperprocess == 0):
            proid = i / linesperprocess + 1
            iterFilePath = os.path.join(outputpath, "iter." + str(proid) + ".madx")
            iterFile = open(iterFilePath, 'w')
            iterFilePaths.append(iterFilePath)
        iterFile.write(lines[i])
        if i == len(lines) - 1 or (i % linesperprocess == linesperprocess - 1):
            iterFile.close()

    # Prepare copies of the job.iterate.madx and all the shell commands
    madxFilePaths = []
    for i in range(1, proid + 1):
        file_in = os.path.join(outputpath, 'job.iterate.madx')
        file_out = os.path.join(outputpath, 'job.iterate.'+str(i)+'.madx')
        # TODO: do not use shell_command
        cmd = "".join(["sed 's/iter.madx/iter.", str(i), ".madx/g' ",
                      file_in, ' > ', file_out])
        _shell_command(cmd)
        madxFilePaths.append(file_out)

    LOG.debug("Running MADX parallel in {:d} processes".format(len(madxFilePaths)))
    process_pool.map(_callMadx, madxFilePaths)
    LOG.debug("Finished MADX parallel")

    LOG.debug("Cleaning temporary MADX files")
    for madxFilePathsItem in madxFilePaths:
        os.remove(madxFilePathsItem)
    for iterFilePathsItem in iterFilePaths:
        os.remove(iterFilePathsItem)

    return process_pool


def _shell_command(cmd):
    ret = os.system(cmd)
    if ret is not 0:
        raise ValueError("COMMAND: %s finished with exit value %i" % (cmd, ret))


def _loadtwiss_beta(varandpath):
    # TODO: Understand what's going on
    (var, path) = varandpath
    twissfile = os.path.join(path, "twiss." + var)
    x = 0
    try:
        x = tfs_pandas.read_tfs(twissfile)
        x = x.set_index('NAME').drop_duplicates()
        x['Q1'] = x.Q1
        x['Q2'] = x.Q2
        os.remove(twissfile)
        LOG.debug("{q1:}, {q2:}".format(q1=x.Q1, q2=x.Q2))
    except IOError as e:
        LOG.error(e.message)
        return []
    return var, x


"""
================================= Main ===========================================
"""

def _start():

    with timeit(lambda t:
                LOG.debug("  Sorted in {:f}s".format(t))):
        generate_fullresponse(_parse_args())


if __name__ == '__main__':
    _start()
