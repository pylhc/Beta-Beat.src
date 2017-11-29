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
                        required=True,
                        dest="accel",
                        type=str)
    parser.add_argument("-o", "--outfile",
                        help="Name of fullresponse file.",
                        dest="outfile",
                        required=True,
                        type=argparse.FileType('w'))
    parser.add_argument("-t", "--temp_dir",
                        help=("Directory to store the temporary MAD-X files, "
                              "it will default to the directory containing 'outfile'."),
                        dest="temp_dir",
                        default=None,
                        type=str)
    parser.add_argument("-v", "--variables",
                        help="Path to file containing variable lists",
                        dest="variable_filepath",
                        required=True,
                        type=argparse.FileType('r'))
    parser.add_argument("-k", "--deltak",
                        help="delta K1L to be applied to quads for sensitivity matrix",
                        default=0.00002,
                        dest="delta_k",
                        type=float)
    parser.add_argument("-p", "--pattern",
                        help="Pattern to be replaced in the MAD-X file by the iterative script calls.",
                        default="%ITER_FILE%",
                        dest="pattern",
                        type=str)
    parser.add_argument("-n", "--num_proc",
                        help="Number of processes to use in parallel.",
                        default=multiprocessing.cpu_count(),
                        dest="num_proc",
                        type=int)

    options = parser.parse_args()
    if not options.temp_dir:
        options.temp_dir = os.path.dirname(options.outfile)

    return options


"""
================================= FullResponse ===========================================
"""


def _generate_fullresponse(options):
    LOG.debug("Generating Fullresponse.")
    variables = _get_variables_from_file(options.variable_filepath)

    process_pool = multiprocessing.Pool(processes=options.num_proc)

    incr_dict = _generate_madx_jobs(variables, options)
    _call_madx(process_pool, options.temp_dir, options.num_proc)
    _clean_up(options.temp_dir, options.num_proc, options.outfile)

    fullresponse = _load_madx_results(variables, process_pool, incr_dict, options.temp_dir)
    _save_fullresponse(options.outfile, fullresponse)


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


def _generate_madx_jobs(variables, options):
    """ Generates madx job-files """
    LOG.debug("Generating MADX jobs.")
    incr_dict = {'0': 0.0}
    vars_per_proc = int(math.ceil(float(len(variables)) / options.num_proc))

    for proc_idx in range(options.num_proc):
        job_filepath, iter_filepath = _get_jobfiles(options.temp_dir, proc_idx)
        _write_jobfile(options.original_job, job_filepath, iter_filepath, options.pattern)
        with open(iter_filepath, "w") as iter_file:
            if proc_idx == 0:
                iter_file.write("twiss, file='{:s}';\n".format(os.path.join(options.temp_dir, "twiss.0")))
            for i in range(vars_per_proc):  # Loop over variables
                var_idx = proc_idx * vars_per_proc + i
                var = variables[var_idx]
                incr_dict[var] = options.delta_k
                iter_file.write("{var:s}={var:s}+({delta:f});\n".format(var=var, delta=options.delta_k))
                iter_file.write("twiss, file='{:s}';\n".format(os.path.join(options.temp_dir, "twiss." + var)))
                iter_file.write("{var:s}={var:s}-({delta:f});\n".format(var=var, delta=options.delta_k))

    return incr_dict


def _get_jobfiles(temp_dir, index):
    job_filepath = os.path.join(temp_dir, "job.iterate.{:d}.madx".format(index))
    iter_filepath = os.path.join(temp_dir, "iter.{:d}.madx".format(index))
    return job_filepath, iter_filepath


def _write_jobfile(original_job, job_filepath, iter_filepath, pattern):
    with open(original_job, "r") as original_file:
        original_str = original_file.read()
    with open(job_filepath, "w") as job_file:
        job_file.write(original_str.replace(
            pattern,
            "call, file={file};".format(iter_filepath),
        ))


def _call_madx(process_pool, temp_dir, num_proc):
    LOG.debug("Starting {} MAD-X jobs...".format(num_proc))
    madx_jobs = [_get_jobfiles(temp_dir, index)[0] for index in range(num_proc)]
    process_pool.map(_launch_single_job, madx_jobs)
    LOG.debug("MAD-X jobs done.")


def _clean_up(temp_dir, num_proc, outfile):
    LOG.debug("Cleaning output and building log...")
    full_log = ""
    for index in range(num_proc):
        job_path, iter_path = _get_jobfiles(temp_dir, index)
        log_path = job_path + ".log"
        with open(log_path, "r") as log_file:
            full_log += log_file.read()
        os.remove(log_path)
        os.remove(job_path)
        os.remove(iter_path)
    with open(outfile + ".log", "w") as full_log_file:
        full_log_file.write(full_log)


def _load_madx_results(variables, process_pool, incr_dict, temp_dir):
    LOG.debug("Loading Madx Results.")
    vars_and_paths = []
    for value in variables + ['0']:
        vars_and_paths.append((value, temp_dir))
    var_to_twiss = {}
    for var, tfs_data in process_pool.map(_load_and_remove_twiss, vars_and_paths):
        tfs_data['incr'] = incr_dict[var]
        var_to_twiss[var] = tfs_data

    resp = pandas.Panel.from_dict(var_to_twiss)
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
    return {'MUX': resp.xs('MUX', axis=0),
            'MUY': resp.xs('MUY', axis=0),
            'BBX': resp.xs('BBX', axis=0),
            'BBY': resp.xs('BBY', axis=0),
            'NDX': resp.xs('NDX', axis=0),
            'Q': resp.loc[['Q1', 'Q2'], :, resp.minor_axis[0]]
           }


def _save_fullresponse(outputfile, fullresponse):
    """ Dumping the FullResponse file """
    LOG.debug("Saving Fullresponse into file '{:s}'".format(outputfile))
    with open(outputfile, 'wb') as dump_file:
        cPickle.Pickler(dump_file, -1).dump(fullresponse)


def _join_with_output(output_root, *path_tokens):
    return os.path.join(output_root, *path_tokens)


def _launch_single_job(path_to_input_file):
    log_file = path_to_input_file + ".log"
    return madx_wrapper.resolve_and_run_file(path_to_input_file, log_file=log_file)


def _load_and_remove_twiss(varandpath):
    (var, path) = varandpath
    twissfile = os.path.join(path, "twiss." + var)
    tfs_data = tfs_pandas.read_tfs(twissfile)
    tfs_data = tfs_data.set_index('NAME').drop_duplicates()
    tfs_data['Q1'] = tfs_data.Q1
    tfs_data['Q2'] = tfs_data.Q2
    os.remove(twissfile)
    return var, tfs_data


if __name__ == '__main__':
    with timeit(lambda t:
                LOG.debug("  Sorted in {:f}s".format(t))):
        _generate_fullresponse(_parse_args())
