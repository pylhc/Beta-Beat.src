"""
Similar to Sequence Parser but with MADX

First: Set all variables to 0
Then: Set one variable at a time to 1

Compare results with case all==0

"""

import madx_wrapper
import math
import multiprocessing
import os

import numpy as np
import pandas

import madx_wrapper
from twiss_optics.optics_class import TwissOptics
from utils import logging_tools
from utils import tfs_pandas
from utils.contexts import timeit, suppress_warnings
from utils.iotools import create_dirs
from correction.fullresponse.response_madx import DEFAULT_PATTERNS

LOG = logging_tools.get_logger(__name__)


# Read Sequence ##############################################################


def evaluate_for_variables(variables, original_jobfile_path, order=4,
                           patterns=DEFAULT_PATTERNS, num_proc=multiprocessing.cpu_count(),
                           temp_dir=None):
    """ Generate a dictionary containing response matrices for
        beta, phase, dispersion, tune and coupling and saves it to a file.

        Args:
            variables (list): List of variables to use.
            original_jobfile_path (str): Name of the original MAD-X job file
                                         defining the sequence file.
            patterns (str): Patterns to be replaced in the MAD-X job file by the iterative
                           script calls. Must contain 'file' and 'twiss_columns'.
            num_proc (int): Number of processes to use in parallel.
            temp_dir (str): temporary directory. If ``None``, uses folder of original_jobfile.
    """
    LOG.debug("Generating Fullresponse via Mad-X.")
    with timeit(lambda t: LOG.debug("  Total time generating fullresponse: {:f}s".format(t))):
        if not temp_dir:
            temp_dir = os.path.dirname(original_jobfile_path)
        create_dirs(temp_dir)

        process_pool = multiprocessing.Pool(processes=num_proc)

        _generate_madx_jobs(variables, original_jobfile_path, order, patterns, num_proc, temp_dir)
        _call_madx(process_pool, temp_dir, num_proc)
        _clean_up(temp_dir, num_proc)

        mapping = _load_madx_results(variables, order, process_pool, temp_dir)
    return mapping


def _generate_madx_jobs(variables, original_jobfile_path, order, patterns, num_proc, temp_dir):
    """ Generates madx job-files """
    def _assign(var, value):
        return "{var:s}={value:d};\n".format(var=var, value=value)

    def _twiss_out(var):
        "twiss, file='{:s}';\n".format(_get_twissfile(temp_dir, var))

    LOG.debug("Generating MADX jobfiles.")
    vars_per_proc = int(math.ceil(float(len(variables)) / num_proc))

    # load template and add columns
    with open(original_jobfile_path, "r") as original_file:
        original_content = original_file.read()
    original_content.replace(
        patterns['twiss_columns'], "NAME,S," + ",".join(_get_orders(order))
    ).replace(
        patterns["element_pattern"], ""  # all elements for beam
    )

    # zero all vars
    all_var_zero = "".join(_assign(var, 0) for var in variables)

    # build content for testing each variable
    for proc_idx in range(num_proc):
        job_content = all_var_zero

        for i in range(vars_per_proc):
            try:
                # var to be tested
                current_var = variables[proc_idx * vars_per_proc + i]
            except IndexError:
                # last thing to do: get baseline
                job_content += _twiss_out(current_var)
                break
            else:
                job_content += _assign(current_var, 1)
                job_content += _twiss_out(current_var)
                job_content += _assign(current_var, 0)

        with open(_get_jobfiles(temp_dir, proc_idx), "w") as job_file:
            job_file.write(original_content.replace(patterns["file"], job_content))


def _call_madx(process_pool, temp_dir, num_proc):
    """ Call madx in parallel """
    LOG.debug("Starting {:d} MAD-X jobs...".format(num_proc))
    madx_jobs = [_get_jobfiles(temp_dir, index) for index in range(num_proc)]
    process_pool.map(_launch_single_job, madx_jobs)
    LOG.debug("MAD-X jobs done.")


def _clean_up(temp_dir, num_proc):
    """ Merge Logfiles and clean temporary outputfiles """
    LOG.debug("Cleaning output and printing log...")
    full_log = ""
    for index in range(num_proc):
        job_path = _get_jobfiles(temp_dir, index)
        log_path = job_path + ".log"
        with open(log_path, "r") as log_file:
            full_log += log_file.read()
        # os.remove(log_path)
        # os.remove(job_path)
    LOG.debug(full_log)


def _load_madx_results(variables, order, process_pool, temp_dir):
    """ Load the madx results in parallel and return var-tfs dictionary """
    LOG.debug("Loading Madx Results.")
    path_and_vars = []
    for value in variables:
        path_and_vars.append((temp_dir, value))

    _, base_tfs = _load_and_remove_twiss((temp_dir, "0"))
    mapping = dict([(o, {}) for o in _get_orders(order)])
    for var, tfs_data in process_pool.map(_load_and_remove_twiss, path_and_vars):
        for order in mapping:
            mapping[order][var] = tfs_data[order] - base_tfs[order]
    return mapping


# Helper #####################################################################


def _get_orders(max_order):
    return ["K{:d}{:s}".format(i, s) for i in range(max_order) for s in ["", "S", "L", "SL"]]


def _get_jobfiles(folder, index):
    """ Return names for jobfile and iterfile according to index """
    jobfile_path = os.path.join(folder, "job.varmap.{:d}.madx".format(index))
    return jobfile_path


def _get_twissfile(folder, var):
    """ Return name of the variable-specific twiss file """
    return os.path.join(folder, "twiss." + var)


def _launch_single_job(inputfile_path):
    """ Function for pool to start a single madx job """
    log_file = inputfile_path + ".log"
    return madx_wrapper.resolve_and_run_file(inputfile_path, log_file=log_file)


def _load_and_remove_twiss(path_and_var):
    """ Function for pool to retrieve results """
    twissfile = _get_twissfile(*path_and_var)
    tfs_data = tfs_pandas.read_tfs(twissfile, index="NAME")
    tfs_data['Q1'] = tfs_data.Q1
    tfs_data['Q2'] = tfs_data.Q2
    # os.remove(twissfile)
    return var, tfs_data


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))