"""
Similar to Sequence Parser but with MADX

First: Set all variables to 0
Then: Set one variable at a time to 1

Compare results with case all==0

"""

import math
import multiprocessing
import os

import madx_wrapper
from utils import logging_tools
from utils import tfs_pandas
from utils.contexts import timeit
from utils.iotools import create_dirs

LOG = logging_tools.get_logger(__name__)

EXT = "varmap"  # Extension Standard


# Read Sequence ##############################################################


def evaluate_for_variables(variables, original_jobfile_path, patterns, step=1e-5, order=2,
                           num_proc=multiprocessing.cpu_count(),
                           temp_dir=None):
    """ Generate a dictionary containing response matrices for
        beta, phase, dispersion, tune and coupling and saves it to a file.

        Args:
            variables (list): List of variables to use.
            original_jobfile_path (str): Name of the original MAD-X job file
                                         defining the sequence file.
            patterns (dict): Patterns to be replaced in the MAD-X job file by the iterative
                             script calls. Must contain 'job_content', 'twiss_columns'
                             and 'element_pattern'.
            num_proc (int): Number of processes to use in parallel.
            temp_dir (str): temporary directory. If ``None``, uses folder of original_jobfile.
    """
    LOG.debug("Generating Fullresponse via Mad-X.")
    with timeit(lambda t: LOG.debug("  Total time generating fullresponse: {:f}s".format(t))):
        if not temp_dir:
            temp_dir = os.path.dirname(original_jobfile_path)
        create_dirs(temp_dir)

        num_proc = num_proc if len(variables) > num_proc else len(variables)
        process_pool = multiprocessing.Pool(processes=num_proc)

        _generate_madx_jobs(variables, original_jobfile_path, step,
                            order, patterns, num_proc, temp_dir)
        _call_madx(process_pool, temp_dir, num_proc)
        _clean_up(temp_dir, num_proc)

        mapping = _load_madx_results(variables, step, order, process_pool, temp_dir)
    return mapping


def _generate_madx_jobs(variables, original_jobfile_path, step,
                        order, patterns, num_proc, temp_dir):
    """ Generates madx job-files """
    def _add_to_var(var, value):
        return "{var:s}={var:s}{value:+e};\n".format(var=var, value=value)

    def _twiss_out(var):
        return "twiss, file='{:s}';\n".format(_get_twissfile(temp_dir, var))

    LOG.debug("Generating MADX jobfiles.")
    vars_per_proc = int(math.ceil(float(len(variables)) / num_proc))

    # load template and add columns
    with open(original_jobfile_path, "r") as original_file:
        original_content = original_file.read()
    original_content = original_content.replace(
        patterns['twiss_columns'], "NAME,S,L," + ",".join(_get_orders(order))
    ).replace(
        patterns["element_pattern"], ""  # all elements for beam
    )

    # zero all vars
    # all_var_zero = "".join([_assign(var, 0) for var in variables])
    all_var_zero = ""

    # build content for testing each variable
    for proc_idx in range(num_proc):
        job_content = all_var_zero

        for i in range(vars_per_proc):
            try:
                # var to be tested
                current_var = variables[proc_idx * vars_per_proc + i]
            except IndexError:
                break
            else:
                job_content += _add_to_var(current_var, step)
                job_content += _twiss_out(current_var)
                job_content += _add_to_var(current_var, -step)

        # last thing to do: get baseline
        if proc_idx+1 == num_proc:
            job_content += _twiss_out("0")

        full_content = original_content.replace(patterns["job_content"], job_content)
        with open(_get_jobfiles(temp_dir, proc_idx), "w") as job_file:
            job_file.write(full_content)


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


def _load_madx_results(variables, step, order, process_pool, temp_dir):
    """ Load the madx results in parallel and return var-tfs dictionary """
    LOG.debug("Loading Madx Results.")
    path_and_vars = []
    for value in variables:
        path_and_vars.append((temp_dir, value))

    _, base_tfs = _load_and_remove_twiss((temp_dir, "0"))
    mapping = dict([(o, {}) for o in _get_orders(order)] +
                   [(o[:-1], {}) for o in _get_orders(order)])
    for var, tfs_data in process_pool.map(_load_and_remove_twiss, path_and_vars):
        for o in _get_orders(order):
            diff = (tfs_data[o] - base_tfs[o]).div(step)
            mask = diff != 0  # drop zeros, maybe abs(diff) < eps ?
            kl_list = diff.loc[mask]
            mapping[o][var] = kl_list
            mapping[o[:-1]][var] = kl_list.div(base_tfs.loc[mask, "L"])
    return mapping


# Helper #####################################################################


def _get_orders(max_order):
    return ["K{:d}{:s}".format(i, s) for i in range(max_order) for s in ["L", "SL"]]


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
    path, var = path_and_var
    twissfile = _get_twissfile(path, var)
    tfs_data = tfs_pandas.read_tfs(twissfile, index="NAME")
    tfs_data['Q1'] = tfs_data.Q1
    tfs_data['Q2'] = tfs_data.Q2
    os.remove(twissfile)
    return var, tfs_data


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
