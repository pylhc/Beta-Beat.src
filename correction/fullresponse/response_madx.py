"""
Provides a function to create the responses of beta, phase, dispersion, tune and coupling
via iterative madx calls.

The variables under investigation need to be provided as a list (which can be gotten from
accelerator class).

For now, the response matrix is stored in a 'pickled' file.

:author: Lukas Malina, Joschua Dilly, Jaime (...) Coello de Portugal
"""
import math
import multiprocessing
import os

import numpy as np
import pandas

import madx_wrapper
from twiss_optics.optics_class import TwissOptics
from utils import logging_tools
from tfs_files import tfs_pandas as tfs
from utils.contexts import timeit, suppress_warnings
from utils.iotools import create_dirs

LOG = logging_tools.get_logger(__name__)


# Full Response Mad-X ##########################################################


def generate_fullresponse(accel_inst, variable_categories,
                          delta_k=0.00002, num_proc=multiprocessing.cpu_count(),
                          temp_dir=None):
    """ Generate a dictionary containing response matrices for
        beta, phase, dispersion, tune and coupling and saves it to a file.

        Args:
            accel_inst : Accelerator Instance.
            variable_categories (list): Categories of the variables/knobs to use. (from .json)
            delta_k (float): delta K1L to be applied to quads for sensitivity matrix
            num_proc (int): Number of processes to use in parallel.
            temp_dir (str): temporary directory. If ``None``, uses folder of original_jobfile.
    """
    LOG.debug("Generating Fullresponse via Mad-X.")
    with timeit(lambda t: LOG.debug("  Total time generating fullresponse: {:f}s".format(t))):
        if not temp_dir:
            temp_dir = accel_inst.model_dir
        create_dirs(temp_dir)

        variables = accel_inst.get_variables(classes=variable_categories)
        if len(variables) == 0:
            raise ValueError("No variables found! Make sure your categories are valid!")

        # try:
        #     variables = variables.tolist()
        # except AttributeError:
        #     pass

        num_proc = num_proc if len(variables) > num_proc else len(variables)
        process_pool = multiprocessing.Pool(processes=num_proc)

        incr_dict = _generate_madx_jobs(accel_inst, variables,
                                        delta_k, num_proc, temp_dir)
        _call_madx(process_pool, temp_dir, num_proc)
        _clean_up(temp_dir, num_proc)

        var_to_twiss = _load_madx_results(variables, process_pool, incr_dict, temp_dir)
        fullresponse = _create_fullresponse_from_dict(var_to_twiss)

    return fullresponse


def _generate_madx_jobs(accel_inst, variables, delta_k, num_proc, temp_dir):
    """ Generates madx job-files """
    LOG.debug("Generating MADX jobfiles.")
    incr_dict = {'0': 0.0}
    vars_per_proc = int(math.ceil(float(len(variables)) / num_proc))

    madx_job = _get_madx_job(accel_inst)

    for proc_idx in range(num_proc):
        jobfile_path = _get_jobfiles(temp_dir, proc_idx)

        current_job = madx_job
        for i in range(vars_per_proc):
            var_idx = proc_idx * vars_per_proc + i
            if var_idx >= len(variables):
                break
            var = variables[var_idx]
            incr_dict[var] = delta_k
            current_job += "{var:s}={var:s}{delta:+f};\n".format(var=var, delta=delta_k)
            current_job += "twiss, file='{:s}';\n".format(os.path.join(temp_dir, "twiss." + var))
            current_job += "{var:s}={var:s}{delta:+f};\n\n".format(var=var, delta=-delta_k)

        if proc_idx == num_proc - 1:
            current_job += "twiss, file='{:s}';\n".format(os.path.join(temp_dir, "twiss.0"))

        with open(jobfile_path, "w") as jobfile:
            jobfile.write(current_job)
    return incr_dict


def _get_madx_job(accel_inst):
    job_content = accel_inst.get_basic_seq_job()
    job_content += (
        "select, flag=twiss, clear;\n"
        "select, flag=twiss, pattern='^BPM.*\.B{beam:d}$', "
        "column=NAME,S,BETX,ALFX,BETY,ALFY,DX,DY,DPX,DPY,X,Y,K1L,MUX,MUY,R11,R12,R21,R22;\n\n"
    ).format(beam=accel_inst.get_beam())
    return job_content


def _call_madx(process_pool, temp_dir, num_proc):
    """ Call madx in parallel """
    LOG.debug("Starting {:d} MAD-X jobs...".format(num_proc))
    madx_jobs = [_get_jobfiles(temp_dir, index) for index in range(num_proc)]
    process_pool.map(_launch_single_job, madx_jobs)
    LOG.debug("MAD-X jobs done.")


def _clean_up(temp_dir, num_proc):
    """ Merge Logfiles and clean temporary outputfiles """
    LOG.debug("Cleaning output and building log...")
    full_log = ""
    for index in range(num_proc):
        job_path = _get_jobfiles(temp_dir, index)
        log_path = job_path + ".log"
        with open(log_path, "r") as log_file:
            full_log += log_file.read()
        os.remove(log_path)
        os.remove(job_path)
    full_log_path = os.path.join(temp_dir, "response_madx_full.log")
    with open(full_log_path, "w") as full_log_file:
        full_log_file.write(full_log)


def _load_madx_results(variables, process_pool, incr_dict, temp_dir):
    """ Load the madx results in parallel and return var-tfs dictionary """
    LOG.debug("Loading Madx Results.")
    vars_and_paths = []
    for value in variables + ['0']:
        vars_and_paths.append((value, temp_dir))
    var_to_twiss = {}
    for var, tfs_data in process_pool.map(_load_and_remove_twiss, vars_and_paths):
        tfs_data['incr'] = incr_dict[var]
        var_to_twiss[var] = tfs_data
    return var_to_twiss


def _create_fullresponse_from_dict(var_to_twiss):
    """ Convert var-tfs dictionary to fullresponse dictionary """
    var_to_twiss = _add_coupling(var_to_twiss)
    resp = pandas.Panel.from_dict(var_to_twiss)
    resp = resp.transpose(2, 1, 0)
    # After transpose e.g: resp[NDX, bpm12l1.b1, kqt3]
    # The magnet called "0" is no change (nominal model)
    resp['NDX'] = resp.xs('DX', axis=0).div(np.sqrt(resp.xs('BETX', axis=0)), axis="index")
    resp['NDY'] = resp.xs('DY', axis=0).div(np.sqrt(resp.xs('BETY', axis=0)), axis="index")
    resp['BBX'] = resp.xs('BETX', axis=0).div(resp.loc['BETX', :, '0'], axis="index")
    resp['BBY'] = resp.xs('BETY', axis=0).div(resp.loc['BETY', :, '0'], axis="index")
    resp = resp.subtract(resp.xs('0', axis=2), axis=2)
    # Remove difference of nominal model with itself (bunch of zeros)
    resp.drop('0', axis=2, inplace=True)
    resp = resp.div(resp.loc['incr', :, :])

    with suppress_warnings(np.ComplexWarning):  # raised as everything is complex-type now
        df = {'MUX': resp.xs('MUX', axis=0).astype(np.float64),
              'MUY': resp.xs('MUY', axis=0).astype(np.float64),
              'BETX': resp.xs('BETX', axis=0).astype(np.float64),
              'BETY': resp.xs('BETY', axis=0).astype(np.float64),
              'BBX': resp.xs('BBX', axis=0).astype(np.float64),
              'BBY': resp.xs('BBY', axis=0).astype(np.float64),
              'DX': resp.xs('DX', axis=0).astype(np.float64),
              'DY': resp.xs('DY', axis=0).astype(np.float64),
              'NDX': resp.xs('NDX', axis=0).astype(np.float64),
              'NDY': resp.xs('NDY', axis=0).astype(np.float64),
              "F1001R": tfs.TfsDataFrame(resp.xs('1001', axis=0).apply(np.real).astype(np.float64)),
              "F1001I": tfs.TfsDataFrame(resp.xs('1001', axis=0).apply(np.imag).astype(np.float64)),
              "F1010R": tfs.TfsDataFrame(resp.xs('1010', axis=0).apply(np.real).astype(np.float64)),
              "F1010I": tfs.TfsDataFrame(resp.xs('1010', axis=0).apply(np.imag).astype(np.float64)),
              'Q': resp.loc[['Q1', 'Q2'], resp.major_axis[0], :].transpose().astype(np.float64),
              }
    return df


def _get_jobfiles(temp_dir, index):
    """ Return names for jobfile and iterfile according to index """
    jobfile_path = os.path.join(temp_dir, "job.iterate.{:d}.madx".format(index))
    return jobfile_path


def _launch_single_job(inputfile_path):
    """ Function for pool to start a single madx job """
    log_file = inputfile_path + ".log"
    madx_wrapper.resolve_and_run_file(inputfile_path, log_file=log_file)


def _load_and_remove_twiss(var_and_path):
    """ Function for pool to retrieve results """
    (var, path) = var_and_path
    twissfile = os.path.join(path, "twiss." + var)
    tfs_data = tfs.read_tfs(twissfile, index="NAME")
    tfs_data['Q1'] = tfs_data.Q1
    tfs_data['Q2'] = tfs_data.Q2
    os.remove(twissfile)
    return var, tfs_data


def _add_coupling(dict_of_tfs):
    """ Adds coupling to the tfs. QUICK FIX VIA LOOP!"""
    with timeit(lambda t: LOG.debug("  Time adding coupling: {:f}s".format(t))):
        for var in dict_of_tfs:
            twopt = TwissOptics(dict_of_tfs[var])
            cpl = twopt.get_coupling("cmatrix")
            dict_of_tfs[var]["1001"] = cpl["F1001"]
            dict_of_tfs[var]["1010"] = cpl["F1010"]
        return dict_of_tfs


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
