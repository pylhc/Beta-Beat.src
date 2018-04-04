"""
Similar to Sequence Parser but with MADX

First: Set all variables to 0
Then: Set one variable at a time to 1

Compare results with case all==0

"""
import cPickle as pickle
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


def evaluate_for_variables(accel_inst, variable_categories, step=1e-5, order=2,
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

        try:
            LOG.debug("Step-size is {:e}.".format(step))
            _generate_madx_jobs(accel_inst, variables, order, num_proc, temp_dir)
            _call_madx(process_pool, temp_dir, num_proc)
            mapping = _load_madx_results(variables, step, order, process_pool, temp_dir)
        except IOError:
                raise IOError("MADX was unable to compute the mapping.")
        finally:
            # _clean_up(variables, temp_dir, num_proc)
            pass
    return mapping


def _get_current_elements(original_content, patterns, folder):
    """ Run madx once to get a list of elements """
    twiss_path = os.path.join(folder, "twiss.get_elements.dat")
    content = original_content.replace(patterns["job_content"],
             "twiss, file = '{:s}';\n".format(twiss_path))
    madx_wrapper.resolve_and_run_string(content)
    twiss = tfs_pandas.read_tfs(twiss_path, index="NAME")
    return twiss.index.values


def _create_macro(name, elements, order, folder):
    """ Create a madx-macro to write out k-values in a twiss-file """
    orders = _get_orders(order)
    macro = "! Macro to build a twiss-like table without calling twiss.\n\n"
    macro += "{:s}(file_out): macro = {{\n".format(name)
    macro += "    create, table=kval, columns=NAME, L, {:s};\n".format(", ".join(orders))
    macro += "    option, -echo;\n\n"
    for elem in elements:
        macro += "    NAME = '{name:s}';\n".format(name=elem)
        macro += "    L = {name:s}->L;\n".format(name=elem)
        for ord_str in orders:
            macro += "    {ord:s} = {name:s}->{ord:s};\n".format(name=elem, ord=ord_str)
        macro += "    fill, table=kval;\n"
        macro += "\n"

    macro += "    option, echo;\n"
    macro += "    select, flag=table, full;\n"
    macro += "    write, table=kval, file=file_out;\n"
    macro += "};\n\n"

    macro_path = _get_macrofile(folder)
    with open(macro_path, "w") as macro_file:
        macro_file.write(macro)

    return "option, -echo;\ncall, file = '{:s}';\noption, echo;\n\n".format(macro_path)


def _generate_madx_jobs(accel_inst, variables, order, num_proc, temp_dir):
    """ Generates madx job-files """
    def _add_to_var(var, value):
        return "{var:s} = {var:s} {value:+e};\n".format(var=var, value=value)

    def _assign(var, value):
        return "{var:s} = {value:d};\n".format(var=var, value=value)

    def _twiss_out(var):
        return "twiss, file = '{:s}';\n".format(_get_twissfile(temp_dir, var))

    def _do_macro(macro, var):
        return "exec, {:s}('{:s}');\n".format(macro, _get_twissfile(temp_dir, var))

    LOG.debug("Generating MADX jobfiles.")
    vars_per_proc = int(math.ceil(float(len(variables)) / num_proc))

    # load template
    madx_script = accel_inst.get_basic_seq_job()



    # create a macro to write out twiss-like table
    macro_name = "write_table_tfs"
    elements = _get_current_elements(original_content, patterns, temp_dir)
    macro_call = _create_macro(macro_name, elements, order, temp_dir)

    # zero all vars
    all_var_zero = "".join([_assign(var, 0) for var in variables])

    # build content for testing each variable
    for proc_idx in range(num_proc):
        job_content = "\n" + macro_call + "\n" + all_var_zero + "\n"

        for i in range(vars_per_proc):
            try:
                # var to be tested
                current_var = variables[proc_idx * vars_per_proc + i]
            except IndexError:
                break
            else:
                job_content += _assign(current_var, 1)
                job_content += _do_macro(macro_name, current_var)
                job_content += _assign(current_var, 0)
                job_content += "\n"

        # last thing to do: get baseline
        if proc_idx+1 == num_proc:
            job_content += _do_macro(macro_name, "0")

        full_content = original_content.replace(patterns["job_content"], job_content)
        with open(_get_jobfiles(temp_dir, proc_idx), "w") as job_file:
            job_file.write(full_content)


def _call_madx(process_pool, temp_dir, num_proc):
    """ Call madx in parallel """
    LOG.debug("Starting {:d} MAD-X jobs...".format(num_proc))
    madx_jobs = [_get_jobfiles(temp_dir, index) for index in range(num_proc)]
    process_pool.map(_launch_single_job, madx_jobs)
    LOG.debug("MAD-X jobs done.")


def _clean_up(variables, temp_dir, num_proc):
    """ Merge Logfiles and clean temporary outputfiles """
    LOG.debug("Cleaning output and printing log...")
    for var in (variables + ["0"]):
        try:
            os.remove(_get_twissfile(temp_dir, var))
        except OSError:
            LOG.info("MADX could not build twiss for '{:s}'".format(var))

    full_log = ""
    for index in range(num_proc):
        job_path = _get_jobfiles(temp_dir, index)
        log_path = job_path + ".log"
        with open(log_path, "r") as log_file:
            full_log += log_file.read()
        os.remove(log_path)
        os.remove(job_path)
    LOG.debug(full_log)

    try:
        os.rmdir(temp_dir)
    except OSError:
        pass


def _load_madx_results(variables, step, order, process_pool, temp_dir):
    """ Load the madx results in parallel and return var-tfs dictionary """
    LOG.debug("Loading Madx Results.")
    path_and_vars = []
    for value in variables:
        path_and_vars.append((temp_dir, value))

    _, base_tfs = _load_and_remove_twiss((temp_dir, "0"))
    mapping = dict([(o, {}) for o in _get_orders(order)] +
                   [(o + "L", {}) for o in _get_orders(order)])
    for var, tfs_data in process_pool.map(_load_and_remove_twiss, path_and_vars):
        for o in _get_orders(order):
            diff = (tfs_data[o] - base_tfs[o])
            mask = diff != 0  # drop zeros, maybe abs(diff) < eps ?
            k_list = diff.loc[mask]
            mapping[o][var] = k_list
            mapping[o + "L"][var] = k_list.mul(base_tfs.loc[mask, "L"])
    return mapping


# Helper #####################################################################


def _get_orders(max_order):
    return ["K{:d}{:s}".format(i, s) for i in range(max_order) for s in ["", "S"]]


def _get_jobfiles(folder, index):
    """ Return names for jobfile and iterfile according to index """
    jobfile_path = os.path.join(folder, "job.varmap.{:d}.madx".format(index))
    return jobfile_path


def _get_twissfile(folder, var):
    """ Return name of the variable-specific twiss file """
    return os.path.join(folder, "twiss." + var)


def _get_macrofile(folder):
    """ Returns the name of the macro """
    return os.path.join(folder, "macro.write_k_table.madx")


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
    return var, tfs_data


# Wrapper ##################################################################


def check_varmap_file(accel_inst, vars_categories):
    """ Checks on varmap file and creates it if not in model folder.
    THIS SHOULD BE REPLACED WITH A CALL TO JAIMES DATABASE, IF IT BECOMES AVAILABLE """
    if accel_inst.optics_file is None:
        raise ValueError("Optics not defined. Please provide modifiers.madx. "
                         "Otherwise MADX evaluation might be unstable.")

    varmapfile_name = "{:s}b{:d}_".format(accel_inst.NAME.lower(), accel_inst.get_beam())
    varmapfile_name += "_".join(sorted(set(vars_categories)))

    varmap_path = os.path.join(accel_inst.model_dir, varmapfile_name + "." + EXT)
    if not os.path.isfile(varmap_path):
        LOG.info("Variable mapping '{:s}' not found. Evaluating it via madx.".format(varmap_path))
        mapping = evaluate_for_variables(accel_inst, vars_categories)
        with open(varmap_path, 'wb') as dump_file:
            pickle.Pickler(dump_file, -1).dump(mapping)

    return varmap_path


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
