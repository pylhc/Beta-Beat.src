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


def evaluate_for_variables(accel_inst, variable_categories, order=2,
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
            _generate_madx_jobs(accel_inst, variables, order, num_proc, temp_dir)
            _call_madx(process_pool, temp_dir, num_proc)
            mapping = _load_madx_results(variables, order, process_pool, temp_dir)
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


def _generate_madx_jobs(accel_inst, variables, order, num_proc, temp_dir):
    """ Generates madx job-files """
    def _assign(var, value):
        return "{var:s} = {value:d};\n".format(var=var, value=value)

    def _do_macro(var):
        return "exec, create_table({table:s}, {f_out:s});\n".format(
            table="table." + var,
            f_out=_get_tablefile(temp_dir, var),
        )

    LOG.debug("Generating MADX jobfiles.")
    vars_per_proc = int(math.ceil(float(len(variables)) / num_proc))

    # load template
    madx_script = _create_basic_job(accel_inst, order, variables)

    # build content for testing each variable
    for proc_idx in range(num_proc):
        job_content = madx_script % {"TEMPFILE": _get_surveyfile(temp_dir, proc_idx)}

        for i in range(vars_per_proc):
            try:
                # var to be tested
                current_var = variables[proc_idx * vars_per_proc + i]
            except IndexError:
                break
            else:
                job_content += _assign(current_var, 1)
                job_content += _do_macro(current_var)
                job_content += _assign(current_var, 0)
                job_content += "\n"

        # last thing to do: get baseline
        if proc_idx+1 == num_proc:
            job_content += _do_macro("0")

        with open(_get_jobfile(temp_dir, proc_idx), "w") as job_file:
            job_file.write(job_content)


def _call_madx(process_pool, temp_dir, num_proc):
    """ Call madx in parallel """
    LOG.debug("Starting {:d} MAD-X jobs...".format(num_proc))
    madx_jobs = [_get_jobfile(temp_dir, index) for index in range(num_proc)]
    process_pool.map(_launch_single_job, madx_jobs)
    LOG.debug("MAD-X jobs done.")


def _clean_up(variables, temp_dir, num_proc):
    """ Merge Logfiles and clean temporary outputfiles """
    LOG.debug("Cleaning output and printing log...")
    for var in (variables + ["0"]):
        try:
            os.remove(_get_tablefile(temp_dir, var))
        except OSError:
            LOG.info("MADX could not build table for '{:s}'".format(var))

    full_log = ""
    for index in range(num_proc):
        survey_path = _get_surveyfile(temp_dir, index)
        job_path = _get_jobfile(temp_dir, index)
        log_path = job_path + ".log"
        with open(log_path, "r") as log_file:
            full_log += log_file.read()
        os.remove(log_path)
        os.remove(job_path)
        os.remove(survey_path)
    LOG.debug(full_log)

    try:
        os.rmdir(temp_dir)
    except OSError:
        pass


def _load_madx_results(variables, order, process_pool, temp_dir):
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


def _get_jobfile(folder, index):
    """ Return names for jobfile and iterfile according to index """
    return os.path.join(folder, "job.varmap.{:d}.madx".format(index))


def _get_tablefile(folder, var):
    """ Return name of the variable-specific table file """
    return os.path.join(folder, "table." + var)


def _get_surveyfile(folder, index):
    """ Returns the name of the macro """
    return os.path.join(folder, "survey.{:d}.tmp".format(index))


def _launch_single_job(inputfile_path):
    """ Function for pool to start a single madx job """
    log_file = inputfile_path + ".log"
    return madx_wrapper.resolve_and_run_file(inputfile_path, log_file=log_file)


def _load_and_remove_twiss(path_and_var):
    """ Function for pool to retrieve results """
    path, var = path_and_var
    twissfile = _get_tablefile(path, var)
    tfs_data = tfs_pandas.read_tfs(twissfile, index="NAME")
    tfs_data['Q1'] = tfs_data.Q1
    tfs_data['Q2'] = tfs_data.Q2
    return var, tfs_data


def _create_basic_job(accel_inst, order, variables):
    """ Create the madx-job basics needed
        TEMPFILE need to be replaced in the returned string.
    """
    all_ks = _get_orders(order)
    # basic sequence creation
    job_content = accel_inst.get_basic_seq_job()

    # create a survey and save it to a temporary file
    job_content += (
        "select, flag=survey, clear;\n"
        "select, flag=survey, pattern='^[MB].*\.B{beam:d}$', COLUMN=NAME, L;\n"
        "survey, file='%(TEMPFILE)s';\n"
        "readmytable, file='%(TEMPFILE)s', table=mytable;\n"
        "n_elem = table(mytable, tablelength);\n"
        "\n"
    ).format(beam=accel_inst.get_beam())

    # create macro for assigning values to the order per element
    job_content += "assign_k_values(element) : macro = {\n"
    for k_val in all_ks:
        job_content += "    {k:s} = element->{k:s};\n".format(k=k_val)
    job_content += "};\n\n"

    # create macro for using the row index as variable (see madx userguide)
    job_content += (
        "create_row(tblname, rowidx) : macro = {\n"
        "    exec, assign_k_values(tabstring(tblname, name, rowidx));\n"
        "};\n\n"
    )

    # create macro to create the full table with loop over elements
    job_content += (
        "create_table(table_id, filename) : macro = {{\n"
        "    create, table=table_id, column=_name, L, {col:s};\n"
        "    i_elem = 0;\n"
        "    while (i_elem < n_elem) {{\n"
        "        i_elem = i_elem + 1;\n"
        "        setvars, table=table_id, row=i_elem;\n"
        "        exec, create_row(table_id, $i_elem);\n"
        "        fill,  table=table_id;\n"
        "    };\n"
        "    write, table=table_id, file='filename';\n"
        "};\n\n"
    ).format(col=",".join(all_ks))

    # set all variables to zero
    for var in variables:
        job_content += "{var:s} = 0;\n".format(var=var)

    job_content += "\n"
    return job_content


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
