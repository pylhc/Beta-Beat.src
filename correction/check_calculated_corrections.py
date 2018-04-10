import os
import re
import sys

sys.path.append(os.path.abspath(os.path.join(__file__, os.pardir)))

import madx_wrapper

from utils.entrypoint import EntryPointParameters, entrypoint
from utils import logging_tools
from utils.plotting import plot_tfs
from model import manager
from utils import iotools
from correction import getdiff

LOG = logging_tools.get_logger(__name__)

# Constants and Parameters #############################################

RESULTS_DIR = "Results"  # name of the (temporary) results folder
DATA_FILES = ['bbx', 'bby', 'phasex', 'phasey', 'dx', 'dy', 'ndx', 'couple']  # getdiff output files


def get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--meas_dir",
        help="Path to the directory containing the measurement files.",
        name="meas_dir",
        required=True,
        type=str,
    )
    params.add_parameter(
        flags="--model_dir",
        help="Path to the model to use.",
        name="model_dir",
        required=True,
        type=str,
    )
    params.add_parameter(
        flags="--corrections_dir",
        help=("Path to the directory containing the correction files. "
              "Defaults to 'measurement_dir'/Corrections if not given."
              ),
        name="corrections_dir",
        type=str,
    )
    params.add_parameter(
        flags="--optics_file",
        help=("Path to the optics file to use. If not present will default to "
              "model_path/modifiers.madx, if such a file exists."),
        name="optics_file",
        type=str,
    )
    params.add_parameter(
        flags="--cleanup",
        help="Clean results-folder from corrections.",
        type=bool,
        default=True,
        name="clean_up",
    )
    params.add_parameter(
        flags="--filepattern",
        help="Filepattern to use to find correction files in subfolders.",
        name="file_pattern",
        type=str,
        default=r"^changeparameters.*?\.madx$",
    )
    params.add_parameter(
        flags="--show",
        help="Show plots.",
        action="store_true",
        name="show_plots",
    )
    params.add_parameter(
        flags="--changemarker",
        help="Changes marker for each line in the plot.",
        action="store_true",
        name="change_marker",
    )
    return params


# Main invocation ############################################################


@entrypoint(get_params())
def main(opt, accel_opt):
    """ Do all corrections in given correction-subfolders.

    Runs Mad-X with corrections applied. The files to use for corrections need to be in subfolders
    and it is assumed, that all corrections of one subfolder are applied at the same time.
    The files are found by the pattern in CHANGEFILE_PATTERN.

    This function also requires accelerator class options.

    Keyword Args:
        Required
        meas_dir (str): Path to the directory containing the measurement files.
                        **Flags**: --meas_dir
        model_dir (str): Path to the model to use.
                         **Flags**: --model_dir

        Optional
        clean_up (bool): Clean results-folder from corrections.
                         **Flags**: --cleanup
                         **Default**: ``True``
        corrections_dir (str): Path to the directory containing the correction files.
                               Defaults to 'measurement_dir'/Corrections if not given.
                               **Flags**: --corrections_dir
        file_pattern (str): Filepattern to use to find correction files in subfolders.
                            **Flags**: --filepattern
                            **Default**: ``^changeparameters.*?\.madx$``
        optics_file (str): Path to the optics file to use.
                           If not present will default to model_path/modifiers.madx,
                           if such a file exists.
                           **Flags**: --optics_file
        change_marker: Changes marker for each line in the plot.
               **Flags**: --changemarker
               **Action**: ``store_true``
        show_plots: Show plots.
            **Flags**: --show
            **Action**: ``store_true``
    """
    # get accelerator class
    accel_cls = manager.get_accel_class(accel_opt)
    accel_inst = accel_cls(model_dir=opt.model_dir)
    if opt.optics_file is not None:
        accel_inst.optics_file = opt.optics_file

    if opt.corrections_dir is None:
        opt.corrections_dir = os.path.join(opt.meas_dir, "Corrections")

    corrections = _get_all_corrections(opt.corrections_dir, opt.file_pattern)
    _call_madx(accel_inst, corrections)
    _get_diffs(corrections, opt.meas_dir)
    _plot(corrections, opt.corrections_dir, opt.show_plots, opt.change_marker)

    if opt.clean_up:
        _clean_up(corrections)


# Private Functions ##########################################################


def _get_all_corrections(source_dir, file_pattern):
    """ Returns a dict of all files found in subfolders of source_dir fitting file_pattern """
    sub_dirs = os.listdir(source_dir)
    corrections = {}
    for sub in sub_dirs:
        fullpath_dir = os.path.join(source_dir, sub)
        if os.path.isdir(fullpath_dir):
            corrections[fullpath_dir] = []

            files = os.listdir(fullpath_dir)
            for file in files:
                fullpath_file = os.path.join(fullpath_dir, file)
                if (os.path.isfile(fullpath_file)
                        and re.search(file_pattern, file)):
                    corrections[fullpath_dir].append(fullpath_file)
    return corrections


def _get_diffs(corrections, meas_dir):
    """ Creates the twissfiles before and after corrections. (Copied into Results folder) """
    for folder in corrections:
        dest = os.path.join(folder, RESULTS_DIR)
        _copy_files(meas_dir, dest)
        getdiff.getdiff(dest)


def _plot(corrections, source_dir, show_plots, change_marker):
    """ Create all plots for the standard parameters """
    legends = [d.replace(source_dir + os.sep, "") for d in corrections.keys()]

    data_files = DATA_FILES

    for data in data_files:
        if data == 'couple':
            meas = 'F1001W_prediction'
            error = ''
        else:
            meas = 'MEA'
            error = 'ERROR'

        data_paths = [os.path.join(folder, "Results", data + ".out") for folder in corrections]
        try:
            plot_tfs.plot(
                files=data_paths,
                optics_params=[meas],
                error_params=[error],
                legends=legends,
                labels=[data],
                output=os.path.join(source_dir, data),
                no_show=not show_plots,
                change_marker=change_marker,
            )
        except IOError:
            LOG.info("Could not plot parameter '{:s}'. ".format(data) +
                     "Probably not calculated by GetLLM.")


def _clean_up(corrections):
    """ Removes the results folders again """
    for folder in corrections:
        iotools.delete_item(os.path.join(folder, RESULTS_DIR))


# MADX-Related ###############################################################


def _call_madx(accel_inst, corrections):
    """ Create and call the madx jobs to apply the corrections """
    original_content = _get_madx_job(accel_inst)
    for dir_correct in corrections:
        dir_out = os.path.join(dir_correct, RESULTS_DIR)
        iotools.create_dirs(dir_out)
        job_content = original_content
        job_content += "twiss, file='{:s}';\n".format(os.path.join(dir_out, 'twiss_no.dat'))
        for file in corrections[dir_correct]:
            job_content += "call, file='{:s}';\n".format(file)
        job_content += "twiss, file='{:s}';\n".format(os.path.join(dir_out, 'twiss_cor.dat'))

        madx_wrapper.resolve_and_run_string(
            job_content,
            output_file=os.path.join(dir_out, "job.corrections.madx"),
            log_file=os.path.join(dir_out, "job.corrections.log"),
        )


def _get_madx_job(accel_inst):
    """ Creates the basic job-string. """
    job_content = accel_inst.get_basic_seq_job()
    job_content += (
        "select, flag=twiss, clear;\n"
        "select, flag=twiss, pattern='^BPM.*\.B{beam:d}$', "
        "column=NAME,S,BETX,ALFX,BETY,ALFY,DX,DY,DPX,DPY,X,Y,K1L,MUX,MUY,R11,R12,R21,R22;\n\n"
    ).format(beam=accel_inst.get_beam())
    return job_content


# Helper #####################################################################

def _copy_files(src, dst):
    """ Copies files only from src to dst directories """
    for item in os.listdir(src):
        src_item = os.path.join(src, item)
        if os.path.isfile(src_item):
            iotools.copy_item(src_item, os.path.join(dst, item))


# Script Mode ################################################################


if __name__ == "__main__":
    main()

