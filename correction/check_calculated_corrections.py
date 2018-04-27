import os
import re
import sys

import numpy as np

sys.path.append(os.path.abspath(os.path.join(__file__, os.pardir)))

import madx_wrapper

from utils.entrypoint import EntryPointParameters, entrypoint
from utils import logging_tools
from utils.plotting import plot_tfs
from model import manager
from utils import iotools
from utils import tfs_pandas
from correction import getdiff

LOG = logging_tools.get_logger(__name__)

# Constants and Parameters #############################################

RESULTS_DIR = "Results"  # name of the (temporary) results folder
BASE_ID = ".tmpbasefile"  # extension of the (temporary) base files


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
    params.add_parameter(
        flags="--autoscale",
        help="Scales the plot, so that this percentage of points is inside the picture.",
        type=float,
        name="auto_scale",
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
    LOG.info("Started 'check_calculated_corrections'.")
    # get accelerator class
    accel_cls = manager.get_accel_class(accel_opt)
    accel_inst = accel_cls(model_dir=opt.model_dir)
    if opt.optics_file is not None:
        accel_inst.optics_file = opt.optics_file

    if opt.corrections_dir is None:
        opt.corrections_dir = os.path.join(opt.meas_dir, "Corrections")
    logging_tools.add_module_handler(
        logging_tools.file_handler(
            os.path.join(opt.corrections_dir, "check_corrections.log")
        )
    )

    corrections = _get_all_corrections(opt.corrections_dir, opt.file_pattern)
    _call_madx(accel_inst, corrections)
    _get_diffs(corrections, opt.meas_dir, opt.file_pattern)
    _plot(corrections, opt.corrections_dir, opt.show_plots, opt.change_marker, opt.auto_scale)

    if opt.clean_up:
        _clean_up(opt.corrections_dir, corrections)


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


def _get_diffs(corrections, meas_dir, file_pattern):
    """ Creates the twissfiles before and after corrections. (Copied into Results folder) """
    for folder in corrections:
        dest = os.path.join(folder, RESULTS_DIR)
        _copy_files(meas_dir, dest, file_pattern)
        getdiff.getdiff(dest)


def _plot(corrections, source_dir, show_plots, change_marker, auto_scale):
    """ Create all plots for the standard parameters """
    data_files = ['bbx', 'bby', 'dx', 'dy', 'ndx']  # 'normal ' getdiff output

    column_map = {  # special cases
        'couple_1001r': {
            'meas': 'F1001re',
            'expect': 'F1001re_prediction',
            'error': '',
            'file': 'couple',
        },
        'couple_1001i': {
            'meas': 'F1001im',
            'expect': 'F1001im_prediction',
            'error': '',
            'file': 'couple',
        },
        'couple_1010r': {
            'meas': 'F1010re',
            'expect': 'F1010re_prediction',
            'error': '',
            'file': 'couple',
        },
        'couple_1010i': {
            'meas': 'F1010im',
            'expect': 'F1010im_prediction',
            'error': '',
            'file': 'couple',
        },
        'chromatic_coupling_r': {
            'meas': 'Cf1001r',
            'expect': 'Cf1001r_prediction',
            'error': 'Cf1001rERR',
            'file': 'chromatic_coupling',
        },
        'chromatic_coupling_i': {
            'meas': 'Cf1001i',
            'expect': 'Cf1001i_prediction',
            'error': 'Cf1001iERR',
            'file': 'chromatic_coupling',
        },
        'phasex': {
            'meas': 'DIFF',
            'expect': 'EXPECT',
            'error': 'ERROR',
            'file': 'phasex',
        },
        'phasey': {
            'meas': 'DIFF',
            'expect': 'EXPECT',
            'error': 'ERROR',
            'file': 'phasey',
        },
    }

    sort_correct = sorted(corrections.keys())
    legends = ["Measurement"] + [d.replace(source_dir + os.sep, "") for d in sort_correct]

    for data in data_files + column_map.keys():
        try:
            meas = column_map[data]['meas']
            expect = column_map[data]['expect']
            error = column_map[data]['error']
            filename = column_map[data]['file']
        except KeyError:
            meas = 'MEA'
            expect = 'EXPECT'
            error = 'ERROR'
            filename = data

        files_c = [os.path.join(folder, RESULTS_DIR, filename + ".out") for folder in sort_correct]

        try:
            file_base = _create_base_file(source_dir, files_c[0], meas, error, expect, data)
            data_paths = [file_base] + files_c

            _log_rms(data_paths, legends, expect)

            plot_tfs.plot(
                files=data_paths,
                y_cols=[expect],
                e_cols=[error],
                legends=legends,
                labels=[data],
                output=os.path.join(source_dir, data),
                no_show=not show_plots,
                change_marker=change_marker,
                auto_scale=auto_scale,
            )

        except IOError:
            LOG.info("Could not plot parameter '{:s}'. ".format(data) +
                     "Probably not calculated by GetLLM.")


def _clean_up(source_dir, corrections):
    """ Removes the results folders again """
    for file in os.listdir(source_dir):
        if file.endswith(BASE_ID):
            iotools.delete_item(os.path.join(source_dir, file))

    for folder in corrections:
        iotools.delete_item(os.path.join(folder, RESULTS_DIR))


# MADX-Related ###############################################################


def _call_madx(accel_inst, corrections):
    """ Create and call the madx jobs to apply the corrections """
    original_content = _get_madx_job(accel_inst)
    for dir_correct in sorted(corrections):
        dir_out = os.path.join(dir_correct, RESULTS_DIR)
        iotools.create_dirs(dir_out)
        job_content = original_content
        job_content += "twiss, file='{:s}';\n".format(os.path.join(dir_out,
                                                                   getdiff.TWISS_NOT_CORRECTED))
        for file in sorted(corrections[dir_correct]):
            job_content += "call, file='{:s}';\n".format(file)
        job_content += "twiss, file='{:s}';\n".format(os.path.join(dir_out,
                                                                   getdiff.TWISS_CORRECTED))

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


def _create_base_file(source_dir, source_file, meas, error, expect, outname):
    """ Copy Measurement into a base-file. """
    data = tfs_pandas.read_tfs(source_file)

    if error == "":
        new_data = data.loc[:, ["S", "NAME", meas]]
        new_data.columns = ["S", "NAME", expect]
    else:
        new_data = data.loc[:, ["S", "NAME", meas, error]]
        new_data.columns = ["S", "NAME", expect, error]

    path_out = os.path.join(source_dir, outname + BASE_ID)
    tfs_pandas.write_tfs(path_out, new_data)
    return path_out


def _copy_files(src, dst, ignore):
    """ Copies files only from src to dst directories """
    old_items = os.listdir(dst) + [getdiff.TWISS_CORRECTED, getdiff.TWISS_NOT_CORRECTED]
    for item in os.listdir(src):
        src_item = os.path.join(src, item)
        if os.path.isfile(src_item) and not re.search(ignore, item) and not item in old_items:
            iotools.copy_item(src_item, os.path.join(dst, item))


def _log_rms(files, legends, column_name):
    """ Calculate and print rms value into log """
    file_name = os.path.splitext(os.path.basename(files[0]))[0]
    LOG.info("Results for '{:s}':".format(file_name))
    LOG.info("  {:<20s} {:11s} {:11s}".format("", " RMS", " Mean"))
    for f, l in zip(files, legends):
        data = tfs_pandas.read_tfs(f)[column_name]
        LOG.info("  {:<20s} {:+11.4e} {:+11.4e}".format(l, _rms(data), np.mean(data)))


def _rms(data):
    return np.sqrt(np.mean(np.square(data)))


# Script Mode ################################################################


if __name__ == "__main__":
    main()

