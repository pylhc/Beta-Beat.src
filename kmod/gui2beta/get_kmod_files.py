"""
Module kmod.get_kmod_files
-----------------------------

Used by the BB-Gui to merge and copy kmod files from subfolders of the chosen folder into
the current optics output folder.
"""
import os
import sys
from os.path import abspath, join, dirname, pardir
new_path = abspath(join(dirname(abspath(__file__)), pardir, pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

from optparse import OptionParser
from tfs_files import tfs_pandas
import re
from utils import logging_tools

from gui2kmod import get_beta_star_filename, get_beta_filename

LOG = logging_tools.get_logger(__name__)


# Main Functions ###############################################################


def merge_and_copy_kmod_output(beam, kmod_dir, res_dir, mod_path):
    """ Merges the needed files into dataframes and saves them.

    Args:
        beam (str): Currently used beam.
        kmod_dir (str): Parent folder of the kmod results
        res_dir (str): Destination folder for the merged output
        mod_path (str): Path to the model.
    """
    LOG.debug("Merging and copying of kmod data started.")
    pattern = re.compile(".*R[0-9]\." + beam.upper())
    model_twiss = tfs_pandas.read_tfs(mod_path, index="NAME")

    ip_dir_names = [d for _, dirs, _ in os.walk(kmod_dir) for d in dirs if pattern.match(d)]

    # copy beta data
    for plane in "xy":
        new_data = tfs_pandas.TfsDataFrame()
        for ip_dir_name in ip_dir_names:
            src = join(kmod_dir, ip_dir_name, get_beta_filename(plane))
            data = _load_source_data(src, "NAME")
            if data is not None:
                new_data = new_data.append(data.loc[data.index.isin(model_twiss.index), :])
        new_data["S"] = model_twiss.loc[new_data.index, "S"]
        dst = join(res_dir, get_beta_merged_filename(plane))
        tfs_pandas.write_tfs(dst, new_data, save_index="NAME")

    # copy beta* data
    new_data = tfs_pandas.TfsDataFrame()
    for ip_dir_name in ip_dir_names:
        src = join(kmod_dir, ip_dir_name, get_beta_star_filename())
        data = _load_source_data(src)
        new_data = new_data.append(data)
    dst = join(res_dir, get_beta_star_merged_filename())
    tfs_pandas.write_tfs(dst, new_data)


def get_beta_star_merged_filename():
    """ Outputfilename of the merged betastar file. """
    return "getkmodbetastar.out"


def get_beta_merged_filename(plane):
    """ Outputfilename of the merged beta file. """
    return get_beta_filename(plane)


# Private Functions ############################################################


def _load_source_data(src, index=None):
    data = None
    try:
        data = tfs_pandas.read_tfs(src, index=index)
    except IOError:
        LOG.warn("Cannot find kmod data in '{:s}', won't copy those files.".format(src))
    else:
        LOG.warn("Loaded kmod data from '{:s}'.".format(src))
    return data


def parse_args():
    usage = 'Usage: %prog -w WORKING_DIR -o OUT_FILE_PATH [options]'
    parser = OptionParser(usage=usage)
    parser.add_option('-k', '--kmod_directory',
                            help='path to kmod directory with stored KMOD measurement files', 
                            action='store', type='string', dest='kmod_dir')
    parser.add_option('-r', '--results_dir',
                            help='Specify results directory of optics analysis', 
                            action='store', type='string', dest='res_dir')
    parser.add_option('-m', '--model_path',
                            help='Specify path to current model',
                            action='store', type='string', dest='mod_path')
    parser.add_option('-b', '--beam',
                            help='define beam used: b1 or b2',
                            action='store', type='string', dest='beam')
    (options, _) = parser.parse_args()
    return options


# Script Mode ##################################################################


if __name__ == '__main__':
    opt = parse_args()
    merge_and_copy_kmod_output(opt.beam, opt.kmod_dir, opt.res_dir, opt.mod_path)
