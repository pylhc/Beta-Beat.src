'''
.. module: Utilities.ndiff

Created on 15 Aug 2013

Utilities.ndiff.py is a wrapper for the program ndiff which is written in c.
This script provides basically the convenience to run ndiff via python methods.

The binaries of ndiff are stored in Beta-Beat.src/binaries/ndiff .

Usage::

    import Utilities.ndiff

    if Utilities.ndiff.compare_files_and_ignore_whitespace("file_a.out", "file_b.out"):
        print "Files are equal"
    else:
        print "Files are not equal"

.. moduleauthor:: vimaier
'''
import sys
import os
import subprocess
import re

import Utilities.iotools

DEFAULT_CFG = "default.cfg"
IGNORE_TFS_HEADER_CFG = "ignore_tfs_header.cfg"



def compare_dirs_with_files_matching_regex_list(dir1, dir2, regex_list=None, file_to_config_file_dict=None, master_config_file="", options_dict=None):
    """
    Compares also subdirectories recursively
    :param string dir1: Path to directory a
    :param string dir2: Path to directory b
    :param list regex_list: List of regular expression patterns. Only files which matches a pattern in the list will be compared.
                            If list is empty, every file will be compared
                            See also http://docs.python.org/2.6/library/re.html
                            Example: ["^.*gitignore$", "^.*out$"]
                            ["(?!^gplot_IP2$)", "(?!^plot_IP2.eps$)", "(?!^var4plot.sh$)", ] # Exclude these three files from comparing
    :param dict file_to_config_file_dict:
                            Keys are filenames(without path) and corresponding values are ndiff config filepaths.
                            If File is not in dict a default/master config file will be used.
    :param string master_config_file:
                            Path to config file which will be used for files which have no entry in
                            file_to_config_file_dict.
                            If not stated a default config file will be used.

    :returns: bool -- True if dirs are equal, otherwise false
    """
    if Utilities.iotools.no_dirs_exist(dir1, dir2):
        print >> sys.stderr, dir1, "or(and)", dir2, "do(es) not exist."
        return False
    if regex_list is None:
        regex_list = []
    if file_to_config_file_dict is None:
        file_to_config_file_dict = {}

    dir1_items = sorted(os.listdir(dir1))

    for item in dir1_items:
        item1 = os.path.join(dir1, item)
        item2 = os.path.join(dir2, item)
        if os.path.isdir(item1):
            if not compare_dirs_with_files_matching_regex_list(item1, item2, regex_list, file_to_config_file_dict, master_config_file, options_dict):
                return False
        else:
            if empty_list_or_str_matches_regex_list(item1, regex_list):
                if item in file_to_config_file_dict:
                    if not compare_files(item1, item2, file_to_config_file_dict[item], options_dict):
                        return False
                else:
                    if not compare_files(item1, item2, master_config_file, options_dict):
                        return False
    return True

def empty_list_or_str_matches_regex_list(file_str, regex_list):
    if 0 == len(regex_list):
        return True
    for re_pattern in regex_list:
        match = re.match(re_pattern, file_str)
        if not match is None:
            return True
    return False

def compare_files_and_ignore_whitespace(file_a, file_b):
    """
    Returns true, if files are equal with the default config. Whitespace will not be compared.
    Otherwise false, and output from ndiff will be printed.
    """
    ignore_whitespace_option = {"--blank":""}
    return compare_files(file_a, file_b, options_dict=ignore_whitespace_option)

def compare_tfs_files_and_ignore_header(file_a, file_b):
    """
    Returns true, if files are equal based on the ignore_tfs_header.cfg file.
    Otherwise false, and output from ndiff will be printed.
    """
    tfs_ignore_header_cfg = os.path.join(get_path_to_root_of_ndiff(), IGNORE_TFS_HEADER_CFG)
    ignore_whitespace_option = {"--blank":""}

    return compare_files(file_a, file_b, tfs_ignore_header_cfg, ignore_whitespace_option)


def compare_files(file_a, file_b, config_file="", options_dict=None ):
    """
    Returns true, if files are equal based on the given config file.
    Otherwise false, and output from ndiff will be printed.
    """
    (exit_code, std_stream, err_stream, call_command) = run_ndiff(file_a, file_b, config_file, options_dict)

    if 0 == exit_code and ndiff_files_were_equal(err_stream):
        return True
    else:
        print "ndiff Command:", call_command
        print std_stream
        print >> sys.stderr, err_stream
        return False


def run_ndiff(file_a, file_b, config_file="", options_dict=None ):
    """
    Runs ndiff in a new subprocess to compare given files (a and b) with the config_file and the
    options_dict.

    :param string file_a: Path to file a
    :param string file_b: Path to file b
    :param string config_file: Path to config file
    :param dict options_dict: option_letter --> argument -- A dictionary with the options.
                                E.g.: a_dict = {"-d":"", "-a":"xyz"}

    :returns: tupel(int, string, string) -- Returns the exit code and the outputstreams of the subprocess.
              (exit_code, std_stream, err_stream)
    """
    options_string = __parse_options_dict(options_dict)

    if "" == config_file:
        config_file = get_path_to_default_config_file()

    call_command = []
    call_command.append(get_os_dependent_path_to_ndiff())
    if "" != options_string:
        call_command.append(options_string)
    call_command += [file_a, file_b, config_file]

    process = subprocess.Popen(call_command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)

    # wait for the process to terminate
    (std_stream, err_stream) = process.communicate()
    exit_code = process.returncode

    return (exit_code, std_stream, err_stream, " ".join(call_command) )




def __parse_options_dict(options_dict):
    if options_dict is None:
        return ""

    result = []
    for key in options_dict:
        result.append(key)
        if "" != options_dict[key]:
            result.append(options_dict[key])
    return " ".join(result)


def get_path_to_default_config_file():
    """ Returns "x/Beta-Beat.src/binaries/ndiff/default.cfg" """
    return os.path.join(get_path_to_root_of_ndiff(), DEFAULT_CFG)


def get_os_dependent_path_to_ndiff():
    path_to_ndiff_root = get_path_to_root_of_ndiff()

    if "posix" == os.name:
        tail = os.path.join("linux","ndiff-linux64")
    elif "nt" == os.name:
        tail = os.path.join("windows","ndiff-win32.exe")
    else:
        raise OSError("ndiff only available for Linux and Windows. Your OS: "+os.name)

    return os.path.join(path_to_ndiff_root, tail)


def get_path_to_root_of_ndiff():
    return os.path.join(Utilities.iotools.get_absolute_path_to_betabeat_root(),
                                      "binaries", "ndiff")


def ndiff_files_were_equal(ndiff_output_string):
    """ Returns True if ndiff_output_string did not print a warning. """
    return "warng" not in ndiff_output_string
