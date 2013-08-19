'''
Created on 15 Aug 2013

@author: vimaier

Utilities.ndiff.py is a wrapper for the program ndiff which is written in c.
This script provides basically the convenience to run ndiff via python methods.

The binaries of ndiff are stored in Beta-Beat.src/binaries/ndiff .

'''
import sys
import os
import subprocess
import re

import Utilities.iotools

DEFAULT_CFG = "default.cfg"
IGNORE_TFS_HEADER_CFG = "ignore_tfs_header.cfg"



def compare_dirs_with_files_mathing_regex_list(dir1, dir2, regex_list=None):
    """
    Compares also subdirectories recursively
    :Parameters:
        'dir1': string
            Path to directory a
        'dir2': string
            Path to directory b
        'regex_list': list
            List of regular expression patterns. Only files which matches a pattern in the list will be compared.
            See also http://docs.python.org/2.6/library/re.html
            If list is empty, every file will be compared
            Example: ["^.*gitignore$", "^.*out$"]
                
        :Return: boolean
            Returns True if dirs are equal, otherwise false
            (exit_code, std_stream, err_stream)
    """
    if Utilities.iotools.no_dirs_exist(dir1, dir2):
        print >> sys.stderr, dir1, "or(and)", dir2, "do(es) not exist."
        return False
    if regex_list is None:
        regex_list = []
        
    dir1_items = sorted(os.listdir(dir1))
    
    for item in dir1_items:
        item1 = os.path.join(dir1, item)
        item2 = os.path.join(dir2, item)
        if os.path.isdir(item1):
            if not compare_dirs_with_files_mathing_regex_list(item1, item2, regex_list):
                return False
        else:
            if empty_list_or_str_matches_regex_list(item1, regex_list):
                if not compare_tfs_files_and_ignore_header(item1, item2):
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

def compare_tfs_files_and_ignore_header(file_a, file_b):
    """
    Returns true, if files are equal based on the ignore_tfs_header.cfg file.
    Otherwise false, and output from ndiff will be printed.
    """
    tfs_ignore_header_cfg = os.path.join(get_path_to_root_of_ndiff(), IGNORE_TFS_HEADER_CFG)
    ignore_whitespace_option = {"--blank":""}
    
    (exit_code, std_stream, err_stream) = run_ndiff(file_a, file_b, tfs_ignore_header_cfg, 
                                                    ignore_whitespace_option)
    
    if 0 == exit_code and ndiff_files_were_equal(err_stream):
        return True
    else:
        print std_stream
        print >> sys.stderr, err_stream
    
    


def run_ndiff(file_a, file_b, config_file="", options_dict=None ):
    """
    Runs ndiff in a new subprocess to compare given files (a and b) with the config_file and the
    options_dict.
    
    :Parameters:
        'file_a': string
            Path to file a
        'file_b': string
            Path to file b
        'config_file': string
            Path to config file
        'options_dict': dictionary option_letter --> argument
            A dictionary with the options.
            E.g.: a_dict = {"-d":"", "-a":"xyz"}
                
        :Return: tupel(int, string, string)
            Returns the exit code and the outputstreams of the subprocess.
            (exit_code, std_stream, err_stream)
    """    
    options_string = __parse_options_dict(options_dict)
    
    if "" == config_file:
        config_file = get_path_to_default_config_file()
    
    
    call_command = [get_os_dependent_path_to_ndiff(), options_string, file_a, file_b, config_file]
    
    process = subprocess.Popen(call_command,
                       stdout=subprocess.PIPE, 
                       stderr=subprocess.PIPE)

    # wait for the process to terminate
    (std_stream, err_stream) = process.communicate()
    exit_code = process.returncode
    
    return (exit_code, std_stream, err_stream)
    
    
    
    
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
    
    return os.path.join(path_to_ndiff_root,tail)


def get_path_to_root_of_ndiff():
    return os.path.join(Utilities.iotools.get_absolute_path_to_betabeat_root(), 
                                      "binaries", "ndiff")
    
    
def ndiff_files_were_equal(ndiff_output_string):
    """ Returns True if ndiff_output_string did not print a warning. """
    return "warng" not in ndiff_output_string