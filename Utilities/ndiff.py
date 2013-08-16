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

import Utilities.iotools

DEFAULT_CFG = "default.cfg"
IGNORE_TFS_HEADER_CFG = "ignore_tfs_header.cfg"



def compare_dirs_with_given_file_endings(dir1, dir2, file_endings_to_compare=None):
    """
    Compares also subdirectories recursively
    :Parameters:
        'dir1': string
            Path to directory a
        'dir2': string
            Path to directory b
        'file_endings_to_compare': list
            If list is empty, every file will be compared
            Example: ["gitignore", "out"]
                
        :Return: boolean
            Returns True if dirs are equal, otherwise false
            (exit_code, std_stream, err_stream)
    """
    if Utilities.iotools.no_dirs_exist(dir1, dir2):
        print >> sys.stderr, dir1, "or(and)", dir2, "do(es) not exist."
        return False
    if file_endings_to_compare is None:
        file_endings_to_compare = []
        
    dir1_items = sorted(os.listdir(dir1))
    dir2_items = sorted(os.listdir(dir2))
    
    if dir1_items != dir2_items:
        print >> sys.stderr, "Items in dirs are not equal:\n",dir1_items, "\n", dir2_items
        return False

    for item in dir1_items:
        item1 = os.path.join(dir1, item)
        item2 = os.path.join(dir2, item)
        if os.path.isdir(item1):
            if not compare_dirs_with_given_file_endings(item1, item2, file_endings_to_compare):
                return False
        else:
            if empty_list_or_ending_matches_list(item1, file_endings_to_compare):
                if not compare_tfs_files_and_ignore_header(item1, item2):
                    return False
    return True

def empty_list_or_ending_matches_list(file_str, file_endings_list):
    if 0 == len(file_endings_list):
        return True
    ending = (file_str.split("."))[-1]
    return ending in file_endings_list

def compare_tfs_files_and_ignore_header(file_a, file_b):
    """
    Returns true, if files are equal based on the ignore_tfs_header.cfg file.
    Otherwise false, and output from ndiff will be printed.
    """
    tfs_ignore_header_cfg = os.path.join(get_path_to_root_of_ndiff(), IGNORE_TFS_HEADER_CFG)
    ignore_whitespace_option = {"--blank":""}
    
    (exit_code, std_stream, err_stream) = run_ndiff(file_a, file_b, tfs_ignore_header_cfg, 
                                                    ignore_whitespace_option)
    
    if 0 == exit_code and ndiff_files_were_equal(std_stream):
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
    return "gwarn" not in ndiff_output_string