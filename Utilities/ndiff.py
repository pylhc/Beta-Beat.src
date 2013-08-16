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