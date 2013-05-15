'''
Created on 19 Mar 2013

@author: vimaier

@version: 1.0.1

This module contains the class TfsFile which handles the output files of GetLLM.

Change history:
 - 1.0.1:
     Suppressed output of the script and saved it in the new attributes _out and _err
'''

import sys
import subprocess

class ScriptRunner(object):
    '''A ScriptRunner runs a python script with optional arguments.
    
    Construct an instance and pass the name of the script with path and optional
    arguments. Use the run() Method to start the script.
    '''


    def __init__(self, script_name, args_dict=None):
        '''Constructor
        
        Keyword arguments:
        script_name -- The path and name to the script which will be executed
        args_dict -- dictionary{string-->string} of arguments whereby the key
            is the prefix/name and the value is the corresponding argument.
            e.g.: {"-f":"/x/y/z.py", "-o":"./output/"}
            
        '''
        if args_dict is None:
            args_dict = {}
        
        self._script_name = script_name
        self._args_dict = args_dict
        
        # Output of the script
        self._out = None
        self._err = None
        
    def get_output(self):
        ''' Returns the output of the run script or None '''
        return self._out
    
    def get_error_output(self):
        ''' Returns the error output of the script or None '''
        return self._err
    
    
    def run_script(self):
        '''Runs the script which was passed to the constructor and returns the
        error code of the script.   
        No output will be printed if error code is 0. 
        Otherwise stdout and stderr will be printed.    
        '''
        
        call_command = [sys.executable, self._script_name]
        
        for key in self._args_dict:
            call_command.append(key+self._args_dict[key])
            
        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE)

        # wait for the process to terminate
        (self._out, self._err) = process.communicate()
        
        
        errcode = process.returncode
        
        if 0 != errcode:
            print "Printing output:-------------------------"
            print self._out
            print "Printing error output:-------------------"
            print self._err
        
        return errcode 
    
    