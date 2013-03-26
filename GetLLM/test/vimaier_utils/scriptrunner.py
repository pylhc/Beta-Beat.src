'''
Created on 19 Mar 2013

@author: vimaier
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
            e.g.: {"-f":"/x/y/z.py", "-o":"./output/"
        '''
        if args_dict is None:
            args_dict = {}
        
        self._script_name = script_name
        self._args_dict = args_dict
    
    
    def run_script(self):
        '''Runs the script which was passed to the constructor and returns the
        return code of the script.        
        '''
        
        call_command = [sys.executable, self._script_name]
        
        for key in self._args_dict:
            call_command.append(key+self._args_dict[key])
        
        return_code = subprocess.call(call_command)
        
        return return_code 