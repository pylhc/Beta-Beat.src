'''
Created on 20 Mar 2013

@author: vimaier

@version: 1.0.0

This module includes the class RunValidator which is used in filecheck.py to check if the 
given data directory for input and output are correct.

Change history:

'''

import os

import metaclass

# Relative paths to the run directory
RELATIVE_PATH_TO_MODEL = "input/model/twiss.dat"
RELATIVE_PATH_TO_SRC_FILES_DIRECTORY = "input/src_files"
RELATIVE_PATH_TO_ACCELERATOR_FILE = "input"
RELATIVE_PATH_TO_VALID_OUTPUT = "output/valid"
RELATIVE_PATH_TO_OUTPUT_TO_CHECK = "output/to_check"

#===================================================================================================
# RunValidator
#===================================================================================================
class RunValidator(object):
    '''
    Validates an input directory for FileCheckGetLLM and provides the 
    arguments for a run.
    '''


    def __init__(self, run_path):
        '''
        Constructor
        
        :Parameters:
            'run_path': string
                The path to the directory which has to be validated.
        '''
        self.__run_path = run_path # Root of the directory with the input and output files
        self.__model_name = ""
        self.__names_of_src_files = ""
        self.__accelerator_type = ""
        self.__valid_output_path = ""
        self.__to_check_output_path = ""
        self.__is_valid = False
        
    def validate(self):
        """ Validates the run_path and its content.
        
        Checks if the run_path is valid and collects thereby the necessary
        information to run GetLLM.py.
        In case of success it returns an empty String otherwise the error
        message.
            
        :Return: string
            an empty string "" if the run_path was validated successfully. 
            Otherwise a not empty string with the reason for not being valid.            
        """
        # Check model file
        self.__model_name = os.path.join(self.__run_path, RELATIVE_PATH_TO_MODEL)
        if not os.path.isfile(self.__model_name):
            return "No model file were found: "+self.__model_name
        
        # Check names of src files
        path_src_file_directory = os.path.join(self.__run_path, 
                                               RELATIVE_PATH_TO_SRC_FILES_DIRECTORY)
        for dirname_dirnames_filenames in os.walk(path_src_file_directory):
            src_filenames = dirname_dirnames_filenames[2]
            
            for filename in src_filenames:
                # Exclude files which ends with the algorithm suffix
                # Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA
                tbt_data_tupel = ("_linx", "_liny", "_svdx", "_svdy", "_hax", "_hay") 
                
                # Downward compatibility for Python interpreter < 2.5
                #if filename.endswith(tbt_data_tupel):
                #    continue
                exclude_file = False
                for tbt_data in tbt_data_tupel:
                    if filename.endswith(tbt_data):
                        exclude_file = True
                        break
                if exclude_file:
                    continue    
                
                
                self.__names_of_src_files += os.path.abspath(
                                          os.path.join(self.__run_path,
                                                       RELATIVE_PATH_TO_SRC_FILES_DIRECTORY,
                                                       filename)
                                                         )+","
            if "" == self.__names_of_src_files:
                return "No files from analysis were found in "+path_src_file_directory
            # Remove last ','   
            self.__names_of_src_files = self.__names_of_src_files[:-1] 
            
            break # There should be no sub-directories
        
        # Check accelerator type
        # The model file contains the type of the accelerator
        model_file = metaclass.twiss(self.__model_name)
        self.__accelerator_type = model_file.SEQUENCE
        if "" == self.__accelerator_type or self.__accelerator_type is None:
            return "No SEQUENCE(accelerator type) found in model: " + self.__model_name
        
        # Check path to valid output
        self.__valid_output_path = os.path.join(self.__run_path, RELATIVE_PATH_TO_VALID_OUTPUT)   
        if not os.path.isdir(self.__valid_output_path):
            return "Valid output path doesn't exist: " + self.__valid_output_path
        
        # Check path to 'to_check' output
        self.__to_check_output_path = os.path.join(self.__run_path, RELATIVE_PATH_TO_OUTPUT_TO_CHECK)
        if not os.path.isdir(self.__to_check_output_path):
            return "'to_check' output path doesn't exist: " + self.__to_check_output_path
            
        self.__is_valid = True
        
        return ""
    # END validate() --------------------------------------------------------
    
    def print_yourself(self):
        """Method for debugging purposes"""
        print "RunValidator\nrun_path: "+self.get_run_path()
        print "model_name: "+self.get_model_name()
        print "src_files: "+self.get_names_of_src_files()
        print "acc-type: "+self.get_accelerator_type()
        print "valid_out: "+self.get_valid_output_path()
        print "to_check_out: "+self.get_to_check_output_path()
        
    
    # GETTERS ---------------------------------------------------------------
    def get_run_path(self):
        """Returns the attribute run_path"""
        return self.__run_path
    
    def get_model_name(self):
        """Returns the attribute model_name"""
        assert self.__is_valid, "No valid run path: "+self.__run_path
        return self.__model_name
    
    def get_names_of_src_files(self):
        """Returns the attribute names_of_src_files"""
        assert self.__is_valid, "No valid run path: "+self.__run_path
        return self.__names_of_src_files
    
    def get_accelerator_type(self):
        """Returns the attribute accelerator_type"""
        assert self.__is_valid, "No valid run path: "+self.__run_path
        return self.__accelerator_type
    
    def get_valid_output_path(self):
        """Returns the attribute run_path"""
        assert self.__is_valid, "No valid run path: "+self.__run_path
        return self.__valid_output_path
    
    def get_to_check_output_path(self):
        """Returns the attribute run_path"""
        assert self.__is_valid, "No valid run path: "+self.__run_path
        return self.__to_check_output_path
        
# END class RunValidator ---------------------------------------------------------------------------