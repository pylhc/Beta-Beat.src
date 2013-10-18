'''
Created on 20 Mar 2013

@author: vimaier

@version: 1.0.1

This module includes the class RunValidator which is used in filecheck.py to check if the
given data directory for input and output are correct.

Change history:
 - 1.0.1: - Added attribute 'run_dir_name' and 'has_valid_output_files'

'''

import os

import Python_Classes4MAD.metaclass as metaclass

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
        self.__run_dir_name = self.determine_dirname()
        self.__model_name = ""
        self.__names_of_src_files = ""
        self.__accelerator_type = ""
        self.__valid_output_path = ""
        self.__has_valid_output_files = False
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
            return "No model file was found: "+self.__model_name

        self.__names_of_src_files = self._validate_src_files()
        if "" == self.__names_of_src_files:
            return "No src files from analysis were found in "+self.__run_path

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

        # Check if output files exist. If not filecheck will run GetLLM_valid.py.
        # Just a simple check. If more as two files(consider .gitignore) in the directory the test will pass.
        if 2 < len( os.listdir(self.__valid_output_path) ):
            self.__has_valid_output_files = True

        # Check path to 'to_check' output
        self.__to_check_output_path = os.path.join(self.__run_path, RELATIVE_PATH_TO_OUTPUT_TO_CHECK)
        if not os.path.isdir(self.__to_check_output_path):
            return "'to_check' output path doesn't exist: " + self.__to_check_output_path

        self.__is_valid = True

        return ""
    # END validate() --------------------------------------------------------


    def _validate_src_files(self):
        """ Valid if source dir contains at least one file with the tbt_dataending_tupel """
        result = set()
        path_src_file_directory = os.path.join(self.__run_path,
                                               RELATIVE_PATH_TO_SRC_FILES_DIRECTORY)
        for dirname_dirnames_filenames in os.walk(path_src_file_directory):
            src_filenames = dirname_dirnames_filenames[2]

            for filename in src_filenames:
                # Exclude files which ends with the algorithm suffix
                # Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA

                filename_without_ending = self._return_name_without_tbt_dataending(filename)
                if "" != filename_without_ending:
                    result.add(
                               os.path.abspath(
                                      os.path.join(path_src_file_directory, filename_without_ending)
                                                         )
                               )
            break # There should be no sub-directories

        return ",".join(result)

    def _return_name_without_tbt_dataending(self, filename):
        if "_linx" in filename:
            return filename.replace("_linx", "")
        elif "_liny" in filename:
            return filename.replace("_liny", "")
        elif "_svdx" in filename:
            return filename.replace("_svdx", "")
        elif "_svdy" in filename:
            return filename.replace("_svdy", "")
        elif "_hax" in filename:
            return filename.replace("_hax", "")
        elif "_hay"  in filename:
            return filename.replace("_hay", "")
        else:
            return ""


    def determine_dirname(self):
        path = self.get_run_path()
        path.replace("/", os.sep)
        path.replace("\\", os.sep)

        splitted_path = self.get_run_path().split(os.sep)

        index = -1

        # For case when the separator is at the end of the string
        while "" == splitted_path[index] :
            index -= 1

        return splitted_path[index]


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

    def get_run_dir_name(self):
        return self.__run_dir_name

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

    def has__no_valid_output_files(self):
        """Returns opposite of has_valid_ouput_files"""
        return not self.__has_valid_output_files

    def get_to_check_output_path(self):
        """Returns the attribute run_path"""
        assert self.__is_valid, "No valid run path: "+self.__run_path
        return self.__to_check_output_path

# END class RunValidator ---------------------------------------------------------------------------