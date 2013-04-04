'''
Created on 19 Mar 2013

@author: vimaier

This module is a test for GetLLM.py. It compares the output files of the GETLLM_SCRIPT with 
expected files from GETLLM_SCRIPT_VALID and prints the result.

If CREATE_VALID_OUTPUT is true then new valid data will be created.
It's important when the test is running on different machines 
because different machines have also a small deviation in the
floating point computation. This deviation should not be a problem
for further processes of the output files.
To summarize it if you're running this script first time on a 
machine set CREATE_VALID_OUTPUT to True otherwise False.

Every directory in PATH_TO_TEST_DATA represents a run and contains the necessary files to run
GetLLM.py. The test contains the execution of each run(directory). In this process will new output 
files be created and compared.
RunValidator checks if a run(directory) is in accordance with the expected structure. 
'''
import os
import sys
import filecmp
import argparse
import unittest
# Add directory to the python search path - needed to run script from command line
# Otherwise ScriptRunner won't be found
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.abspath("../../../Python_Classes4MAD"))
import vimaier_utils.scriptrunner
import runvalidator

'''
# Decides whether the valid script will be run and produce output or not
CREATE_VALID_OUTPUT = True
# Path to original/valid GetLLM.py script
GETLLM_SCRIPT_VALID = "./GetLLM_valid.py"
# Path to GetLLM.py script
GETLLM_SCRIPT = "../GetLLM.py"
PATH_TO_TEST_DATA = "./data"
'''

#===================================================================================================
# Parse argument
#===================================================================================================
description = ("     filecheck.py [options]\n"
               "Examples:"
               "\filecheck.py -o 1"
               "\filecheck.py -o 0 -v ../main/GetLLM_valid.py -m ./GetLLM_valid.py -p ./data"
               "\n\nThis module is a test for GetLLM.py. \n\nIt compares the output files of " 
               "the GETLLM_SCRIPT with expected files from GETLLM_SCRIPT_VALID and prints "
               "the result.\n\n")

parser = argparse.ArgumentParser(description)
parser.add_argument("-o","--valid_output", dest="CREATE_VALID_OUTPUT",
                    type=int, default=0,
                    help="Decides whether the valid script will be run(==0) and produce output or not(!=0)")
parser.add_argument("-v","--valid_getllm_script", dest="GETLLM_SCRIPT_VALID",
                    default="./GetLLM_valid.py",
                    help="Path to original/valid GetLLM.py script")
parser.add_argument("-m","--modified_getllm_script", dest="GETLLM_SCRIPT",
                    default="../GetLLM.py",
                    help="Path to the modified GetLLM.py script")
parser.add_argument("-p","--path_to_test_data", dest="PATH_TO_TEST_DATA",
                    default="./data",
                    help="Path to the root of the test data directory")

args = parser.parse_args()




#===================================================================================================
# # TestFileOutputGetLLM
#===================================================================================================
class TestFileOutputGetLLM(unittest.TestCase):
    """TestCase which compares the output of GetLLM.py"""

    def setUp(self):
        pass


    def tearDown(self):
        pass

    
    def test_file_output(self):
        """ Tests the output files of the GetLLM.py script."""
        all_tests_valid = True
        num_of_runs = 0
        num_of_valid_runs = 0
        # Run test for every directory in PATH_TO_TEST_DATA
        for element in os.listdir(args.PATH_TO_TEST_DATA):
            run_path = os.path.join(args.PATH_TO_TEST_DATA, element)
            if os.path.isdir(run_path):
                run_validator = runvalidator.RunValidator(run_path)
                validation_msg = run_validator.validate()
                if not "" == validation_msg:
                    print "Run cancelled for: "+run_path+"\tReason: "+validation_msg
                    continue
                # Valid directory structure to run GetLLM.py
                if not self.run_single_test(run_validator):
                    all_tests_valid = False
                    num_of_valid_runs -= 1
                    
                num_of_valid_runs += 1
                num_of_runs += 1
                
        print "Valid runs: %s/%s \n" % (str(num_of_valid_runs), str(num_of_runs))        
        if not all_tests_valid:
            raise AssertionError()
    

    #====================================================================
    # Helper methods
    #====================================================================

    def run_single_test(self, run_validator):
        """Runs the script with one single run directory
        
        :Parameters:
            'run_validator': RunValidator
                A valid RunValidator object.
                
        :Return: boolean
            True if the test run successfully otherwise False.     
        """
        # Run script with corresponding arguments
        dict_args = {"-f":run_validator.get_names_of_src_files(), 
                     "-m":run_validator.get_model_name(), 
                     "-a":run_validator.get_accelerator_type(),
                     "-o":run_validator.get_valid_output_path()}
        
        if args.CREATE_VALID_OUTPUT:
            # Run original/valid script
            valid_script_runner = vimaier_utils.scriptrunner.ScriptRunner(
                                                args.GETLLM_SCRIPT_VALID, dict_args)
            print "Starting GetLLM_valid.py("+run_validator.get_run_path()+")..."
            print 
            valid_script_runner.run_script()
            print 
            print "GetLLM_valid.py("+run_validator.get_run_path()+") finished...\n"
        
        dict_args["-o"] = run_validator.get_to_check_output_path()
        
        script_runner = vimaier_utils.scriptrunner.ScriptRunner(args.GETLLM_SCRIPT, dict_args)
        print "Starting GetLLM.py("+run_validator.get_run_path()+")..."
        print 
        script_runner.run_script()
        print 
        print "GetLLM.py finished("+run_validator.get_run_path()+")..."
        
        # Check output of the directory. Therefore:
        # Read filenames of to_check_output_path
        to_check_filenames = []
        for dirname_dirnames_outfilenames in os.walk(run_validator.get_to_check_output_path()):
            to_check_filenames = dirname_dirnames_outfilenames[2]
            break #there are no subdirectories
        
        # Downward compatibility for Python interpreter < 2.7
        #self.assertGreater(len(to_check_filenames), 0, 
        #                   "No to_check output files in "+run_validator.get_to_check_output_path())
        self.assertNotEqual(len(to_check_filenames), 0, 
                           "No to_check output files in "+run_validator.get_to_check_output_path())
        
        # Read filenames of valid_output_path
        valid_filenames = []
        for dirname_dirnames_validfilenames in os.walk(run_validator.get_valid_output_path()):
            valid_filenames = dirname_dirnames_validfilenames[2]
            break #there are no subdirectories
        
        # Downward compatibility for Python interpreter < 2.7
        #self.assertGreater(len(valid_filenames), 0, 
        #                   "No valid output files in "+run_validator.get_valid_output_path())
        self.assertNotEqual(len(valid_filenames), 0, 
                           "No valid output files in "+run_validator.get_valid_output_path())
        
        # Compare each valid file with corresponding to_check file
        print "Checking output files for run: "+ run_validator.get_run_path()
        
        match_mismatch_error = filecmp.cmpfiles(run_validator.get_valid_output_path(),
                                                run_validator.get_to_check_output_path(),
                                                valid_filenames, 
                                                shallow=False)
        
        equal_files = len(match_mismatch_error[0])
        overall_files = len(valid_filenames)
        
        print str(equal_files)+" of "+str(overall_files)+ " are equal\n"
        
        if equal_files == overall_files:
            return True
        else:
            print "Following files are not matching([mismatches], [errors]):"
            print match_mismatch_error[1:], "\n"
            return False
    # END run_single_test() --------------------------------------------
        
# END class TestFileOutPutGetLLM -------------------------------------------------------------------

def main():
    # Remove arguments from sys.argv
    # Needed to avoid conflicts with the options from unittest module
    arguments_tpl = ('-o', "--create_valid_output",
                   "-v","--valid_getllm_script",
                   "-m","--modified_getllm_script",
                   "-p","--path_to_test_data")
    del_lst = []
    for i,option in enumerate(sys.argv):
        if option in arguments_tpl:
            del_lst.append(i)
            del_lst.append(i+1)

    del_lst.reverse()
    for i in del_lst:
        del sys.argv[i]
        
    
    unittest.main()
    


if __name__ == "__main__":
    main()