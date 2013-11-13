''' Created on 19 Mar 2013

:author: vimaier


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

Returns as exit code the number of failed runs.

Change history:
 - 1.0.1:
    Added option "special_output". It is possible to choose a separate output directory. The
    produced output will be deleted again.
 - 1.0.2 29th April 2013:
    Changed compare algorithm. Created compare_tfs_files(...) to check content independent from
    whitespace. The lines with the descriptors "GetLLMVersion", "MAD_FILE" and "FILES" will be
    ignored.
Use git to view history:
    git log <file>
'''

import os
import sys
import optparse
import unittest
import shutil

import __init__ # @UnusedImport init will include paths

import vimaier_utils.scriptrunner
import Utilities.iotools
import vimaier_utils.compare_utils
import runvalidator

#===================================================================================================
# Parse argument
#===================================================================================================
DESCRIPTION = ("     filecheck.py [options]\n"
               "Examples:"
               "\filecheck.py -o 1"
               "\filecheck.py -o 0 -v ../main/GetLLM_valid.py -m ./GetLLM_valid.py -p ./data"
               "\n\nThis module is a test for GetLLM.py. \n\nIt compares the output files of "
               "the GETLLM_SCRIPT with expected files from GETLLM_SCRIPT_VALID and prints "
               "the result.\n"
               "Use special_output if you want to create output files in a different directory. The"
               " produced files in this directory will be deleted again if  the test passes.\n\n")

parser = optparse.OptionParser(usage=DESCRIPTION)
parser.add_option("-o", "--valid_output",
                    help="Decides whether the valid script will be run(!=0) and produce output or not(==0)",
                    default=1, dest="CREATE_VALID_OUTPUT")
parser.add_option("-v", "--valid_getllm_script", dest="GETLLM_SCRIPT_VALID",
                    default="./GetLLM_valid.py",
                    help="Path to original/valid GetLLM.py script")
parser.add_option("-m", "--modified_getllm_script", dest="GETLLM_SCRIPT",
                    default="../GetLLM.py",
                    help="Path to the modified GetLLM.py script")
parser.add_option("-p", "--path_to_test_data", dest="PATH_TO_TEST_DATA",
                    default="./data",
                    help="Path to the root of the test data directory")
parser.add_option("-s", "--special_output", dest="SPECIAL_OUTPUT",
                    default="",
                    help="If special_output is given the output will be produced into this directory."+
                            " Valid output will also be produced.")
parser.add_option("-t", "--test_name", dest="TEST_NAME",
                    default=None,
                    help="If defined, only run the test with this name")

(options, args) = parser.parse_args()

#===================================================================================================
# Checking arguments
#===================================================================================================
# Decides whether the valid script will be run and produce output or not
try:
    CREATE_VALID_OUTPUT = bool( int(options.CREATE_VALID_OUTPUT) )
except ValueError:
    print "Wrong option 'valid_output': ", options.CREATE_VALID_OUTPUT
    print "Will produce valid output"
    CREATE_VALID_OUTPUT = True
# Path to original/valid GetLLM.py script
GETLLM_SCRIPT_VALID = options.GETLLM_SCRIPT_VALID
if Utilities.iotools.not_exists_directory(GETLLM_SCRIPT_VALID):
    GETLLM_SCRIPT_VALID = os.path.join(
                                     Utilities.iotools.get_absolute_path_to_betabeat_root(),
                                     "GetLLM", "test", "./GetLLM_valid.py"
                                     )
# Path to GetLLM.py script
GETLLM_SCRIPT = options.GETLLM_SCRIPT
if Utilities.iotools.not_exists_directory(GETLLM_SCRIPT):
    GETLLM_SCRIPT = os.path.join(
                                     Utilities.iotools.get_absolute_path_to_betabeat_root(),
                                     "GetLLM", "GetLLM.py"
                                     )
# Path to special output
SPECIAL_OUTPUT = options.SPECIAL_OUTPUT
if "" != SPECIAL_OUTPUT :
    CREATE_VALID_OUTPUT = True
    if not os.path.isdir(SPECIAL_OUTPUT):
        print "special_output is not a directory: ", SPECIAL_OUTPUT
        SPECIAL_OUTPUT = ""
PATH_TO_TEST_DATA = options.PATH_TO_TEST_DATA
if Utilities.iotools.not_exists_directory(PATH_TO_TEST_DATA):
    PATH_TO_TEST_DATA = os.path.join(
                                     Utilities.iotools.get_absolute_path_to_betabeat_root(),
                                     "GetLLM", "test", "data"
                                     )

# If defined, run only this test, otherwise run all tests
TEST_NAME = options.TEST_NAME

#===================================================================================================
# # TestFileOutputGetLLM
#===================================================================================================
class TestFileOutputGetLLM(unittest.TestCase):
    """TestCase which compares the output of GetLLM.py"""


    num_of_failed_tests = 0
    """Static variable which is used as exit code."""

    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_file_output(self):
        """ Tests the output files of the GetLLM.py script."""
        all_tests_valid = True
        num_of_runs = 0
        num_of_valid_runs = 0
        # if test name is given, only this test, otherwise all
        if options.TEST_NAME != None:
            test_dirs_to_run = [options.TEST_NAME]
        else:
            test_dirs_to_run = os.listdir(PATH_TO_TEST_DATA)
        # Run test for every directory in PATH_TO_TEST_DATA
        for element in test_dirs_to_run:
            run_path = os.path.join(PATH_TO_TEST_DATA, element)
            if os.path.isdir(run_path):
                run_validator = runvalidator.RunValidator(run_path)
                validation_msg = run_validator.validate()
                if not "" == validation_msg:
                    print "Run cancelled for: "+run_path+"\tReason: "+validation_msg
                    continue
                # Valid directory structure to run GetLLM.py

                print "==============================================================="
                print "Run test for", run_validator.get_run_dir_name()
                print "==============================================================="
                if not self.run_single_test(run_validator):
                    all_tests_valid = False
                    num_of_valid_runs -= 1

                num_of_valid_runs += 1
                num_of_runs += 1


        TestFileOutputGetLLM.num_of_failed_tests = num_of_runs - num_of_valid_runs

        print "==============================================================="
        print "Valid runs: %s/%s " % (str(num_of_valid_runs), str(num_of_runs))
        print "==============================================================="
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

        :Return: bool
            True if the test run successfully otherwise False.
        """
        valid_output_path = os.path.abspath(run_validator.get_valid_output_path())
        to_check_output_path = os.path.abspath(run_validator.get_to_check_output_path())

        # Create separate output folder if wanted
        if "" != SPECIAL_OUTPUT:
            valid_output_path = os.path.join(SPECIAL_OUTPUT,
                                             run_validator.get_run_dir_name(),
                                             "valid")
            to_check_output_path = os.path.join(SPECIAL_OUTPUT,
                                                 run_validator.get_run_dir_name(),
                                                "to_check" )
            if not os.path.exists(valid_output_path):
                os.makedirs(valid_output_path)
            if not os.path.exists(to_check_output_path):
                os.makedirs(to_check_output_path)

        # Run script with corresponding arguments
        dict_args = {"-f":run_validator.get_names_of_src_files(),
                     "-m":run_validator.get_model_name(),
                     "-a":run_validator.get_accelerator_type(),
                     "-o":valid_output_path}

        if CREATE_VALID_OUTPUT or run_validator.has__no_valid_output_files():
            # Run original/valid script
            print "Deleting old output in "+valid_output_path
            if not Utilities.iotools.deleteFilesWithoutGitignore(valid_output_path):
                print >> sys.stderr,"Could not delete old output files. Will run anyway..."
            valid_script_runner = vimaier_utils.scriptrunner.ScriptRunner(
                                                GETLLM_SCRIPT_VALID, dict_args)
            print "Starting GetLLM_valid.py("+run_validator.get_run_path()+")..."
            errorcode = valid_script_runner.run_script()
            print "GetLLM_valid.py("+run_validator.get_run_path()+") finished...\n"
            if errorcode != 0:
                sys.exit(errorcode)

        dict_args["-o"] = to_check_output_path

        script_runner = vimaier_utils.scriptrunner.ScriptRunner(GETLLM_SCRIPT, dict_args)
        print "Deleting old output in "+to_check_output_path
        if not Utilities.iotools.deleteFilesWithoutGitignore(to_check_output_path):
            print >> sys.stderr,"Could not delete old output files. Will run anyway..."
        print "Starting GetLLM.py("+run_validator.get_run_path()+")..."
        errorcode = script_runner.run_script()
        print "GetLLM.py finished("+run_validator.get_run_path()+")..."
        if errorcode != 0:
            sys.exit(errorcode)

        # Check output of the directory. Therefore:
        # Read filenames of to_check_output_path
        to_check_filenames = []
        for dirname_dirnames_outfilenames in os.walk(to_check_output_path):
            to_check_filenames = dirname_dirnames_outfilenames[2]
            break #there are no subdirectories

        # Downward compatibility for Python interpreter < 2.7
        #self.assertGreater(len(to_check_filenames), 0,
        #                   "No to_check output files in "+to_check_output_path)
        self.assertNotEqual(len(to_check_filenames), 0,
                           "No to_check output files in "+to_check_output_path)

        # Read filenames of valid_output_path
        valid_filenames = []
        for dirname_dirnames_validfilenames in os.walk(valid_output_path):
            valid_filenames = dirname_dirnames_validfilenames[2]
            break #there are no subdirectories

        # Downward compatibility for Python interpreter < 2.7
        #self.assertGreater(len(valid_filenames), 0,
        #                   "No valid output files in "+valid_output_path)
        self.assertNotEqual(len(valid_filenames), 0,
                           "No valid output files in "+valid_output_path)

        # Compare each valid file with corresponding to_check file
        print "Checking output files for run: "+ run_validator.get_run_path()

        num_equal_files = 0
        num_overall_files = len(valid_filenames)
        for name in valid_filenames:
#             err_msg = vimaier_utils.compare_utils.compare_tfs_files(os.path.join(valid_output_path, name),
#                                                                    os.path.join(to_check_output_path, name) )
#             if "" != err_msg:
#                 print name, "are not equal:", err_msg
            succesful = vimaier_utils.compare_utils.compare_tfs_files(os.path.join(valid_output_path, name),
                                                                          os.path.join(to_check_output_path, name) )
            if succesful:
                num_equal_files += 1


        print str(num_equal_files)+" of "+str(num_overall_files)+ " are equal"

        if num_equal_files == num_overall_files:
            # Delete created special folders if files are equal
            if "" != SPECIAL_OUTPUT:
                path_to_created_folder = os.path.join(SPECIAL_OUTPUT, run_validator.get_run_dir_name())
                print "Deleting output folder: ", path_to_created_folder,"\n"
                shutil.rmtree(path_to_created_folder)
            return True
        else:
            return False
    # END run_single_test() --------------------------------------------

# END class TestFileOutPutGetLLM -------------------------------------------------------------------

def main():
    # Remove arguments from sys.argv
    # Needed to avoid conflicts with the options from unittest module
    arguments_tpl = ('-o', "--create_valid_output",
                   "-v","--valid_getllm_script",
                   "-m","--modified_getllm_script",
                   "-p","--path_to_test_data"
                   "-t","--test_name")
    del_lst = []
    for i, option in enumerate(sys.argv):
        if option in arguments_tpl:
            del_lst.append(i)
            del_lst.append(i+1)

    del_lst.reverse()
    for i in del_lst:
        del sys.argv[i]

    # Run the test
    text_test_runner = unittest.TextTestRunner().run(unittest.TestLoader().loadTestsFromTestCase(TestFileOutputGetLLM))

    if 0 != len(text_test_runner.errors):
        sys.exit(len(text_test_runner.errors))
    elif 0 != len(text_test_runner.failures):
        sys.exit(len(text_test_runner.failures))

if __name__ == "__main__":
    main()

