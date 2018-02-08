'''
Created on  5 Nov 2013

:author: vimaier

This module tests the output files of correct_ChromCoup.py . First it runs the 'valid' version in the dir test/chrom_coup
and then it runs the modified version in the dir Correction. After that it compares the content of the dirs
test.chrom_coup.data.output.valid and to_check.
'''
import unittest
import os
import sys
import subprocess
import shutil


import __init__ # @UnusedImport init will include paths
import Utilities.iotools
import Utilities.ndiff

CURRENT_PATH = os.path.dirname(__file__)
# Path 'x/Beta-Beat.src/Correction/test/chrom_coup'

_SHORT_RUN = False # If True, correct_coupleDy will only run on first dir
_ARGUMENTS_FILE_NAME = "arguments.txt" # Optional file which is located inside the run dir to manipulate the arguments for the script
_DELETE_VALID_OUTPUT = True

class TestCorrectChromCoup(unittest.TestCase):

    PATH_TO_MODIFIED_SCRIPT = os.path.join(CURRENT_PATH, "..", "..", "correct_ChromCoup.py")
    PATH_TO_VALID_SCRIPT = os.path.join(CURRENT_PATH, "valid_correct_ChromCoup.py")

    PATH_TO_VALID = os.path.join(CURRENT_PATH, "data", "output", "valid")
    PATH_TO_TO_CHECK = os.path.join(CURRENT_PATH, "data", "output", "to_check")
    PATH_TO_INPUT = os.path.join(CURRENT_PATH, "data", "input")

    BETABEAT_ROOT = Utilities.iotools.get_absolute_path_to_betabeat_root()

    __successful = False
    __have_to_run_valid_file = False

    def setUp(self):
        self._run_dirs = Utilities.iotools.get_all_dir_names_in_dir(TestCorrectChromCoup.PATH_TO_INPUT)
        self._check_if_valid_output_exist_and_set_attribute()
        self._delete_modified_output()

    def tearDown(self):
        if TestCorrectChromCoup.__successful:
            self._delete_output()

    def testOutput(self):
        print "Start TestCorrectChromCoup"
        for run_dir in self._run_dirs:
            self._run_valid_file_if_necesaary(run_dir)
            self._run_modified_file(run_dir)
            self._compare_output_dirs(run_dir)
            if _SHORT_RUN:
                break
        TestCorrectChromCoup.__successful = True
        print "End TestCorrectChromCoup"

    #===============================================================================================
    # helper
    #===============================================================================================
    def _check_if_valid_output_exist_and_set_attribute(self):
        if self._no_valid_output_exists():
            print "No valid output. Have to run valid file."
            TestCorrectChromCoup.__have_to_run_valid_file = True


    def _no_valid_output_exists(self):
        """
        Returns true if cannot find valid output.
        """
        return not self._valid_output_exists()
    def _valid_output_exists(self):
        all_dirs = Utilities.iotools.get_all_dir_names_in_dir(TestCorrectChromCoup.PATH_TO_VALID)
        if 0 == len(all_dirs) or len(self._run_dirs) != len(all_dirs):
            return False
        is_valid = False
        for dir_in_valid_output in all_dirs:
            is_valid = self._dir_has_valid_files(
                            os.path.join(TestCorrectChromCoup.PATH_TO_VALID, dir_in_valid_output)
                            )
            if not is_valid:
                return False
        return True
    def _dir_has_valid_files(self, path_to_dir):
        """ Returns True if dir contains x, x.knob, x.madx and x.tfs.(x:=changeparameters_chromcouple). """
        has_all_files = True
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "changeparameters_chromcouple"))
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "changeparameters_chromcouple.knob"))
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "changeparameters_chromcouple.madx"))
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "changeparameters_chromcouple.tfs"))
        return has_all_files

    def _delete_output(self):
        if _DELETE_VALID_OUTPUT:
            Utilities.iotools.delete_content_of_dir(TestCorrectChromCoup.PATH_TO_VALID)
        self._delete_modified_output()

    def _delete_modified_output(self):
        ''' Deletes content in PATH_TO_TO_CHECK. '''
        Utilities.iotools.delete_content_of_dir(TestCorrectChromCoup.PATH_TO_TO_CHECK)

    def _run_valid_file_if_necesaary(self, run_dir):
        ''' Runs valid script for given dir '''
        if self.__have_to_run_valid_file:
            print "  Run valid script for "+run_dir
            self._run_in_src_dir(TestCorrectChromCoup.PATH_TO_VALID_SCRIPT, run_dir, os.path.join(TestCorrectChromCoup.PATH_TO_VALID, run_dir))

    def _run_in_src_dir(self, path_to_script, run_dir_name, output_dst_path):
        """
        correct_ChromCoup.py produces files which have the output dir as path in it. Thus two produced files in different
        output dirs produce different results due to the saved path.
        Solution: Both scripts will be run in the input folder and the output files will be moved to the output_dst_path
        for comparing.
        """
        output_path = os.path.join(TestCorrectChromCoup.PATH_TO_INPUT, run_dir_name)
        self._run_script(path_to_script, output_path)
        self._move_produced_files(output_path, output_dst_path)


    def _run_modified_file(self, run_dir):
        print "  Run modified script for "+run_dir
        self._run_in_src_dir(TestCorrectChromCoup.PATH_TO_MODIFIED_SCRIPT, run_dir, os.path.join(TestCorrectChromCoup.PATH_TO_TO_CHECK, run_dir))


    def _run_script(self, path_to_script, path_to_out_run_dir_root):
        ''' Runs given script with path_to_dir.
        path_to_dir has to contain 'FullResponse_chromcouple' and 'chromcoupling.out'(produced by getsuper.py)
        and if needed the file ___ARGUMENTS_FILE_NAME.
        '''
        arguments_list = self._prepare_and_get_arguments_list(path_to_out_run_dir_root)
        call_command = sys.executable + " " + os.path.abspath(path_to_script) + " " + " ".join(arguments_list)

        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)

        # wait for the process to terminate
        (out_stream, err_stream) = process.communicate()

        errcode = process.returncode

        if 0 != errcode:
            print "Error running command:", call_command
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream
        self.assertEqual(0, errcode, "Execution failed")

    def _move_produced_files(self, src_path, dst_path):
        Utilities.iotools.create_dirs(dst_path)
        Utilities.iotools.delete_content_of_dir(dst_path) # Otherwise shutil.move() would raise an error
        for abs_filename in Utilities.iotools.get_all_absolute_filenames_in_dir_and_subdirs(src_path):
            if self._is_produced_output_file(abs_filename):
                shutil.move(abs_filename, dst_path)
    def _is_produced_output_file(self, file_path):
        return "changeparameters_chromcouple" in file_path

    def _prepare_and_get_arguments_list(self, path_to_out_run_dir_root):
        # args_dict contains all arguments of the script. Changeable args were set to default values.
        args_dict = {}
        args_dict["--accel"] = "LHCB1"
        args_dict["--path"] = path_to_out_run_dir_root
        args_dict["--cut"] = "0.1"
        args_dict["--errorcut"] = "0.1,0.1"
        args_dict["--modelcut"] = "0.1,0.1"
        args_dict["--rpath"] = TestCorrectChromCoup.BETABEAT_ROOT
        args_dict["--MinStr"] = "0.000001"
        args_dict["--Dy"] = "1,1,0,0,0"
        args_dict["--opt"] = path_to_out_run_dir_root
        args_dict["--Variables"] = "kss"

        self._set_args_from_file_in_run_dir(path_to_out_run_dir_root, args_dict)

        args_list = []
        for key in args_dict:
            args_list.append(key+"="+args_dict[key])

        return args_list


    def _set_args_from_file_in_run_dir(self, path_to_run, args_dict):
        """ Reads the arguments from file _ARGUMENTS_FILE_NAME
        File should have following syntax:
        --accel=LHCB1
        --cut=0.01
        --errorcut=0.02,0.02
        --modelcut=0.0,0.01
        --MinStr=0.000001
        --Dy=1,1,0,0,1
        # path, rpath and opt should not be set.
        """
        arguments_file_path = os.path.join(path_to_run, _ARGUMENTS_FILE_NAME)
        if os.path.exists(arguments_file_path):
            with open(arguments_file_path) as arg_file:
                for arg_line in arg_file:
                    key = arg_line.split("=", 1)[0].strip()
                    value = arg_line.split("=", 1)[1].strip()
                    args_dict[key] = value


    def _compare_output_dirs(self, run_dir):
        ''' Compares output by using filecmp '''
        print "  Comparing output files"
        valid_dir = os.path.join(self.PATH_TO_VALID, run_dir)
        to_check_dir = os.path.join(self.PATH_TO_TO_CHECK, run_dir)
        dict_file_to_config_file = {
                        "changeparameters_chromcouple.knob" : os.path.join(TestCorrectChromCoup.PATH_TO_INPUT, "ignore_date_in_line_3.cfg"),
                        "changeparameters_chromcouple.tfs" : os.path.join(TestCorrectChromCoup.PATH_TO_INPUT, "ignore_date_in_line_3.cfg")
                                    }
        self.assertTrue(
                        Utilities.ndiff.compare_dirs_with_files_matching_regex_list(valid_dir, to_check_dir, file_to_config_file_dict=dict_file_to_config_file),
                        "Output files are not equal"
                        )

# END TestCorrectChromCoup -----------------------------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()
