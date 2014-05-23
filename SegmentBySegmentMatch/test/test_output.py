import __init__  # @UnusedImport
import unittest
import os
from Utilities import iotools, ndiff
import sys
import subprocess


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class TestOutput(unittest.TestCase):

    path_to_valid = os.path.join(CURRENT_PATH, "data", "valid")
    path_to_check = os.path.join(CURRENT_PATH, "data", "to_check")
    path_to_input = os.path.join(CURRENT_PATH, "data", "input")
    path_to_modified_sbs_match = os.path.join(CURRENT_PATH, "..", "SegmentBySegmentMatch.py")

    def tearDown(self):
        input_files = iotools.get_all_filenames_in_dir(self.path_to_input)
        for input_file_name in input_files:
            directory_to_check = os.path.join(self.path_to_check, input_file_name)
            input_arguments = " clean --temp=" + directory_to_check
            call_command = os.path.abspath(self.path_to_modified_sbs_match) + " " + input_arguments
            call_command = sys.executable + " " + call_command
            process = subprocess.Popen(call_command,
                           shell=True)
            process.communicate()
            iotools.delete_item(directory_to_check)

    def testOutput(self):
        input_files = iotools.get_all_filenames_in_dir(self.path_to_input)
        for input_file_name in input_files:
            directory_valid = os.path.join(self.path_to_valid, input_file_name)
            self.assertTrue(iotools.exists_directory(directory_valid),
                            "Valid directory for input file " + input_file_name + " doesn't exists")
            print "Starting run for " + input_file_name
            self._run_sbsmatch_for_input(input_file_name)
            print "Checking results..."
            self._check_results_for_input(input_file_name)
            print "Done"

    def _run_sbsmatch_for_input(self, input_file_name):
        directory_to_check = os.path.join(self.path_to_check, input_file_name)
        file_input = os.path.join(self.path_to_input, input_file_name)

        input_arguments_file = open(file_input, "r")
        input_arguments = input_arguments_file.readline()
        self.assertFalse(input_arguments is None, "Can't read input file " + input_file_name)
        input_arguments += " --temp=" + directory_to_check

        call_command = os.path.abspath(self.path_to_modified_sbs_match) + " " + input_arguments
        call_command = sys.executable + " " + call_command

        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)
        (out_stream, err_stream) = process.communicate()

        errcode = process.returncode

        if 0 != errcode:
            print "Error running command:", call_command
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream

    def _check_results_for_input(self, input_file_name):
        directory_valid = os.path.join(self.path_to_valid, input_file_name)
        directory_to_check = os.path.join(self.path_to_check, input_file_name, "match")
        self._compare_directories_with_ndiff(directory_valid, directory_to_check)

    def _compare_directories_with_ndiff(self, directory_valid, directory_to_check):
        valid_item_list = os.listdir(directory_valid)
        for valid_item in valid_item_list:
            valid_item_path = os.path.join(directory_valid, valid_item)
            to_check_item_path = os.path.join(directory_to_check, valid_item)
            self.assertTrue(os.path.exists(to_check_item_path), "Unable to find item " + to_check_item_path)
            if(os.path.isdir(valid_item_path)):
                self._compare_directories_with_ndiff(valid_item_path, to_check_item_path)
            else:
                self._compare_files_with_ndiff(valid_item_path, to_check_item_path)

    def _compare_files_with_ndiff(self, valid_file_path, to_check_file_path):
        self.assertTrue(ndiff.compare_tfs_files_and_ignore_header(valid_file_path, to_check_file_path),
                        "Differences found in file: " + to_check_file_path)


if __name__ == "__main__":
    unittest.main()
