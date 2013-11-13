'''
Created on 16 Oct 2013

@author: vimaier
'''
import os
import sys
import unittest
import subprocess

import __init__ # @UnusedImport init will include paths
import Utilities.iotools
import Utilities.ndiff

class TestSvdClean(unittest.TestCase):


    def setUp(self):
        self.__root = os.path.dirname(os.path.abspath(__file__))
        self.__path_svd_clean = os.path.join(Utilities.iotools.get_absolute_path_to_betabeat_root(),
                                             "Analyses", "svd_clean.py")
        self.__path_to_input = os.path.join(self.__root, "data", "input")
        self.__path_to_check = os.path.join(self.__root, "data", "output", "to_check")
        self.__path_expected = os.path.join(self.__root, "data", "output", "expected")
        Utilities.iotools.create_dirs(self.__path_to_check)


    def tearDown(self):
        self.__delete_created_files()


    def testSvdClean(self):
        self.__run_svd_clean_for_every_file_in_input_dir()
        self.__compare_expected_and_to_check_dirs()

    # helper ---------------------------------------------------------------------------------------
    def __delete_created_files(self):
        Utilities.iotools.delete_content_of_dir(self.__path_to_check)
        self.__delete_bad_bpms_files_in_input_dir()

    def __delete_bad_bpms_files_in_input_dir(self):
        file_paths = Utilities.iotools.get_all_absolute_filenames_in_dir_and_subdirs(self.__path_to_input)
        for file_path in file_paths:
            if file_path.endswith(".bad"):
                Utilities.iotools.delete_item(file_path)


    def __run_svd_clean_for_every_file_in_input_dir(self):
        print "Run svd_clean..."
        input_file_names = Utilities.iotools.get_all_filenames_in_dir_and_subdirs(self.__path_to_input)
        for file_name in input_file_names:
            self.__run_svd_clean_for_file(file_name)

    def __run_svd_clean_for_file(self, file_name):
        debug_option_for_python = "-d" # Needed to run script in DEBUG mode
        call_command = [sys.executable, debug_option_for_python]
        call_command.append(self.__path_svd_clean)
        call_command.append("--file")
        call_command.append(os.path.join(self.__path_to_input, file_name))
        call_command.append("--newfile")
        call_command.append(os.path.join(self.__path_to_check, file_name))
        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

        # wait for the process to terminate
        (out_stream, err_stream) = process.communicate()

        errcode = process.returncode

        if 0 != errcode:
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream

    def __compare_expected_and_to_check_dirs(self):
        print "Compare dirs..."
        ndiff_cfg_file = os.path.join(self.__root, "data", "ignore_line_modified_in_header.cfg")
        self.assertTrue(
                Utilities.ndiff.compare_dirs_with_files_matching_regex_list(
                            self.__path_expected, self.__path_to_check, master_config_file=ndiff_cfg_file
                            ),
                "Created files differ from expected files." )


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()