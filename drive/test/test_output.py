'''
Created on 1 Jul 2013

@author: vimaier

This module tests the output files of drive. First it runs Drive_God_lin in the dir drive.test.valid.
Then it runs the modified version in the dir drive. After that it compares the content of the dirs
drive.test.data.valid and to_check.
'''
import unittest
import os
import sys
import subprocess
import filecmp

CURRENT_PATH = os.path.dirname(__file__)
# Path 'x/Beta-Beat.src/drive/test'

import __init__ # @UnusedImport init will include paths
import Utilities.iotools
import Utilities.ndiff

_SHORT_RUN = False # If True, Drive will only run on first dir

class TestOutput(unittest.TestCase):

    path_to_valid_drive = os.path.join(CURRENT_PATH, "valid", "Drive_God_lin")

    path_to_valid = os.path.join(CURRENT_PATH, "data", "valid")
    path_to_to_check = os.path.join(CURRENT_PATH, "data", "to_check")
    path_to_input = os.path.join(CURRENT_PATH, "data", "input")

    __have_to_run_valid_file = False
    successful = False

    @staticmethod
    def get_os_dependent_path_to_modified_drive():
        if sys.platform == "linux" or sys.platform == "linux2":
            return os.path.join(CURRENT_PATH, "..", "Drive_God_lin")
        elif sys.platform == "win32":
            return os.path.join(CURRENT_PATH, "..", "Drive_God_lin_win.exe")
            # Windows...
        else:
            raise OSError("No drive version for given os: "+sys.platform)

    def setUp(self):
        self._check_if_valid_output_exist_and_set_attribute()
        self._delete_modified_output()
        self._copied_files = self._copy_input_files()


    def tearDown(self):
        if TestOutput.successful:
            self._delete_modified_output()


    def testOutput(self):
        print "Start TestOutput of drive"
        self._run_modified_file()
        self._compare_output_dirs()
        print "End TestOutput of drive"

    #===============================================================================================
    # helper
    #===============================================================================================
    def _check_if_valid_output_exist_and_set_attribute(self):
        if self._no_valid_output_exists():
            print "No valid output. Have to run valid file."
            self.__have_to_run_valid_file = True


    def _no_valid_output_exists(self):
        """
        Returns true if cannot find valid output. Assuming that every subdir of path_to_valid has
        a linx and liny file(valid output).
        """
        return not self._valid_output_exists()

    def _valid_output_exists(self):
        """ True, if valid output available.
            Assuming that every subdir of path_to_valid has a linx and liny file(valid output)."""
        if Utilities.iotools.not_exists_directory(TestOutput.path_to_valid):
            return False
        is_valid = False
        for item in os.listdir(TestOutput.path_to_valid):
            abs_item = os.path.join(TestOutput.path_to_valid, item)
            if os.path.isdir(abs_item):
                is_valid = self._dir_has_linx_and_liny(abs_item)
                if not is_valid:
                    return False
        return True



    def _dir_has_linx_and_liny(self, path_to_dir):
        """ Returns True if dir contains a *linx and a *liny file. """
        has_linx = False
        has_liny = False
        for item in os.listdir(path_to_dir):
            if item.endswith("linx"):
                has_linx = True
            if item.endswith("liny"):
                has_liny = True
        return has_linx and has_liny


    def _delete_modified_output(self):
        ''' Deletes content in path_to_valid and path_to_to_check. '''
        Utilities.iotools.delete_content_of_dir(TestOutput.path_to_to_check)

    def _copy_input_files(self):
        ''' Copies input data for drive from path_to_input. Return all copied files '''
        print "Copying input files"
        # Only Drive.inp and DrivingTerms are needed. DrivingTerms holds abs path to sdds file (vimaier)
        all_run_dirs = Utilities.iotools.get_all_dir_names_in_dir(TestOutput.path_to_input)
        for run_dir in all_run_dirs:
            src = os.path.join(TestOutput.path_to_input, run_dir, "Drive.inp")
            dst = os.path.join(TestOutput.path_to_to_check, run_dir)
            Utilities.iotools.create_dirs(dst)
            Utilities.iotools.copy_item(src, dst)
            src = os.path.join(TestOutput.path_to_input, run_dir, "DrivingTerms")
            dst = os.path.join(TestOutput.path_to_to_check, run_dir)
            Utilities.iotools.copy_item(src, dst)
        return Utilities.iotools.get_all_absolute_filenames_in_dir_and_subdirs(TestOutput.path_to_to_check)


    def _run_dir(self, path_to_run_dir, path_to_drive):
        """ Runs drive from path_to_drive for all subdirs in path_to_run_dir """
        for index, directory in enumerate(os.listdir(path_to_run_dir)):
            if self._break_after_first_run(index):
                break
            single_dir_path = os.path.join(path_to_run_dir, directory)
            print "    Run:", single_dir_path
            self._run_drive(path_to_drive, single_dir_path)

    def _break_after_first_run(self, index):
        return 0 < index and _SHORT_RUN


    def _run_modified_file(self):
        ''' Runs drive.Drive_God_lin for all directories in drive.test.data.to_check. '''
        print "  Run to_check files"
        self._run_dir(self.path_to_to_check, TestOutput.get_os_dependent_path_to_modified_drive())


    def _run_drive(self, path_to_drive, path_to_dir):
        ''' Runs given drive with path_to_dir.
        path_to_dir has to contain Drive.inp and DrivingTerms.
        '''
        call_command = os.path.abspath(path_to_drive) + " " + os.path.abspath(path_to_dir)

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


    def _delete_copied_files(self):
        for copied_file in self._copied_files:
            Utilities.iotools.delete_item(copied_file)


    def _compare_output_dirs(self):
        ''' Compares output by using filecmp '''
        print "  Comparing output files"
        self._delete_copied_files()
        for index, directory in enumerate(os.listdir(self.path_to_valid)):
            if self._break_after_first_run(index):
                break
            valid_dir = os.path.join(self.path_to_valid, directory)
            to_check_dir = os.path.join(self.path_to_to_check, directory)
#            self._compare_dir_with_filecmp(valid_dir, to_check_dir)
#             self._compare_dir_almost_equal_doubles(valid_dir, to_check_dir)
            self._compare_dirs_with_ndiff(valid_dir, to_check_dir)
        TestOutput.successful = True


    def _compare_dir_with_filecmp(self, valid_dir, to_check_dir):
        dir_compare = filecmp.dircmp(valid_dir, to_check_dir)
        dir_compare.report()
        self.assertEqual(0, len(dir_compare.diff_files), "Files are not equal.")
        self.assertEqual(0, len(dir_compare.left_only), "Files existing in only one dir")
        self.assertEqual(0, len(dir_compare.right_only), "Files existing in only one dir")


    def _compare_dir_almost_equal_doubles(self, valid_dir, to_check_dir):
        self.assertTrue(
                        Utilities.compare.equal_dirs_with_double_epsilon_comparing(valid_dir, to_check_dir, except_files=self._copied_files),
                        "Directories not equal: "+valid_dir+" and "+to_check_dir
                        )


    def _compare_dirs_with_ndiff(self, valid_dir, to_check_dir):
        regex_list = [r"^.*\.lin(x|y)$", r"^.*\.(x|y)$"] # *.linx; *.liny; *.x; *.y  ==> Output of drive
        self.assertTrue(
                        Utilities.ndiff.compare_dirs_with_files_matching_regex_list(valid_dir, to_check_dir, regex_list),
                        "Directories not equal: "+valid_dir+" and "+to_check_dir
                        )

# END TestOutput -----------------------------------------------------------------------------------


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestOutput.testOutput']
    unittest.main()
