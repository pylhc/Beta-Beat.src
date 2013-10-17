'''
Created on 4 Jun 2013

@author: vimaier

@version: 0.0.1

MODEL.LHCB.model.Corrections.test.filecheck
Compares output files of
    MODEL.LHCB.model.Corrections.getdiff
and
    MODEL.LHCB.model.Corrections.test.getdiff_valid
to validate the correctness of getdiff.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import os
import sys
import unittest
import filecmp
import subprocess

import __init__  # @UnusedImport
import Utilities.iotools



GETDIFF_PATH = os.path.dirname( os.path.dirname(os.path.abspath(__file__)) )

#===================================================================================================
# # TestFileOutputGetdiff
#===================================================================================================
class TestFileOutputGetdiff(unittest.TestCase):
    """TestCase which compares the output of getdiff.py"""


    def setUp(self):
        self.path_to_valid = os.path.join(GETDIFF_PATH, "test", "data", "output", "valid")
        self.path_to_to_check = os.path.join(GETDIFF_PATH, "test", "data", "output", "to_check")
        self.path_to_input_files = os.path.join(GETDIFF_PATH, "test", "data", "input")
        """ List which holds all copied filenames from input."""

        self._create_output_dirs_if_not_exist()
        self._delete_content_of_output_dirs()
        self._copy_input_files()


    def tearDown(self):
        self._delete_content_of_output_dirs()


    def test_file_output(self):
        """ Tests the output files of the getdiff.py script."""
        command = "python " + os.path.join(GETDIFF_PATH,"getdiff.py ") + self.path_to_valid
        ret = subprocess.call(command, shell=True)
        self.assertEqual(0, ret, "Execution of getdiff.py failed: "+command)
        command = "python "+os.path.join(GETDIFF_PATH,"test","getdiff_valid.py ") + self.path_to_to_check
        ret = subprocess.call(command, shell=True)
        self.assertEqual(0, ret, "Execution of getdiff.py failed: "+command)

        dir_compare = filecmp.dircmp(self.path_to_valid, self.path_to_to_check)
        dir_compare.report()
        self.assertEqual(0, len(dir_compare.diff_files), "Files are not equal.")

    #===============================================================================================
    # helper
    #===============================================================================================
    def _create_output_dirs_if_not_exist(self):
        Utilities.iotools.create_dirs(self.path_to_valid)
        Utilities.iotools.create_dirs(self.path_to_to_check)

    def _delete_content_of_output_dirs(self):
        ''' Deletes content in path_to_valid and path_to_to_check. '''
        Utilities.iotools.delete_content_of_dir(self.path_to_valid)
        Utilities.iotools.delete_content_of_dir(self.path_to_to_check)


    def _copy_input_files(self):
        ''' Copies input data for drive from path_to_input_files. '''
        Utilities.iotools.copy_content_of_dir(self.path_to_input_files, self.path_to_valid)
        Utilities.iotools.copy_content_of_dir(self.path_to_input_files, self.path_to_to_check)

# END class TestFileOutPutGetLLM -------------------------------------------------------------------



def main():
    # Run the test
    text_test_runner = unittest.TextTestRunner().run(unittest.TestLoader().loadTestsFromTestCase(TestFileOutputGetdiff))

    if 0 != len(text_test_runner.errors):
        sys.exit(len(text_test_runner.errors))
    elif 0 != len(text_test_runner.failures):
        sys.exit(len(text_test_runner.failures))

if __name__ == "__main__":
    main()
