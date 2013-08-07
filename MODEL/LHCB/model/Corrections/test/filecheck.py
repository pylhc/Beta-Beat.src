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
import shutil
import filecmp
import subprocess

import __init__
import Utilities.iotools



GETDIFF_PATH = os.path.dirname( os.path.dirname(os.path.abspath(__file__)) )

#===================================================================================================
# # TestFileOutputGetdiff
#===================================================================================================
class TestFileOutputGetdiff(unittest.TestCase):
    """TestCase which compares the output of getdiff.py"""

    
    num_of_failed_tests = 0
    """Static variable which is used as exit code."""
    
    list_of_in_files = []
    """ List which holds all copied filenames from input."""
    path_to_valid = os.path.join(GETDIFF_PATH,"test","data","output","valid")
    path_to_to_check = os.path.join(GETDIFF_PATH,"test","data","output","to_check")
    path_to_input_files = os.path.join(GETDIFF_PATH,"test","data","input")
    
    def setUp(self):
        self._create_output_dirs_if_not_exist()
        self._delete_output_dirs()
        self._copy_input_files()


    def tearDown(self):
        self._delete_copied_input_files()

    
    def test_file_output(self):
        """ Tests the output files of the getdiff.py script."""
        command ="python "+os.path.join(GETDIFF_PATH,"getdiff.py ") + self.path_to_valid
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
    
    def _delete_output_dirs(self):
        ''' Deletes content in path_to_valid and path_to_to_check. '''
        Utilities.iotools.delete_content_of_dir(self.path_to_valid)
        Utilities.iotools.delete_content_of_dir(self.path_to_to_check)


    def _copy_input_files(self):
        ''' Copies input data for drive from path_to_input_files. '''
        Utilities.iotools.copy_content_of_dir(self.path_to_input_files, self.path_to_valid)
        Utilities.iotools.copy_content_of_dir(self.path_to_input_files, self.path_to_to_check)
        
    def _delete_copied_input_files(self):
        for file_name in self.list_of_in_files:
            os.remove(os.path.join(self.path_to_valid, file_name))
            os.remove(os.path.join(self.path_to_to_check, file_name))

        
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
