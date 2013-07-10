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

# Add directory to the python search path - needed to run script from command line
# Otherwise Utilities.iotools.py won't be found
sys.path.append(os.path.join(CURRENT_PATH, "..", ".."))

import Utilities.iotools


class TestOutput(unittest.TestCase):
    
    path_to_modified_drive = os.path.join(CURRENT_PATH, "..", "Drive_God_lin")
    path_to_valid_drive = os.path.join(CURRENT_PATH, "valid", "Drive_God_lin")
    
    path_to_valid = os.path.join(CURRENT_PATH, "data", "valid")
    path_to_to_check = os.path.join(CURRENT_PATH, "data", "to_check")
    path_to_input = os.path.join(CURRENT_PATH, "data", "input")
    
    successful = False

    def setUp(self):
        self._delete_output_dirs()
        self._copy_input_files()


    def tearDown(self):
        if TestOutput.successful:
            self._delete_output_dirs()


    def testOutput(self):
        print "Start TestOutput of drive"
        self._run_valid_file()
        self._run_modified_file()
        self._compare_output_dirs()
        print "End TestOutput of drive"
    
    #===============================================================================================
    # helper
    #===============================================================================================
    def _delete_output_dirs(self):
        ''' Deletes content in path_to_valid and path_to_to_check. '''
        Utilities.iotools.delete_content_of_dir(TestOutput.path_to_valid)
        Utilities.iotools.delete_content_of_dir(TestOutput.path_to_to_check)


    def _copy_input_files(self):
        ''' Copies input data for drive from path_to_input. '''
        Utilities.iotools.copy_content_of_dir(TestOutput.path_to_input, TestOutput.path_to_valid)
        Utilities.iotools.copy_content_of_dir(TestOutput.path_to_input, TestOutput.path_to_to_check)
        
    
    def _run_valid_file(self):
        ''' Runs drive.test.valid.Drive_God_lin for all directories in drive.test.data.valid '''
        print "  Run valid files"
        for directory in os.listdir(self.path_to_valid):
            valid_dir_path = os.path.join(self.path_to_valid, directory)
            print "    Run:", valid_dir_path
            self._run_drive(self.path_to_valid_drive, valid_dir_path)
            
            
    def _run_modified_file(self):
        ''' Runs drive.Drive_God_lin for all directories in drive.test.data.to_check. '''
        print "  Run to_check files"
        for directory in os.listdir(self.path_to_to_check):
            to_check_dir_path = os.path.join(self.path_to_to_check, directory)
            print "    Run:", to_check_dir_path
            self._run_drive(self.path_to_modified_drive, to_check_dir_path)
            
    
    def _run_drive(self, path_to_drive, path_to_dir):
        ''' Runs given drive with path_to_dir. 
        path_to_dir has to contain Drive.inp and DrivingTerms and the sdds file.
        '''
        call_command = [path_to_drive, path_to_dir]
        
        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE)

        # wait for the process to terminate
        (out_stream, err_stream) = process.communicate()
        
        errcode = process.returncode
        
        if 0 != errcode:
            print "Error running command:", " ".join(call_command)
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream
            

    def _compare_output_dirs(self):
        ''' Compares output by using filecmp '''
        print "  Comparing output files"

        for directory in os.listdir(self.path_to_valid):
            valid_dir = os.path.join(self.path_to_valid, directory)
            to_check_dir = os.path.join(self.path_to_to_check, directory)
            dir_compare = filecmp.dircmp(valid_dir, to_check_dir)
            dir_compare.report()
            self.assertEqual(0, len(dir_compare.diff_files), "Files are not equal.")
        
        TestOutput.successful = True
        
# END TestOutput -----------------------------------------------------------------------------------


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestOutput.testOutput']
    unittest.main()