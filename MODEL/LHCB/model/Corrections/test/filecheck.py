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
# Add directory to the python search path - needed to run script from command line
# Otherwise ScriptRunner won't be found
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join( os.path.dirname(os.path.abspath(__file__)),"../../../Python_Classes4MAD" ))



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
    
    def setUp(self):
        # Copy needed input files
        p_from = os.path.join(GETDIFF_PATH,"test","data","input")
        
        self.list_of_in_files = os.listdir(p_from)
        for file_name in self.list_of_in_files:
            shutil.copy(os.path.join(p_from,file_name), self.path_to_valid)
            shutil.copy(os.path.join(p_from,file_name), self.path_to_to_check)


    def tearDown(self):
        # Delete copied input files again
        for file_name in self.list_of_in_files:
            os.remove(os.path.join(self.path_to_valid, file_name))
            os.remove(os.path.join(self.path_to_to_check, file_name))

    
    def test_file_output(self):
        """ Tests the output files of the getdiff.py script."""
        command ="python "+os.path.join(GETDIFF_PATH,"getdiff.py ") + self.path_to_valid
        print subprocess.call(command)
        command = "python "+os.path.join(GETDIFF_PATH,"test","getdiff_valid.py ") + self.path_to_to_check
        print subprocess.call(command)
        
        dir_compare = filecmp.dircmp(self.path_to_valid, self.path_to_to_check)
        
        dir_compare.report()

        
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
