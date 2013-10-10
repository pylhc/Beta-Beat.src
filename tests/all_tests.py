'''
Created on 6 Aug 2013

@author: vimaier

This module shall contain all PyUnit tests from the Beta-Beat.src directory.
'''
import sys
import unittest

import __init__ # @UnusedImport init will include paths
#import drive.test.test_output
import MODEL.LHCB.model.Corrections.test.filecheck
#import GetLLM.test.filecheck
import GetLLM.test.vimaier_utils.test_compare_utils
import Utilities.test.tfs_file_writer_test

def suite():
    """
    Creates a TestSuit, adds all available testcases and returns the suite.
    """
    test_loader = unittest.TestLoader()
    test_suite = unittest.TestSuite()
    tests = [
             #drive.test.test_output.TestOutput, # Commented out because this test has his own test on CDash
             MODEL.LHCB.model.Corrections.test.filecheck.TestFileOutputGetdiff,
             #GetLLM.test.filecheck.TestFileOutputGetLLM, # Commented out because this test has his own test on CDash
             GetLLM.test.vimaier_utils.test_compare_utils.TestCompareUtils,
             Utilities.test.tfs_file_writer_test.TestTfsFileWriter
             ]

    for t in tests:
        test_suite.addTest(
                        test_loader.loadTestsFromTestCase(t)
                           )
    return test_suite


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    sys.exit(not result.wasSuccessful())
