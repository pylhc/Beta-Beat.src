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
import Utilities.test.tfs_file_writer_test
import Utilities.test.iotools_test
import Analyses.test.svd_clean_test
#import MODEL.LHCB.fullresponse.test.test_fullresponse_parallel
import Correction.test.test_correct_coupleDy

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
             Utilities.test.tfs_file_writer_test.TestTfsFileWriter,
             Utilities.test.iotools_test.TestReplacingKeywords,
             Utilities.test.iotools_test.TestGetFilenamesInDir,
             #MODEL.LHCB.fullresponse.test.test_fullresponse_parallel.TestGenFullRespParallel, # Commented out because this test has his own test on CDash
             Analyses.test.svd_clean_test.TestSvdClean,
             Correction.test.test_correct_coupleDy.TestCorrectCoupleDy
             ]

    for t in tests:
        test_suite.addTest(
                        test_loader.loadTestsFromTestCase(t)
                           )
    return test_suite


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    sys.exit(not result.wasSuccessful())
