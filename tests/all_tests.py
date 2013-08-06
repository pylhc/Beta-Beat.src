'''
Created on 6 Aug 2013

@author: vimaier

This module shall contain all PyUnit tests from the Beta-Beat.src directory.
'''
import sys
import unittest

import drive.test.test_output
import MODEL.LHCB.model.Corrections.test.filecheck
import GetLLM.test.filecheck


def suite():
    """
    Creates a TestSuit, adds all available testcases and returns the suite.
    """
    test_suite = unittest.TestSuite()
    test_suite.addTest(drive.test.test_output.TestOutput)
    test_suite.addTest(MODEL.LHCB.model.Corrections.test.filecheck.TestFileOutputGetdiff)
    test_suite.addTest(GetLLM.test.filecheck.TestFileOutputGetLLM)
    return test_suite


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    sys.exit(not result.wasSuccessful())