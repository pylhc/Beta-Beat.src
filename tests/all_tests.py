'''
Created on 6 Aug 2013

@author: vimaier

This module shall contain all PyUnit tests from the Beta-Beat.src directory.
'''
import sys
import unittest

import __init__ # @UnusedImport init will include paths
import drive.test.test_output
import MODEL.LHCB.model.Corrections.test.filecheck
import GetLLM.test.filecheck


def suite():
    """
    Creates a TestSuit, adds all available testcases and returns the suite.
    """
    test_suite = unittest.TestSuite()
    test_suite.addTest(
                       get_all_test_as_suite(drive.test.test_output.TestOutput)
                       )
    test_suite.addTest(
                       get_all_test_as_suite(MODEL.LHCB.model.Corrections.test.filecheck.TestFileOutputGetdiff)
                       )
    test_suite.addTest(
                       get_all_test_as_suite(GetLLM.test.filecheck.TestFileOutputGetLLM)
                       )
    return test_suite

    
def get_all_test_as_suite(TestClass):
    suite_of_testclass = unittest.TestSuite()
    for method in dir(TestClass):
        if method.startswith("test"):
            suite_of_testclass.addTest(TestClass(method))
    return suite_of_testclass


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=2).run(suite())
    sys.exit(not result.wasSuccessful())
    