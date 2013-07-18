'''
Created on 18 Jul 2013

@author: vimaier
'''
import unittest

import compare_utils

class TestCompareUtils(unittest.TestCase):


    def test_is_whitespace(self):
        self.assertTrue(compare_utils.is_whitespace(" "), "is_whitespace cannot detect ' '")
        self.assertTrue(compare_utils.is_whitespace("\t"), "is_whitespace cannot detect '\t'")
        
        
    def test_parse_descriptor_line(self):
        line = '@ ABC %s "A test"'
        split_line = compare_utils.parse_descriptor_line(line)
        self.assertEqual(['@','ABC','%s','A test'], split_line, "Split line is wrong:"+(" ".join(split_line)))

        line = '@ABC %s test'
        split_line = compare_utils.parse_descriptor_line(line)
        self.assertEqual(['@','ABC','%s','test'], split_line, "Split line is wrong:"+(" ".join(split_line)))

        line = '@     ABC         %s           "test            asf            sdfsf    "         '
        split_line = compare_utils.parse_descriptor_line(line)
        self.assertEqual(['@','ABC','%s','test            asf            sdfsf    '], split_line, "Split line is wrong:"+(" ".join(split_line)))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_is_whitespace']
    unittest.main()