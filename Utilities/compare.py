'''
Created on 8 Aug 2013

@author: vimaier

Hold tools for comparing.

'''
import sys
import os

import Utilities

def equal_dirs_with_double_epsilon_comparing(dir1, dir2):
    """ 
    Compares file in given dirs by going through each line, spliting line and trying to parse to 
    double and comparing double with epsilon value. 
    In case of a subdir recursion will be used.
    """
    if Utilities.iotools.no_dirs_exist(dir1, dir2):
        print >> sys.stderr, dir1, "or(and)", dir2, "do(es) not exist."
        return False
    
    dir1_items = os.listdir(dir1)
    dir2_items = os.listdir(dir2)
    
    if dir1_items != dir2_items:
        print >> sys.stderr, "Items in dirs are not equal:\n",dir1_items, "\n", dir2_items
        return False

    for item in dir1:
        item1 = os.path.join(dir1, item)
        item2 = os.path.join(dir2, item)
        if os.path.isdir(item1):
            equal_dirs_with_double_epsilon_comparing(item1, item2)
        else:
            if not equal_files_with_double_epsilon_comparing(open(item1), open(item2)):
                return False
    return True
    
    
def equal_files_with_double_epsilon_comparing(file1, file2):
    for line in file1:
        split_line1 = line.split()
        split_line2 = file2.readline()
        if not _equal_splitted_lines_with_double_epsilon(split_line1, split_line2):
            return False
    return True

def _equal_splitted_lines_with_double_epsilon(str_list1, str_list2):
    if len(str_list1) != len(str_list2):
        print >> sys.stderr, "Line lengths not equal:\n",str_list1, "\n", str_list2
        return False
    for i in xrange(len(str_list1)):
        token1 = str_list1[i]
        token2 = str_list2[i]
        
        if not equal_strings_with_double_epsilon(token1, token2):
            return False
    return True


def equal_strings_with_double_epsilon(str1, str2):    
    equal = True
    float_or_false = is_float(str1)
    
    if False != float_or_false:
        str1 = float_or_false
        float_or_false = is_float(str2)
        if False == float_or_false:
            equal = False
        else:
            equal = almost_equal_double(str1, float_or_false)
    else:
        equal = str1 == str2
    
    if not equal:
        print >> sys.stderr, "Tokens not Equal:", str1, ";", str2
    
    return equal
    
        
def is_float(s):
    try:
        result = float(s)
        return result
    except ValueError:
        return False
    
def almost_equal_double(d1, d2, eps=0.000001):
    return abs(d1-d2) < abs(d1)*eps

    