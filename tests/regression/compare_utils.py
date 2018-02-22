from __future__ import print_function
import os
import re


# Directories comparators #####################################################

def compare_dirs(dir1, dir2, ignore=None, function=None):
    """Test if directories are equal.

        Arguments:
            dir1: First directory to check.
            dir2: Second directory to check.
            ignore: If not None, it must be a list of regular expressions.
                Whatever file or directory in dir1 or dir2 that matches one of
                the patterns will be ignored.
            function: Function to use to compare individual files. If None it
                defaults to compare_text_files_exact. It must take two paths as
                input and return True only if the files are equal.
    """
    if function is None:
        function = compare_text_files_exact
    try:
        childs1 = os.listdir(dir1)
    except os.error:
        print("Cannot read dir {}".format(dir1))
        return False
    try:
        childs2 = os.listdir(dir2)
    except os.error:
        print("Cannot read dir {}".format(dir1))
        return False
    if ignore:
        childs1 = _regexs_filter(childs1, ignore)
        childs2 = _regexs_filter(childs2, ignore)
    in_1_not_in_2 = [name for name in childs1 if name not in childs2]
    in_2_not_in_1 = [name for name in childs2 if name not in childs1]
    if in_1_not_in_2 or in_2_not_in_1:
        print(
            ("Non-common elements found: In {dir1}: {in_1_not_in_2}, "
             "in {dir2}: {in_2_not_in_1}").format(
                 dir1=dir1, dir2=dir2,
                 in_1_not_in_2=in_1_not_in_2, in_2_not_in_1=in_2_not_in_1)
        )
        return False
    result = True
    common_names = childs1
    for name in common_names:
        full_path1 = os.path.join(dir1, name)
        full_path2 = os.path.join(dir2, name)
        if os.path.isfile(full_path1) and os.path.isfile(full_path2):
            result &= function(full_path1, full_path2)
        elif os.path.isfile(full_path1) and os.path.isfile(full_path2):
            result &= compare_dirs(full_path1, full_path2)
        else:
            print("{fp1} and {fp2} are of different type!"
                  .format(fp1=full_path1, fp2=full_path2))
            result = False
    return result


# Single file comparators #####################################################

def compare_text_files_exact(file1, file2):
    """Reads the two files and returns True if they are exactly equal.
    """
    with open(file1, "r") as f1_data, open(file2, "r") as f2_data:
        str1, str2 = f1_data.read(), f2_data.read()
        if str1 != str2:
            print(
                "Files {f1} and {f2} differ."
                .format(f1=file1, f2=file2)
            )
            return False
        return True


# Helpers #####################################################################

def _regexs_filter(names, regexs):
    new_names = []
    matchers = [re.compile(regex) for regex in regexs]
    for name in names:
        if not any([matcher.findall(name) for matcher in matchers]):
            new_names.append(name)
    return new_names
