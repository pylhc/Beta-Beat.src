from __future__ import print_function
import os
import filecmp


_IGNORE_DIRS = (".git", )


def compare_dirs_with(valid_dir, test_dir, function=None):
    try:
        _rec_compare_dirs_with(function, valid_dir, test_dir)
    except _DifferFound as diff:
        print(diff.message)
        return False
    return True


def compare_tfs_files(valid_file, test_file, margin_dict=None):
    pass


def _report_on_failure(dircpm_inst, valid_dir, test_dir):
    if dircpm_inst.left_only or dircpm_inst.right_only:
        raise _DifferFound(
            "Found non-common files:\n" +
            "-Only in {}: {}\n".format(valid_dir, dircpm_inst.left_only) +
            "-Only in {}: {}\n".format(test_dir, dircpm_inst.right_only)
        )
    if dircpm_inst.funny_files:
        raise _DifferFound(
            "Couldn't read files {} from {} and {}"
            .format(dircpm_inst.funny_files, valid_dir, test_dir)
        )
    if dircpm_inst.diff_files:
        raise _DifferFound(
            "Files differ: {} from {} and {}"
            .format(dircpm_inst.diff_files, valid_dir, test_dir)
        )


def _rec_compare_dirs_with(function, valid_dir, test_dir):
    my_dircmp = _MyDircmp(valid_dir, test_dir, cmp_funct=function)
    _report_on_failure(my_dircmp, valid_dir, test_dir)
    for subdir in my_dircmp.common_dirs:
        new_difcmp = _rec_compare_dirs_with(function,
                                            os.path.join(valid_dir, subdir),
                                            os.path.join(test_dir, subdir))
        _report_on_failure(new_difcmp, valid_dir, test_dir)
    return my_dircmp


class _MyDircmp(filecmp.dircmp):
    """
    """

    def __init__(self, a, b, ignore=None, hide=None,
                 cmp_funct=filecmp.cmpfiles):
        filecmp.dircmp.__init__(self, a, b, ignore, hide)
        self.cmp_funct = cmp_funct

    def phase3(self):
        # Find out differences between common files
        xx = self.cmp_funct(self.left, self.right, self.common_files)
        self.same_files, self.diff_files, self.funny_files = xx

    def phase4(self):
        # Find out differences between common subdirectories
        # A new dircmp object is created for each common subdirectory,
        # these are stored in a dictionary indexed by filename.
        # The hide and ignore properties are inherited from the parent
        self.subdirs = {}
        for x in self.common_dirs:
            a_x = os.path.join(self.left, x)
            b_x = os.path.join(self.right, x)
            self.subdirs[x] = _MyDircmp(a_x, b_x, self.ignore, self.hide)


class _DifferFound(Exception):
    pass
