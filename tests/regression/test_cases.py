from __future__ import print_function
import os
from regression import TestCase, launch_test_set, find_regressions
import compare_utils


def get_test_cases():
    return (
        TestCase(name="should_error",
                 script="tests/regression/something.py",
                 arguments="5",
                 output="tests/regression/_test_out/",
                 test_function=compare_utils.compare_dirs_with),
        TestCase(name="should_fail",
                 script="tests/regression/something.py",
                 arguments="5",
                 output="tests/regression/_test_out/",
                 test_function=lambda x, y: False),
        TestCase(name="should_succeed",
                 script="tests/regression/something.py",
                 arguments="5",
                 output="tests/regression/_test_out/",
                 test_function=lambda x, y: True),
    )


def test_function(valid_path, test_path):
    print("Compare valid: {} and test: {}.".format(valid_path, test_path))
    raise NotImplementedError()


if __name__ == "__main__":
    launch_test_set(get_test_cases(), os.path.abspath("../.."))
    # find_regressions(get_test_cases(), "./_test_repo", os.path.abspath("../.."))
