import os
from regression import TestCase, launch_test_set


TEST_CASES = (
    TestCase("test_the_test",
             "tests/regression/something.py", "5", "./_test_out/",
             lambda dir1, dir2: False),
)


if __name__ == "__main__":
    launch_test_set(TEST_CASES, os.path.abspath("../.."))
