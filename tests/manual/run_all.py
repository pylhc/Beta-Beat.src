from __future__ import print_function
import sys
import pytest
from os.path import abspath, join, dirname, pardir
sys.path.append(abspath(join(dirname(__file__), pardir, pardir)))


def launch_tests():
    """Launches all the available tests.
    """
    print("Launching unit tests...")
    unit_result = pytest.main([join(dirname(__file__), "unit",)])
    if unit_result:  # Non 0 return means failed tests.
        raise UnitTestsFailed()


class UnitTestsFailed(Exception):
    pass


if __name__ == "__main__":
    launch_tests()
