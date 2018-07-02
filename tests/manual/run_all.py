from __future__ import print_function
import os
import sys
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)))

from tests.alltests import UnitTestsFailed


def launch_tests():
    """Launches all the available tests.
    """
    print("Launching unit tests...")
    unit_result = pytest.main([os.path.join(
        os.path.dirname(__file__),
        "unit",
    )])
    if unit_result:  # Non 0 return means failed tests.
        raise UnitTestsFailed()


if __name__ == "__main__":
    launch_tests()