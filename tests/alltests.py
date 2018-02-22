from __future__ import print_function
import os
import pytest
import regression.test_cases


def launch_tests(do_unit=True, do_regression=True):
    """Launches all the available tests.
    """
    if do_unit:
        print("Launching unit tests...")
        unit_result = pytest.main([os.path.join(
            os.path.dirname(__file__),
            "unit",
        )])
        if unit_result:  # Non 0 return means failed tests.
            raise UnitTestsFailed()

    if do_regression:
        print("\nLaunching regression tests...")
        regression.test_cases.run_tests()

    # TODO: Precision checks (circular tests)


class UnitTestsFailed(Exception):
    """Raised when the unit tests fail.
    """
    pass


if __name__ == "__main__":
    launch_tests()
