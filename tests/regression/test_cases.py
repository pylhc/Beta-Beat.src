import sys
import os
from os.path import join, abspath, dirname
import compare_utils
import regression

ABS_ROOT = abspath(join(dirname(__file__), "..", ".."))

REGR_DIR = join("tests", "regression")
TBTS = join("tests", "inputs", "tbt_files")
MODELS = join("tests", "inputs", "models")


TEST_CASES_HOLE_IN_ONE = (
    regression.TestCase(
        name="hole_in_one_test_flat_3dkick",
        script=join("hole_in_one", "hole_in_one.py"),
        arguments=("--file={file} --model={model} --output={output} clean "
                   "harpy --tunex 0.27 --tuney 0.322 --tunez 4.5e-4 "
                   "--nattunex 0.28 --nattuney 0.31".format(
                       file=join(TBTS, "flat_beam1_3d.sdds"),
                       model=join(MODELS, "flat_beam1", "twiss.dat"),
                       output=join(REGR_DIR, "_out_hole_in_one_test_flat_3dkick"))),
        output=join(REGR_DIR, "_out_hole_in_one_test_flat_3dkick"),
        test_function=lambda d1, d2: compare_utils.compare_dirs(d1, d2, ignore=[r".*\.log"]),
        pre_hook=lambda dir: os.makedirs(join(dir, REGR_DIR, "_out_hole_in_one_test_flat_3dkick")),
    ),
)


def run_tests():
    result = regression.launch_test_set(TEST_CASES_HOLE_IN_ONE, ABS_ROOT)
    if not result:
        # Test failed
        sys.exit(-1)


if __name__ == "__main__":
    run_tests()
