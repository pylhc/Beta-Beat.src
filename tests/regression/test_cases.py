import sys
import os
from os.path import join, abspath, dirname
import compare_utils
import regression

ABS_ROOT = abspath(join(dirname(__file__), "..", ".."))

REGR_DIR = join("tests", "regression")
TBTS = join("tests", "inputs", "tbt_files")
MODELS = join("tests", "inputs", "models")
HARM_FILES = join("tests", "inputs", "harmonic_results")


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


TEST_CASES_GETLLM = (
    regression.TestCase(
        name="getllm_test_flat_disp",
        script=join("GetLLM", "GetLLM.py"),
        arguments=("--accel=LHCB1 "
                   "--model={model} "
                   "--files={files_dir}/on_mom_file1.sdds,{files_dir}/on_mom_file2.sdds,{files_dir}/neg_mom_file1.sdds,{files_dir}/pos_mom_file1.sdds "
                   "--output={output} "
                   "--tbtana=SUSSIX --bpmu=mm --lhcphase=1 "
                   "--errordefs={errordefs}")
                   .format(model=join(MODELS, "flat_beam1", "twiss.dat"),
                           files_dir=join(HARM_FILES, "flat_60_15cm_b1"),
                           output=join(REGR_DIR, "_out_getllm_test_flat_disp"),
                           errordefs=join(MODELS, "flat_beam1", "error_deff.txt")),
        output=join(REGR_DIR, "_out_getllm_test_flat_disp"),
        test_function=compare_utils.compare_dirs,
        pre_hook=lambda dir: None,
    ),
)


def run_tests():
    """Run the test cases and raise RegressionTestFailed on failure.
    """
    alltests = list(TEST_CASES_HOLE_IN_ONE) + list(TEST_CASES_GETLLM)
    regression.launch_test_set(alltests, ABS_ROOT, tag_regexp="^BBGUI_.*$")


if __name__ == "__main__":
    run_tests()
