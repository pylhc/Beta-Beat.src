import sys
import os
import regression
import compare_utils

ABS_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

MYSELF = os.path.join("tests", "regression", "example.py")
CURRENT_DIR = os.path.dirname(MYSELF)

TEST_CASES = (
    regression.TestCase(name="should_succeed",
                        script=MYSELF,
                        arguments=os.path.join(CURRENT_DIR, "_test_out_success"),
                        output=os.path.join(CURRENT_DIR, "_test_out_success"),
                        test_function=compare_utils.compare_dirs_with),
    regression.TestCase(name="should_fail",
                        script=MYSELF,
                        arguments=os.path.join(CURRENT_DIR, "_test_out_fail"),
                        output=os.path.join(CURRENT_DIR, "_test_out_fail"),
                        test_function=lambda dir1, dir2: False),
    regression.TestCase(name="should_raise",
                        script=MYSELF,
                        arguments="valid_string make_it_raise",
                        output=os.path.join("tests", "regression", "_test_out_success"),
                        test_function=compare_utils.compare_dirs_with),
)


def launch_examples():
    regression.launch_test_set(TEST_CASES, ABS_ROOT)


def _fake_test(args):
    # This should raise if len(args) > 2
    _, output = args
    print(ABS_ROOT, output)
    output = os.path.join(ABS_ROOT, output)
    if not os.path.isdir(output):
        os.makedirs(output)
    with open(os.path.join(output, "test_file.txt"), "w") as test_file:
        test_file.write("This is a file!")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        launch_examples()
    else:
        _fake_test(sys.argv)
