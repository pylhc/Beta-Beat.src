from __future__ import print_function
import sys
import os
import subprocess
import traceback
import contextlib
import tempfile
import shutil
from collections import namedtuple
import re
import git


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_TEST_REGEXP = "^before.*$"


class TestCase(namedtuple(
        "TestCase",
        ("name", "script", "arguments", "output", "test_function")
)):
    """Data class to hold information about the test to run.

    Attributes:
        script: Path, relative to the root of the repository, of the script to
            run.
        arguments: Command line arguments to pass to the script.
        output: Path, relative to the root of the repository, where to expect
            the output of the script.
        test_function: Function of signature (path, path) -> bool to compare
            the results. It must return true only if the test passed without
            issues.
    """
    __slots__ = ()  # This makes the class lightweight and immutable.


def launch_test_set(test_cases, repo_path, tag_regexp=_TEST_REGEXP):
    """
    """
    test_tag = find_tag(repo_path, tag_regexp)
    print("Testing against tag \"{tag}\":\n"
          "    -> {tag.commit.summary} ({tag.commit.hexsha})"
          .format(tag=test_tag))
    with _temporary_dir() as new_repo_path:
        print("Cloning repository into {}".format(new_repo_path))
        clone_revision(repo_path, test_tag.commit.hexsha, new_repo_path)
        print("\n################# Start test cases #################")
        find_regressions(test_cases, new_repo_path, repo_path)


def find_regressions(test_cases, valid_path, test_path):
    """Test the directory test_path for regressions against valid_path.

    It will print the results of each regression search test case in test_path,
    firstly showing a summary like .FE, where '.' means valid, 'F' failed and
    'E' errored. Afterwards it will show a more detailed result, showing
    exception thrown and differences found.

    Attributes:
        test_cases: an iterable of TestCase to test the test_path.
        valid_path: Path to the reference directory to test against.
        test_path: Path to the directory to be tested for regressions.
    """
    results = []
    for test_case in test_cases:
        try:
            result = TestResult(
                result=run_test_case(test_case, valid_path, test_path),
                is_error=False
            )
        # Sorry pylint, but I have to capture every exception:
        except:  # pylint: disable=W0702
            result = TestResult(
                result=traceback.format_exc(),
                is_error=True,
            )
        results.append(result)
        print(result.get_char(), end="")
    print("\n\n################### Test results ###################")
    report = ""
    for result, test_case in zip(results, test_cases):
        report += "Test {}: ".format(test_case.name)
        if result.is_error:
            report += "errored ->\n"
            report += "{}".format(result.result)
        else:
            report += "succeeded" if result else "failed"
        report += "\n\n"
    print(report)


def run_test_case(test_case, valid_path, test_path):
    """Launch single test_case and compare the results.

    Runs the given test, comparing the repositories at valid_path against
    test_path.

    Arguments:
        test_case: The test case to run.
        valid_path: The reference directory to compare against.
        test_path: The test directory to check for regressions.
    Returns:
        Whatever 'test_case.test_function' returns, usually a boolean.
    """
    valid_script = os.path.join(valid_path, test_case.script)
    test_script = os.path.join(test_path, test_case.script)
    _launch_command(valid_script, test_case.arguments)
    _launch_command(test_script, test_case.arguments)
    valid_outpath = os.path.join(valid_path, test_case.output)
    test_outpath = os.path.join(test_path, test_case.output)
    return test_case.test_function(valid_outpath, test_outpath)


def clone_revision(source, revision, to_path):
    """Clone the repository source to certain path.

    Clones the repository present in source (can be a local path or an remote
    source) at given revision in the path given by to_path.

    Arguments:
        source: Local path, remote repository url, or git.Repo instance to
            clone.
        revision: Commit hash or tag to clone.
        to_path: Output path to write the repository.

    Returns:
        The new repository instance.
    """
    try:
        os.makedirs(to_path)
    except os.error:
        pass  # The directory is already there.
    my_repo = _force_repo(source)
    new_repo = my_repo.clone(to_path)
    new_repo.head.reference = new_repo.commit(revision)
    assert new_repo.head.is_detached
    new_repo.head.reset(index=True, working_tree=True)
    return new_repo


def find_tag(repo, tag_regexp):
    """Returns the most recent git tag that matches a regular expression.

    Arguments:
        repo: The Repo object, url or path to use.
        tag_regexp: Regular expression to try and match in the tag names.
    Returns:
        The most recent tag in 'repo' that matches tag_regexp.
    Raises:
        TestError: if 'tag_regexp' doesn't match any tag in 'repo'.
    """
    repo = _force_repo(repo)
    matcher = re.compile(tag_regexp)
    ordered_tags = sorted(repo.tags,
                          key=lambda tag: tag.commit.committed_datetime,
                          reverse=True)
    for tag in ordered_tags:
        if matcher.findall(tag.name):
            return tag
    raise TestError("Can't find a tag that matches {}".format(tag_regexp))


class TestError(Exception):
    """Raised when an error in the test configuration happens.
    """
    pass


class TestResult(object):
    """Hold the results of the tests.

    Attributes:
        result: The result of the test, either a boolean for success or fail,
            or an string containing the stack trace of an Exception.
        is_error: True only if the test errored.
    """
    __slots__ = ("result", "is_error", "get_char")

    def __init__(self, result, is_error):
        self.result = result
        self.is_error = is_error

    def get_char(self):
        if self.is_error:
            return "E"
        return "." if self.result else "F"


def _force_repo(repo):
    if isinstance(repo, git.Repo):
        return repo
    return git.Repo(repo)


def _launch_command(script, args):
    # TODO: Not sure about this... but easiest way of capturing exceptions.
    with _temporary_args(args.split(" ")):
        execfile(script)


@contextlib.contextmanager
def _temporary_dir():
    try:
        dir_path = tempfile.mkdtemp()
        yield dir_path
    finally:
        shutil.rmtree(dir_path)


@contextlib.contextmanager
def _temporary_args(args):
    old_args = sys.argv
    try:
        sys.argv = args
        yield
    finally:
        sys.argv = old_args
