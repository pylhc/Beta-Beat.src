from __future__ import print_function
import sys
import os
import subprocess
import time
import traceback
import contextlib
import tempfile
import shutil
from cStringIO import StringIO
from collections import namedtuple
import re
import git
import yaml


_PYTHON = sys.executable


class TestCase(namedtuple(
        "TestCase",
        ("name", "script", "arguments", "output", "test_function", "pre_hook")
)):
    """Data class to hold information about the test to run.

    Attributes:
        name: A string to identify the test case.
        script: Path, relative to the root of the repository, of the script to
            run.
        arguments: Command line arguments to pass to the script.
        output: Path, relative to the root of the repository, where to expect
            the output of the script.
        test_function: Function of signature (path, path) -> bool to compare
            the results. It must return true only if the test passed without
            issues. The script will be called with the repository root as
            working directory, so all path must take this into account.
        pre_hook: If not None, this function will be called receiving the root
            of each repository as parameter. Useful to perform pre-test
            operations, like creating the output dirs.
    """
    __slots__ = ()  # This makes the class lightweight and immutable.


def launch_test_set(test_cases, repo_path,
                    yaml_conf=None, tag_regexp=None):
    """
    """
    _print_sep("Test session starts")
    commit_hexsha = _find_commit(repo_path, yaml_conf, tag_regexp)
    result = False
    with _temporary_dir() as new_repo_path, Timer() as regression_timer:
        print("Cloning repository into {}".format(new_repo_path))
        clone_revision(repo_path, commit_hexsha, new_repo_path)
        result = find_regressions(test_cases, new_repo_path, repo_path)
    _print_sep("Regression finished in {:.2f}s"
               .format(regression_timer.get_duration()))
    if not result:
        raise RegressionTestFailed()


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
    summary = "Testing...: "
    sys.stdout.flush()
    for test_case in test_cases:
        print("\r{} (running: '{}')".format(summary, test_case.name), end="")
        sys.stdout.flush()
        result = run_test_case(test_case, valid_path, test_path)
        results.append(result)
        summary += result.get_microsummary()
    print("\r{}\n".format(summary))
    _print_sep("Test results")
    report = ""
    for result in results:
        report += result.get_full_summary()
    print(report)
    return all([result.is_success for result in results])


def run_test_case(test_case, valid_path, test_path):
    """Launch single test_case and compare the results.

    Runs the given test, comparing the repositories at valid_path against
    test_path. After the test is run the output directories will be
    deleted.

    Arguments:
        test_case: The test case to run.
        valid_path: The reference directory to compare against.
        test_path: The test directory to check for regressions.
    Returns:
        TODO
    """
    valid_script = os.path.join(valid_path, test_case.script)
    test_script = os.path.join(test_path, test_case.script)
    valid_outpath = os.path.join(valid_path, test_case.output)
    test_outpath = os.path.join(test_path, test_case.output)
    try:
        with Timer() as test_timer:
            if test_case.pre_hook:
                test_case.pre_hook(valid_path)
                test_case.pre_hook(test_path)
            valid_result = _launch_command(valid_path, valid_script, test_case.arguments)
            test_result = _launch_command(test_path, test_script, test_case.arguments)
        result = TestResult(test_case, valid_result, test_result,
                            valid_outpath, test_outpath, test_timer.get_duration())
    finally:
        _remove_if_exists(valid_outpath)
        _remove_if_exists(test_outpath)
    return result


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


def _find_commit(repo, yaml_conf, tag_regexp):
    repo = _force_repo(repo)
    if yaml_conf and tag_regexp:
        raise TestError("Only one of yaml_conf or tag_regexp is allowed.")
    if yaml_conf:
        commit = commit_from_yaml(repo, yaml_conf)
        print("Testing against commit: {commit.summary} ({commit.hexsha})"
              .format(commit=commit))
        return commit.hexsha
    if tag_regexp:
        test_tag = find_tag(repo, tag_regexp)
        print("Testing against tag \"{tag}\":\n"
              "    -> {tag.commit.summary} ({tag.commit.hexsha})"
              .format(tag=test_tag))
        return test_tag.commit.hexsha
    raise TestError("No yaml_conf or tag_regexp given.")


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


def commit_from_yaml(repo, yaml_conf):
    """Returns the commit object configured in yaml_conf.

    Reads the YAML file at 'yaml_conf', searches for a 'regression' category
    and a ref_commit under it, containing the hash of a commit in 'repo' and
    return the associated commit object.

    Arguments:
        repo: The Repo object, url or path to use.
        yaml_conf: The path to the configuration YAML file.
    Returns:
        The commit associated with the hash given in yaml_conf.
    Raises:
        TestError: If no 'regression' category or 'ref_commit' is found in the
            YAML file.
        ValueError: If the given commit hash is not in the repository.
    """
    repo = _force_repo(repo)
    with open(yaml_conf, "r") as yaml_file:
        yaml_data = yaml.load(yaml_file)
    try:
        strhash = yaml_data["regression"]["ref_commit"]
    except KeyError:
        raise TestError(
            "Unable to find 'ref_commit' under 'regression' in file {}"
            .format(yaml_conf)
        )
    commit = repo.commit(strhash)
    return commit


class TestError(Exception):
    """Raised when an error in the test configuration happens.
    """
    pass


class TestResult(object):
    """
    TODO
    """

    RunResult = namedtuple("RunResult", ("stdout", "stderr", "raised"))

    def __init__(self, test_case,
                 valid_result, test_result, valid_output, test_output, duration=None):
        self.test_case = test_case
        self.valid_result = valid_result
        self.test_result = test_result
        self.valid_output = valid_output
        self.test_output = test_output
        self.duration = duration
        self._compare_stdout = None
        self._compare_stderr = None
        self._compare_error = None
        self.is_exception = self.valid_result.raised or self.test_result.raised
        self.is_regression = self._check_regression()
        self.is_success = not self.is_exception and not self.is_regression

    def get_name(self):
        """Returns the name of the test case.
        """
        return self.test_case.name

    def get_microsummary(self):
        """Returns a small summary: . -> success, F -> failed, E -> error.
        """
        if self.is_exception:
            return "E"
        return "." if not self.is_regression else "F"

    def get_full_summary(self):
        """Returns a full summary of this tests result.
        """
        report = "-> Test {}: ".format(self.test_case.name)
        if self.is_exception:
            report += "errored\n"
            if self._compare_error:
                report += self._compare_error
            else:
                for name, run_res in (("Valid", self.valid_result),
                                      ("Test", self.test_result)):
                    if run_res.raised:
                        report += "{} case execution failed.\n".format(name)
                        if run_res.stdout != "":
                            report += "Standard output:\n{}\n".format(run_res.stdout)
                        if run_res.stderr != "":
                            report += "Error output:\n{}\n".format(run_res.stderr)
                        break
        elif self.is_regression:
            report += "failed, differences found.\n"
            if self._compare_stdout != "":
                report += "Standard output:\n{}\n".format(self._compare_stdout)
            if self._compare_stderr != "":
                report += "Error output:\n{}\n".format(self._compare_stderr)
        else:
            report += "succeeded"
            if self.duration is not None:
                report += " ({:.2f}s)".format(self.duration)
            report += "\n"
        report += "\n"
        return report

    def _check_regression(self):
        if self.is_exception:
            return True
        old_stdout, old_stderr = sys.stdout, sys.stderr
        try:
            sys.stdout = mystdout = StringIO()
            sys.stderr = mystderr = StringIO()
            compare_res = not self.test_case.test_function(
                self.valid_output,
                self.test_output,
            )
            self._compare_stdout = mystdout.getvalue()
            self._compare_stderr = mystderr.getvalue()
        # User provided function, I can't know what will raise...
        except:  # pylint: disable=W0702
            self._compare_error = (
                "An exception happened while comparing the directories:\n" +
                traceback.format_exc()
            )
            self.is_exception = True
            return True
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
        return compare_res


class RegressionTestFailed(Exception):
    """Raised when the regression test fails
    """
    pass


def _force_repo(repo):
    if isinstance(repo, git.Repo):
        return repo
    return git.Repo(repo)


def _launch_command(repo_root, script, args):
    if isinstance(args, str):
        args = args.split(" ")
    comm = subprocess.Popen([_PYTHON] + [script] + args,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            cwd=repo_root)
    stdout, stderr = comm.communicate()
    raised = comm.returncode != 0
    return TestResult.RunResult(stdout=stdout, stderr=stderr, raised=raised)


def _remove_if_exists(dir_path):
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)


def _print_sep(text=None, length=80):
    if text is None:
        print("=" * length)
    else:
        left_len = (length - len(text) - 2)/2
        right_len = length - left_len - len(text) - 2
        print("=" * left_len + " " + text + " " + "=" * right_len)


# Contexts ####################################################################


@contextlib.contextmanager
def _temporary_dir():
    try:
        dir_path = tempfile.mkdtemp()
        yield dir_path
    finally:
        shutil.rmtree(dir_path)


class Timer(object):
    def __enter__(self):
        self.start = time.time()
        self.finish = None
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish = time.time()

    def get_duration(self):
        if self.finish is None:
            return time.time() - self.start
        return self.finish - self.start
