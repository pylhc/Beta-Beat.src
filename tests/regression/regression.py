from __future__ import print_function
import sys
import os
import traceback
from collections import namedtuple
import re
import git


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_TEST_DIRNAME = "_test_repo"
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


def launch_test_set(test_cases, repo_path):
    """
    """
    test_rev = find_test_revision(repo_path, _TEST_REGEXP)
    print("Testing against revision: {}".format(test_rev))
    new_repo_path = os.path.join(_THIS_DIR, _TEST_DIRNAME)
    print("Cloning repository at revision...")
    clone_revision(repo_path, test_rev, new_repo_path)
    print("Start test cases...")
    results = []
    for test_case in test_cases:
        print("Test case: {}".format(test_case.name))
        try:
            result = run_test_case(test_case, repo_path, new_repo_path)
        except Exception as err:
            result = err
            print("E")
        else:
            print("." if result else "F")
        results.append(result)
    for idx, (result, test_case) in enumerate(zip(results, test_cases)):
        report = "Test #{}: ".format(idx)
        if isinstance(result, Exception):
            report += "errored ->"
            print(err)
            continue
        report += "succeeded" if result else "failed"


def run_test_case(test_case, orig_repo_path, test_repo_path):
    """Launch test_case and compare the results.
    """
    orig_path = os.path.join(orig_repo_path, test_case.script)
    test_path = os.path.join(test_repo_path, test_case.script)
    _launch_command(orig_path, test_case.arguments)
    _launch_command(test_path, test_case.arguments)
    orig_outpath = os.path.join(orig_repo_path, test_case.output)
    test_outpath = os.path.join(test_repo_path, test_case.output)
    return test_case.test_function(orig_outpath, test_outpath)


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
    return new_repo


def find_test_revision(repo, tag_regexp):
    """
    TODO
    """
    repo = _force_repo(repo)
    matcher = re.compile(tag_regexp)
    ordered_tags = sorted(repo.tags,
                          key=lambda tag: tag.commit.committed_datetime,
                          reverse=True)
    for tag in ordered_tags:
        if matcher.findall(tag.name):
            return tag.commit.hexsha
    raise TestError("Can't find a tag that matches {}".format(tag_regexp))


def _force_repo(repo):
    if isinstance(repo, git.Repo):
        return repo
    return git.Repo(repo)


def _launch_command(script, args):
    python = sys.executable
    os.system("{python} {script} {args}"
              .format(python=python, script=script, args=args))


class TestError(Exception):
    """Raised when an error in the test configuration happens.
    """
    pass
