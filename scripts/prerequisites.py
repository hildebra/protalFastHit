"""prerequisites.py - skips of the Python test suites for what a machine lacks, failures under PROTAL_TESTS_REQUIRED=1.

A test whose prerequisite is missing (scikit-learn for the trainer, $PROTAL, art_illumina, git, ...) is skipped, so
that each suite runs anywhere. With PROTAL_TESTS_REQUIRED=1 in the environment (CI sets it, as for tests/e2e) a missing
prerequisite is a failure instead: a test that should run cannot pass unseen as a skip.

  import prerequisites
  @prerequisites.requires(HAVE_SKLEARN, "needs scikit-learn")      # a class or a test method
  class T(unittest.TestCase): ...
  prerequisites.missing("no art_illumina")                         # in a test or setUpClass: skip, or fail
"""

import functools
import os
import shutil
import unittest

REQUIRED = os.environ.get("PROTAL_TESTS_REQUIRED", "") not in ("", "0")


def missing(reason):
    """Skip the test (or, in setUpClass, its class) for `reason`; fail with it under PROTAL_TESTS_REQUIRED."""
    if REQUIRED:
        raise AssertionError(f"missing prerequisite: {reason} (PROTAL_TESTS_REQUIRED is set)")
    raise unittest.SkipTest(reason)


def requires(condition, reason):
    """unittest.skipUnless(condition, reason) for a test class or method; under PROTAL_TESTS_REQUIRED a missing
    prerequisite fails the class (its setUpClass) or the test instead."""
    if condition:
        return lambda item: item
    if not REQUIRED:
        return unittest.skip(reason)

    def failing(item):
        if isinstance(item, type):
            item.setUpClass = classmethod(lambda cls: missing(reason))
            return item

        @functools.wraps(item)
        def test(*args, **kwargs):
            missing(reason)
        return test
    return failing


def executable(path):
    """Whether `path` (e.g. $PROTAL) names an executable file."""
    return bool(path) and os.path.isfile(path) and os.access(path, os.X_OK)


def on_path(*commands):
    """Whether every command is on the PATH."""
    return all(shutil.which(c) for c in commands)
