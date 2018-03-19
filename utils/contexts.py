import sys
import os
import time
import warnings
import pandas as pd
from logging_tools import NEWLINE
from inspect import currentframe
from contextlib import contextmanager


@contextmanager
def log_out(stdout=sys.stdout, stderr=sys.stderr):
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = stdout
    sys.stderr = stderr
    try:
        yield
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr


@contextmanager
def silence():
    devnull = open(os.devnull, "w")
    with log_out(stdout=devnull, stderr=devnull):
        try:
            yield
        finally:
            devnull.close()


@contextmanager
def timeit(function):
    start_time = time.time()
    try:
        yield
    finally:
        time_used = time.time() - start_time
        function(time_used)


@contextmanager
def suppress_warnings(warning_classes):
    with warnings.catch_warnings(record=True) as warn_list:
        yield
    for w in warn_list:
        if not issubclass(w.category, warning_classes):
            warnings.warn(w)


@contextmanager
def log_pandas_settings_with_copy(log_func):
    caller_line = currentframe().f_back.f_back.f_lineno  # one frame for contextmanager
    old_mode = pd.options.mode.chained_assignment
    pd.options.mode.chained_assignment = 'warn'
    try:
        with warnings.catch_warnings(record=True) as warn_list:
            yield
        for w in warn_list:
            if not issubclass(w.category, pd.core.common.SettingWithCopyWarning):
                warnings.warn(w)
            else:
                message = w.message.args[0].split("\n")
                log_func("{:s} (l. {:d})".format(message[1], caller_line))
    finally:
        pd.options.mode.chained_assignment = old_mode
