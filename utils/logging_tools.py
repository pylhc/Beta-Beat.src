import logging
import inspect
import sys
import os
from utils import iotools
import datetime

DIVIDER = "|"
BASIC_FORMAT = '%(levelname)7s {div:s} %(message)s {div:s} %(name)s'.format(div=DIVIDER)
COLOR_LEVEL = '\33[38;2;150;150;255m'
COLOR_MESSAGE = '\33[38;2;255;255;180m'
COLOR_WARN = '\33[38;2;255;161;53m'
COLOR_ERROR = '\33[38;2;216;31;42m'
COLOR_NAME = '\33[38;2;150;150;150m'
COLOR_DIVIDER = '\33[38;2;150;150;150m'
COLOR_RESET = '\33[0m'


NOTSET = logging.NOTSET
DEBUG = logging.DEBUG
INFO = logging.INFO
WARN = logging.WARN
WARNING = logging.WARNING
ERROR = logging.ERROR
CRITICAL = logging.CRITICAL
FATAL = logging.FATAL


# Classes and Contexts #########################################################


class MaxFilter(object):
    """ To get messages only up to a certain level """
    def __init__(self, level):
        self.__level = level

    def filter(self, log_record):
        return log_record.levelno <= self.__level


class TempFile(object):
    """ Context Manager.
    Lets another function write into a temporary file and logs its contents.

    It won't open the file though, so only the files path is returned.

    Args:
        file_path (str): Place to write the tempfile to.
        log_func (func): The function with which the content should be logged (e.g. LOG.info)
    """
    def __init__(self, file_path, log_func):
        self.path = file_path
        self.log_func = log_func

    def __enter__(self):
        return self.path

    def __exit__(self, type, value, traceback):
        try:
            with open(self.path, "r") as f:
                content = f.read()
            self.log_func("{:s}:\n".format(self.path) + content)
        except IOError:
            self.log_func("{:s}: -file does not exist-".format(self.path))
        else:
            os.remove(self.path)


# Public Methods ###############################################################


def get_logger(name, level_root=DEBUG, level_console=INFO, fmt=BASIC_FORMAT):
    """
    Sets up logger if name is __main__. Returns logger based on module name)

    Args:
        name: only used to check if __name__ is __main__
        level_root: main logging level, default DEBUG
        level_console: console logging level, default INFO
        fmt: Format of the logging. For default see BASIC_FORMAT

    Returns:
        Logger instance.
    """
    # get current module
    caller_file = _get_caller()
    current_module = _get_current_module(caller_file)

    if name == "__main__":
        # set up root logger
        root_logger = logging.getLogger("")
        root_logger.setLevel(level_root)

        # print logs to the console
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level_console)
        console_formatter = logging.Formatter(_bring_color(fmt))
        console_handler.setFormatter(console_formatter)
        console_handler.addFilter(MaxFilter(INFO))
        root_logger.addHandler(console_handler)

        # print console warnings
        console_warn_handler = logging.StreamHandler(sys.stdout)
        console_warn_handler.setLevel(max(WARNING, level_console))
        logging.Formatter(_bring_color(fmt, WARNING))
        console_warn_handler.setFormatter(console_formatter)
        console_warn_handler.addFilter(MaxFilter(WARNING))
        root_logger.addHandler(console_warn_handler)

        # print errors to error-stream
        error_handler = logging.StreamHandler(sys.stderr)
        error_handler.setLevel(ERROR)
        error_formatter = logging.Formatter(_bring_color(fmt, ERROR))
        error_handler.setFormatter(error_formatter)
        root_logger.addHandler(error_handler)

    # logger for the current file
    return logging.getLogger(".".join([current_module, os.path.basename(caller_file)]))


def file_handler(logfile, level=DEBUG, format=BASIC_FORMAT):
    """ Convenience function so the caller does not have to import logging """
    handler = logging.FileHandler(logfile, mode='w', )
    handler.setLevel(level)
    formatter = logging.Formatter(format)
    handler.setFormatter(formatter)
    return handler


def add_module_handler(handler):
    """ Add handler at current module level """
    current_module = _get_current_module()
    logging.getLogger(current_module).addHandler(handler)


def add_root_handler(handler):
    """ Add handler at root level """
    logging.getLogger("").addHandler(handler)


def getLogger(name):
    """ Convenience function so the caller does not have to import logging """
    return logging.getLogger(name)


def start_debug_mode():
        """ Setup Logger for debugging mode """
        # get current module
        caller_file = _get_caller()
        current_module = _get_current_module(caller_file)

        logger = logging.getLogger(".".join([current_module, os.path.basename(caller_file)]))

        logger.setLevel(DEBUG)
        logger.debug("Running in Debug-Mode.")

        now = str(datetime.datetime.now().isoformat())
        log_file = os.path.abspath(caller_file).replace(".pyc", "").replace(".py", "") + ".log"
        log_file = os.path.join(os.path.dirname(log_file), now + os.path.basename(log_file))
        logger.debug("Writing log to file '{:s}'.".format(log_file))
        logging.getLogger(current_module).addHandler(file_handler(log_file))


# Private Methods ##############################################################


def _get_caller():
    """ Find the caller of the current log-function """
    this_file, _ = os.path.splitext(__file__)
    caller_file = this_file
    caller_frame = inspect.currentframe()
    while this_file == caller_file:
        caller_frame = caller_frame.f_back
        (caller_file_full, _, _, _, _) = inspect.getframeinfo(caller_frame)
        caller_file, _ = os.path.splitext(caller_file_full)
    return caller_file


def _get_current_module(current_file=None):
    """ Find the name of the current module """
    if not current_file:
        current_file = _get_caller()
    path_parts = os.path.abspath(current_file).split(os.path.sep)
    repo_parts = iotools.get_absolute_path_to_betabeat_root().split(os.path.sep)
    current_module = '.'.join(path_parts[len(repo_parts):-1])
    return current_module


def _bring_color(format_string, colorlevel=INFO):
    """ Adds color to the logs (can only be used in a terminal) """
    if not sys.stdout.isatty():
        # Not a tty. You're being piped or redirected
        return format_string

    level = "%(levelname)"
    message = "%(message)"
    name = "%(name)"
    format_string = format_string.replace(level, COLOR_LEVEL + level)
    
    if colorlevel <= INFO:
        format_string = format_string.replace(message, COLOR_MESSAGE + message)
    elif colorlevel <= WARNING:
        format_string = format_string.replace(message, COLOR_WARN + message)
    else:
        format_string = format_string.replace(message, COLOR_ERROR + message)

    format_string = format_string.replace(name, COLOR_NAME + name)
    format_string = format_string.replace(DIVIDER, COLOR_DIVIDER + DIVIDER)
    format_string = format_string + COLOR_RESET

    return format_string
