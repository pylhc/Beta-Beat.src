import logging
import inspect
import sys
import os


def get_logger(name, logfile=None, level_main=logging.DEBUG, level_console=logging.INFO):
    """
    Sets up logger if name is __main__. Returns logger based on module name)
    :param name: only used to check if __name__ is __main__
    :param logfile: log output file ('__file__.log' of caller if 'None')
    :param level_main: main logging level, default DEBUG
    :param level_console: console logging level, default INFO
    :return:
    """

    # get current module
    caller_file = _get_caller()
    current_module = _get_current_module(caller_file)

    if name == "__main__":
        if not logfile:
            logfile = caller_file + '.log'

        # set up root logger
        root_logger = logging.getLogger("")
        root_logger.setLevel(level_main)

        # print all logs to the console
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level_console)
        console_formatter = logging.Formatter("%(name)s: %(message)s")
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

        # set up log file
        file_handler = logging.FileHandler(logfile, mode="w", )
        file_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

    # logger for the current file
    return logging.getLogger(".".join([current_module, os.path.basename(caller_file)]))


def file_handler(logfile, level=logging.DEBUG, format="%(name)s - %(levelname)s - %(message)s"):
    """ Convenience function so the caller does not have to import logging """
    handler = logging.FileHandler(logfile, mode='w', )
    handler.setLevel(level)
    formatter = logging.Formatter(format)
    handler.setFormatter(formatter)
    return handler


def add_module_handler(handler):
    """ Add handler at current module level """
    current_module = _get_current_module()
    logger = logging.getLogger(current_module)
    logger.addHandler(handler)
    # logger.disabled = True


def getLogger(name):
    """ Convenience function so the caller does not have to import logging """
    return logging.getLogger(name)


def _get_caller():
    """ Find the caller of the current log-function """
    caller_frame = inspect.currentframe().f_back
    this_file, _ = os.path.splitext(__file__)

    caller_file = this_file
    while this_file == caller_file:
        (caller_file_full, _, _, _, _) = inspect.getframeinfo(caller_frame)
        caller_file, _ = os.path.splitext(caller_file_full)
        caller_frame = caller_frame.f_back

    return caller_file


def _get_current_module(current_file=None):
    """ Find the name of the current module """
    if not current_file:
        current_file = _get_caller()

    path_parts = current_file.split(os.path.sep)
    current_module = '.'.join(path_parts[(path_parts.index('Beta-Beat.src') + 1):-1])
    return current_module
