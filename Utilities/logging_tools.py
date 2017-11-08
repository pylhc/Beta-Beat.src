import logging
import inspect
import sys
import os


def get_logger(name, filename=None, level_main=logging.DEBUG, level_console=logging.INFO):
    """
    Sets up logger if name is __main__, otherwise just returns the logger like logging.getLogger()
    :param name: __name__ of caller
    :param filename: __file__ of caller
    :param level_main: main logging level, default DEBUG
    :param level_console: console logging level, default INFO
    :return:
    """
    if name == "__main__":
        if not filename:
            previous_frame = inspect.currentframe().f_back
            (filename, line_number,
             function_name, lines, index) = inspect.getframeinfo(previous_frame)

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
        file_cropped = os.path.basename(filename)[:-3]
        logfile = os.path.join(os.path.dirname(filename), file_cropped + '.log')
        file_handler = logging.FileHandler(logfile, mode="w", )
        file_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

        return logging.getLogger(file_cropped)
    else:
        # logger for the current file
        return logging.getLogger(name)


def getLogger(name):
    """ Convenience function so the caller does not have to import logging """
    return logging.getLogger(name)
