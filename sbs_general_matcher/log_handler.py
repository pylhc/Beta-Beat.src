import sys
import os
import logging
from logging import StreamHandler
from logging.handlers import RotatingFileHandler


def set_up_console_logger(logger):
    logger.setLevel(logging.DEBUG)

    console_handler = StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(_get_formatter())

    logger.addHandler(console_handler)
    return logger


def add_file_handler(logger, log_dir):
    log_file = os.path.join(log_dir, "log.txt")
    file_handler = RotatingFileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(_get_formatter())
    logger.addHandler(file_handler)
    logger.debug("Set up file handler")


def _get_formatter():
    return logging.Formatter("%(name)s - %(levelname)s - %(message)s")


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
