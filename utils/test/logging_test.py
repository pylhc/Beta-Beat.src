import sys
import os
from os.path import abspath, join, dirname
new_path = abspath(join(dirname(abspath(__file__)), os.pardir, os.pardir))
if new_path not in sys.path:
    sys.path.append(new_path)
from utils import logging_tools


def _write_messages(log):
    log.debug("I am a debug message.")
    log.debug("Another debug message.")
    log.info("And I am an info message.")
    log.info("Another info message following.")
    log.warning("This is a warning!")
    log.warning("You have been warned.")
    log.error("Dont be afraid, this is only a test for an error.")
    log.error("But this might be a real one! Kidding, it isn't")


def main(level_console=0):
    log = logging_tools.get_logger(__name__, level_console=level_console)
    _write_messages(log)


if __name__ == "__main__":
    main()
