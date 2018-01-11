""" Entry Point Decorator

Allows a function to be decorated as entry point.
This function will then automatically accept console parameters, config files, json files,
kwargs and dictionaries as input and will parse this input according to the parameters
given to the EntryPoint-Decorator.


FAR FROM FINISHED!!!!
DO NOT USE YET!!!
I WILL REMOVE THIS DISCLAIMER ONCE IT IS FINISHED, I JUST WANTED TO BACKUP MY PROGRESS HERE!!
cheers, Joschua
"""

import ConfigParser
from argparse import ArgumentParser
from Utilities import logging_tools as logtools
from Utilities.dict_tools import DictParser
from Utilities.dict_tools import DotDict

LOG = logtools.get_logger(__name__)


ID_CONFIG = "entry_cfg"
ID_CONFIG_SECTION = "section"
ID_DICT = "entry_dict"
ID_JSON = "entry_json"


# TODO:
# set up dict_parser (only accept kw arguments)
# set up commandline parser
# if kw_argument = entry_cfg = -> parse config file

class EntryError(Exception):
    pass


class EntryPoint(object):
    def __init__(self, param):
        """ Initialize decoration: Handle the desired input parameters. """
        self.argparse = self._create_argument_parser(param)
        self.dictparse = self._create_dict_parser(param)
        self.configparse = ConfigParser.ConfigParser()

    def __call__(self, func):
        """ Builds wrapper around the function 'func' (called on decoration)

        Whenever the decorated function is called, actually this wrapper is called
        """
        def wrapper(*args, **kwargs):
            """ Parse the arguments to the decorated function accordingly """
            options = DotDict({})
            if len(args) > 0 and len(kwargs) > 0:
                raise EntryError("Cannot combine positional arguments with keyword arguments.")

            if len(args) > 1:
                raise EntryError("Only one positional argument allowed (dict or config file)")

            dict_p = self.dictparse
            cfg_p = self.configparse
            arg_p = self.argparse

            if args:
                if isinstance(args[0], basestring):
                    # assume config file
                    pass
                elif isinstance(args[0], dict):
                    # dictionary
                    pass
                else:
                    raise EntryError("Only dictionary or configfiles"
                                     "are allowed as positional arguments")

            elif len(kwargs) == 0:
                # assume commandline arguments
                arg_p.parse_args()
            else:
                if ID_CONFIG in kwargs:
                    with open(kwargs[ID_CONFIG]) as config_file:
                        cfg_p.readfp(config_file)

                    if ID_CONFIG_SECTION in kwargs:
                        cfg_items = cfg_p.items(kwargs[ID_CONFIG_SECTION])
                    else:
                        cfg_items = []
                        for section in cfg_p.sections():
                            cfg_items += cfg_p.items(section)
                dict_p.parse_config_items(cfg_items)


                dict_p.parse_options(kwargs)

            for p in self.param:
                options[p] = kwargs.pop(p, self.param[p])
                return func(options)
        return wrapper

    @staticmethod
    def _create_argument_parser(param):
        parser = ArgumentParser()
        return parser

    @staticmethod
    def _create_dict_parser(param):
        parser = DictParser()
        return parser


"""
======================== Example Usage ========================
"""


def _param_definitions():
    return {"x": 7}


@EntryPoint(_param_definitions())
def some_function(options):
    # options can be handled like an object
    return options.x * 8


if __name__ == "__main__":
    print some_function(x=6)
