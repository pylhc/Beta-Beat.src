""" Entry Point Decorator

Allows a function to be decorated as entry point.
This function will then automatically accept console parameters, config files, json files,
kwargs and dictionaries as input and will parse this input according to the parameters
given to the EntryPoint-Decorator.


kwarguments to use:
        name: name of the variable (e.g. later use options.NAME)
        required: True or False
        default: Default value, if variable not present
        help: String
        type: Value type (if nargs is given, set to list for dicts!)
        choices: choices to choose from (choices need to be of type, if given)

        argparse only:
                flags: commandline flag, e.g. "--file"
                nargs: number of arguments to consume (if given convert entries to type)





FAR FROM FINISHED!!!!
DO NOT USE YET!!!
I WILL REMOVE THIS DISCLAIMER ONCE IT IS FINISHED, I JUST WANTED TO BACKUP MY PROGRESS HERE!!
cheers, Joschua
"""

import ConfigParser
import json
from argparse import ArgumentParser
from Utilities import logging_tools as logtools
from Utilities.dict_tools import DictParser
from Utilities.dict_tools import DotDict
from Utilities.dict_tools import Argument
from Utilities.dict_tools import ArgumentError

LOG = logtools.get_logger(__name__)


ID_CONFIG = "entry_cfg"
ID_CONFIG_SECTION = "section"
ID_DICT = "entry_dict"
ID_JSON = "entry_json"


# TODO:
# set up dict_parser (only accept kw arguments)
# set up commandline parser
# if kw_argument = entry_cfg = -> parse config file

"""
======================== EntryPoint Errors ========================
"""


class EntryError(Exception):
    pass


"""
======================== EntryPoint Wrapper ========================
"""


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
            if len(args) > 0 and len(kwargs) > 0:
                raise EntryError("Cannot combine positional arguments with keyword arguments.")

            if len(args) > 1:
                raise EntryError("Only one positional argument allowed (dict or config file)")

            dict_p = self.dictparse
            arg_p = self.argparse

            if args:
                if isinstance(args[0], basestring):
                    # assume config file
                    options = dict_p.parse_config_items(self._read_config(args[0]))
                elif isinstance(args[0], dict):
                    # dictionary
                    options = dict_p.parse_options(args[0])
                else:
                    raise EntryError("Only dictionary or configfiles"
                                     "are allowed as positional arguments")
            elif len(kwargs) == 0:
                # assume commandline arguments
                options = arg_p.parse_args()
            else:
                if ID_CONFIG in kwargs:
                    if len(kwargs) > 2 or (len(kwargs) == 2 and ID_CONFIG_SECTION not in kwargs):
                        raise EntryError(
                            "Only '{:s}' and '{:s}'".format(ID_CONFIG, ID_CONFIG_SECTION) +
                            " arguments are allowed, when using a config file.")
                    options = self._read_config(kwargs[ID_CONFIG],
                                                kwargs.get(ID_CONFIG_SECTION, None))
                    options = dict_p.parse_config_items(options)

                elif ID_DICT in kwargs:
                    if len(kwargs) > 1:
                        raise EntryError("Only one argument allowed when using a dictionary")
                    options = dict_p.parse_options(kwargs[ID_DICT])

                elif ID_JSON in kwargs:
                    if len(kwargs) > 1:
                        raise EntryError("Only one argument allowed when using a json-file")
                    with open(kwargs[ID_JSON], 'r') as json_file:
                        json_dict = json.load(json_file)
                    options = dict_p.parse_options(json_dict)

                else:
                    options = dict_p.parse_options(kwargs)
            return func(options)
        return wrapper

    @staticmethod
    def _create_argument_parser(param):
        parser = ArgumentParser()
        #TODO: add arguments
        return parser

    @staticmethod
    def _create_dict_parser(param):
        if isinstance(param, dict):
            parser = DictParser()
        # TODO: add arguments
        return parser

    def _read_config(self, cfgfile_path, section=None):
        """ Get content from config file"""
        cfgparse = self.configparse
        dictparse = self.dictparse

        # config File
        with open(cfgfile_path) as config_file:
            cfgparse.readfp(config_file)

        cfg_items = {}
        if section:
            cfg_items[section] = cfgparse.items(section)
        else:
            for section in cfgparse.sections():
                cfg_items[section] = cfgparse.items(section)

        return dictparse.parse_config_items(cfg_items)


"""
======================== EntryPoint Arguments ========================
"""


class EntryPointArguments(object):
    """ Helps to build a simple dictionary structure via add_argument functions.

    You really don't need that, but old habits die hard."""
    def __init__(self):
        pass  #TODO

    def add_argument(self, **kwargs):
        """ Add arguments """
        pass  #TODO


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
