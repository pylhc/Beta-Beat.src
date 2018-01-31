""" Entry Point Decorator

Allows a function to be decorated as entry point.
This function will then automatically accept console parameters, config files, json files,
kwargs and dictionaries as input and will parse this input according to the parameters
given to the EntryPoint-Decorator.


Parameters need to be a list or a dictionary of dictionaries with the following keys:
        name (required): Name of the variable (e.g. later use options.NAME).
                         If 'params' is a dictionary, the key will be used as name.
        flags (required): Commandline flag, e.g. "--file"
        required (optional): Boolean
        default (optional): Default value, if variable not present
        help (optional): String
        type (optional): Value type (if nargs is given, set to list for dicts!)
        choices (optional): choices to choose from (choices need to be of type, if given)
        nargs (optional): number of arguments to consume (commandline only)

The 'strict' option changes the behaviour for unknown parameters:
    strict=True raises exceptions, strict=False loggs debug messages.
"""

import ConfigParser
import copy
import json
import argparse
from argparse import ArgumentParser
from Utilities import logging_tools as logtools
from Utilities.dict_tools import DictParser
from Utilities.dict_tools import DotDict
from Utilities.dict_tools import REMAINDER

LOG = logtools.get_logger(__name__)


ID_CONFIG = "entry_cfg"
ID_DICT = "entry_dict"
ID_JSON = "entry_json"
ID_SECTION = "section"

ARGPARSE_ONLY = ["nargs", "flags"]

"""
======================== EntryPoint Errors ========================
"""


class EntryError(Exception):
    pass


"""
======================== EntryPoint Wrapper ========================
"""


class EntryPoint(object):
    def __init__(self, arguments, strict=True):
        """ Initialize decoration: Handle the desired input parameters. """
        self.strict = strict

        # add argument dictionary to EntryPoint
        self.remainder = None
        self.arguments = EntryPoint._dict2list_args(arguments)
        self._check_arguments()

        # create parsers from arguments
        self.argparse = self._create_argument_parser()
        self.dictparse = self._create_dict_parser()
        self.configparse = self._create_config_parser()

    def __call__(self, func):
        """ Builds wrapper around the function 'func' (called on decoration)

        Whenever the decorated function is called, actually this wrapper is called
        """
        def wrapper(*args, **kwargs):
            return func(self.parse(*args, **kwargs))
        return wrapper

    def parse(self, *args, **kwargs):
        """ Parse whatever input arguments come """
        if len(args) > 0 and len(kwargs) > 0:
            raise EntryError("Cannot combine positional arguments with keyword arguments.")

        if len(args) > 1:
            raise EntryError("Only one positional argument allowed (dict or config file).")

        if args:
            options = self._handle_arg(args[0])
        elif len(kwargs) > 0:
            options = self._handle_kwargs(kwargs)
        else:
            options = self._handle_commandline()
        return options

    #########################
    # Create Parsers
    #########################

    def _create_argument_parser(self):
        """ Creates the ArgumentParser from arguments. """
        parser = ArgumentParser()
        parser = add_arguments_to_generic(parser, self.arguments)
        return parser

    def _create_dict_parser(self):
        """ Creates the DictParser from arguments. """
        parser = DictParser(strict=self.strict)
        parser = add_arguments_to_generic(parser, self.arguments)
        return parser

    def _create_config_parser(self):
        """ Creates the config parser. Maybe more to do here later with arguments. """
        parser = ConfigParser.ConfigParser()
        return parser

    #########################
    # Handlers
    #########################

    def _handle_commandline(self, args=None):
        """ No input to function """
        options, unknown_args = self.argparse.parse_known_args(args)
        options = DotDict(vars(options))
        if unknown_args:
            if self.remainder:
                options[self.remainder] = unknown_args
            elif self.strict:
                raise EntryError("Unknown arguments: {:s}".format(unknown_args))
            else:
                LOG.debug("Unknown arguments skipped: {:s}".format(unknown_args))
        return options

    def _handle_arg(self, arg):
        """ *args has been input """
        if isinstance(arg, basestring):
            # assume config file
            options = self.dictparse.parse_config_items(self._read_config(arg))
        elif isinstance(arg, dict):
            # dictionary
            options = self.dictparse.parse_options(arg)
        elif isinstance(arg, list):
            # list of commandline arguments
            options = self._handle_commandline(arg)
        else:
            raise EntryError("Only dictionary or configfiles"
                             "are allowed as positional arguments")
        return options

    def _handle_kwargs(self, kwargs):
        """ **kwargs been input """
        if ID_CONFIG in kwargs:
            if len(kwargs) > 2 or (len(kwargs) == 2 and ID_SECTION not in kwargs):
                raise EntryError(
                    "Only '{:s}' and '{:s}'".format(ID_CONFIG, ID_SECTION) +
                    " arguments are allowed, when using a config file.")
            options = self._read_config(kwargs[ID_CONFIG],
                                        kwargs.get(ID_SECTION, None))
            options = self.dictparse.parse_config_items(options)

        elif ID_DICT in kwargs:
            if len(kwargs) > 1:
                raise EntryError("Only one argument allowed when using a dictionary")
            options = self.dictparse.parse_options(kwargs[ID_DICT])

        elif ID_JSON in kwargs:
            if len(kwargs) > 2 or (len(kwargs) == 2 and ID_SECTION not in kwargs):
                raise EntryError(
                    "Only '{:s}' and '{:s}'".format(ID_JSON, ID_SECTION) +
                    " arguments are allowed, when using a json file.")
            with open(kwargs[ID_JSON], 'r') as json_file:
                json_dict = json.load(json_file)

            if ID_SECTION in kwargs:
                json_dict = json_dict[kwargs[ID_SECTION]]

            options = self.dictparse.parse_options(json_dict)

        else:
            options = self.dictparse.parse_options(kwargs)

        return options

    #########################
    # Helpers
    #########################

    def _check_arguments(self):
        """ EntryPoint specific checks for arguments """
        for arg in self.arguments:
            arg_name = arg.get("name", None)
            if arg_name is None:
                raise AttributeError("An Argument needs a Name!")

            if arg.get("nargs", None) == argparse.REMAINDER:
                raise AttributeError("Argument '{:s}' is set as remainder.".format(arg_name) +
                                     "This method is really buggy," +
                                     "use 'type=entrypoint.REMAINDER' instead.")

            if arg.get("type", None) == REMAINDER:
                if self.remainder is not None:
                    raise AttributeError("More than one remainder argument found!")
                self.remainder = arg_name

                if arg.get("flags", None) is not None:
                    raise AttributeError("Argument '{:s}'".format(arg_name) +
                                         "is set to collect remaining arguments." +
                                         "It is not supposed to have flags.")

            elif arg.get("flags", None) is None:
                raise AttributeError("Argument '{:s}'".format(arg_name) +
                                     "does not have flags. " +
                                     "Only remainder arguments are allowed that privilege")

    def _read_config(self, cfgfile_path, section=None):
        """ Get content from config file"""
        cfgparse = self.configparse

        with open(cfgfile_path) as config_file:
            cfgparse.readfp(config_file)

        sections = cfgparse.sections()
        if not section and len(sections) == 1:
            section = sections[0]
        elif not section:
            raise EntryError("'{:s}' contains multiple sections. Please specify one!")

        return cfgparse.items(section)

    @staticmethod
    def _dict2list_args(param):
        """ Convert dictionary to list and add name by key """
        if isinstance(param, dict):
            out = []
            for key in param:
                item = param[key]
                item["name"] = key
                out.append(item)
            return out
        else:
            return param


"""
======================== EntryPoint Arguments ========================
"""


class EntryPointArguments(DotDict):
    """ Helps to build a simple dictionary structure via add_argument functions.

    You really don't need that, but old habits die hard."""
    def add_argument(self, **kwargs):
        """ Add arguments """
        name = kwargs.pop("name")
        if name in self:
            raise ValueError("'{:s}' is already an argument.".format(name))
        else:
            self[name] = kwargs


"""
======================== Public Helpers ========================
"""


def add_arguments_to_generic(parser, args):
    """ Adds entry-point style arguments to either
    ArgumentParser, DictParser or EntryPointArguments
    """
    args = copy.deepcopy(args)
    if isinstance(parser, EntryPointArguments):
        for arg in args:
            parser.add_argument(arg)

    elif isinstance(parser, ArgumentParser):
        for arg in args:
            arg["dest"] = arg.pop("name", None)
            flags = arg.pop("flags", None)
            if flags is None:
                parser.add_argument(**arg)
            else:
                if isinstance(flags, basestring):
                    flags = [flags]
                parser.add_argument(*flags, **arg)

    elif isinstance(parser, DictParser):
        for arg in args:
            if "nargs" in arg and arg["nargs"] != "?":
                arg["subtype"] = arg.get("type", None)
                arg["type"] = list

            for label in ARGPARSE_ONLY:
                arg.pop(label, None)

            name = arg.pop("name")
            parser.add_argument(name, **arg)
    else:
        raise TypeError("Parser not recognised.")
    return parser

"""
======================== Script Mode ========================
"""


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))


