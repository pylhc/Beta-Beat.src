""" Entry Point Decorator

Allows a function to be decorated as entry point.
This function will then automatically accept console arguments, config files, json files,
kwargs and dictionaries as input and will parse it according to the parameters
given to the EntryPoint-Decorator.

Terminology:
    Parameter - Items containing info on how to parse Arguments
    Arguments - The input to the wrapped-function
    Options - The parsed arguments and hence the options of the function

Hence, an ArgumentError will be raised in case of something going wrong during parsing,
while ParameterErrors will be raised if something goes wrong when adding parameters to the list.

It is also possible to use the EntryPoint Class similar to a normal parser:
    ep_parser = EntryPoint(parameters)
    options, unknown = ep_parser.parse(arguments)

    NOT TO DO: "ep_parser(arguments)"
    This invokes the wrapping behaviour

Parameters need to be a list or a dictionary of dictionaries with the following keys:
        name (required): Name of the variable (e.g. later use options.NAME).
                         If 'params' is a dictionary, the key will be used as name.
        flags (required): Commandline flag, e.g. "--file"
        required (optional): Boolean
        default (optional): Default value, if variable not present
        help (optional): String
        type (optional): Value type (if nargs is given, set to list for dicts!)
        choices (optional): choices to choose from (choices need to be of type, if given)
        nargs (optional): number of arguments to consume (commandline only, do not use REMAINDER!)


The 'strict' option changes the behaviour for unknown parameters:
    strict=True raises exceptions, strict=False loggs debug messages and returns the options
    Hence a wrapped function with strict=True must accept one input, with strict=False two.
    Default: False


"""

import ConfigParser
import copy
import json
import argparse
from argparse import ArgumentParser
from Utilities import logging_tools as logtools
from Utilities.dict_tools import DictParser
from Utilities.dict_tools import DotDict
from Utilities.dict_tools import ArgumentError
from Utilities.dict_tools import ParameterError

LOG = logtools.get_logger(__name__)


ID_CONFIG = "entry_cfg"
ID_DICT = "entry_dict"
ID_JSON = "entry_json"
ID_SECTION = "section"

ARGPARSE_ONLY = ["nargs", "flags"]


"""
======================== EntryPoint Wrapper ========================
"""


class EntryPoint(object):
    def __init__(self, parameter, strict=False):
        """ Initialize decoration: Handle the desired input parameter. """
        self.strict = strict

        # add argument dictionary to EntryPoint
        self.remainder = None
        self.parameter = EntryPoint._dict2list_param(parameter)
        self._check_parameter()

        # create parsers from parameter
        self.argparse = self._create_argument_parser()
        self.dictparse = self._create_dict_parser()
        self.configparse = self._create_config_parser()

    def __call__(self, func):
        """ Builds wrapper around the function 'func' (called on decoration)

        Whenever the decorated function is called, actually this wrapper is called
        """
        if self.strict:
            def wrapper(*args, **kwargs):
                return func(self.parse(*args, **kwargs))
        else:
            def wrapper(*args, **kwargs):
                options, unknown_options = self.parse(*args, **kwargs)
                return func(options, unknown_options)
        return wrapper

    def parse(self, *args, **kwargs):
        """ Parse whatever input parameter come.

            This is the heart of EntryPoint and will recognize the input and parse it
            accordingly.
            Allowed inputs are:
                - Dictionary with arguments as key-values
                - Key-Value arguments
                - Path to a one-section config file
                - Commandline Arguments
                - Commandline Arguments in string-form (as list)
                - Special Key-Value arguments are:
                    entry_dict: Value is a dict of arguments
                    entry_cfg: Path to config file
                    entry_json: Path to json file
                    section: Section to use in config file, or subdirectory to use in json file.
                             Only works with the key-value version of config file.
                             If not given only one-section config files are allowed.
         """
        if len(args) > 0 and len(kwargs) > 0:
            raise ArgumentError("Cannot combine positional parameter with keyword parameter.")

        if len(args) > 1:
            raise ArgumentError("Only one positional argument allowed (dict or config file).")

        if args:
            options = self._handle_arg(args[0])
        elif len(kwargs) > 0:
            options = self._handle_kwargs(kwargs)
        else:
            options = self._handle_commandline()

        return options  # options might include known and unknown options

    #########################
    # Create Parsers
    #########################

    def _create_argument_parser(self):
        """ Creates the ArgumentParser from parameter. """
        parser = ArgumentParser()
        parser = add_params_to_generic(parser, self.parameter)
        return parser

    def _create_dict_parser(self):
        """ Creates the DictParser from parameter. """
        parser = DictParser(strict=self.strict)
        parser = add_params_to_generic(parser, self.parameter)
        return parser

    def _create_config_parser(self):
        """ Creates the config parser. Maybe more to do here later with parameter. """
        parser = ConfigParser.ConfigParser()
        return parser

    #########################
    # Handlers
    #########################

    def _handle_commandline(self, args=None):
        """ No input to function """
        options, unknown_opts = self.argparse.parse_known_args(args)
        options = DotDict(vars(options))
        if self.strict:
            if unknown_opts:
                raise ArgumentError("Unknown options: {:s}".format(unknown_opts))
            return options
        else:
            if unknown_opts:
                LOG.debug("Unknown options: {:s}".format(unknown_opts))
            return options, unknown_opts

    def _handle_arg(self, arg):
        """ *args has been input """
        if isinstance(arg, basestring):
            # assume config file
            options = self.dictparse.parse_config_items(self._read_config(arg))
        elif isinstance(arg, dict):
            # dictionary
            options = self.dictparse.parse_arguments(arg)
        elif isinstance(arg, list):
            # list of commandline parameter
            options = self._handle_commandline(arg)
        else:
            raise ArgumentError("Only dictionary or configfiles"
                                "are allowed as positional arguments")
        return options  # options might include known and unknown options

    def _handle_kwargs(self, kwargs):
        """ **kwargs been input """
        if ID_CONFIG in kwargs:
            if len(kwargs) > 2 or (len(kwargs) == 2 and ID_SECTION not in kwargs):
                raise ArgumentError(
                    "Only '{:s}' and '{:s}'".format(ID_CONFIG, ID_SECTION) +
                    " arguments are allowed, when using a config file.")
            options = self._read_config(kwargs[ID_CONFIG],
                                        kwargs.get(ID_SECTION, None))
            options = self.dictparse.parse_config_items(options)

        elif ID_DICT in kwargs:
            if len(kwargs) > 1:
                raise ArgumentError("Only one argument allowed when using a dictionary")
            options = self.dictparse.parse_arguments(kwargs[ID_DICT])

        elif ID_JSON in kwargs:
            if len(kwargs) > 2 or (len(kwargs) == 2 and ID_SECTION not in kwargs):
                raise ArgumentError(
                    "Only '{:s}' and '{:s}'".format(ID_JSON, ID_SECTION) +
                    " arguments are allowed, when using a json file.")
            with open(kwargs[ID_JSON], 'r') as json_file:
                json_dict = json.load(json_file)

            if ID_SECTION in kwargs:
                json_dict = json_dict[kwargs[ID_SECTION]]

            options = self.dictparse.parse_arguments(json_dict)

        else:
            options = self.dictparse.parse_arguments(kwargs)

        return options   # options might include known and unknown options

    #########################
    # Helpers
    #########################

    def _check_parameter(self):
        """ EntryPoint specific checks for parameter """
        for param in self.parameter:
            arg_name = param.get("name", None)
            if arg_name is None:
                raise ParameterError("A Parameter needs a Name!")

            if param.get("nargs", None) == argparse.REMAINDER:
                raise ParameterError("Parameter '{:s}' is set as remainder.".format(arg_name) +
                                     "This method is really buggy, hence it is forbidden.")

            if param.get("flags", None) is None:
                raise ParameterError("Parameter '{:s}'".format(arg_name) +
                                     "does not have flags.")

    def _read_config(self, cfgfile_path, section=None):
        """ Get content from config file"""
        cfgparse = self.configparse

        with open(cfgfile_path) as config_file:
            cfgparse.readfp(config_file)

        sections = cfgparse.sections()
        if not section and len(sections) == 1:
            section = sections[0]
        elif not section:
            raise ArgumentError("'{:s}' contains multiple sections. Please specify one!")

        return cfgparse.items(section)

    @staticmethod
    def _dict2list_param(param):
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


class EntryPointParameters(DotDict):
    """ Helps to build a simple dictionary structure via add_argument functions.

    You really don't need that, but old habits die hard."""
    def add_parameter(self, **kwargs):
        """ Add parameter """
        name = kwargs.pop("name")
        if name in self:
            raise ParameterError("'{:s}' is already a parameter.".format(name))
        else:
            self[name] = kwargs


"""
======================== Public Helpers ========================
"""


def add_params_to_generic(parser, params):
    """ Adds entry-point style parameter to either
    ArgumentParser, DictParser or EntryPointArguments
    """
    params = copy.deepcopy(params)
    if isinstance(parser, EntryPointParameters):
        for param in params:
            parser.add_parameter(param)

    elif isinstance(parser, ArgumentParser):
        for param in params:
            param["dest"] = param.pop("name", None)
            flags = param.pop("flags", None)
            if flags is None:
                parser.add_argument(**param)
            else:
                if isinstance(flags, basestring):
                    flags = [flags]
                parser.add_argument(*flags, **param)

    elif isinstance(parser, DictParser):
        for param in params:
            if "nargs" in param and param["nargs"] != "?":
                param["subtype"] = param.get("type", None)
                param["type"] = list

            for label in ARGPARSE_ONLY:
                param.pop(label, None)

            name = param.pop("name")
            parser.add_parameter(name, **param)
    else:
        raise TypeError("Parser not recognised.")
    return parser

"""
======================== Script Mode ========================
"""


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))


