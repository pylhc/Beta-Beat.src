"""

"""
import json
from Utilities import logging_tools
LOG = logging_tools.get_logger(__name__)


_TC = {  # Tree Characters
    '|': u'\u2502',  # Horizontal
    '-': u'\u2500',  # Vertical
    'L': u'\u2514',  # L-Shape
    'S': u'\u251C',  # Split
}

"""
======================== Additional Dictionary Classes and Functions ========================
"""


class DotDict(dict):
    """ Make dict fields accessible by . """
    def __init__(self, *args, **kwargs):
        super(DotDict, self).__init__(*args, **kwargs)
        for key in self:
            if isinstance(self[key], dict):
                self[key] = DotDict(self[key])

    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def print_dict_tree(dictionary, name='Dictionary'):
    """ Prints a dictionary as a tree """
    def print_tree(tree, level_char):
        for i, key in enumerate(sorted(tree.keys())):
            if i == len(tree) - 1:
                node_char = _TC['L'] + _TC['-']
                level_char_pp = level_char + '   '
            else:
                node_char = _TC['S'] + _TC['-']
                level_char_pp = level_char + _TC['|'] + '  '

            if isinstance(tree[key], dict):
                LOG.info(u"{:s}{:s} {:s}"
                         .format(level_char, node_char, str(key)))
                print_tree(tree[key], level_char_pp)
            else:
                LOG.info(u"{:s}{:s} {:s}: {:s}"
                         .format(level_char, node_char, str(key), str(tree[key])))

    LOG.info('{:s}:'.format(name))
    print_tree(dictionary, '')


"""
======================== Dict Parser ========================
"""


class ArgumentError(Exception):
    pass


class OptionError(Exception):
    pass


class Argument(object):
    """ Helper Class for DictParser """
    def __init__(self, name, **kwargs):
        self.name = name
        self.required = kwargs.pop('required', False)
        self.default = kwargs.pop('default', None)
        self.help = kwargs.pop('help', '')
        self.type = kwargs.pop('type', None)
        self.subtype = kwargs.pop('subtype', None)
        self.choices = kwargs.pop('choices', None)

        if len(kwargs) > 0:
            ArgumentError("'{:s}' are not valid parameters for Argument.".format(kwargs.keys()))

        self._validate()

    def _validate(self):
        if not isinstance(self.name, basestring):
            raise ArgumentError("Argument '{:s}': ".format(str(self.name)) +
                                "Name is not a valid string.")

        if self.type and not isinstance(self.default, self.type):
            raise ArgumentError("Argument '{:s}': ".format(self.name) +
                                "Default value not of specified type.")

        if self.choices and not isinstance(self.choices, list):
                raise ArgumentError("Argument '{:s}': ".format(self.name) +
                                    "Choices need to be a list.")

        if self.choices and self.default not in self.choices:
            raise ArgumentError("Argument '{:s}': ".format(self.name) +
                                "Default value not found in choices.")

        if self.choices and self.type:
            for choice in self.choices:
                if not isinstance(choice, self.type):
                    raise ArgumentError("Choice '{:s}'".format(choice) +
                                        "of Argument '{:s}': ".format(self.name) +
                                        "is not of type '{:s}'.".format(self.type))

        if self.subtype and not (self.type or self.type == list):
            raise ArgumentError("Argument '{:s}': ".format(self.name) +
                                "Parameter 'subtype' is only accepted if 'type' is list.")

        if self.required and self.default is not None:
            LOG.warn("Argument '{:s}': ".format(self.name) +
                     "Value is required but default value is given. The latter will be ignored.")


class DictParser(object):
    """ Provides functions to parse a dictionary.

    First build a dictionary structure with Arguments as leafs via add_argument or on init.
    A similar structured option dictionary with the values as leafs can then be parsed.
    """

    def __init__(self, dictionary=None, strict=True):
        """
        Initialize Class either empty or with preconfigured dictionary

        Args:
            dictionary: Preconfigured Dictionary for parsing
            strict: Strict Parsers don't accept unknown options. If False, it just prints warnings.
        """
        self.strict = strict

        if dictionary:
            self._validate_arguments(dictionary)
            self.dictionary = dictionary
        else:
            self.dictionary = {}

    #########################
    # Static Methods (private)
    #########################

    @staticmethod
    def _validate_arguments(dictionary):
        """ Validates an input dictionary that can be used as arguments.

        Args:
            dictionary: Dictionary to validate
        """
        for key in dictionary:
            arg = dictionary[key]
            if isinstance(arg, dict):
                try:
                    DictParser._validate_arguments(arg)
                except ArgumentError as e:
                    e.message = "'{:s}.{:s}".format(key, e.message[1:])
                    e.args = (e.message,)
                    raise
            elif not isinstance(arg, Argument):
                raise ArgumentError("'{:s}' is not a valid entry.".format(key))
            else:
                if key != arg.name:
                    raise ArgumentError("'{:s}': Key and name need to be the same.".format(key))

    @staticmethod
    def _check_value(key, opt_dict, arg_dict):
        """ Checks if in opt_dict[key] satisfies arg_dict[key]

        Args:
            key: key to check
            opt_dict: Options-structure. Can be None or empty.
            arg_dict: Argument-structure. Needs to contain 'key'

        Returns:
            The appropriate value for opt_dict[key]
        """
        arg = arg_dict[key]
        if not opt_dict or key not in opt_dict:
            if arg.required:
                raise OptionError("'{:s}' required in options.\n{:s}".format(
                    key, arg.help)
                )
            else:
                return arg.default

        opt = opt_dict[key]
        if arg.type and not isinstance(opt, arg.type):
            raise OptionError("'{:s}' is not of type {:s}.\n{:s}".format(
                key, arg.type.__name__, arg.help)
            )
        if arg.choices and opt not in arg.choices:
            raise OptionError("'{:s}' needs to be one of {:s}.\n{:s}".format(
                key, arg.choices, arg.help)
            )
        return opt

    @staticmethod
    def _parse_options(opt_dict, arg_dict, strict):
        """ Use parse_options()!

        This is a helper Function for parsing options. It does all the work. Called recursively.
        Needs to be static, because of the recursion! Sorry about that.

        Args:
            opt_dict: Dictionary with the options
            arg_dict: Dictionary with the arguments to check the options against
            strict: True: raises error on unknown options
                    False: Just adds them to default dict (prints warning)

        Returns:
            Dictionary with parsed options
        """
        checked_dict = {}
        for key in arg_dict:
            if isinstance(arg_dict[key], Argument):
                checked_dict[key] = DictParser._check_value(key, opt_dict, arg_dict)
                if checked_dict[key] is None and key not in opt_dict:
                    checked_dict.pop(key)  # CAREFUL HERE
            elif isinstance(arg_dict[key], dict):
                try:
                    if not opt_dict or not (key in opt_dict):
                        checked_dict[key] = DictParser._parse_options(None,
                                                                      arg_dict[key], strict)
                    else:
                        checked_dict[key] = DictParser._parse_options(opt_dict[key],
                                                                      arg_dict[key], strict)
                except OptionError as e:
                    e.message = "'{:s}.{:s}".format(key, e.message[1:])
                    e.args = (e.message,)
                    raise

            opt_dict.pop(key, None)  # Default value avoids KeyError

        if len(opt_dict) > 0 and strict:
            raise OptionError("Unknown Options: '{:s}'.".format(opt_dict.keys()))

        for key in opt_dict:
            LOG.warn("Unknown Option: '{:s}'.".format(key))
            checked_dict[key] = opt_dict[key]
        return checked_dict

    #########################
    # Public Methods
    #########################

    def parse_options(self, options):
        """ Parse a given option dictionary and return parsed options.

        Args:
            options: Options to parse

        Return:
            Parsed options
        """
        return self._parse_options(options.copy(), self.dictionary, self.strict)

    def parse_config_items(self, items):
        """ Parse a list of (name, value) items, where the values are all strings.

        Args:
            items: list of (name, value) items.

        Returns:
            Parsed options
        """
        options = self._convert_config_items(items)
        return self._parse_options(options, self.dictionary, self.strict)

    def add_argument(self, arg, **kwargs):
        """ Adds an argument to the parser.

        If you want it to be an argument of a sub-dictionary add
        the 'loc=subdict.subdict' keyword to the input.

        Args:
            arg: Argument to add (either of object of class argument or string defining the name)
            kwargs: Any of the argument-fields (apart from 'name') and/or 'loc'

        Returns:
            This object
        """
        loc = kwargs.pop('loc', None)
        if not isinstance(arg, Argument):
            arg = Argument(arg, **kwargs)
        self._add_arg_to_dict(arg, loc)
        return self

    def add_argument_dict(self, dictionary, loc):
        """ Appends a complete subdictionary to existing argument structure at node 'loc'.

        Args:
            loc: locination of the node to append the sub-dictionary
            dictionary: The dictionary to append

        Returns:
            This object
        """
        fields = loc.split('.')
        name = fields[-1]
        sub_dict = self._traverse_dict('.'.join(fields[:-1]))

        if name in sub_dict:
            raise ArgumentError("'{:s}' already exists in parser!".format(name))

        self._validate_arguments(dictionary)
        sub_dict[name] = dictionary
        return self

    def help(self):
        # TODO: Print Help-Message
        pass

    def tree(self):
        """ Prints the current Argument Tree (I made dis :) ) """
        def print_tree(tree, level_char):
            for i, key in enumerate(sorted(tree.keys())):
                if i == len(tree) - 1:
                    node_char = _TC['L'] + _TC['-']
                    level_char_pp = level_char + '   '
                else:
                    node_char = _TC['S'] + _TC['-']
                    level_char_pp = level_char + _TC['|'] + '  '

                LOG.info(u"{:s}{:s} {:s}".format(level_char, node_char, key))
                if isinstance(tree[key], dict):

                    print_tree(tree[key], level_char_pp)
                else:
                    leaf = tree[key]
                    LOG.info(u"{:s}{:s} {:s}: {:s}".format(
                             level_char_pp, _TC['S'] + _TC['-'],
                             'Required', str(leaf.required)))
                    LOG.info(u"{:s}{:s} {:s}: {:s}".format(
                             level_char_pp, _TC['S'] + _TC['-'],
                             'Default', str(leaf.default)))
                    LOG.info(u"{:s}{:s} {:s}: {:s}".format(
                             level_char_pp, _TC['S'] + _TC['-'],
                             'Type', leaf.type.__name__ if leaf.type else 'None'))
                    LOG.info(u"{:s}{:s} {:s}: {:s}".format(
                             level_char_pp, _TC['S'] + _TC['-'],
                             'Choices', str(leaf.choices)))
                    LOG.info(u"{:s}{:s} {:s}: {:s}".format(
                             level_char_pp, _TC['L'] + _TC['-'],
                             'Help', leaf.help))

        LOG.info('Argument Dictionary')
        print_tree(self.dictionary, '')

    #########################
    # Private Methods
    #########################

    def _add_arg_to_dict(self, arg, loc=None):
        """ Adds and argument to the argument dictionary.

        These will be used to parse an incoming option structure.

        Args:
            arg: Argument to add
            loc: Path to sub-dictionary as string (e.g. subdict.subdict.loc[.arg])

        Returns:
            This object
        """
        sub_dict = self._traverse_dict(loc)
        if arg.name in sub_dict:
            raise ArgumentError("'{:s}' already exists in parser!".format(arg.name))
        sub_dict[arg.name] = arg
        return self

    def _traverse_dict(self, loc=None):
        """ Traverses the dictionary to the subdict defined by loc.

        Adds non-existing substructures automatically.

        Args:
            loc: Path to sub-dictionary as string (e.g. argument.subargument.locination)

        Returns:
            Sub-dictionary
        """
        d = self.dictionary
        if loc:
            traverse = loc.split('.')
            for i, t in enumerate(traverse):
                try:
                    d = d[t]
                except KeyError:
                    d[t] = {}
                    d = d[t]
                if isinstance(d, Argument):
                    raise ArgumentError(
                        "'{:s}' is already an argument and hence cannot be a subdict.".format(
                            '.'.join(traverse[:i] + [t])))
        return d

    def _convert_config_items(self, items):
        """ Converts items list to a dictionary with types already in place """
        def evaluate(item):
            try:
                return eval(value)  # sorry for using that
            except NameError:
                return value  # might be a string in the first place.

        out = {}
        for name, value in items:
            if name in self.dictionary:
                arg = self.dictionary[name]
                if arg.type == list:
                    if not value.startswith("["):
                        value = "[" + value + "]"
                    value = evaluate(value)
                    if arg.subtype:
                        for idx, entry in enumerate(value):
                            value[idx] = arg.subtype(value)
                out[name] = value
            else:
                # could check self.strict here, but result is passed to get checked anyway
                out[name] = evaluate(value)
        return out



"""
======================== Script Testing Mode ========================
"""

if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
