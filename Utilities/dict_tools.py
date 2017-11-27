"""

"""

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
                print(u"{:s}{:s} {:s}"
                      .format(level_char, node_char, str(key)))
                print_tree(tree[key], level_char_pp)
            else:
                print(u"{:s}{:s} {:s}: {:s}"
                      .format(level_char, node_char, str(key), str(tree[key])))

    print('{:s}:'.format(name))
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
        self.required = kwargs.get('required', False)
        self.default = kwargs.get('default', None)
        self.help = kwargs.get('help', '')
        self.type = kwargs.get('type', None)
        self.choices = kwargs.get('choices', None)

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

        if self.required and self.default is not None:
            Warning("Argument '{:s}': ".format(self.name) +
                    "Value is required but default value is given. The latter will be ignored.")


class DictParser(object):
    """
    Provides functions to parse a dictionary.
    First build a dictionary structure with Arguments as leafs via add_argument or on init.
    A similar structured option dictionary with the values as leafs can then be parsed.
    """

    def __init__(self, dictionary=None):
        """
        Initialize Class either empty or with preconfigured dictionary
        :param dictionary: Preconfigured Dictionary for parsing
        """
        if dictionary:
            DictParser._validate_arguments(dictionary)
            self.dictionary = dictionary
        else:
            self.dictionary = {}

    #########################
    # Static Methods (private)
    #########################

    @staticmethod
    def _validate_arguments(dictionary):
        """
        Validates an input dictionary that can be used as arguments.
        :param dictionary: Dictionary to validate
        :return: Nothing
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
        """
        Checks if in opt_dict[key] satisfies arg_dict[key]
        :param key: key to check
        :param opt_dict: Options-structure. Can be None or empty.
        :param arg_dict: Argument-structure. Needs to contain 'key'
        :return: The approriate value for opt_dict[key]
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
    def _parse_options(opt_dict, arg_dict):
        """
        Helper Function for parsing options. Does all the work. Called recursively.
        :param opt_dict: Dictionary with the options
        :param arg_dict: Dictionary with the arguments to check the options against
        :return: Dictionary with parsed options
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
                        checked_dict[key] = DictParser._parse_options(None, arg_dict[key])
                    else:
                        checked_dict[key] = DictParser._parse_options(opt_dict[key], arg_dict[key])
                except OptionError as e:
                    e.message = "'{:s}.{:s}".format(key, e.message[1:])
                    e.args = (e.message,)
                    raise
            try:
                opt_dict.pop(key)
            except KeyError:
                pass

        for key in opt_dict:
            print("WARNING: '{:s}' is not a valid option.".format(key))
            checked_dict[key] = opt_dict[key]
        return checked_dict

    #########################
    # Public Methods
    #########################

    def parse_options(self, options):
        """
        Parse a given option dictionary. Returns parsed options.
        :param options: Options to parse
        :return: Parsed options
        """
        return DictParser._parse_options(options.copy(), self.dictionary)

    def add_argument(self, arg, **kwargs):
        """
        Adds an argument to the parser. If you want it to be an argument of a sub-dictionary add
        the 'dest=subdict.subdict' keyword to the input.
        :param arg: Argument to add (either of object of class argument or string defining the name)
        :param kwargs: Any of the argument-fields (apart from 'name') and/or 'dest'
        :return: This object
        """
        dest = kwargs.pop('dest', None)
        if not isinstance(arg, Argument):
            arg = Argument(arg, **kwargs)
        self._add_arg_to_dict(arg, dest)
        return self

    def add_argument_dict(self, dictionary, dest):
        """
        Appends a complete subdictionary to existing argument structure at node 'dest'.
        :param dest: Destination of the node to append the sub-dictionary
        :param dictionary: The dictionary to append
        :return: This object
        """
        fields = dest.split('.')
        name = fields[-1]
        sub_dict = self._traverse_dict('.'.join(fields[:-1]))

        if name in sub_dict:
            raise ArgumentError("'{:s}' already exists in parser!".format(name))

        DictParser._validate_arguments(dictionary)
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

                print(u"{:s}{:s} {:s}".format(level_char, node_char, key))
                if isinstance(tree[key], dict):

                    print_tree(tree[key], level_char_pp)
                else:
                    leaf = tree[key]
                    print(u"{:s}{:s} {:s}: {:s}".format(
                        level_char_pp, _TC['S'] + _TC['-'],
                        'Required', str(leaf.required)))
                    print(u"{:s}{:s} {:s}: {:s}".format(
                        level_char_pp, _TC['S'] + _TC['-'],
                        'Default', str(leaf.default)))
                    print(u"{:s}{:s} {:s}: {:s}".format(
                        level_char_pp, _TC['S'] + _TC['-'],
                        'Type', leaf.type.__name__ if leaf.type else 'None'))
                    print(u"{:s}{:s} {:s}: {:s}".format(
                        level_char_pp, _TC['S'] + _TC['-'],
                        'Choices', str(leaf.choices)))
                    print(u"{:s}{:s} {:s}: {:s}".format(
                        level_char_pp, _TC['L'] + _TC['-'],
                        'Help', leaf.help))

        print('Argument Dictionary')
        print_tree(self.dictionary, '')

    #########################
    # Private Methods
    #########################

    def _add_arg_to_dict(self, arg, dest=None):
        """
        Adds and argument to the argument dictionary.
        These will be used to parse an incoming option structure.
        :param arg: Argument to add
        :param dest: Path to sub-dictionary as string (e.g. subdict.subdict.dest[.arg])
        :return: This object
        """
        sub_dict = self._traverse_dict(dest)
        if arg.name in sub_dict:
            raise ArgumentError("'{:s}' already exists in parser!".format(arg.name))
        sub_dict[arg.name] = arg
        return self

    def _traverse_dict(self, dest=None):
        """
        Traverses the dictionary to the subdict defined by dest.
        Adds non-existing substructures automatically.
        :param dest: Path to sub-dictionary as string (e.g. argument.subargument.destination)
        :return: sub-dictionary
        """
        d = self.dictionary
        if dest:
            traverse = dest.split('.')
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


"""
======================== Script Testing Mode ========================
"""

if __name__ == '__main__':
    parser = DictParser({
        'test': Argument(
            name='test',
            default='test',
            required=False,
        ),
        'test1': Argument(
            name='test1'
        ),
        'sub': {
            'test': Argument(
                'test',
                type=int,
                default=2
            )
        },
        'sub2': {
            'subsub': {
                'buub': Argument(
                    'buub',
                    type=int,
                    default=2
                )
            }

        }
    }
    )

    # parser.tree()
    opt = {
        'sub': {
            'test': 2
        },
        'sub2': {
            'subsub': {
                'buub': 2
            }
        }
    }

    opt = parser.parse_options(opt)
    print_dict_tree(opt, "Options")
