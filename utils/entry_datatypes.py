""" Advanced Datatypes to add as type to entrypoint. Or any parser, really. """
from __future__ import print_function

import abc
import sys
import six

# for testing purposes:
from os.path import abspath, join, dirname, pardir
sys.path.append(abspath(join(dirname(__file__), pardir)))

from utils.entrypoint import entrypoint, EntryPointParameters


# Meta Class Helper ############################################################


def get_instance_faker_meta(*classes):
    """ Returns the metaclass that fakes the isinstance() checks. """
    class FakeMeta(abc.ABCMeta):
        def __instancecheck__(cls, inst):
            return isinstance(inst, classes)
    return FakeMeta


# 'Merge' Classes ##############################################################


def get_multi_class(*classes):
    """ Create a class 'behaving' like all classes in `classes`.

    In case a value needs to be converted to a class in this list,
    it is attempted to cast the input to the classes in the given order
    (i.e. string-classes need to go to the end, as they 'always' succeed).
    """
    class MultiClass(object):
        __metaclass__ = get_instance_faker_meta(*classes)

        @classmethod
        def _convert_to_a_type(cls, value):
            for c in classes:
                try:
                    return c.__new__(c, value)
                except ValueError:
                    pass
            else:
                raise ValueError(
                    "The value '{}' cant be converted to any of the classes '{}' ".format(
                        value,
                        ",".join([c.__name__ for c in classes]),
                    )
                )

        def __new__(cls, value):
            if isinstance(value, six.string_types) or not isinstance(value, classes):
                return cls._convert_to_a_type(value)
            return value

    return MultiClass


# More Fake Classes ############################################################


class DictAsString(object):
    """ Use dicts in command line like {"key":value} """
    __metaclass__ = get_instance_faker_meta(str, dict)

    def __new__(cls, s):
        if isinstance(s, dict):
            return s

        d = eval(s)
        if not isinstance(d, dict):
            raise ValueError("'{}' can't be converted to a dictionary.".format(s))
        return d


class BoolOrString(object):
    """ A class that behaves like a boolean when possible, otherwise like a string."""
    __metaclass__ = get_instance_faker_meta(bool, str)

    def __new__(cls, value):
        if value in ["True", "1", True, 1]:
            return bool.__new__(bool, True)

        elif value in ["False", "0", False, 0]:
            return bool.__new__(bool, False)

        else:
            return str.__new__(str, value)


# Test #########################################################################


def get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--dict",
        name="dict",
        type=DictAsString,
        help="Use a dictionary!",
    )
    params.add_parameter(
        flags="--int_or_str",
        name="int_or_str",
        type=get_multi_class(int, str),
        help="either int or string",
    )
    params.add_parameter(
        flags="--bool_or_str",
        name="bool_or_str",
        type=BoolOrString,
        help="either bool or string",
    )
    return params


@entrypoint(get_params(), strict=False)
def tester(opt, other):
    print("Opt:")
    print(opt)
    print("\n")
    print("Other:")
    print(other)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        tester(dict={"test": 5}, bool_or_str="test", int_or_str="test")
        tester(dict={"test": 5}, bool_or_str=False, int_or_str=10)
    else:
        tester()
        # e.g.
        # --dict '{"hello":4, "there":{"General":"Kenobi"}}' --int_or_str hello --bool_or_str there
        # --dict '{"hello":4, "there":{"General":"Kenobi"}}' --int_or_str 5 --bool_or_str 0
        # --dict '{"hello":4, "there":{"General":"Kenobi"}}' --int_or_str 0 --bool_or_str True
