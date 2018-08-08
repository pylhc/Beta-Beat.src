""" Tester for EntryPoint

Tests some very basic stuff, including errors, so some ERROR-LOGs are expected.
"""
from __future__ import print_function

import os

# noinspection PyUnresolvedReferences
import sys
import os
from os.path import abspath, join, dirname
new_path = abspath(join(dirname(abspath(__file__)), os.pardir, os.pardir))
if new_path not in sys.path:
    sys.path.append(new_path)
from utils import logging_tools
from utils.dict_tools import print_dict_tree
from utils.entrypoint import ArgumentError
from utils.entrypoint import EntryPoint
from utils.entrypoint import EntryPointParameters
from utils.entrypoint import entrypoint

from model import manager
manager._get_params()

LOG = logging_tools.get_logger(__name__, level_console=0, fmt="%(levelname)7s | %(message)s")
THISDIR = os.path.dirname(os.path.abspath(__file__))


# Example Parameter Definitions ################################################


def _get_params_arguments():
    """ Parameters defined with EntryPointArguments (which is a dict *cough*) """
    args = EntryPointParameters()
    args.add_parameter(name="accel",
                       flags=["-a", "--accel"],
                       help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
                       choices=["LHCB1","LHCB2","LHCB5"],
                       default="LHCB1")
    args.add_parameter(name="dict",
                       flags=["-d", "--dictionary"],
                       help="File with the BPM dictionary",
                       default="/test.notafile",
                       type=str)
    args.add_parameter(name="anumber",
                       flags=["-num", "--anum"],
                       help="Just a number.",
                       type=float,
                       default=19.,
                       )
    args.add_parameter(name="anint",
                       flags=["-i", "--int"],
                       help="Just a number.",
                       type=int,
                       required=True,
                       )
    args.add_parameter(name="alist",
                       flags=["-l", "--lint"],
                       help="Just a number.",
                       type=int,
                       nargs="+",
                       required=True,
                       )
    args.add_parameter(name="anotherlist",
                       flags=["-k", "--alint"],
                       help="list.",
                       type=str,
                       nargs=3,
                       ),
    return args


def _get_params_dict():
    """ Parameters defined as a dictionary of dictionaries directly """
    args = {
        "accel": dict(
            flags=["-a", "--accel"],
            help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
            choices=["LHCB1", "LHCB2", "LHCB5"],
            default="LHCB1"),
        "dict": dict(
            flags=["-d", "--dictionary"],
            help="File with the BPM dictionary",
            default="/test.notafile",
            type=str),
        "anumber": dict(
            flags=["-num", "--anum"],
            help="Just a number.",
            type=float,
            default=19.,
        ),
        "anint": dict(
            flags=["-i", "--int"],
            help="Just a number.",
            type=int,
            required=True,
        ),
        "alist": dict(
            flags=["-l", "--lint"],
            help="Just a number.",
            type=int,
            nargs="+",
            required=True,
        ),
        "anotherlist": dict(
            flags=["-k", "--alint"],
            help="list.",
            type=str,
            nargs=3,
        ),
    }
    return args


def _get_params_list():
    """ Parameters defined as list """
    return [
        dict(name="accel",
             flags=["-a", "--accel"],
             help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
             choices=["LHCB1","LHCB2","LHCB5"],
             default="LHCB1"),
        dict(name="dict",
             flags=["-d", "--dictionary"],
             help="File with the BPM dictionary",
             default="/test.notafile",
             type=str),
        dict(name="anumber",
             flags=["-n", "--anum"],
             help="Just a number.",
             type=float,
             default=19.,
             ),
        dict(name="anint",
             flags=["-i", "--int"],
             help="Just a number.",
             type=int,
             required=True,
             ),
        dict(name="alist",
             flags=["-l", "--lint"],
             help="Just a number.",
             type=int,
             nargs="+",
             required=True,
             ),
        dict(name="anotherlist",
             flags=["-k", "--alint"],
             help="list.",
             type=str,
             nargs=3,
             ),
    ]


# Example Wrapped Functions ####################################################


@entrypoint(_get_params_arguments())
def some_function(options, unknown_options):
    LOG.info("Some Function")
    print_dict_tree(options)
    LOG.info("Unknown Options: \n {:s}".format(unknown_options))
    LOG.info("\n")


@entrypoint(_get_params_arguments(), strict=True)
def strict_function(options):
    LOG.info("Strict Function")
    print_dict_tree(options)
    LOG.info("\n")


@entrypoint(_get_params_dict())
def some_other_function(options, unknown_options):
    LOG.info("Some Other Function")
    print_dict_tree(options)
    LOG.info("Unknown Options: \n {:s}".format(unknown_options))
    LOG.info("\n")


@entrypoint(_get_params_list())
def some_function_list_param(options, unknown_options):
    LOG.info("Some Function with list params")
    print_dict_tree(options)
    LOG.info("Unknown Options: \n {:s}".format(unknown_options))
    LOG.info("\n")


class TestClass(object):
    @entrypoint(_get_params_dict())
    def instance_function(self, options, unknowns):
        LOG.info("Class Function with dict params")
        print_dict_tree(options)
        LOG.info("Unknown Options: \n {:s}".format(unknowns))
        LOG.info("\n")

    @classmethod
    @entrypoint(_get_params_dict())
    def class_function(cls, options, unknowns):
        LOG.info("Class Function with dict params")
        print_dict_tree(options)
        LOG.info("Unknown Options: \n {:s}".format(unknowns))
        LOG.info("\n")


    @classmethod
    @entrypoint(_get_params_dict(), strict=True)
    def class_function_strict(cls, options):
        LOG.info("Class Function with dict params")
        print_dict_tree(options)
        LOG.info("\n")


# Tests ########################################################################


def some_function_test():
    arg_list = "-a LHCB5 -n 5.6 -i 10 -l 1 2 3 -k hubba dubba subba"
    arg_list_unknowns = arg_list + " --xx was -j ist das"
    kw_dict = dict(accel="LHCB5", anumber=5.6, anint=10,
                   alist=[1, 2, 3], anotherlist=["hubba", "dubba", "subba"])
    kw_w_unknowns = kw_dict.copy()
    kw_w_unknowns.update(un="what", known="is that?")

    LOG.info("# KW-Args ########################")
    LOG.info("KW-Arguments")
    some_function(**kw_dict)
    some_other_function(**kw_dict)
    some_function_list_param(**kw_dict)
    LOG.info("\n\n")

    LOG.info("KW-Arguments, with unknowns")
    some_function(**kw_w_unknowns)
    some_other_function(**kw_w_unknowns)
    some_function_list_param(**kw_w_unknowns)
    LOG.info("\n\n")

    LOG.info("# KW-Args entry dict ########################")
    LOG.info("KW-Arguments, entry_dict")
    some_function(entry_dict=kw_dict)
    LOG.info("KW-Arguments, entry_dict w/ unknowns")
    some_function(entry_dict=kw_w_unknowns)
    LOG.info("\n\n")

    LOG.info("# Positional Dict ########################")
    LOG.info("Positional argument, dict")
    some_function(kw_dict)
    LOG.info("Positional argument, dict w/ unknowns")
    some_function(kw_w_unknowns)
    LOG.info("\n\n")

    LOG.info("# Positional CFG ########################")
    LOG.info("Positional argument, configfile")
    some_function(os.path.join(THISDIR, "entrypoint_config_test.cfg"))
    LOG.info("\n\n")

    LOG.info("# Positional List ########################")
    LOG.info("Positional argument, list")
    some_function(arg_list.split(" "))
    LOG.info("Positional argument, list w/ unknowns")
    some_function(arg_list_unknowns.split(" "))
    LOG.info("\n\n")

    LOG.info("# Strict Function ########################")
    LOG.info("KW-Arguments")
    strict_function(**kw_dict)
    LOG.info("KW-Arguments, entry_dict")
    strict_function(entry_dict=kw_dict)
    LOG.info("Positional argument, dict")
    strict_function(kw_dict)
    LOG.info("Positional argument, list")
    strict_function(arg_list.split(" "))
    LOG.info("\n\n")

    LOG.info("# Strict Functions Unknowns ########################")
    LOG.info("KW-Arguments")
    try:
        strict_function(**kw_w_unknowns)
    except ArgumentError as e:
        LOG.error(e.message)
        LOG.error("\n")
    LOG.info("KW-Arguments, entry_dict")
    try:
        strict_function(entry_dict=kw_w_unknowns)
    except ArgumentError as e:
        LOG.error(e.message)
        LOG.error("\n")
    LOG.info("Positional argument, dict")
    try:
        strict_function(kw_w_unknowns)
    except ArgumentError as e:
        LOG.error(e.message)
        LOG.error("\n")
    LOG.info("Positional argument, list")
    try:
        strict_function(arg_list_unknowns.split(" "))
    except ArgumentError as e:
        LOG.error(e.message)
        LOG.error("\n")
    LOG.info("\n\n")


def class_function_test():
    kw_dict = dict(accel="LHCB5", anumber=5.6, anint=10,
                   alist=[1, 2, 3], anotherlist=["hubba", "dubba", "subba"])
    kw_w_unknowns = kw_dict.copy()
    kw_w_unknowns.update(un="what", known="is that?")

    LOG.info("# Class Functions ########################")

    LOG.info("# KW-Args ########################")
    LOG.info("KW-Arguments")
    c = TestClass()
    c.instance_function(**kw_w_unknowns)
    TestClass.class_function(**kw_w_unknowns)
    TestClass.class_function_strict(**kw_dict)
    LOG.info("\n\n")


def as_parser_test():
    kw_dict = dict(accel="LHCB5", anumber=5.6, anint=10, alist=[1, 2, 3])
    p = EntryPoint(_get_params_list())
    opt, uopt = p.parse(**kw_dict)
    print_dict_tree(opt)


# Script Mode ##################################################################


if __name__ == "__main__":
    as_parser_test()
    some_function_test()
    class_function_test()


