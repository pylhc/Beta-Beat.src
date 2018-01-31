from __future__ import print_function
import __init__
import os
from Utilities.entrypoint import EntryPoint
from Utilities.entrypoint import EntryPointArguments
from Utilities.dict_tools import print_dict_tree
from Utilities import logging_tools
LOG = logging_tools.get_logger(__name__, level_console=0, fmt="%(message)s")
THISDIR = os.path.dirname(os.path.abspath(__file__))

"""
======================== Example Parameter Definitions ========================
"""


def _get_params_arguments():
    """ Parameters defined with EntryPointArguments (which is a dict *cough*) """
    args = EntryPointArguments()
    args.add_argument(name="accel",
                      flags=["-a", "--accel"],
                      help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
                      choices=["LHCB1","LHCB2","LHCB5"],
                      default="LHCB1")
    args.add_argument(name="dict",
                      flags=["-d", "--dictionary"],
                      help="File with the BPM dictionary",
                      default="/test.notafile",
                      type=str)
    args.add_argument(name="anumber",
                      flags=["-num", "--anum"],
                      help="Just a number.",
                      type=float,
                      default=19.,
                      )
    args.add_argument(name="anint",
                      flags=["-i", "--int"],
                      help="Just a number.",
                      type=int,
                      required=True,
                      )
    args.add_argument(name="alist",
                      flags=["-l", "--lint"],
                      help="Just a number.",
                      type=int,
                      nargs="+",
                      required=True,
                      )
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
        "alist": dict(flags=["-l", "--lint"],
                      help="Just a number.",
                      type=int,
                      nargs="+",
                      required=True,
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
             flags=["-num", "--anum"],
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
    ]



@EntryPoint(_get_params_arguments())
def some_function(options):
    LOG.info("Some Function")
    print_dict_tree(options)
    LOG.info("\n\n")


@EntryPoint(_get_params_dict())
def some_other_function(options):
    LOG.info("Some Other Function")
    print_dict_tree(options)
    LOG.info("\n\n")


@EntryPoint(_get_params_list())
def some_function_list_param(options):
    LOG.info("Some Function with list params")
    print_dict_tree(options)
    LOG.info("\n\n")


def some_function_test():
    kw_dict = dict(accel="LHCB5", anumber=5.6, anint=10, alist=[1,2,3])

    LOG.info("KW-Arguments")
    some_function(**kw_dict)
    some_other_function(**kw_dict)
    some_function_list_param(**kw_dict)

    LOG.info("KW-Arguments, entry_dict")
    some_function(entry_dict=kw_dict)

    LOG.info("Positional argument, dict")
    some_function(kw_dict)

    LOG.info("Positional argument, configfile")
    some_function(os.path.join(THISDIR, "entrypoint_config_test.cfg"))


if __name__ == "__main__":
    some_function_test()

