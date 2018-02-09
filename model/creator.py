import __init__
import sys
import os
import argparse
from Utilities.entrypoint import split_arguments

import manager  # noqa
from model_creators import model_creator  # noqa
from model_creators.lhc_model_creator import (  # noqa
    LhcModelCreator,
    LhcBestKnowledgeCreator,
    LhcSegmentCreator,
    LhcCouplingCreator,
)
from model_creators.psbooster_model_creator import PsboosterModelCreator


CREATORS = {
    "lhc": {"nominal": LhcModelCreator,
            "best_knowledge": LhcBestKnowledgeCreator,
            "segment": LhcSegmentCreator,
            "coupling_correction": LhcCouplingCreator},
    "psbooster": {"nominal": PsboosterModelCreator},
}


def create_model(accel_inst, model_type, output_path, **kwargs):
    CREATORS[accel_inst.NAME][model_type].create_model(
        accel_inst,
        output_path,
        **kwargs
    )


def _i_am_main():
    args = sys.argv[1:]
    options, accel_args = _parse_creator_args(args)
    accel_inst = manager.get_accel_instance(accel_args)
    create_model(
        accel_inst,
        options.type,
        options.output,
        writeto=options.writeto,
        logfile=options.logfile,
    )


def _parse_creator_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "type",
        help="Type of model to create, either nominal or best_knowledge",
        choices=("nominal", "best_knowledge", "coupling_correction"),
    )
    parser.add_argument(
        "--output",
        help="Output path for model, twiss files will be writen here.",
        dest="output",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--writeto",
        help="Path to the file where to write the resulting MAD-X script. ",
        dest="writeto",
        type=str,
    )
    parser.add_argument(
        "--logfile",
        help=("Path to the file where to write the MAD-X script output."
              "If not provided it will be written to sys.stdout."),
        dest="logfile",
        type=str,
    )
    options, unknown_args = parser.parse_known_args(args)
    return options, unknown_args


if __name__ == "__main__":
    _i_am_main()
