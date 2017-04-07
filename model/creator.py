import sys
import os
import argparse

CURRENT_DIR = os.path.dirname(__file__)
sys.path.append(os.path.abspath(os.path.join(CURRENT_DIR, "..")))

import manager  # noqa
from model_creators import model_creator  # noqa
from model_creators.lhc_model_creator import (  # noqa
    LhcModelCreator,
    LhcBestKnowledgeCreator,
)


CREATORS = {
    "lhc": {"nominal": LhcModelCreator,
            "best_knowledge": LhcBestKnowledgeCreator},
}


def create_model(accel_name, accel_inst, model_type, output_path):
    CREATORS[accel_name][model_type].create_model(accel_inst, output_path)


def _i_am_main():
    rest_args = sys.argv[1:]

    accel_name, rest_args = manager.parse_accel_name(rest_args)
    accel_cls, rest_args = manager.get_accel_class_from_args(
        accel_name,
        rest_args
    )
    accel_inst, rest_args = accel_cls.init_from_args(rest_args)
    model_type, output_path = _parse_rest_args(rest_args)
    create_model(accel_name, accel_inst, model_type, output_path)


def _parse_rest_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "type",
        help="Type of model to create, either nominal or best_knowledge",
        choices=("nominal", "best_knowledge"),
    )
    parser.add_argument(
        "--output",
        help="Output path for model, twiss files will be writen here.",
        dest="output",
        required=True,
        type=str,
    )
    options = parser.parse_args(args)
    return options.type, options.output


if __name__ == "__main__":
    _i_am_main()
