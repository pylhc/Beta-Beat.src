import argparse
from accelerators import lhc


ACCELS = {
    "lhc": lhc.Lhc,
}


def get_accel_class(name, *args, **kwargs):
    try:
        accel = ACCELS[name]
    except KeyError:
        raise ValueError(
            "name should be one of: " +
            str(ACCELS.keys())
        )
    accel_cls, rest_args = accel.get_class(*args, **kwargs)
    return accel_cls


def get_accel_class_from_args(name, args):
    try:
        accel = ACCELS[name]
    except KeyError:
        raise ValueError(
            "name should be one of: " +
            str(ACCELS.keys())
        )
    accel_cls, rest_args = accel.get_class_from_args(args)
    return accel_cls, rest_args


def parse_accel_name(args):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--accel",
        help=(
            "Choose the accelerator to use"
        ),
        dest="accel",
        required=True,
        choices=ACCELS.keys(),
    )
    options, rest_args = parser.parse_known_args(args)
    return options.accel, rest_args
