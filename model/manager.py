import argparse
from model.accelerators import lhc, ps


ACCELS = {
    lhc.Lhc.NAME: lhc.Lhc,
    ps.Cps.NAME: ps.Cps
}


def get_accel_class(name, *args, **kwargs):
    accel = _try_to_get_class(name, ACCELS)
    accel_cls = accel.get_class(*args, **kwargs)
    return accel_cls


def get_accel_class_from_args(args):
    name, args = _parse_accel_name(args)
    accel = _try_to_get_class(name, ACCELS)
    accel_cls, rest_args = accel.get_class_from_args(args)
    return accel_cls, rest_args


def _try_to_get_class(name, cls_dict):
    try:
        return cls_dict[name]
    except KeyError:
        raise ValueError(
            "name should be one of: " +
            str(ACCELS.keys())
        )


def _parse_accel_name(args):
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
