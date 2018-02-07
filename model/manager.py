from Utilities.entrypoint import entrypoint, EntryPoint, EntryPointParameters, split_arguments
from model.accelerators import lhc, esrf, psbooster

ACCELS = {
    lhc.Lhc.NAME: lhc.Lhc,
    esrf.Esrf.NAME: esrf.Esrf,
    psbooster.Psbooster.NAME: psbooster.Psbooster,
}


def _get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags=["--accel"],
        help=("Choose the accelerator to use."
              "Can be the class already, which is then returned."
              ),
        name="accel",
        required=True,
        choices=ACCELS.keys(),
    )
    return params


@entrypoint(_get_params())
def get_accel_class(opt, cls_opt):
    """Returns accelerator class

    Keyword Args:
        accel: Choose the accelerator to use. Can be the class already, which is then returned.
    """
    if not isinstance(opt.accel, str):
        # assume it's the class
        return opt.accel

    accel = _get_parent_class(opt.accel)
    accel_cls = accel.get_class(cls_opt)
    return accel_cls


def get_accel_class_from_args(args=None):
    """ LEGACY-FUNCTION SHOULD BE REPLACED BY USING get_accel_class """
    parser = EntryPoint(_get_params())
    name, args = parser.parse(args)

    accel = _get_parent_class(name)

    accel_args, rest_args = split_arguments(args, accel.get_class_params())
    accel_cls, rest_args = accel.get_class(accel_args)
    return accel_cls, rest_args


def _get_parent_class(name):
    try:
        return ACCELS[name]
    except KeyError:
        raise ValueError(
            "name should be one of: " +
            str(ACCELS.keys())
        )


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
