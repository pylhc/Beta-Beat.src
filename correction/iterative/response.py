import os

from correction.iterative import response_madx
from correction.iterative.correct_iterative import DEFAULTS
from model import manager
from twiss_optics.response_class import TwissResponse
from utils.entrypoint import EntryPointParameters, entrypoint


def get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--mode",
        name="mode",
        type=str,
        choices=("madx", "analytical"),
        default="madx",
    )
    params.add_parameter(
        flags="--variables",
        help="Comma separated names of the variables classes to use.",
        name="variable_categories",
        default=DEFAULTS["variables"],
    )
    params.add_parameter(
        flags=["--model_dir"],
        help="Path to the model directory.",
        name="model_dir",
    )


@entrypoint(get_params())
def main(opt, other_opt):
    accel_cls, other_opt = manager.get_accel_class_and_unkown(other_opt)
    variables = accel_cls.get_variables(classes=opt.variable_categrories)

    if opt.mode == "madx":
        try:
            other_opt += ["--variables", str(variables)]
        except TypeError:
            other_opt.variables = variables

        response_madx.generate_fullresponse(other_opt)
    elif opt.mode == "analytical":
        if not opt.model_dir:
            raise ValueError("Please provide a path to the model directory.")

        accel_inst = accel_cls(model_dir=opt.model_dir)


        if not os.path.isfile(os.path.join(opt.model_dir,
                accel_inst.name.lower() + "b" + str(accel_inst.get_beam()))):

        tr = TwissResponse(accel_cls.get)


if __name__ == "__main__":
    pass