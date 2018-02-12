import __init__
import manager  # noqa
from Utilities.entrypoint import EntryPointParameters, entrypoint
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


def _get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags=["--type"],
        name="type",
        help="Type of model to create, either nominal or best_knowledge",
        choices=("nominal", "best_knowledge", "coupling_correction"),
    )
    params.add_parameter(
        flags=["--output"],
        help="Output path for model, twiss files will be writen here.",
        name="output",
        required=True,
        type=str,
    )
    params.add_parameter(
        flags=["--writeto"],
        help="Path to the file where to write the resulting MAD-X script. ",
        name="writeto",
        type=str,
    )
    params.add_parameter(
        flags=["--logfile"],
        help=("Path to the file where to write the MAD-X script output."
              "If not provided it will be written to sys.stdout."),
        name="logfile",
        type=str,
    )
    return params


# Main functions ###############################################################


@entrypoint(_get_params())
def create_instance_and_model(opt, accel_opt):
    accel_inst = manager.get_accel_instance(accel_opt)
    create_model(
        accel_inst,
        opt.type,
        opt.output,
        writeto=opt.writeto,
        logfile=opt.logfile,
    )


def create_model(accel_inst, model_type, output_path, **kwargs):
    CREATORS[accel_inst.NAME][model_type].create_model(
        accel_inst,
        output_path,
        **kwargs
    )


# Script Mode ##################################################################


if __name__ == "__main__":
    create_instance_and_model()
