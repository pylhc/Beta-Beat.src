import cPickle as pickle
import os

from correction.iterative import response_madx
from correction.iterative.correct_iterative import DEFAULTS
from madx import madx_wrapper as madx
from model import manager
from twiss_optics.response_class import TwissResponse
from twiss_optics.sequence_parser import EXT as VARMAP_EXT
from utils import logging_tools
from utils.contexts import timeit
from utils.entrypoint import EntryPointParameters, entrypoint

LOG = logging_tools.get_logger(__name__)


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
    params.add_parameter(
        flags=["-o", "--outfile"],
        help="Name of fullresponse file.",
        name="outfile_path",
        required=True,
        type=str
    )


@entrypoint(get_params())
def create_response(opt, other_opt):
    LOG.info("Creating response.")
    accel_cls, other_opt = manager.get_accel_class_and_unkown(other_opt)
    variables = accel_cls.get_variables(classes=opt.variable_categrories)

    if opt.mode == "madx":
        try:
            other_opt += ["--variables", str(variables)]
            other_opt += ["--outfile", opt.outfile_path]
        except TypeError:
            other_opt.variables = variables
            other_opt.outfile_path = opt.outfile_path
        response_madx.generate_fullresponse(other_opt)

    elif opt.mode == "analytical":
        if not opt.model_dir:
            raise ValueError("Please provide a path to the model directory.")

        accel_inst = accel_cls(model_dir=opt.model_dir)

        file_name = accel_inst.name.lower() + "b" + str(accel_inst.get_beam())
        varmap_path = os.path.join(opt.model_dir, file_name + "." + VARMAP_EXT)
        if not os.path.isfile(varmap_path):
            with logging_tools.TempFile("correct_iter_madxout.tmp", LOG.debug) as log_file:
                madx.resolve_and_run_string(
                    os.path.join(opt.model_dir, "job.save_sequence.madx"),
                    log_file=log_file.path,
                )
            varmap_path = varmap_path.replace("." + VARMAP_EXT, ".seq")

        LOG.debug("Creating Fullresponse via TwissResponse")
        with timeit(lambda t: LOG.debug("  Total time getting TwissResponse: {:f}s".format(t))):
            tr = TwissResponse(varmap_path, accel_inst.model, variables)
            fullresponse = tr.get_fullresponse()

        LOG.debug("Saving Response into file '{:s}'".format(opt.outfile_path))
        with open(opt.outfile_path, 'wb') as dump_file:
            pickle.Pickler(dump_file, -1).dump(fullresponse)


if __name__ == "__main__":
    create_response()
