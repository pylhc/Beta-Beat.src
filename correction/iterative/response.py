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
        flags="--creator",
        name="creator",
        type=str,
        choices=("madx", "twiss"),
        default="madx",
    )
    params.add_parameter(
        flags="--variables",
        help="Comma separated names of the variables classes to use.",
        name="variable_categories",
        default=DEFAULTS["variables"],
    )
    params.add_parameter(
        flags=["-m", "--model_dir"],
        help="Path to the model directory.",
        name="model_dir",
        required=True,
        type=str
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
    """ Entry point for creating pandas-based response matrices.

    The response matrices can be either created by response_madx or TwissResponse.
    """
    LOG.info("Creating response.")
    accel_cls, other_opt = manager.get_accel_class_and_unkown(other_opt)
    variables = accel_cls.get_variables(classes=opt.variable_categrories)

    if opt.creator == "madx":
        # add some more arguments to the given ones
        jobfile_path = os.path.join(opt.model_dir, "job.iterpandas.madx")
        try:
            other_opt += ["--variables", str(variables)]
            other_opt += ["--jobfile", jobfile_path]
            other_opt += ["--outfile", opt.outfile_path]
        except TypeError:
            other_opt.variables = variables
            other_opt.original_jobfile_path = jobfile_path
            other_opt.outfile_path = opt.outfile_path

        # generate the response (saved into file)
        response_madx.generate_fullresponse(other_opt)

    elif opt.creator == "twiss":
        accel_inst = accel_cls(model_dir=opt.model_dir)

        # handle the varmap/sequence file
        varmapfile_name = accel_inst.name.lower() + "b" + str(accel_inst.get_beam())
        varmap_path = os.path.join(opt.model_dir, varmapfile_name + "." + VARMAP_EXT)
        if not os.path.isfile(varmap_path):
            with logging_tools.TempFile("save_sequence_madxout.tmp", LOG.debug) as log_file:
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
