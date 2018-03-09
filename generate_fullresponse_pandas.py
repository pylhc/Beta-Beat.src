import cPickle as pickle
import os

import madx_wrapper as madx
from correction.fullresponse import response_madx
from correction.fullresponse.response_analytical import TwissResponse
from global_correct_iterative import DEFAULT_ARGS
from model import manager
from twiss_optics.sequence_parser import EXT as VARMAP_EXT
from utils import logging_tools
from utils.contexts import timeit
from utils.entrypoint import EntryPointParameters, entrypoint

LOG = logging_tools.get_logger(__name__)


def get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--creator",
        help="Create either with madx or analytically from twiss file.",
        name="creator",
        type=str,
        choices=("madx", "twiss"),
        default="madx",
    )
    params.add_parameter(
        flags="--variables",
        help="List of the variables classes to use.",
        name="variable_categories",
        default=DEFAULT_ARGS["variables"],
        nargs="+",
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
    params.add_parameter(
        flags=["-k", "--deltak"],
        help="Delta K1L to be applied to quads for sensitivity matrix (madx-only).",
        default=0.00002,
        name="delta_k",
        type=float
    )
    params.add_parameter(
        flags="--debug",
        help="Print debug information.",
        name="debug",
        action="store_true",
    )
    return params


@entrypoint(get_params())
def create_response(opt, other_opt):
    """ Entry point for creating pandas-based response matrices.

    The response matrices can be either created by response_madx or TwissResponse.

    Keyword Args:
        Required
        model_dir (str): Path to the model directory.
                         **Flags**: ['-m', '--model_dir']
        outfile_path (str): Name of fullresponse file.
                            **Flags**: ['-o', '--outfile']
        Optional
        creator (str): Create either with madx or analytically from twiss file.
                       **Flags**: --creator
                       **Choices**: ('madx', 'twiss')
                       **Default**: ``madx``
        debug: Print debug information.
               **Flags**: --debug
               **Action**: ``store_true``
        delta_k (float): Delta K1L to be applied to quads for sensitivity matrix (madx-only).
                         **Flags**: ['-k', '--deltak']
                         **Default**: ``2e-05``
        variable_categories: List of the variables classes to use.
                             **Flags**: --variables
                             **Default**: ``['MQM', 'MQT', 'MQTL', 'MQY']``
    """
    with logging_tools.DebugMode(active=opt.debug):
        LOG.info("Creating response.")
        accel_cls, other_opt = manager.get_accel_class_and_unkown(other_opt)
        variables = accel_cls.get_variables(classes=opt.variable_categories)

        if opt.creator == "madx":
            jobfile_path = os.path.join(opt.model_dir, "job.iterpandas.madx")

            # generate the response (saved into file)
            response_madx.generate_fullresponse(opt.outfile_path, variables, jobfile_path,
                                                delta_k=opt.delta_k)

        elif opt.creator == "twiss":
            accel_inst = accel_cls(model_dir=opt.model_dir)

            # handle the varmap/sequence file
            varmapfile_name = accel_inst.NAME.lower() + "b" + str(accel_inst.get_beam())
            varmap_path = os.path.join(opt.model_dir, varmapfile_name + "." + VARMAP_EXT)
            if not os.path.isfile(varmap_path):
                with logging_tools.TempFile("save_sequence_madxout.tmp", LOG.debug) as log_file:
                    madx.resolve_and_run_file(
                        os.path.join(opt.model_dir, "job.save_sequence.madx"),
                        log_file=log_file,
                    )
                varmap_path = varmap_path.replace("." + VARMAP_EXT, ".seq")

            LOG.debug("Creating Fullresponse via TwissResponse")
            with timeit(lambda t: LOG.debug("  Total time getting TwissResponse: {:f}s".format(t))):
                model = accel_inst.get_elements_tfs()
                tr = TwissResponse(varmap_path, model, variables)
                fullresponse = tr.get_fullresponse()

            LOG.debug("Saving Response into file '{:s}'".format(opt.outfile_path))
            with open(opt.outfile_path, 'wb') as dump_file:
                pickle.Pickler(dump_file, -1).dump(fullresponse)


if __name__ == "__main__":
    create_response()
