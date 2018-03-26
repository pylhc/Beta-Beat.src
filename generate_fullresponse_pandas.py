"""
Provides a response generation wrapper.
The response matrices can be either created by response_madx or analytically via TwissResponse.

:author: Joschua Dillys
"""
import cPickle as pickle
import os

from correction.fullresponse import response_madx
from correction.fullresponse.response_twiss import TwissResponse
from global_correct_iterative import DEFAULT_ARGS, check_varmap_file
from model import manager
from utils import logging_tools
from utils.contexts import timeit
from utils.entrypoint import EntryPointParameters, entrypoint

LOG = logging_tools.get_logger(__name__)

DEFAULT_PATTERNS = {
    "job_content": "%JOB_CONTENT%",  # used in lhc_model_creator, sequence_evaluation
    "twiss_columns": "%TWISS_COLUMNS%",  # used in lhc_model_creator, sequence_evaluation
    "element_pattern": "%ELEMENT_PATTERN%",  # used in lhc_model_creator, sequence_evaluation
}

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
        flags="--optics_params",
        help="List of parameters to correct upon (e.g. BBX BBY; twiss-only).",
        name="optics_params",
        type=str,
        nargs="+",
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
        optics_params (str): List of parameters to correct upon (e.g. BBX BBY; twiss-only).
                             **Flags**: --optics_params
        variable_categories: List of the variables classes to use.
                             **Flags**: --variables
                             **Default**: ``['MQM', 'MQT', 'MQTL', 'MQY']``
    """
    with logging_tools.DebugMode(active=opt.debug,
                                 log_file=os.path.join(opt.model_dir, "generate_fullresponse.log")):
        LOG.info("Creating response.")
        accel_cls, other_opt = manager.get_accel_class_and_unkown(other_opt)
        accel_inst = accel_cls(model_dir=opt.model_dir)
        variables = accel_inst.get_variables(classes=opt.variable_categories)

        if opt.creator == "madx":
            jobfile_path = os.path.join(opt.model_dir, "tmpl.generate_fullresponse.madx")
            patterns = {
                "job_content": "%JOB_CONTENT%",
                "twiss_columns": "%TWISS_COLUMNS%",
                "element_pattern": "%ELEMENT_PATTERN%",
            }
            madx_script = accel_inst.get_basic_twiss_job(patterns["job_content"],
                                                         patterns["twiss_columns"],
                                                         patterns["element_pattern"])
            with open(jobfile_path, "w") as f:
                f.write(madx_script)

            fullresponse = response_madx.generate_fullresponse(variables, jobfile_path,
                                                               patterns=patterns,
                                                               delta_k=opt.delta_k)
        elif opt.creator == "twiss":

            varmap_path = check_varmap_file(accel_inst, variables)

            LOG.debug("Creating response via TwissResponse.")
            with timeit(lambda t:
                        LOG.debug("Total time getting TwissResponse: {:f}s".format(t))):
                tr = TwissResponse(varmap_path, accel_inst.get_elements_tfs(), variables)
                fullresponse = tr.get_response_for(opt.optics_params)

        LOG.debug("Saving Response into file '{:s}'".format(opt.outfile_path))
        with open(opt.outfile_path, 'wb') as dump_file:
            pickle.Pickler(dump_file, -1).dump(fullresponse)


if __name__ == "__main__":
    create_response()
