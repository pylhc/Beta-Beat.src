"""
Entrypoint Extract Online Model
---------------------

Entrypoint for the online model extractor wrapper.
"""
from online_model import extractor_wrapper
from utils import logging_tools, iotools
from utils.entrypoint import entrypoint, EntryPointParameters, ArgumentError

LOG = logging_tools.get_logger(__name__)

FUNCTION_OVERVIEW = 'overview'
FUNCTION_DEFINITION = 'definition'


def get_params():
    return EntryPointParameters({
        "function": dict(
            flags=["-f", "--functionality"],
            type=str,
            required=True,
            choices=["overview", "definition"],
            help="Which functionality to run."
        ),
        "knob_names": dict(
            flags=["-k", "--knobs", "--knobnames"],
            type=str,
            nargs="+",
            default=[],
            help="Names of the knobs to show."
        ),
        "time": dict(
            flags=["-t", "--time"],
            type=str,
        ),
        "cwd": dict(
            flags=["-c", "--cwd"],
            type=str,
            default="./",
            help="Path of the current working directory.",
        ),
        "server": dict(
            flags=["-s", "--server"],
            type=str,
            help="Server to use."
        ),
        "show_plot": dict(
            flags=["--showplots"],
            action="store_true",
            help="Whether to show plots or not."
            "(Only for {:s} functionality.)".format(FUNCTION_OVERVIEW)
        )
    })


@entrypoint(get_params(), strict=True)
def main(opt):
    """ Entrypoint for the online model extractor python wrapper.

    Creates either the overview or knob definitions, depending on the functionality chosen.

    Keyword Args:
        Required
        function (str): Which functionality to run.
                        **Flags**: ['-f', '--functionality']
                        **Choices**: ['overview', 'definition']

        Optional
        cwd (str): Path of the current working directory.
                   **Flags**: ['-c', '--cwd']
                   **Default**: ``./``
        knob_names (str): Names of the knobs to show.
                          **Flags**: ['-k', '--knobs', '--knobnames']
                          **Default**: __See Source__
        server (str): Server to use.
                      **Flags**: ['-s', '--server']
        show_plot: Whether to show plots or not.(Only for overview functionality.)
                   **Flags**: ['--showplots']
                   **Action**: ``store_true``
        time (str): -Help not available-
                    **Flags**: ['-t', '--time']

    """
    iotools.create_dirs(opt.cwd)
    if opt.function == FUNCTION_OVERVIEW:
        extractor_wrapper.extract_overview(opt.knob_names, opt.time, opt.cwd,
                                           server=opt.server, show_plot=opt.show_plot)
    elif opt.function == FUNCTION_DEFINITION:
        if not opt.knob_names:
            raise ArgumentError(
                "Argument 'knob_names' required for function '{:s}'.".format(FUNCTION_DEFINITION)
            )
        if opt.time is None:
            raise ArgumentError(
                "Argument 'time' required for function '{:s}'.".format(FUNCTION_DEFINITION)
            )
        if opt.show_plot:
            LOG.warn(
                "Argument 'show_plot' has no effect in function '{:s}'".format(FUNCTION_DEFINITION)
            )
        extractor_wrapper.extract_knob_value_and_definition(opt.knob_names, opt.time, opt.cwd,
                                                            server=opt.server)


def _create_help():
    """ Print help for parameters. """
    import sys
    from utils.entrypoint import create_parameter_help
    create_parameter_help(sys.modules[__name__])


if __name__ == '__main__':
    # _create_help()
    main()
