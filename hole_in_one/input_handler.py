import argparse


def parse_args():
    main_parser = _get_main_parser()
    main_options, rest = main_parser.parse_known_args()

    subparsers = {
        "clean": _get_clean_parser,
        "harpy": _get_harpy_parser,
    }
    suboptions = {}
    while len(rest) > 0:
        next_option = rest.pop(0)
        try:
            subparse_funct = subparsers[next_option]
        except KeyError:
            raise ArgumentError(
                "Too many arguments: " + str([next_option] + rest)
            )
        subparser = subparse_funct()
        del subparsers[next_option]
        suboption, rest = subparser.parse_known_args(rest)
        suboptions[next_option] = suboption

    main_input = MainInput.init_from_options(main_options)
    clean_input = None
    harpy_input = None
    try:
        clean_input = CleanInput.init_from_options(suboptions["clean"])
    except KeyError:
        pass
    try:
        harpy_input = HarpyInput.init_from_options(suboptions["harpy"])
    except KeyError:
        pass
    return main_input, clean_input, harpy_input


class ArgumentError(Exception):
    pass


class MainInput(object):
    DEFAULTS = {
        "write_raw": False,
    }

    def __init__(self):
        self.file = None
        self.model = None
        self.outputdir = None
        self.write_raw = MainInput.DEFAULTS["write_raw"]

    @staticmethod
    def init_from_options(options):
        self = MainInput()
        self.file = options.file
        self.model = options.model
        self.outputdir = options.outputdir
        self.write_raw = options.write_raw
        return self


def _get_main_parser():
    parser = argparse.ArgumentParser()

    # Obligatory arguments ###########
    parser.add_argument(
        "--file",
        help="Binary File to clean",
        required=True,
        dest="file"
    )
    parser.add_argument(
        "--model",
        help="Model for BPM locations",
        required=True,
        dest="model"
    )
    ################################

    # Optional arguments ###########
    parser.add_argument(
        "--outputdir",
        help="""Output directory.
                If not present it will use the input file directory.""",
        dest="outputdir"
    )
    parser.add_argument(
        "--write_raw",
        help=("If present, it will write the raw" +
              "ASCII file in the --outputdir."),
        action="store_true",
    )
    return parser
    ################################


class CleanInput(object):
    DEFAULTS = {
        "svd_mode": "numpy",
        "startturn": 0,
        "endturn": 50000,
        "sing_val": 12,
        "peak_to_peak": 0.00001,
        "max_peak": 20.,
        "single_svd_bpm_threshold": 0.925,
        "noresync": False,
        "write_clean": False,
    }

    def __init__(self):
        self.svd_mode = CleanInput.DEFAULTS["svd_mode"]
        self.startturn = CleanInput.DEFAULTS["startturn"]
        self.endturn = CleanInput.DEFAULTS["endturn"]
        self.sing_val = CleanInput.DEFAULTS["sing_val"]
        self.peak_to_peak = CleanInput.DEFAULTS["peak_to_peak"]
        self.max_peak = CleanInput.DEFAULTS["max_peak"]
        self.single_svd_bpm_threshold = CleanInput.DEFAULTS[
            "single_svd_bpm_threshold"
        ]
        self.noresync = CleanInput.DEFAULTS["noresync"]
        self.write_clean = CleanInput.DEFAULTS["write_clean"]

    @staticmethod
    def init_from_options(options):
        self = CleanInput()
        self.svd_mode = options.svd_mode
        self.startturn = options.startturn
        self.endturn = options.endturn
        self.sing_val = options.sing_val
        self.peak_to_peak = options.peak_to_peak
        self.max_peak = options.max_peak
        self.single_svd_bpm_threshold = options.single_svd_bpm_threshold
        self.noresync = options.noresync
        self.write_clean = options.write_clean
        return self


def _get_clean_parser():
    parser = argparse.ArgumentParser()

    # Optional arguments ###########
    parser.add_argument(
        "--svd_mode",
        help="""Defines modes of the SVD funtion itself. Should be one of:
                - numpy: numpy.linalg.svd
                - sparse: scipy.sparse.linalg.svds
                - random: Randomized SVD
        """,
        default=CleanInput.DEFAULTS["svd_mode"],
        dest="svd_mode", type=str,
        choices=("numpy", "sparse", "random"),
    )
    parser.add_argument(
        "--startturn",
        help="Turn index to start. Default is first turn: %(default)s",
        default=CleanInput.DEFAULTS["startturn"],
        dest="startturn", type=int
    )
    parser.add_argument(
        "--endturn",
        help="""First turn index to be ignored.
                Default is a number that is lower than the maximum which
                can be handled by drive: %(default)s""",
        default=CleanInput.DEFAULTS["endturn"],
        dest="endturn", type=int
    )
    parser.add_argument(
        "--sing_val",
        help="""Keep this amount of singular values in decreasing order,
                rest will be cut (set to 0).
                Default is a large number: %(default)s""",
        default=CleanInput.DEFAULTS["sing_val"],
        dest="sing_val", type=int
    )
    parser.add_argument(
        "--pk-2-pk",
        help="""Peak to peak amplitude cut. This removes BPMs where
                abs(max(turn values) - min(turn values)) <= threshold.
                Default: %(default)s""",
        default=CleanInput.DEFAULTS["peak_to_peak"],
        dest="peak_to_peak", type=float)
    parser.add_argument(
        "--max-peak-cut",
        help="""Maximum peak tolerance in mm. This removes BPMs where the
                maximum measured oscillation > threshold.
             Default: %(default)s""",
        default=CleanInput.DEFAULTS["max_peak"],
        dest="max_peak", type=float)
    parser.add_argument(
        "--single_svd_bpm_threshold",
        help="""Threshold for single BPM dominating a mode.
                Should be > 0.9 for LHC. Default: %(default)s""",
        default=CleanInput.DEFAULTS["single_svd_bpm_threshold"],
        dest="single_svd_bpm_threshold", type=float)
    parser.add_argument(
        "--noresync",
        help="""Ignore the synchronization error of the BPMs.""",
        dest="noresync", action="store_true"
    )
    parser.add_argument(
        "--write_clean",
        help="""If present, it will write the cleaned ASCII file
                in the --outputdir.""",
        dest="write_clean", action="store_true"
    )
    ################################
    return parser


class HarpyInput(object):
    DEFAULTS = {
        "tunez": 0.0,
        "nattunez": 0.0,
        "tolerance": 0.01,
        "harpy_mode": "bpm",
        "sequential": False,
    }

    def __init__(self):
        self.tunex = None
        self.tuney = None
        self.nattunex = None
        self.nattuney = None
        self.tunez = HarpyInput.DEFAULTS["tunez"]
        self.nattunez = HarpyInput.DEFAULTS["nattunez"]
        self.tolerance = HarpyInput.DEFAULTS["tolerance"]
        self.harpy_mode = HarpyInput.DEFAULTS["harpy_mode"]
        self.sequential = HarpyInput.DEFAULTS["sequential"]

    @staticmethod
    def init_from_options(options):
        self = HarpyInput()
        self.tunex = options.tunex
        self.tuney = options.tuney
        self.nattunex = options.nattunex
        self.nattuney = options.nattuney
        self.tunez = options.tunez
        self.nattunez = options.nattunez
        self.tolerance = options.tolerance
        self.harpy_mode = options.harpy_mode
        self.sequential = options.sequential
        return self


def _get_harpy_parser():
    parser = argparse.ArgumentParser()

    # Obligatory arguments ###########
    parser.add_argument(
        "--tunex", help="Guess for the main horizontal tune.",
        dest="tunex", type=float,
    )
    parser.add_argument(
        "--tuney", help="Guess for the main vertical tune.",
        dest="tuney", type=float,
    )
    ################################

    # Optional arguments ###########
    parser.add_argument(
        "--nattunex", help="Guess for the natural horizontal tune.",
        dest="nattunex", type=float,
    )
    parser.add_argument(
        "--nattuney", help="Guess for the natural vertical tune.",
        dest="nattuney", type=float,
    )
    parser.add_argument(
        "--tunez", help="Guess for the main synchrotron tune.",
        default=HarpyInput.DEFAULTS["tunez"],
        dest="tunez", type=float,
    )
    parser.add_argument(
        "--nattunez", help="Guess for the natural synchrotron tune.",
        default=HarpyInput.DEFAULTS["nattunez"],
        dest="nattunez", type=float,
    )
    parser.add_argument(
        "--tolerance", help="Tolerance on the guess for the tunes.",
        default=HarpyInput.DEFAULTS["tolerance"],
        dest="tolerance", type=float,
    )
    parser.add_argument(
        "--harpy_mode", help="Harpy resonance computation mode.",
        dest="harpy_mode", type=str,
        choices=("bpm", "svd"),
        default=HarpyInput.DEFAULTS["harpy_mode"],
    )
    parser.add_argument(
        "--sequential", help="If set, it will run in only one process.",
        dest="sequential", action="store_true",
    )
    ################################
    return parser
