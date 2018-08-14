import sys
import os
from os.path import abspath, dirname
import argparse
from optics_measurements.optics_input import _get_optics_parser, OpticsInput
from model import manager


def parse_args(args=None):
    main_parser = _get_main_parser()
    main_options, rest = main_parser.parse_known_args(args)
    if "optics" in rest:
        rest.append('--files={}'.format(main_options.file))
        rest.append('--outputdir={}'.format(main_options.outputdir))
        rest.append('--model_dir={}'.format(dirname(abspath(main_options.model))))
    subparsers = {
        "clean": _get_clean_parser,
        "harpy": _get_harpy_parser,
        "optics": _get_optics_parser,
    }
    suboptions = {}
    for subparser_name in subparsers.keys():
        if subparser_name in rest:
            rest.remove(subparser_name)
            suboption, rest = subparsers[subparser_name]().parse_known_args(rest)
            suboptions[subparser_name] = suboption
    main_input = MainInput.init_from_options(main_options)
    clean_input = None
    harpy_input = None
    optics_input = None

    if "clean" in suboptions:
        clean_input = CleanInput.init_from_options(suboptions["clean"])

    if "harpy" in suboptions:
        harpy_input = HarpyInput.init_from_options(suboptions["harpy"])

    if "optics" in suboptions:
        accelerator = manager.get_accel_instance(rest)
        suboptions["optics"].accelerator = accelerator
        optics_input = OpticsInput.init_from_options(suboptions["optics"])



    python_path = "/afs/cern.ch/work/o/omc/anaconda/bin/python "
    return main_input, clean_input, harpy_input, optics_input, (python_path + " ".join(sys.argv))


class ArgumentError(Exception):
    pass


class MainInput(object):
    DEFAULTS = {
        "write_raw": False,
        "startturn": 0,
        "endturn": 50000,
        "skip_files": False,
    }

    def __init__(self):
        self.file = None
        self.model = None
        self.outputdir = None
        self.write_raw = MainInput.DEFAULTS["write_raw"]
        self.startturn = MainInput.DEFAULTS["startturn"]
        self.endturn = MainInput.DEFAULTS["endturn"]
        self.skip_files = MainInput.DEFAULTS["skip_files"]

    @staticmethod
    def init_from_options(options):
        self = MainInput()
        self.file = options.file#[file_name.strip() for file_name in options.file.strip("\"").split(",")]
        self.model = options.model
        self.outputdir = options.outputdir
        self.write_raw = options.write_raw
        self.skip_files = options.skip_files
        if self.outputdir is None:
            outdir = os.path.dirname(self.file)
            self.outputdir = outdir
        self.startturn = options.startturn
        self.endturn = options.endturn
        return self


def _get_main_parser():
    parser = argparse.ArgumentParser()

    # Obligatory arguments ###########
    parser.add_argument(
        "--file", "--files",
        help="Comma separated binary files to clean",
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
    parser.add_argument(
        "--skip_files",
        help=("If present, it will not write the .amps and the .freqs files" +
              "ASCII file in the --outputdir."),
        action="store_true",    
    )
    parser.add_argument(
        "--startturn",
        help="Turn index to start. Default is first turn: %(default)s",
        default=MainInput.DEFAULTS["startturn"],
        dest="startturn", type=int
    )
    parser.add_argument(
        "--endturn",
        help="""First turn index to be ignored.
                Default is a number that is lower than the maximum which
                can be handled by drive: %(default)s""",
        default=MainInput.DEFAULTS["endturn"],
        dest="endturn", type=int
    )
    return parser
    ################################


class CleanInput(object):
    DEFAULTS = {
        "svd_mode": "numpy",
        "sing_val": 12,
        "peak_to_peak": 0.00001,
        "max_peak": 20.,
        "single_svd_bpm_threshold": 0.925,
        "bad_bpms": "",
        "wrong_polarity_bpms": "",
        "noresync": False,
        "write_clean": False,
        "no_exact_zeros": False,
    }

    def __init__(self):
        self.svd_mode = CleanInput.DEFAULTS["svd_mode"]
        self.sing_val = CleanInput.DEFAULTS["sing_val"]
        self.peak_to_peak = CleanInput.DEFAULTS["peak_to_peak"]
        self.max_peak = CleanInput.DEFAULTS["max_peak"]
        self.single_svd_bpm_threshold = CleanInput.DEFAULTS[
            "single_svd_bpm_threshold"
        ]
        self.bad_bpms = CleanInput.DEFAULTS["bad_bpms"]
        self.wrong_polarity_bpms = CleanInput.DEFAULTS["wrong_polarity_bpms"]
        self.noresync = CleanInput.DEFAULTS["noresync"]
        self.write_clean = CleanInput.DEFAULTS["write_clean"]
        self.no_exact_zeros = CleanInput.DEFAULTS["no_exact_zeros"]

    @staticmethod
    def init_from_options(options):
        self = CleanInput()
        self.svd_mode = options.svd_mode
        self.sing_val = options.sing_val
        self.peak_to_peak = options.peak_to_peak
        self.max_peak = options.max_peak
        self.single_svd_bpm_threshold = options.single_svd_bpm_threshold
        if options.bad_bpms == "":
            self.bad_bpms = []
        else:
            self.bad_bpms = [bad_bpm.strip() for bad_bpm in options.bad_bpms.strip("\"").split(",")]
        if options.wrong_polarity_bpms == "":
            self.wrong_polarity_bpms = []
        else:
            self.wrong_polarity_bpms = [wrong_polarity_bpm.strip() for wrong_polarity_bpm in options.wrong_polarity_bpms.strip("\"").split(",")]
        self.noresync = options.noresync
        self.write_clean = options.write_clean
        self.no_exact_zeros = options.no_exact_zeros
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
        "--bad_bpms",
        help=""""Comma separated bad BPMs to clean. Default: %(default)s""",
        default=CleanInput.DEFAULTS["bad_bpms"],
        dest="bad_bpms", type=str)
    parser.add_argument(
        "--wrong_polarity_bpms",
        help=""""Comma separated BPMs with swapped polarity in both planes. Default: %(default)s""",
        default=CleanInput.DEFAULTS["wrong_polarity_bpms"],
        dest="wrong_polarity_bpms", type=str)
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
    parser.add_argument(
        "--no_exact_zeros",
        help="""If present, will not remove files with a single zero .""",
        dest="no_exact_zeros", action="store_true"
    )
    ################################
    return parser


class HarpyInput(object):
    DEFAULTS = {
        "tunez": 0.0,
        "tolerance": 0.01,
        "harpy_mode": "svd",
        "sequential": False,
        "no_tune_clean": False,
        "tune_clean_limit": 1e-5,
        "is_free_kick": False,
    }

    def __init__(self):
        self.tunex = None
        self.tuney = None
        self.tunez = HarpyInput.DEFAULTS["tunez"]
        self.nattunex = None
        self.nattuney = None
        self.nattunez = None
        self.tolerance = HarpyInput.DEFAULTS["tolerance"]
        self.harpy_mode = HarpyInput.DEFAULTS["harpy_mode"]
        self.sequential = HarpyInput.DEFAULTS["sequential"]
        self.no_tune_clean = HarpyInput.DEFAULTS["no_tune_clean"]
        self.tune_clean_limit = HarpyInput.DEFAULTS["tune_clean_limit"]
        self.is_free_kick = HarpyInput.DEFAULTS["is_free_kick"]

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
        self.no_tune_clean = options.no_tune_clean
        self.tune_clean_limit = options.tune_clean_limit
        self.is_free_kick = options.is_free_kick
        return self


def _get_harpy_parser():
    parser = argparse.ArgumentParser()

    # Obligatory arguments ###########
    parser.add_argument(
        "--tunex", help="Guess for the main horizontal tune.",
        dest="tunex", type=float,
        required=True,
    )
    parser.add_argument(
        "--tuney", help="Guess for the main vertical tune.",
        dest="tuney", type=float,
        required=True,
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
        default=None,
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
        choices=("bpm", "svd", "fast"),
        default=HarpyInput.DEFAULTS["harpy_mode"],
    )
    parser.add_argument(
        "--sequential", help="If set, it will run in only one process.",
        dest="sequential", action="store_true",
    )
    parser.add_argument(
        "--no_tune_clean",
        help="If present, the automatic tune cleaning will be deactivated.",
        dest="no_tune_clean", action="store_true",
    )
    parser.add_argument(
        "--tune_clean_limit",
        help="The autoclean wont remove tune deviation lower than this limit.",
        dest="tune_clean_limit", type=float,
        default=HarpyInput.DEFAULTS["tune_clean_limit"]
    )
    parser.add_argument(
        "--free_kick",
        help="If present, it will perform the free kick phase correction",
        dest="is_free_kick", action="store_true",
    )
    ################################
    return parser
