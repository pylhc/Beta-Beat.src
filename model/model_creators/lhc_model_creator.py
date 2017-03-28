from __future__ import print_function
import os
import sys
import logging
from accelerators.lhc import LhcExcitationMode
from accelerators import lhc
import argparse
import model_creator

AFS_ROOT = "/afs"
if "win" in sys.platform and sys.platform != "darwin":
    AFS_ROOT = "\\AFS"

LOGGER = logging.getLogger(__name__)


class LhcModelCreator(model_creator.ModelCreator):

    LHC_MODES = {
        "lhc_runI": lhc.LhcRunI,
        "lhc_runII": lhc.LhcRunII2015,
        "lhc_runII_2016": lhc.LhcRunII2016,
        "lhc_runII_2016_ats": lhc.LhcRunII2016Ats,
        "hllhc": lhc.HlLhc,
    }
    ERR_DEF_PATH = os.path.join(AFS_ROOT, "cern.ch", "work", "o", "omc",
                                "Error_definition_files")
    ERR_DEF_FILES = {
        "4.5": "0450GeV", "1.0": "1000GeV",
        "1.5": "1500GeV", "2.0": "2000GeV",
        "2.5": "2500GeV", "3.0": "3000GeV",
        "3.5": "3500GeV", "4.0": "4000GeV",
        "4.5": "4500GeV", "5.0": "5000GeV",
        "5.5": "5500GeV", "6.0": "6000GeV",
        "6.5": "6500GeV",
    }

    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        madx_template = None
        with open(lhc_instance.get_nominal_template()) as textfile:
            madx_template = textfile.read()
        iqx, iqy = cls._get_full_tunes(lhc_instance)
        use_acd = "1" if (lhc_instance.excitation ==
                          LhcExcitationMode.ACD) else "0"
        use_adt = "1" if (lhc_instance.excitation ==
                          LhcExcitationMode.ADT) else "0"
        replace_dict = {
            "RUN": lhc_instance.MACROS_NAME,
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.beam,
            "PATH": output_path,
            "DPP": lhc_instance.dpp,
            "QMX": iqx,
            "QMY": iqy,
            "USE_ACD": use_acd,
            "USE_ADT": use_adt,
            "QX": "", "QY": "", "QDX": "", "QDY": "",
        }
        if (lhc_instance.excitation in
                (LhcExcitationMode.ACD, LhcExcitationMode.ADT)):
            replace_dict["QX"] = lhc_instance.nat_tune_x
            replace_dict["QY"] = lhc_instance.nat_tune_y
            replace_dict["QDX"] = lhc_instance.drv_tune_x
            replace_dict["QDY"] = lhc_instance.drv_tune_y
        madx_script = madx_template % replace_dict
        return madx_script

    @classmethod
    def prepare_run(cls, lhc_instance, output_path):
        file_name = cls.ERR_DEF_FILES[str(lhc_instance.energy)]
        file_path = os.path.join(cls.ERR_DEF_PATH, file_name)
        # TODO: Windows?
        link_path = os.path.join(output_path, "error_deff.txt")
        if os.path.isfile(link_path):
            os.unlink(link_path)
        os.symlink(file_path, link_path)

    @classmethod
    def start_from_terminal(cls):
        parser = cls._get_arg_parser()
        options = parser.parse_args()
        lhc_mode = options.lhc_mode
        instance = cls.LHC_MODES[lhc_mode]()
        instance.beam = options.beam
        instance.nat_tune_x = options.nat_tune_x
        instance.nat_tune_y = options.nat_tune_y
        if options.acd and options.adt:
            raise model_creator.ModelCreationError(
                "Select only one excitation type."
            )
        if options.acd:
            instance.excitation = LhcExcitationMode.ACD
        elif options.adt:
            instance.excitation = LhcExcitationMode.ADT
        else:
            instance.excitation = LhcExcitationMode.FREE
        if options.acd or options.adt:
            instance.drv_tune_x = options.drv_tune_x
            instance.drv_tune_y = options.drv_tune_y
        instance.dpp = options.dpp
        instance.energy = options.energy
        instance.optics_file = options.optics
        cls.create_model(instance, options.output)

    @classmethod
    def _get_full_tunes(cls, lhc_instance):
        iqx, iqy = (lhc_instance.INT_TUNE_X +
                    lhc_instance.nat_tune_x,
                    lhc_instance.INT_TUNE_Y +
                    lhc_instance.nat_tune_y)
        return iqx, iqy

    @classmethod
    def _get_arg_parser(cls):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "lhc",
            help=(
                "This should always be the first" +
                "argument in LHC model creator."
            ),
            type=str,
        )
        parser.add_argument(
            "--lhcmode",
            help=("LHC mode to use. Should be one of: " +
                  str(cls.LHC_MODES.keys())),
            required=True,
            dest="lhc_mode",
            type=str,
        )
        parser.add_argument(
            "--beam",
            help="Beam to use.",
            required=True,
            dest="beam",
            type=int,
        )
        parser.add_argument(
            "--nattunex",
            help="Natural tune X without integer part.",
            required=True,
            dest="nat_tune_x",
            type=float,
        )
        parser.add_argument(
            "--nattuney",
            help="Natural tune Y without integer part.",
            required=True,
            dest="nat_tune_y",
            type=float,
        )
        parser.add_argument(
            "--acd",
            help="Activate excitation with ACD.",
            dest="acd",
            action="store_true",
        )
        parser.add_argument(
            "--adt",
            help="Activate excitation with ADT.",
            dest="adt",
            action="store_true",
        )
        parser.add_argument(
            "--drvtunex",
            help="Driven tune X without integer part.",
            dest="drv_tune_x",
            type=float,
        )
        parser.add_argument(
            "--drvtuney",
            help="Driven tune Y without integer part.",
            dest="drv_tune_y",
            type=float,
        )
        parser.add_argument(
            "--dpp",
            help="Delta p/p to use.",
            dest="dpp",
            default=0.0,
            type=float,
        )
        parser.add_argument(
            "--energy",
            help="Energy in Tev.",
            dest="energy",
            type=float,
        )
        parser.add_argument(
            "--optics",
            help="Path to the optics file to use (modifiers file).",
            dest="optics",
            type=str,
        )
        parser.add_argument(
            "--output",
            help="Output path for model, twiss files will be writen here.",
            dest="output",
            type=str,
        )
        return parser


class LhcBestKnowledgeCreator(LhcModelCreator):

    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        madx_template = None
        with open(lhc_instance.get_best_knowledge_tmpl()) as textfile:
            madx_template = textfile.read()
        iqx, iqy = LhcBestKnowledgeCreator._get_full_tunes(lhc_instance)
        replace_dict = {
            "RUN": lhc_instance.MACROS_NAME,
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.beam,
            "PATH": output_path,
            "DPP": lhc_instance.dpp,
            "QMX": iqx,
            "QMY": iqy,
            "ENERGY": lhc_instance.energy,
        }
        madx_script = madx_template % replace_dict
        return madx_script

    @classmethod
    def _get_arg_parser(cls):
        parser = LhcModelCreator._get_arg_parser()
        options = parser.parse_args()
        if options.acd or options.adt:
            raise model_creator.ModelCreationError(
                "Don't set ACD or ADT for best knowledge model."
            )
        return parser
