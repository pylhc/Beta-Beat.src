from __future__ import print_function
import os
import logging
from accelerators.lhc import LhcExcitationMode
from accelerators import lhc
import argparse
import model_creator

ACCEL_DIR = os.path.join(model_creator.CURRENT_DIR, "accelerators")

LOGGER = logging.getLogger(__name__)


class LhcModelCreator(model_creator.ModelCreator):

    LHC_DIR = os.path.join(ACCEL_DIR, "lhc")
    NOMINAL_TEMPL = os.path.join(LHC_DIR, "lhc_nominal.madx")
    LHC_MODES = {
        "lhc_runI": lhc.LhcRunI,
        "lhc_runII": lhc.LhcRunII2015,
        "lhc_runII_2016": lhc.LhcRunII2016,
        "lhc_runII_2016_ats": lhc.LhcRunII2016Ats,
        "hllhc": lhc.HlLhc,
    }

    @staticmethod
    def get_madx_script(lhc_instance, output_path):
        madx_template = None
        with open(LhcModelCreator.NOMINAL_TEMPL) as textfile:
            madx_template = textfile.read()
        iqx, iqy = (lhc_instance.INT_TUNE_X +
                    lhc_instance.nat_tune_x,
                    lhc_instance.INT_TUNE_Y +
                    lhc_instance.nat_tune_y)
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

    @staticmethod
    def start_from_terminal():
        options = LhcModelCreator._parse_args()
        lhc_mode = options.lhc_mode
        instance = LhcModelCreator.LHC_MODES[lhc_mode]()
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
        instance.optics_file = options.optics
        LhcModelCreator.create_model(instance, options.output)

    @staticmethod
    def _parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "lhc",
            help=(
                "This should always be the first" +
                "argument in lhc model creator."
            ),
            type=str,
        )
        parser.add_argument(
            "--lhcmode",
            help=("LHC mode to use. Should be one of: " +
                  str(LhcModelCreator.LHC_MODES.keys())),
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
        options = parser.parse_args()
        return options
