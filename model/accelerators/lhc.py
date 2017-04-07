from __future__ import print_function
import sys
import os
import argparse
from model_creators import model_creator
from accelerator import Accelerator


CURRENT_DIR = os.path.dirname(__file__)


def get_lhc_modes():
    return {
        "lhc_runI": LhcRunI,
        "lhc_runII": LhcRunII2015,
        "lhc_runII_2016": LhcRunII2016,
        "lhc_runII_2016_ats": LhcRunII2016Ats,
        "hllhc": HlLhc,
    }


class LhcExcitationMode(object):
    FREE, ACD, ADT = range(3)


class Lhc(Accelerator):
    INT_TUNE_X = 64.
    INT_TUNE_Y = 59.

    def __init__(self):
        self.optics_file = None
        self.nat_tune_x = None
        self.nat_tune_y = None
        self._excitation = None
        self.drv_tune_x = None
        self.drv_tune_y = None
        self.energy = None
        self.dpp = 0.0

    @classmethod
    def init_from_args(cls, args):
        parser = cls._get_arg_parser()
        options, rest_args = parser.parse_known_args(args)
        instance = cls()
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
        instance.fullresponse = options.fullresponse
        instance.verify_object()
        return instance, rest_args

    @classmethod
    def get_class(cls, lhc_mode, beam):
        return Lhc._get_beamed_class(get_lhc_modes()[lhc_mode], beam)

    @classmethod
    def get_class_from_args(cls, args):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--lhcmode",
            help=("LHC mode to use. Should be one of: " +
                  str(get_lhc_modes().keys())),
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
        options, rest_args = parser.parse_known_args(args)
        lhc_mode = options.lhc_mode
        beam = options.beam
        return Lhc.get_class(lhc_mode, beam), rest_args

    @classmethod
    def _get_beamed_class(this_class, new_class, beam):
        beam_mixin = _LhcB1Mixin if beam == 1 else _LhcB2Mixin
        beamed_class = type(new_class.__name__ + "B" + str(beam),
                            (new_class, beam_mixin),
                            {})
        return beamed_class

    @classmethod
    def _get_arg_parser(cls):
        parser = argparse.ArgumentParser()
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
            required=True,
            type=str,
        )
        parser.add_argument(
            "--fullresponse",
            help=("If present, fullresponse template will" +
                  "be filled and put in the output directory."),
            dest="fullresponse",
            action="store_true",
        )
        return parser

    def verify_object(self):  # TODO: Maybe more checks?
        if self.excitation is None:
            raise LhcError("Excitation mode not set.")
        if (self.excitation == LhcExcitationMode.ACD or
                self.excitation == LhcExcitationMode.ADT):
            if self.drv_tune_x is None or self.drv_tune_y is None:
                raise LhcError("Driven tunes not set.")

    def get_nominal_tmpl(self):
        return self.get_file("nominal.madx")

    def get_best_knowledge_tmpl(self):
        return self.get_file("best_knowledge.madx")

    def get_file(self, filename, beam=None):
        if beam is None:
            return os.path.join(CURRENT_DIR, "lhc", filename)
        elif beam == 1:
            return os.path.join(CURRENT_DIR, "lhc", "beam1", filename)
        elif beam == 2:
            return os.path.join(CURRENT_DIR, "lhc", "beam2", filename)
        else:
            raise model_creator.ModelCreationError("Beam should be 1 or 2.")

    @property
    def excitation(self):
        return self._excitation

    @excitation.setter
    def excitation(self, excitation_mode):
        if excitation_mode not in (LhcExcitationMode.FREE,
                                   LhcExcitationMode.ACD,
                                   LhcExcitationMode.ADT):
            raise ValueError("Wrong excitation mode.")
        self._excitation = excitation_mode


class _LhcB1Mixin(object):
    @classmethod
    def get_beam(cls):
        return 1


class _LhcB2Mixin(object):
    @classmethod
    def get_beam(cls):
        return 2


class LhcAts(Lhc):
    INT_TUNE_X = 62.
    INT_TUNE_Y = 60.


class LhcRunI(Lhc):
    MACROS_NAME = "lhc_runI"


class LhcRunII2015(Lhc):
    MACROS_NAME = "lhc_runII"


class LhcRunII2016(Lhc):
    MACROS_NAME = "lhc_runII_2016"


class LhcRunII2016Ats(LhcAts):
    MACROS_NAME = "lhc_runII_2016_ats"


class HlLhc(LhcAts):
    MACROS_NAME = "hllhc"


class LhcError(Exception):
    """
    Raised when an LHC creation error occurs.
    """
    pass


if __name__ == "__main__":
    print("Import this module.", file=sys.stderr)
