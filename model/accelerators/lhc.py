from __future__ import print_function
import os
import argparse
import re
import numpy as np
from accelerator import Accelerator, AcceleratorDefinitionError


CURRENT_DIR = os.path.dirname(__file__)
LHC_DIR = os.path.join(CURRENT_DIR, "lhc")


def get_lhc_modes():
    return {
        "lhc_runI": LhcRunI,
        "lhc_runII": LhcRunII2015,
        "lhc_runII_2016": LhcRunII2016,
        "lhc_runII_2016_ats": LhcRunII2016Ats,
        "lhc_runII_2017": LhcRunII2017,
        "hllhc10": HlLhc10,
        "hllhc12": HlLhc12,
    }


class LhcExcitationMode(object):
    FREE, ACD, ADT = range(3)


class Lhc(Accelerator):
    NAME = "lhc"
    MACROS_NAME = "lhc"

    def __init__(self):
        self.optics_file = None
        self.nat_tune_x = None
        self.nat_tune_y = None
        self._excitation = None
        self.drv_tune_x = None
        self.drv_tune_y = None
        self.energy = None
        self.dpp = 0.0
        self.xing = None

    @classmethod
    def init_from_args(cls, args):
        parser = cls._get_arg_parser()
        options, rest_args = parser.parse_known_args(args)
        instance = cls()
        instance.nat_tune_x = options.nat_tune_x
        instance.nat_tune_y = options.nat_tune_y
        if options.acd and options.adt:
            raise AcceleratorDefinitionError(
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
        instance.xing = options.xing
        instance.verify_object()
        return instance, rest_args

    @classmethod
    def get_class(cls, lhc_mode=None, beam=None):
        new_class = cls
        if lhc_mode is not None:
            new_class = get_lhc_modes()[lhc_mode]
        if beam is not None:
            new_class = cls._get_beamed_class(new_class, beam)
        return new_class

    @classmethod
    def get_class_from_args(cls, args):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--lhcmode",
            help=("LHC mode to use. Should be one of: " +
                  str(get_lhc_modes().keys())),
            dest="lhc_mode",
            type=str,
        )
        parser.add_argument(
            "--beam",
            help="Beam to use.",
            dest="beam",
            type=int,
        )
        options, rest_args = parser.parse_known_args(args)
        lhc_mode = options.lhc_mode
        beam = options.beam
        return cls.get_class(lhc_mode, beam), rest_args

    @classmethod
    def _get_beamed_class(cls, new_class, beam):
        beam_mixin = _LhcB1Mixin if beam == 1 else _LhcB2Mixin
        beamed_class = type(new_class.__name__ + "B" + str(beam),
                            (new_class, beam_mixin),
                            {})
        return beamed_class

    @classmethod
    def get_arc_bpms_mask(cls, list_of_elements):
        mask = []
        pattern = re.compile("BPM.*\.([0-9]+)[RL].\..*", re.IGNORECASE)
        for element in list_of_elements:
            match = pattern.match(element)
            # The arc bpms are from BPM.14... and up
            if match and int(match.group(1)) > 14:
                mask.append(True)
            else:
                mask.append(False)
        return np.array(mask)

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
        parser.add_argument(
            "--xing",
            help=("If present, x-ing  angles will be applied to model"),
            dest="xing",
            action="store_true",
        )
        return parser

    def verify_object(self):  # TODO: Maybe more checks?
        try:
            self.get_beam()
        except AttributeError:
            raise AcceleratorDefinitionError(
                "The accelerator definition is incomplete, beam " +
                "has to be specified (--beam option missing?)."
            )
        if self.optics_file is None:
            raise AcceleratorDefinitionError(
                "The accelerator definition is incomplete, optics "
                "file has not been specified."
            )
        if self.excitation is None:
            raise AcceleratorDefinitionError("Excitation mode not set.")
        if self.xing is None:
            raise AcceleratorDefinitionError("Crossing on or off not set.")
        if (self.excitation == LhcExcitationMode.ACD or
                self.excitation == LhcExcitationMode.ADT):
            if self.drv_tune_x is None or self.drv_tune_y is None:
                raise AcceleratorDefinitionError("Driven tunes not set.")

    @classmethod
    def get_nominal_tmpl(cls):
        return cls.get_file("nominal.madx")

    @classmethod
    def get_best_knowledge_tmpl(cls):
        return cls.get_file("best_knowledge.madx")

    @classmethod
    def get_segment_tmpl(cls):
        return cls.get_file("segment.madx")

    @classmethod
    def get_file(cls, filename, beam=None):
        if beam is None:
            return os.path.join(CURRENT_DIR, "lhc", filename)
        elif beam == 1:
            return os.path.join(CURRENT_DIR, "lhc", "beam1", filename)
        elif beam == 2:
            return os.path.join(CURRENT_DIR, "lhc", "beam2", filename)
        else:
            raise AcceleratorDefinitionError("Beam should be 1 or 2.")

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


class LhcSegment(Lhc):

    def __init__(self):
        self.optics_file = None
        self.label = None
        self.start = None
        self.end = None
        self.xing = None

    @classmethod
    def init_from_args(cls, args):
        raise NotImplementedError(
            "LHC segments can only be instantiated by class definition."
        )

    @classmethod
    def get_class(cls, lhc_mode=None, beam=None):
        new_class = cls
        if lhc_mode is not None:
            new_class = cls._get_specific_class(get_lhc_modes()[lhc_mode])
        if beam is not None:
            new_class = cls._get_beamed_class(new_class, beam)
        return new_class

    @classmethod
    def _get_specific_class(cls, lhc_cls):
        specific_cls = type(lhc_cls.__name__ + "Segment",
                            (cls, lhc_cls),
                            {})
        return specific_cls

    def verify_object(self):
        try:
            self.get_beam()
        except AttributeError:
            raise AcceleratorDefinitionError(
                "The accelerator definition is incomplete, beam "
                "has to be specified (--beam option missing?)."
            )
        if self.optics_file is None:
            raise AcceleratorDefinitionError(
                "The accelerator definition is incomplete, optics "
                "file has not been specified."
            )
        if self.xing is None:
            raise AcceleratorDefinitionError("Crossing on or off not set.")
        if self.label is None:
            raise AcceleratorDefinitionError("Segment label not set.")
        if self.start is None:
            raise AcceleratorDefinitionError("Segment start not set.")
        if self.end is None:
            raise AcceleratorDefinitionError("Segment end not set.")


class _LhcB1Mixin(object):
    @classmethod
    def get_beam(cls):
        return 1


class _LhcB2Mixin(object):
    @classmethod
    def get_beam(cls):
        return 2


class LhcAts(Lhc):
    MACROS_NAME = "lhc_runII_ats"


# Specific accelerator definitions ###########################################

class LhcRunI(Lhc):

    @classmethod
    def load_main_seq_madx(cls):
        load_main_seq = _get_call_main_for_year("2012")
        load_main_seq += _get_madx_call_command(
            os.path.join(LHC_DIR, "2012", "install_additional_elements.madx")
        )
        return load_main_seq


class LhcRunII2015(Lhc):

    @classmethod
    def load_main_seq_madx(cls):
        return _get_call_main_for_year("2015")


class LhcRunII2016(Lhc):

    @classmethod
    def load_main_seq_madx(cls):
        return _get_call_main_for_year("2016")


class LhcRunII2016Ats(LhcAts, LhcRunII2016):
    pass


class LhcRunII2017(LhcAts):

    @classmethod
    def load_main_seq_madx(cls):
        return _get_call_main_for_year("2017")


class HlLhc10(LhcAts):
    MACROS_NAME = "hllhc"

    @classmethod
    def load_main_seq_madx(cls):
        load_main_seq = _get_call_main_for_year("2015")
        load_main_seq += _get_call_main_for_year("hllhc1.0")
        return load_main_seq


class HlLhc12(LhcAts):
    MACROS_NAME = "hllhc"

    @classmethod
    def load_main_seq_madx(cls):
        load_main_seq = _get_call_main_for_year("2015")
        load_main_seq += _get_call_main_for_year("hllhc1.2")
        return load_main_seq

##############################################################################


# General functions ##########################################################

def _get_madx_call_command(path_to_call):
    command = "call, file = \""
    command += path_to_call
    command += "\";\n"
    return command


def _get_call_main_for_year(year):
    call_main = _get_madx_call_command(
        os.path.join(LHC_DIR, year, "main.seq")
    )
    return call_main

##############################################################################
