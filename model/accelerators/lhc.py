from __future__ import print_function
import os
import argparse
import re
import sys
import json
from collections import OrderedDict
import numpy as np
import pandas as pd
from Utilities import tfs_pandas
from accelerator import Accelerator, AcceleratorDefinitionError, Element, get_commonbpm, AccExcitationMode
from time import time

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
        
        # for GetLLM
        self.model = None
        self.model_driven = None
        self.model_best_knowledge = None
        self.elements = None
        self.elements_centre = None

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
            instance.excitation = AccExcitationMode.ACD
        elif options.adt:
            instance.excitation = AccExcitationMode.ADT
        else:
            instance.excitation = AccExcitationMode.FREE
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
    def init_from_model_dir(cls, model_dir):  # prints only for debugging
        
        print("Creating accelerator instance from model dir")
        instance = cls()
        
        instance.model_tfs = tfs_pandas.read_tfs(os.path.join(model_dir, "twiss.dat")).set_index("NAME")
        print("model path =", os.path.join(model_dir, "twiss.dat"))
            
        instance._excitation = AccExcitationMode.FREE
        ac_filename = os.path.join(model_dir, "twiss_ac.dat")
        adt_filename = os.path.join(model_dir, "twiss_adt.dat")
        
        if os.path.isfile(ac_filename):
            instance.model_driven = tfs_pandas.read_tfs(ac_filename).set_index("NAME")
            instance.excitation = AccExcitationMode.ACD
            driven_filename = ac_filename
            
        if os.path.isfile(adt_filename):
            if instance.excitation == AccExcitationMode.ACD:
                raise AcceleratorDefinitionError("ADT as well as ACD models provided. What do you want? Please come back to me once you have made up your mind.")

            instance.model_driven = tfs_pandas.read_tfs(adt_filename).set_index("NAME")
            instance.excitation = AccExcitationMode.ADT
            driven_filename = adt_filename
        
        try:
            model_best_knowledge_path = os.path.join(model_dir, "twiss_best_knowledge.dat")
            if os.path.isfile(model_best_knowledge_path):
                instance.model_best_knowledge = tfs_pandas.read_tfs(model_best_knowledge_path).set_index("NAME")
        except IOError:
            instance.model_best_knowledge = None
            
        elements_path = os.path.join(model_dir, "twiss_elements.dat")
        if os.path.isfile(elements_path):
            instance.elements = tfs_pandas.read_tfs(elements_path).set_index("NAME")
        else:
            raise AcceleratorDefinitionError("Elements twiss not found")
        elements_path = os.path.join(model_dir, "twiss_elements_centre.dat")
        if os.path.isfile(elements_path):
            instance.elements_centre = tfs_pandas.read_tfs(elements_path).set_index("NAME")
        else:
            instance.elements_centre = instance.elements
        
        instance.drv_tune_x = float(instance.get_driven_tfs().headers["Q1"])
        instance.drv_tune_y = float(instance.get_driven_tfs().headers["Q2"])
        instance.nat_tune_x = float(instance.model_tfs.headers["Q1"])
        instance.nat_tune_y = float(instance.model_tfs.headers["Q2"])
        
        return instance

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
    def get_segment(cls, label, first_elem, last_elem, optics_file):
        segment_cls = type(cls.__name__ + "Segment",
                           (_LhcSegmentMixin, cls),
                           {})
        segment_inst = segment_cls()
        beam = cls.get_beam()
        bpms_file_name = "beam1bpms.tfs" if beam == 1 else "beam2bpms.tfs"
        bpms_file = _get_file_for_year(cls.YEAR, bpms_file_name)
        bpms_file_data = tfs_pandas.read_tfs(bpms_file).set_index("NAME")
        first_elem_s = bpms_file_data.loc[first_elem, "S"]
        last_elem_s = bpms_file_data.loc[last_elem, "S"]
        segment_inst.label = label
        segment_inst.start = Element(first_elem, first_elem_s)
        segment_inst.end = Element(last_elem, last_elem_s)
        segment_inst.optics_file = optics_file
        segment_inst.xing = False
        segment_inst.verify_object()
        return segment_inst

    @classmethod
    def _get_beamed_class(cls, new_class, beam):
        beam_mixin = _LhcB1Mixin if beam == 1 else _LhcB2Mixin
        beamed_class = type(new_class.__name__ + "B" + str(beam),
                            (beam_mixin, new_class),
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
        if (self.excitation == AccExcitationMode.ACD or
                self.excitation == AccExcitationMode.ADT):
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
    def get_file(cls, filename):
        return os.path.join(CURRENT_DIR, "lhc", filename)

    @classmethod
    def get_variables(cls, frm=None, to=None, classes=None):
        correctors_dir = os.path.join(LHC_DIR, "2012", "correctors")
        all_corrs = _merge_jsons(
            os.path.join(correctors_dir, "correctors_b" + str(cls.get_beam()),
                         "beta_correctors.json"),
            os.path.join(correctors_dir, "correctors_b" + str(cls.get_beam()),
                         "coupling_correctors.json"),
            cls._get_triplet_correctors_file(),
        )
        my_classes = classes
        if my_classes is None:
            my_classes = all_corrs.keys()
        vars_by_class = set(_flatten_list(
            [all_corrs[corr_cls] for corr_cls in my_classes if corr_cls in all_corrs])
        )
        elems_matrix = tfs_pandas.read_tfs(
            cls._get_corrector_elems()
        ).sort_values("S").set_index("S").loc[frm:to, :]
        vars_by_position = _remove_dups_keep_order(_flatten_list(
            [raw_vars.split(",") for raw_vars in elems_matrix.loc[:, "VARS"]]
        ))
        return _list_intersect_keep_order(vars_by_position, vars_by_class)

    @classmethod
    def _get_triplet_correctors_file(cls):
        correctors_dir = os.path.join(LHC_DIR, "2012", "correctors")
        return os.path.join(correctors_dir, "triplet_correctors.json")

    @classmethod
    def _get_corrector_elems(cls):
        correctors_dir = os.path.join(LHC_DIR, "2012", "correctors")
        return os.path.join(correctors_dir,
                            "corrector_elems_b" + str(cls.get_beam()) + ".tfs")

    @property
    def excitation(self):
        return self._excitation

    @excitation.setter
    def excitation(self, excitation_mode):
        if excitation_mode not in (AccExcitationMode.FREE,
                                   AccExcitationMode.ACD,
                                   AccExcitationMode.ADT):
            raise ValueError("Wrong excitation mode.")
        self._excitation = excitation_mode
        
        
    # For GetLLM --------------------------------------------------------------
     
    def get_exciter_bpm(self, plane, commonbpms):
        
        if self.get_beam() == 1:
            if self.excitation == AccExcitationMode.ACD:
                return get_commonbpm("BPMYA.5L4.B1", "BPMYB.6L4.B1", commonbpms)
               
            elif self.excitation == AccExcitationMode.ADT:
                if plane == "H":
                    return get_commonbpm("BPMWA.B5L4.B1", "BPMWA.A5L4.B1", commonbpms)
                elif plane == "V":
                    return get_commonbpm("BPMWA.B5R4.B1", "BPMWA.A5R4.B1", commonbpms)
        elif self.get_beam() == 2:
            if self.excitation == AccExcitationMode.ACD:
                return get_commonbpm("BPMYB.5L4.B2", "BPMYA.6L4.B2", commonbpms)
            elif self.excitation == AccExcitationMode.ADT:
                if plane == "H":
                    return get_commonbpm("BPMWA.B5R4.B2", "BPMWA.A5R4.B2", commonbpms)
                elif plane == "V":
                    return get_commonbpm("BPMWA.B5L4.B2", "BPMWA.A5L4.B2", commonbpms)
        return None
    
    def get_exciter_name(self, plane):
        if self.beam() == 1:
            if self.excitation == AccExcitationMode.ACD:
                if plane == "H":
                    return 'MKQA.6L4.B1'
                elif plane == "V":
                    return 'MKQA.6L4.B1'
            elif self.excitation == AccExcitationMode.ADT:
                if plane == "H":
                    return "ADTKH.C5L4.B1"
                elif plane == "V":
                    return "ADTKV.B5R4.B1"
        elif self.beam() == 2:
            if self.excitation == AccExcitationMode.ACD:
                if plane == "H":
                    return 'MKQA.6L4.B2'
                elif plane == "V":
                    return 'MKQA.6L4.B2'
            elif self.excitation == AccExcitationMode.ADT:
                if plane == "H":
                    return "ADTKH.B5R4.B2"
                elif plane == "V":
                    return "ADTKV.C5L4.B2"
        return None
    
    def get_s_first_BPM(self):
        if self.get_beam() == 1:
            return self.model_tfs["S"][self.model_tfs.indx["BPMSW.1L2.B1"]]
        elif self.get_beam() == 2:
            return self.model_tfs["S"][self.model_tfs.indx["BPMSW.1L8.B2"]]
        return None
        
    def get_model_tfs(self):
        return self.model_tfs
        
    def get_driven_tfs(self):
        if self.model_driven is None:
            return self.model_tfs
        return self.model_driven

    def get_best_knowledge_model_tfs(self):
        if self.model_best_knowledge is None:
            return self.model_tfs
        return self.model_best_knowledge
    
    def get_elements_tfs(self):
        return self.elements

    def get_elements_centre_tfs(self):
        return self.elements_centre

class _LhcSegmentMixin(object):

    def __init__(self):
        self._start = None
        self._end = None

    def get_segment_vars(self, classes=None):
        return self.get_variables(frm=self.start.s,
                                  to=self.end.s,
                                  classes=classes)

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
    
    @classmethod
    def get_beam_direction(cls):
        return 1


class _LhcB2Mixin(object):
    @classmethod
    def get_beam(cls):
        return 2
    
    @classmethod
    def get_beam_direction(cls):
        return -1

class LhcAts(Lhc):
    MACROS_NAME = "lhc_runII_ats"


# Specific accelerator definitions ###########################################

class LhcRunI(Lhc):
    YEAR = "2012"

    @classmethod
    def load_main_seq_madx(cls):
        load_main_seq = _get_call_main_for_year("2012")
        load_main_seq += _get_madx_call_command(
            os.path.join(LHC_DIR, "2012", "install_additional_elements.madx")
        )
        return load_main_seq


class LhcRunII2015(Lhc):
    YEAR = "2015"

    @classmethod
    def load_main_seq_madx(cls):
        return _get_call_main_for_year("2015")


class LhcRunII2016(Lhc):
    YEAR = "2016"

    @classmethod
    def load_main_seq_madx(cls):
        return _get_call_main_for_year("2016")


class LhcRunII2016Ats(LhcAts, LhcRunII2016):
    pass


class LhcRunII2017(LhcAts):
    YEAR = "2017"

    @classmethod
    def load_main_seq_madx(cls):
        return _get_call_main_for_year("2017")


class HlLhc10(LhcAts):
    MACROS_NAME = "hllhc"
    YEAR = "hllhc1.0"

    @classmethod
    def load_main_seq_madx(cls):
        load_main_seq = _get_call_main_for_year("2015")
        load_main_seq += _get_call_main_for_year("hllhc1.0")
        return load_main_seq


class HlLhc12(LhcAts):
    MACROS_NAME = "hllhc"
    YEAR = "hllhc1.2"

    @classmethod
    def load_main_seq_madx(cls):
        load_main_seq = _get_call_main_for_year("2015")
        load_main_seq += _get_call_main_for_year("hllhc1.2")
        return load_main_seq

    @classmethod
    def _get_triplet_correctors_file(cls):
        correctors_dir = os.path.join(LHC_DIR, "hllhc1.2", "correctors")
        return os.path.join(correctors_dir, "triplet_correctors.json")

    @classmethod
    def _get_corrector_elems(cls):
        correctors_dir = os.path.join(LHC_DIR, "hllhc1.2", "correctors")
        return os.path.join(correctors_dir,
                            "corrector_elems_b" + str(cls.get_beam()) + ".tfs")


class HlLhc12NewCircuit(LhcAts):
    MACROS_NAME = "hllhc"
    YEAR = "hllhc12"

    

class HlLhc12NoQ2Trim(HlLhc12):
    MACROS_NAME = "hllhc"
    YEAR = "hllhc12"

##############################################################################


# General functions ##########################################################

def _get_call_main_for_year(year):
    call_main = _get_madx_call_command(
        _get_file_for_year(year, "main.seq")
    )
    return call_main


def _get_madx_call_command(path_to_call):
    command = "call, file = \""
    command += path_to_call
    command += "\";\n"
    return command


def _get_file_for_year(year, filename):
    return os.path.join(LHC_DIR, year, filename)


def _merge_jsons(*files):
    full_dict = {}
    for json_file in files:
        with open(json_file, "r") as json_data:
            json_dict = json.load(json_data)
            for key, value in json_dict.iteritems():
                full_dict[key] = value
    return full_dict


def _flatten_list(my_list):
    return [item for sublist in my_list for item in sublist]


def _remove_dups_keep_order(my_list):
    return list(OrderedDict.fromkeys(my_list))


def _list_intersect_keep_order(primary_list, secondary_list):
    return [elem for elem in primary_list if elem in secondary_list]


##############################################################################
