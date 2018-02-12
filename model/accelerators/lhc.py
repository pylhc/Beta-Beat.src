from __future__ import print_function
import os
import re
import json
from collections import OrderedDict
import numpy as np
from Utilities import tfs_pandas
from accelerator import Accelerator, AcceleratorDefinitionError, Element
from Utilities.entrypoint import EntryPoint, EntryPointParameters, split_arguments

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
    """ Parent Class for Lhc-Types.

    Keyword Args:
        Required
        nat_tune_x (float): Natural tune X without integer part.
                            **Flags**: ['--nattunex']
        nat_tune_y (float): Natural tune Y without integer part.
                            **Flags**: ['--nattuney']
        optics (str): Path to the optics file to use (modifiers file).
                      **Flags**: ['--optics']

        Optional
        acd (bool): Activate excitation with ACD.
                    **Flags**: ['--acd']
                    **Default**: ``False``
        adt (bool): Activate excitation with ADT.
                    **Flags**: ['--adt']
                    **Default**: ``False``
        dpp (float): Delta p/p to use.
                     **Flags**: ['--dpp']
                     **Default**: ``0.0``
        drv_tune_x (float): Driven tune X without integer part.
                            **Flags**: ['--drvtunex']
        drv_tune_y (float): Driven tune Y without integer part.
                            **Flags**: ['--drvtuney']
        energy (float): Energy in Tev.
                        **Flags**: ['--energy']
        fullresponse (bool): If True, fullresponse template will be filled
        and put in the output directory.
                             **Flags**: ['--fullresponse']
                             **Default**: ``False``
        xing (bool): If True, x-ing  angles will be applied to model
                     **Flags**: ['--xing']
                     **Default**: ``False``
    """
    NAME = "lhc"
    MACROS_NAME = "lhc"

    @staticmethod
    def get_class_parameters():
        params = EntryPointParameters()
        params.add_parameter(
            flags=["--lhcmode"],
            help=("LHC mode to use. Should be one of: " +
                  str(get_lhc_modes().keys())),
            name="lhc_mode",
            type=str,
            choices=get_lhc_modes().keys()
        )
        params.add_parameter(
            flags=["--beam"],
            help="Beam to use.",
            name="beam",
            type=int,
        )
        return params

    @staticmethod
    def get_instance_parameters():
        params = EntryPointParameters()
        params.add_parameter(
            flags=["--nattunex"],
            help="Natural tune X without integer part.",
            required=True,
            name="nat_tune_x",
            type=float,
        )
        params.add_parameter(
            flags=["--nattuney"],
            help="Natural tune Y without integer part.",
            required=True,
            name="nat_tune_y",
            type=float,
        )
        params.add_parameter(
            flags=["--acd"],
            help="Activate excitation with ACD.",
            name="acd",
            action="store_true"
        )
        params.add_parameter(
            flags=["--adt"],
            help="Activate excitation with ADT.",
            name="adt",
            action="store_true",
        )
        params.add_parameter(
            flags=["--drvtunex"],
            help="Driven tune X without integer part.",
            name="drv_tune_x",
            type=float,
        )
        params.add_parameter(
            flags=["--drvtuney"],
            help="Driven tune Y without integer part.",
            name="drv_tune_y",
            type=float,
        )
        params.add_parameter(
            flags=["--dpp"],
            help="Delta p/p to use.",
            name="dpp",
            default=0.0,
            type=float,
        )
        params.add_parameter(
            flags=["--energy"],
            help="Energy in Tev.",
            name="energy",
            type=float,
        )
        params.add_parameter(
            flags=["--optics"],
            help="Path to the optics file to use (modifiers file).",
            name="optics",
            required=True,
            type=str,
        )
        params.add_parameter(
            flags=["--fullresponse"],
            help=("If True, fullresponse template will "
                  "be filled and put in the output directory."),
            name="fullresponse",
            action="store_true",
        )
        params.add_parameter(
            flags=["--xing"],
            help=("If True, x-ing  angles will be applied to model"),
            name="xing",
            action="store_true",
        )
        return params

    # Entry-Point Wrappers #####################################################

    def __init__(self, *args, **kwargs):
        # for reasons of import-order and class creation, decoration was not possible
        parser = EntryPoint(self.get_instance_parameters(), strict=True)
        opt = parser.parse(*args, **kwargs)
        self.nat_tune_x = opt.nat_tune_x
        self.nat_tune_y = opt.nat_tune_y
        if opt.acd and opt.adt:
            raise AcceleratorDefinitionError(
                "Select only one excitation type."
            )
        if opt.acd:
            self.excitation = LhcExcitationMode.ACD
        elif opt.adt:
            self.excitation = LhcExcitationMode.ADT
        else:
            self.excitation = LhcExcitationMode.FREE

        if opt.acd or opt.adt:
            # "required"
            self.drv_tune_x = opt.drv_tune_x
            self.drv_tune_y = opt.drv_tune_y

        # required
        self.optics_file = opt.optics

        # optional with default
        self.dpp = opt.dpp
        self.fullresponse = opt.fullresponse

        # optional no default
        self.energy = opt.get("energy", None)
        self.xing = opt.get("xing", None)
        self.verify_object()

    @classmethod
    def init_and_get_unknowns(cls, args=None):
        """ Initializes but also returns unknowns.

         For the desired philosophy of returning parameters all the time,
         try to avoid this function, e.g. parse outside parameters first.
         """
        opt, rest_args = split_arguments(args, cls.get_instance_parameters())
        return cls(opt), rest_args

    @classmethod
    def get_class(cls, *args, **kwargs):
        """ Returns LHC subclass .

        Keyword Args:
            Optional
            beam (int): Beam to use.
                        **Flags**: ['--beam']
            lhc_mode (str): LHC mode to use.
                            **Flags**: ['--lhcmode']
                            **Choices**: ['lhc_runII_2016_ats', 'hllhc12', 'hllhc10', 'lhc_runI',
                            'lhc_runII', 'lhc_runII_2016', 'lhc_runII_2017']

        Returns:
            Lhc subclass.
        """
        parser = EntryPoint(cls.get_class_parameters(), strict=True)
        opt = parser.parse(*args, **kwargs)
        return cls._get_class(opt)

    @classmethod
    def get_class_and_unknown(cls, *args, **kwargs):
        """ Returns LHC subclass and unkown args .

        For the desired philosophy of returning parameters all the time,
        try to avoid this function, e.g. parse outside parameters first.
        """
        parser = EntryPoint(cls.get_class_parameters(), strict=False)
        opt, unknown_opt = parser.parse(*args, **kwargs)
        return cls._get_class(opt), unknown_opt

    @classmethod
    def _get_class(cls, opt):
        """ Actual get_class function """
        new_class = cls
        if opt.lhc_mode is not None:
            new_class = get_lhc_modes()[opt.lhc_mode]
        if opt.beam is not None:
            new_class = cls._get_beamed_class(new_class, opt.beam)
        return new_class

    # Public Methods ##########################################################

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
    def get_nominal_multidpp_tmpl(cls):
        return cls.get_file("nominal_multidpp.madx")
    
    @classmethod
    def get_coupling_tmpl(cls):
        return cls.get_file("coupling_correct.madx")

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

    # Private Methods ##########################################################

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
        if excitation_mode not in (LhcExcitationMode.FREE,
                                   LhcExcitationMode.ACD,
                                   LhcExcitationMode.ADT):
            raise ValueError("Wrong excitation mode.")
        self._excitation = excitation_mode


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


class _LhcB2Mixin(object):
    @classmethod
    def get_beam(cls):
        return 2


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


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
