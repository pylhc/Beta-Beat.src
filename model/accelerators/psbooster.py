import os
import argparse
from accelerator import Accelerator


CURRENT_DIR = os.path.dirname(__file__)
PSB_DIR = os.path.join(CURRENT_DIR, "psbooster")


class Psbooster(Accelerator):
    
    NAME = "psbooster"

    def __init__(self):
        self.nat_tune_x = None
        self.nat_tune_y = None
        self.acd = None
        self.drv_tune_x = None
        self.drv_tune_y = None
        self.energy = None
        self.fullresponse = False

    @classmethod
    def get_class_from_args(cls, args):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--ring",
            help="Ring to use.",
            dest="ring",
            type=int,
        )
        options, rest_args = parser.parse_known_args(args)
        ring = options.ring
        return cls.get_class(ring), rest_args

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
            "--energy",
            help="Energy in Tev.",
            dest="energy",
            type=float,
        )
        parser.add_argument(
            "--fullresponse",
            help=("If present, fullresponse template will" +
                  "be filled and put in the output directory."),
            dest="fullresponse",
            action="store_true",
        )
        return parser

    @classmethod
    def init_from_args(cls, args):
        parser = cls._get_arg_parser()
        options, rest_args = parser.parse_known_args(args)
        instance = cls()
        instance.nat_tune_x = options.nat_tune_x
        instance.nat_tune_y = options.nat_tune_y
        instance.acd = options.acd
        if options.acd:
            instance.drv_tune_x = options.drv_tune_x
            instance.drv_tune_y = options.drv_tune_y
        instance.energy = options.energy
        instance.fullresponse = options.fullresponse
        return instance, rest_args

    def verify_object(self):
	pass

    @classmethod
    def get_nominal_tmpl(cls):
        return os.path.join(PSB_DIR, "nominal.madx")

    @classmethod
    def get_psb_dir(cls):
        return PSB_DIR

    @classmethod
    def get_arc_bpms_mask(cls, list_of_elements):
        # TODO: Anaaaaaa
        pass

    @classmethod
    def get_class(cls, ring=None):
        new_class = cls
        if ring not in (1, 2, 3, 4):
            raise ValueError("Ring must be 1, 2, 3 or 4.")
        if ring is not None:
            new_class = type(
                new_class.__name__ + "Ring{}".format(ring),
                (new_class, ),
                {"get_ring": classmethod(lambda cls: ring)}
            )
        return new_class
