import os
from accelerator import Accelerator
from utils.entrypoint import EntryPoint, EntryPointParameters, split_arguments

CURRENT_DIR = os.path.dirname(__file__)
PSB_DIR = os.path.join(CURRENT_DIR, "psbooster")


class Psbooster(Accelerator):
    """ Parent Class for Psbooster-Types.

    Keyword Args:
        Required
        nat_tune_x (float): Natural tune X without integer part.
                            **Flags**: ['--nattunex']
        nat_tune_y (float): Natural tune Y without integer part.
                            **Flags**: ['--nattuney']

        Optional
        acd (bool): Activate excitation with ACD.
                    **Flags**: ['--acd']
                    **Default**: ``False``
        drv_tune_x (float): Driven tune X without integer part.
                            **Flags**: ['--drvtunex']
        drv_tune_y (float): Driven tune Y without integer part.
                            **Flags**: ['--drvtuney']
        energy (float): Energy in Tev.
                        **Flags**: ['--energy']
        fullresponse (bool): If present, fullresponse template willbe filled and put
                             in the output directory.
                             **Flags**: ['--fullresponse']
                             **Default**: ``False``

    """
    NAME = "psbooster"

    @staticmethod
    def get_class_parameters():
        params = EntryPointParameters()
        params.add_parameter(
            flags=["--ring"],
            help="Ring to use.",
            name="ring",
            type=int,
            choices=[1, 2, 3, 4]
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
            flags=["--energy"],
            help="Energy in Tev.",
            name="energy",
            type=float,
        )
        params.add_parameter(
            flags=["--fullresponse"],
            help=("If present, fullresponse template will" +
                  "be filled and put in the output directory."),
            name="fullresponse",
            action="store_true",
        )
        return params

    # Entry-Point Wrappers #####################################################

    def __init__(self, *args, **kwargs):
        # for reasons of import-order and class creation, decoration was not possible
        parser = EntryPoint(self.get_instance_parameters(), strict=True)
        opt = parser.parse(*args, **kwargs)

        # required
        self.nat_tune_x = opt.nat_tune_x
        self.nat_tune_y = opt.nat_tune_y
        self.acd = opt.acd
        if self.acd:
            self.drv_tune_x = opt.drv_tune_x
            self.drv_tune_y = opt.drv_tune_y

        # optional with default
        self.fullresponse = opt.fullresponse

        # optional w/o default
        self.energy = opt.get("energy", None)

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
        """ Returns Psbooster class.

        Keyword Args:
            Optional
            ring (int): Ring to use.
                        **Flags**: ['--ring']
                        **Choices**: [1, 2, 3, 4]

        Returns:
            Psbooster class.
        """
        parser = EntryPoint(cls.get_class_parameters(), strict=True)
        opt = parser.parse(*args, **kwargs)
        return cls._get_class(opt)

    @classmethod
    def get_class_and_unknown(cls, *args, **kwargs):
        """ Returns Psbooster subclass and unkown args .

        For the desired philosophy of returning parameters all the time,
        try to avoid this function, e.g. parse outside parameters first.
        """
        parser = EntryPoint(cls.get_class_parameters(), strict=False)
        opt, unknown_opt = parser.parse(*args, **kwargs)
        return cls._get_class(opt), unknown_opt

    @classmethod
    def _get_class(cls, opt):
        new_class = cls
        if opt.ring is not None:
            new_class = type(
                new_class.__name__ + "Ring{}".format(opt.ring),
                (new_class,),
                {"get_ring": classmethod(lambda cls: opt.ring)}
            )
        return new_class

    # Public Methods ##########################################################

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

    # Private Methods ##########################################################


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
