import os
import pandas as pd
from model.accelerators.accelerator import Accelerator,Element
from utils.entrypoint import EntryPoint, EntryPointParameters, split_arguments
from utils import logging_tools
from tfs_files import tfs_pandas
from Python_Classes4MAD.metaclass import twiss
import logging

LOGGER = logging.getLogger(__name__)

CURRENT_DIR = os.path.dirname(__file__)
PSB_DIR = os.path.join(CURRENT_DIR, "psbooster")


class Psbooster(Accelerator):
    """ Parent Class for Psbooster-Types.

    Keyword Args:
        Required
        ring (int): Ring number.
                            **Flags**: ['--ring']
                            
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
            flags=["--dpp"],
            help="Delta p/p to use.",
            name="dpp",
            default=0.0,
            type=float,
        )

        params.add_parameter(
            flags=["--optics"],
            help="Path to the optics file to use (modifiers file).",
            name="optics",
            type=str,
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

        self.dpp = opt.get("dpp", 0.0)
        
        self.optics_file = opt.get("optics", None)

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
        else:
            print("No ring info in options")
        return new_class

    # Public Methods ##########################################################
    @classmethod
    def get_segment(cls, label, first_elem, last_elem, optics_file, twiss_file):
        segment_cls = type(cls.__name__ + "Segment",
                          (_PsboosterSegmentMixin,cls),
                          {})


        LOGGER.debug('twiss_file is <%s>',twiss_file)
        tw = twiss(twiss_file)
        
        LOGGER.debug('twiss_file has tunes %f %f ',tw.Q1,tw.Q2)
        ring = _get_ring_from_seqname(tw.SEQUENCE)

        #ring = cls.get_ring()

        segment_inst = segment_cls()
        

        bpms_file = _get_file_for_ring(ring)
        bpms_file_data = tfs_pandas.read_tfs(bpms_file).set_index("NAME")
        first_elem_s = bpms_file_data.loc[first_elem, "S"]
        last_elem_s = bpms_file_data.loc[last_elem, "S"]
        segment_inst.label = label
        segment_inst.start = Element(first_elem, first_elem_s)
        segment_inst.end = Element(last_elem, last_elem_s)
        segment_inst.optics_file = optics_file
        segment_inst.fullresponse = None
        

        segment_inst.nat_tune_x = tw.Q1
        segment_inst.nat_tune_y = tw.Q2
        segment_inst.energy = tw.ENERGY
        segment_inst.sequence = tw.SEQUENCE
        segment_inst.ring = ring
        segment_inst.kind = '' # '' means beta from phase, can be 'betaamp', in the future 'betakmod'
        
        return segment_inst    


    def verify_object(self):
        pass

    @classmethod
    def get_nominal_tmpl(cls):
        return os.path.join(PSB_DIR, "nominal.madx")

    @classmethod
    def get_segment_tmpl(cls):
        return cls.get_file("segment.madx")

    @classmethod
    def get_iteration_tmpl(cls):
        return cls.get_file("template.iterate.madx")

    @classmethod
    def get_corrtest_tmpl(cls):
        return cls.get_file("correction_test.madx")

    @classmethod
    def get_psb_dir(cls):
        return PSB_DIR

    @classmethod
    def get_element_types_mask(cls, list_of_elements, types):
        # TODO: Anaaaaaa
        raise NotImplementedError("Anaaaaa!")
        re_dict = {
            "bpm": r".*",
            "magnet": r".*",
            "arc_bpm": r".*",
        }

        unknown_elements = [ty for ty in types if ty not in re_dict]
        if len(unknown_elements):
            raise TypeError("Unknown element(s): '{:s}'".format(str(unknown_elements)))

        series = pd.Series(list_of_elements)

        mask = series.str.match(re_dict[types[0]], case=False)
        for ty in types[1:]:
            mask = mask | series.str.match(re_dict[ty], case=False)
        return mask.values

    @classmethod
    def get_file(cls, filename):
        return os.path.join(CURRENT_DIR, "psbooster", filename)


class _PsboosterSegmentMixin(object):

   def __init__(self):
       self._start = None
       self._end = None


    # Private Methods ##########################################################


def _get_file_for_ring(ring):
    return os.path.join(PSB_DIR, "twiss_ring" + str(ring) + ".dat")

def _get_ring_from_seqname(seq):
    ring = None
    
    if (seq.upper() == "PSB1"):
        ring = 1
    if (seq.upper() == "PSB2"):
        ring = 2
    if (seq.upper() == "PSB3"):
        ring = 3
    if (seq.upper() == "PSB4"):
        ring = 4
    
    if (ring == None):
        LOGGER.error("Sequence name is none of the expected ones (PSB1,PSB2,PSB3,PSB4)")

    return ring


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
