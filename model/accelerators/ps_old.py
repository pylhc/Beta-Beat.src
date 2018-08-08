from __future__ import print_function
import os
import argparse
from tfs_files import tfs_pandas
from accelerator import Accelerator, AcceleratorDefinitionError, AccExcitationMode

CURRENT_DIR = os.path.dirname(__file__)
#LHC_DIR = os.path.join(CURRENT_DIR, "lhc")



class Cps(Accelerator):
    NAME = "CPS"
    MACROS_NAME = "CPS"

    def __init__(self):
        self.optics_file = None
        self.nat_tune_x = None
        self.nat_tune_y = None
        self.energy = None
        self.dpp = 0.0
        self.xing = None
        
        # for GetLLM
        self.model = None
        self.model_driven = None
        self.model_best_knowledge = None
        self.elements = None
        self.elements_centre = None
        self.excitation = AccExcitationMode.FREE
        self.modelpath = None

    @classmethod
    def init_from_args(cls, args):
        parser = cls._get_arg_parser()
        options, rest_args = parser.parse_known_args(args)
        instance = cls()
        instance.nat_tune_x = options.nat_tune_x
        instance.nat_tune_y = options.nat_tune_y
        
        
        instance.dpp = options.dpp
        instance.energy = options.energy
        instance.optics_file = options.optics
        instance.fullresponse = options.fullresponse
        instance.verify_object()
        
        return instance, rest_args
        
    @classmethod
    def init_from_model_dir(cls, model_dir):  # prints only for debugging
        
        print("Creating accelerator instance from model dir")
        instance = cls()
        
        instance.modelpath = model_dir

        if os.path.isfile(model_dir):
            model_dir = os.path.dirname(model_dir)
        instance.model_tfs = tfs_pandas.read_tfs(os.path.join(model_dir, "twiss.dat")).set_index("NAME")
        print("model path =", os.path.join(model_dir, "twiss.dat"))
            
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
        
        instance.nat_tune_x = float(instance.model_tfs.headers["Q1"])
        instance.nat_tune_y = float(instance.model_tfs.headers["Q2"])
        
        return instance
    
    @classmethod
    def get_class(cls):
        return Ps()
    
    @classmethod
    def get_class_from_args(cls, args):
        parser = argparse.ArgumentParser()
       
        options, rest_args = parser.parse_known_args(args)
        return cls, rest_args

  
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
        if self.optics_file is None:
            raise AcceleratorDefinitionError( 
                "The accelerator definition is incomplete, optics "
                "file has not been specified."
            )
            
    def get_arc_bpms_mask(cls, list_of_elements):
        return [True] * len(list_of_elements)

    
    
    def get_s_first_BPM(self):
        return 0

    def get_errordefspath(self):
        return os.path.join(self.modelpath, "errordefs")
    
    def get_k_first_BPM(self, list_of_bpms):
        return len(list_of_bpms)
        
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

