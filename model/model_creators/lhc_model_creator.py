from __future__ import print_function
import os
import sys
import logging
from accelerators.lhc import LhcExcitationMode
import model_creator

AFS_ROOT = "/afs"
if "win" in sys.platform and sys.platform != "darwin":
    AFS_ROOT = "\\AFS"

LOGGER = logging.getLogger(__name__)


class LhcModelCreator(model_creator.ModelCreator):
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
        with open(lhc_instance.get_nominal_tmpl()) as textfile:
            madx_template = textfile.read()
        iqx, iqy = cls._get_full_tunes(lhc_instance)
        use_acd = "1" if (lhc_instance.excitation ==
                          LhcExcitationMode.ACD) else "0"
        use_adt = "1" if (lhc_instance.excitation ==
                          LhcExcitationMode.ADT) else "0"
        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
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
        if lhc_instance.fullresponse:
            cls._prepare_fullresponse(lhc_instance, output_path)
        if lhc_instance.energy is not None:
            file_name = cls.ERR_DEF_FILES[str(lhc_instance.energy)]
            file_path = os.path.join(cls.ERR_DEF_PATH, file_name)
            # TODO: Windows?
            link_path = os.path.join(output_path, "error_deff.txt")
            if os.path.isfile(link_path):
                os.unlink(link_path)
            os.symlink(file_path, link_path)

    @classmethod
    def _prepare_fullresponse(cls, lhc_instance, output_path):
        with open(lhc_instance.get_file("fullresponse.madx")) as textfile:
            fullresponse_template = textfile.read()
        iqx, iqy = cls._get_full_tunes(lhc_instance)
        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "QMX": iqx,
            "QMY": iqy,
        }
        fullresponse_script = fullresponse_template % replace_dict
        with open(os.path.join(output_path,
                               "job.iterate.madx"), "w") as textfile:
            textfile.write(fullresponse_script)

    @classmethod
    def _get_full_tunes(cls, lhc_instance):
        iqx, iqy = (lhc_instance.INT_TUNE_X +
                    lhc_instance.nat_tune_x,
                    lhc_instance.INT_TUNE_Y +
                    lhc_instance.nat_tune_y)
        return iqx, iqy


class LhcBestKnowledgeCreator(LhcModelCreator):

    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        if lhc_instance.excitation is not LhcExcitationMode.FREE:
            raise model_creator.ModelCreationError(
                "Don't set ACD or ADT for best knowledge model."
            )
        if lhc_instance.energy is None:
            raise model_creator.ModelCreationError(
                "Best knowledge model requires energy."
            )
        with open(lhc_instance.get_best_knowledge_tmpl()) as textfile:
            madx_template = textfile.read()
        iqx, iqy = cls._get_full_tunes(lhc_instance)
        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "DPP": lhc_instance.dpp,
            "QMX": iqx,
            "QMY": iqy,
            "ENERGY": lhc_instance.energy,
        }
        madx_script = madx_template % replace_dict
        return madx_script
