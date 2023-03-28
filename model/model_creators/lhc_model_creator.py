from __future__ import print_function

import logging
import os
import shutil
import sys

from model.accelerators.accelerator import AccExcitationMode

import model_creator

AFS_ROOT = "/afs"
if "win" in sys.platform and sys.platform != "darwin":
    AFS_ROOT = "\\AFS"

LOGGER = logging.getLogger(__name__)

ERROR_TABLES_ROOT = "/afs/cern.ch/eng/sl/lintrack/error_tables"

class LhcModelCreator(model_creator.ModelCreator):
    ERR_DEF_PATH = os.path.join(AFS_ROOT, "cern.ch", "eng", "sl", "lintrack",
                                "omc_repositories", "omc3", "omc3", "model", 
                                "accelerators", "lhc", "systematic_errors")
   
    ERR_DEF_FILES = {
        "0.45": "0450GeV.tfs",
        "1.0": "1000GeV.tfs",
        "1.5": "1500GeV.tfs",
        "2.0": "2000GeV.tfs",
        "2.5": "2500GeV.tfs",
        "3.0": "3000GeV.tfs",
        "3.5": "3500GeV.tfs",
        "4.0": "4000GeV.tfs",
        "4.5": "4500GeV.tfs",
        "5.0": "5000GeV.tfs",
        "5.5": "5500GeV.tfs",
        "6.0": "6000GeV.tfs",
        "6.5": "6500GeV.tfs",
        "6.8": "6800GeV.tfs",
    }

    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        """Returns the MAD-X script.

        In case of 2022 afs is assumed at the moment and a symblink is created.
        If not this has to be handled by the user.
        """
        use_acd = "1" if (lhc_instance.excitation == AccExcitationMode.ACD) else "0"
        use_adt = "1" if (lhc_instance.excitation == AccExcitationMode.ADT) else "0"
        crossing_on = "1" if lhc_instance.xing else "0"
        beam = lhc_instance.get_beam()
        if lhc_instance.YEAR in ["2022", "2023"]:
            src = "/afs/cern.ch/eng/acc-models/lhc/" + lhc_instance.YEAR
            dst = output_path + "/acc-models-lhc"
            if os.path.exists(dst):
                LOGGER.warn("model exists, overwriting")
                os.unlink(dst)
            os.symlink(src, dst)

        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": beam,
            "PATH": output_path,
            "QMX": lhc_instance.nat_tune_x,
            "QMY": lhc_instance.nat_tune_y,
            "USE_ACD": use_acd,
            "USE_ADT": use_adt,
            "DPP": lhc_instance.dpp,
            "CROSSING_ON": crossing_on,
            "QX": "",
            "QY": "",
            "QDX": "",
            "QDY": "",
        }
        if lhc_instance.excitation in (AccExcitationMode.ACD, AccExcitationMode.ADT):
            replace_dict["QX"] = lhc_instance.nat_tune_x
            replace_dict["QY"] = lhc_instance.nat_tune_y
            replace_dict["QDX"] = lhc_instance.drv_tune_x
            replace_dict["QDY"] = lhc_instance.drv_tune_y

        with open(lhc_instance.get_nominal_tmpl()) as textfile:
            madx_template = textfile.read()

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
            try:
                os.unlink(link_path)
            except OSError:
                pass
            os.symlink(file_path, link_path)
            # for omc3
            error_deffs_path = os.path.join(output_path, "error_deffs.txt")
            try:
                os.unlink(error_deffs_path)
            except OSError:
                pass

            if os.path.exists(file_path):
                shutil.copy(file_path, error_deffs_path)
            else:
                os.symlink(file_path, error_deffs_path)

    @classmethod
    def _prepare_fullresponse(cls, lhc_instance, output_path):
        with open(lhc_instance.get_iteration_tmpl()) as textfile:
            iterate_template = textfile.read()

        crossing_on = "1" if lhc_instance.xing else "0"
        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "QMX": lhc_instance.nat_tune_x,
            "QMY": lhc_instance.nat_tune_y,
            "CROSSING_ON": crossing_on,
        }

        with open(os.path.join(output_path, "job.iterate.madx"), "w") as textfile:
            textfile.write(iterate_template % replace_dict)


class LhcBestKnowledgeCreator(LhcModelCreator):
    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        if lhc_instance.excitation is not AccExcitationMode.FREE:
            raise model_creator.ModelCreationError(
                "Don't set ACD or ADT for best knowledge model."
            )
        if lhc_instance.energy is None:
            raise model_creator.ModelCreationError(
                "Best knowledge model requires energy."
            )
        with open(lhc_instance.get_best_knowledge_tmpl()) as textfile:
            madx_template = textfile.read()
        crossing_on = "1" if lhc_instance.xing else "0"

        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "DPP": lhc_instance.dpp,
            "QMX": lhc_instance.nat_tune_x,
            "QMY": lhc_instance.nat_tune_y,
            "ENERGY": lhc_instance.energy,
            "CROSSING_ON": crossing_on,
        }
        if lhc_instance.error_table is not None:
            replace_dict["ERROR_TABLE"] = """
readmytable, file = "/afs/cern.ch/eng/sl/lintrack/error_tables/Beam{}/{}.errors", table=errtab;
seterr, table=errtab;
call, file = "/afs/cern.ch/eng/sl/lintrack/error_tables/Beam{}/{}.madx";
""".format(lhc_instance.get_beam(), lhc_instance.error_table,
           lhc_instance.get_beam(), lhc_instance.error_table)
        else:
            replace_dict["ERROR_TABLE"] = ""

        madx_script = madx_template % replace_dict
        return madx_script


class LhcSegmentCreator(model_creator.ModelCreator):
    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        with open(lhc_instance.get_segment_tmpl()) as textfile:
            madx_template = textfile.read()
        if lhc_instance.YEAR in ["2022", "2023"] and not os.path.exists(
            output_path + "/acc-models-lhc"
        ):
            os.symlink(
                "/afs/cern.ch/eng/acc-models/lhc/" + lhc_instance.YEAR, output_path + "/acc-models-lhc"
            )
        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "LABEL": lhc_instance.label,
            "BETAKIND": lhc_instance.kind,
            "STARTFROM": lhc_instance.start.name,
            "ENDAT": lhc_instance.end.name,
        }
        madx_script = madx_template % replace_dict
        return madx_script


class LhcCouplingCreator(model_creator.ModelCreator):
    @classmethod
    def get_madx_script(cls, lhc_instance, output_path):
        with open(lhc_instance.get_coupling_tmpl()) as textfile:
            madx_template = textfile.read()
            print(madx_template)
        crossing_on = "1" if lhc_instance.xing else "0"
        replace_dict = {
            "LIB": lhc_instance.MACROS_NAME,
            "MAIN_SEQ": lhc_instance.load_main_seq_madx(),
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "QMX": lhc_instance.nat_tune_x,
            "QMY": lhc_instance.nat_tune_y,
            "CROSSING_ON": crossing_on,
        }
        madx_script = madx_template % replace_dict
        return madx_script
