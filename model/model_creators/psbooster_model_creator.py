import model_creator
import os
import shutil


class PsboosterModelCreator(model_creator.ModelCreator):

    @classmethod
    def get_madx_script(cls, instance, output_path):
        replace_dict = {
            "FILES_DIR": instance.get_psb_dir(),
            "RING": instance.get_ring(),
            "USE_ACD": 1 if instance.acd else 0,
            "NAT_TUNE_X": instance.nat_tune_x,
            "NAT_TUNE_Y": instance.nat_tune_y,
            "KINETICENERGY": instance.energy,
            "DPP": instance.dpp,
            "OUTPUT": output_path,
            "DRV_TUNE_X": "", "DRV_TUNE_Y": "",
        }
        if (instance.acd):
            replace_dict["DRV_TUNE_X"] = instance.drv_tune_x
            replace_dict["DRV_TUNE_Y"] = instance.drv_tune_y

        with open(instance.get_nominal_tmpl()) as textfile:
            madx_template = textfile.read()
        
        return madx_template % replace_dict

    @classmethod
    def _prepare_fullresponse(cls, instance, output_path):
        with open(instance.get_iteration_tmpl()) as textfile:
            iterate_template = textfile.read()

       
        replace_dict = {
            "FILES_DIR": instance.get_psb_dir(),
            "RING": instance.get_ring(),
            "LIB": instance.NAME, # "psbooster"
            "OPTICS_PATH": instance.optics_file,
            "PATH": output_path,
            "KINETICENERGY": instance.energy,
            "NAT_TUNE_X": instance.nat_tune_x,
            "NAT_TUNE_Y": instance.nat_tune_y,
            "DRV_TUNE_X": "", 
            "DRV_TUNE_Y": "",
            "DPP": instance.dpp,
            "OUTPUT": output_path,
        }

        with open(os.path.join(output_path,
                               "job.iterate.madx"), "w") as textfile:
            textfile.write(iterate_template % replace_dict)

    @classmethod
    def prepare_run(cls, instance, output_path):
        if instance.fullresponse:
            cls._prepare_fullresponse(instance, output_path)
        
        file_name = "error_deff_ring" + str(instance.get_ring()) +".txt"
        file_path = instance.get_psb_dir()
        src_path  = os.path.join(file_path, file_name)
        dest_path = os.path.join(output_path, "error_deff.txt")
        
        shutil.copy(src_path, dest_path)
        
        #os.link(src, dst) (file_path, link_path)
            
