import model_creator
import os

class PsModelCreator(model_creator.ModelCreator):

    @classmethod
    def get_madx_script(cls, instance, output_path):
        replace_dict = {
            "FILES_DIR": instance.get_ps_dir(),
            "USE_ACD": 1 if instance.acd else 0,
            "NAT_TUNE_X": instance.nat_tune_x,
            "NAT_TUNE_Y": instance.nat_tune_y,
            "KINETICENERGY": instance.energy,
            "DPP": instance.dpp,
            "OUTPUT": output_path,
            "DRV_TUNE_X": "", 
            "DRV_TUNE_Y": "",
            "OPTICS_PATH": instance.optics_file,
        }
        
        if (instance.acd):
            replace_dict["DRV_TUNE_X"] = instance.drv_tune_x
            replace_dict["DRV_TUNE_Y"] = instance.drv_tune_y

        with open(instance.get_nominal_tmpl()) as textfile:
            madx_template = textfile.read()
        
        #print(replace_dict)
        #print(madx_template)
        
        out = madx_template % replace_dict
        
        #print(out)
        
        return out
    
    @classmethod
    def _prepare_fullresponse(cls, instance, output_path):
        with open(instance.get_iteration_tmpl()) as textfile:
            iterate_template = textfile.read()

       
        replace_dict = {
            "FILES_DIR": instance.get_ps_dir(),
            "LIB": instance.MACROS_NAME,
            "OPTICS_PATH": instance.optics_file,
            "PATH": output_path,
            "KINETICENERGY": instance.energy,
            "NAT_TUNE_X": instance.nat_tune_x,
            "NAT_TUNE_Y": instance.nat_tune_y,
            "DRV_TUNE_X": "", 
            "DRV_TUNE_Y": "",
        }

        with open(os.path.join(output_path,
                               "job.iterate.madx"), "w") as textfile:
            textfile.write(iterate_template % replace_dict)
