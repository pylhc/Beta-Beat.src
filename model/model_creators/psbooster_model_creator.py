import model_creator


class PsboosterModelCreator(model_creator.ModelCreator):

    @classmethod
    def get_madx_script(cls, instance, output_path):
        replace_dict = {
            "FILES_DIR": instance.get_psb_dir(),
            "RING": instance.get_ring(),
            "USE_ACD": 1 if instance.acd else 0,
            "NAT_TUNE_X": instance.nat_tune_x,
            "NAT_TUNE_Y": instance.nat_tune_y,
            "OUTPUT": output_path,
            "DRV_TUNE_X": "", "DRV_TUNE_Y": "",
        }
        if (instance.acd):
            replace_dict["DRV_TUNE_X"] = instance.drv_tune_x
            replace_dict["DRV_TUNE_Y"] = instance.drv_tune_y

        with open(instance.get_nominal_tmpl()) as textfile:
            madx_template = textfile.read()
	return madx_template % replace_dict
