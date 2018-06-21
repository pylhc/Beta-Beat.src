import pandas as pd

from accelerator import Accelerator


class Esrf(Accelerator):

    NAME = "esrf"

    @classmethod
    def get_arc_bpms_mask(cls, list_of_elements):
        """ Chooses the bpms with large dispersion.
        Which are:
            bpms 1-5 in even cells.
            bpms 3-7 in odd cells. """
        return pd.Series(list_of_elements).str.match(
            r"BPM\.(\d*[02468]\.[1-5]|\d*[13579]\.[3-7])",
            case=False).values

    @classmethod
    def get_class(cls):
        new_class = cls
        return new_class
