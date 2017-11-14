import re
import numpy as np
from accelerator import Accelerator


class Esrf(Accelerator):

    NAME = "esrf"

    @classmethod
    def get_arc_bpms_mask(cls, list_of_elements):
        mask = []
        pattern = re.compile("BPM\.([0-9]+)\.([1-7])", re.IGNORECASE)
        for element in list_of_elements:
            match = pattern.match(element)
            # The arc bpms are from BPM.14... and up
            if match:
                cell = int(match.group(1))
                bpm_number = int(match.group(2))
                mask.append(not ((cell % 2 == 0 and bpm_number in [6, 7]) or
                    (cell % 2 != 0 and bpm_number in [1, 2])))
            else:
                mask.append(False)
        return np.array(mask)

    @classmethod
    def get_class(cls):
        new_class = cls
        return new_class
