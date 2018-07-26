"""
.. module: tune

Created on 17/07/18

:author: Lukas Malina

It computes betatron tunes and provides structures to store them.
"""
import numpy as np
from model.accelerators.accelerator import AccExcitationMode
from utils import stats
PLANES = ('X', 'Y')


def calculate_tunes(measure_input, input_files):
    tune_d = TuneDict()
    accelerator = measure_input.accelerator
    tune_d["X"]["QM"] = accelerator.get_model_tfs().headers["Q1"]
    tune_d["Y"]["QM"] = accelerator.get_model_tfs().headers["Q2"]
    for plane in PLANES:
        char = {"X": "1", "Y": "2"}
        tunes = np.array([df.headers["Q" + char[plane]] for df in input_files._get_zero_dpp_frames(plane)])
        tunes_rms = np.array(
            [df.headers["Q" + char[plane] + "RMS"] for df in input_files._get_zero_dpp_frames(plane)])
        tune = stats.weighted_mean(tunes, errors=tunes_rms)
        tune_d[plane]["Q"] = tune
        tune_d[plane]["QF"] = tune
        tune_d[plane]["QFM"] = accelerator.nat_tune_x if plane is "X" else accelerator.nat_tune_y

        if measure_input.accelerator.excitation != AccExcitationMode.FREE:
            tune_d[plane]["QM"] = accelerator.drv_tune_x if plane is "X" else accelerator.drv_tune_y
            tune_d[plane]["QF"] = tune_d[plane]["Q"] - tune_d[plane]["QM"] + tune_d[plane]["QFM"]
    return tune_d


class TuneDict(dict):
    """
    Used as data structure to hold tunes
    """
    def __init__(self):
        super(TuneDict, self).__init__(zip(PLANES, ({"Q": 0.0, "QF": 0.0, "QM": 0.0, "QFM": 0.0},
                                                    {"Q": 0.0, "QF": 0.0, "QM": 0.0, "QFM": 0.0})))


class _TuneData(object):

    def __init__(self, tune_dict):
        self.q1 = tune_dict["X"]["Q"]  # Driven horizontal tune
        self.q2 = tune_dict["Y"]["Q"]  # Driven vertical tune
        # Free is from analytic equation
        self.q1f = tune_dict["X"]["QF"]  # Free horizontal tune
        self.q2f = tune_dict["Y"]["QF"]  # Free vertical tune
        self.q1mdl = tune_dict["X"]["QM"]
        self.q2mdl = tune_dict["Y"]["QM"]
        self.q1mdlf = tune_dict["X"]["QFM"]
        self.q2mdlf = tune_dict["Y"]["QFM"]
