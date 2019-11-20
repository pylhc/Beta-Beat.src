import sys
import math
import copy
import numpy as np
from numpy import sin, cos, tan
from collections import OrderedDict
import pandas as pd
import utils.bpm
import phase
from constants import PI, TWOPI, kEPSILON
from tfs_files import tfs_pandas
import numpy as np
import matplotlib.pyplot as plt
DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#def get_amplitude_calibration_and_beta(Files, plane, bd,all_bpms):
#    amp,calibration,calibration_error,psid = get_amplitude_calibration_and_phase(Files,all_bpms,plane,bd)
#    return amp,calibration,calibration_error,psid


def get_amplitude_calibration_and_phase(Files,plane,bd,all_bpms):
    amp = []
    calibration = []
    calibration_error = []
    psid = []
    for i in range(len(Files)):
        if plane == 'H':
            amp.append(np.array(
                    [2 * Files[i].AMPX[Files[i].indx[b]] for b in all_bpms]
            ))
            try:
                 calibration.append(np.array(
                        [Files[i].CALIBRATION[Files[i].indx[b]] for b in all_bpms]
                ))
                 calibration_error.append(np.array(
                        [Files[i].ERROR_CALIBRATION[Files[i].indx[b]] for b in all_bpms]
                ))
            except AttributeError:
                  calibration.append(np.ones(len(np.array([2 * Files[i].AMPX[Files[i].indx[b]] for b in all_bpms]))))
                  calibration_error.append(np.zeros(len(np.array([2 * Files[i].AMPX[Files[i].indx[b]] for b in all_bpms]))))
            psid.append(bd * 2 * np.pi * np.array([Files[i].MUX[Files[i].indx[b]] for b in all_bpms]))
 #           print "all bpms"
 #           print all_bpms
        if plane == 'V':
 #           print "all bpms"
 #           print len(all_bpms)
            amp.append(np.array(
                    [2 * Files[i].AMPY[Files[i].indx[b]] for b in all_bpms]
            ))
            try:
                calibration.append(np.array(
                        [Files[i].CALIBRATION[Files[i].indx[b]] for b in all_bpms]
                ))
                calibration_error.append(np.array(
                        [Files[i].ERROR_CALIBRATION[Files[i].indx[b]] for b in all_bpms]
                ))
            except AttributeError:
                  calibration.append(np.ones(len(np.array([2 * Files[i].AMPY[Files[i].indx[b]] for b in all_bpms]))))
                  calibration_error.append(np.zeros(len(np.array([2 * Files[i].AMPY[Files[i].indx[b]] for b in all_bpms]))))
            psid.append(bd * 2 * np.pi * np.array([Files[i].MUY[Files[i].indx[b]] for b in all_bpms]))

    calibration = np.array(calibration)
    calibration_error = np.array(calibration_error)
    amp = np.array(amp)
    amp_transpose = amp.transpose()
    calibration_transpose = calibration_error.transpose()
    calibration_error_transpose = calibration_error.transpose()
    psid_array = np.array(psid)
    psid_transpose = psid_array.transpose()
#    print "amplitude"
#    print amp
    return amp,calibration,calibration_error,psid_transpose


