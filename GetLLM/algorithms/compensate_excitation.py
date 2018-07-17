'''
Created on 27 May 2013

@author: ?

@version: 0.0.1

GetLLM.algorithms.compensate_ac_effect.py stores helper functions to compensate the AC dipole effect
based on analytic formulae (by R. Miyamoto).
This module is not intended to be executed. It stores only functions for GetLLM.


Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import numpy as np
from numpy import sin, cos, tan
from constants import PI, TWOPI, kEPSILON
from utils import logging_tools, stats

LOGGER = logging_tools.get_logger(__name__)

DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


def GetACPhase_AC2BPMAC(commonbpms, Qd, Q, plane, acc):
    """Returns the necessary values for the exciter compensation.

    Args:
        model (pandas.DataFrame): model twiss as DataFrame.
        commonbpms (pandas.DataFrame): commonbpms (see GetLLM._get_commonbpms)
        Qd, Q (float): Driven and natural fractional tunes.
        plane (char): H, V
        accelerator: accelerator class instance.

    Returns tupel(a,b,c,d):
        a (string): name of the nearest BPM.
        b (float): compensated phase advance between the exciter and the nearest BPM.
        c (int): k of the nearest BPM.
        d (string): name of the exciter element.
    """
    bd = acc.get_beam_direction()
    r = sin(PI * (Qd - Q)) / sin(PI * (Qd + Q))
    plane_mu = "MUX" if plane == "H" else "MUY"
    [k, bpmac1], exciter = acc.get_exciter_bpm(plane, commonbpms)
    model_driven = acc.get_driven_tfs()
    ##return bpmac1, np.arctan((1 + r) / (1 - r) * tan(TWOPI * model.loc[bpmac1, plane_mu] + PI * Q)) % PI - PI * Qd, k
    try:
        psi = (
            np.arctan((1 + r) / (1 - r) *
                      tan(TWOPI * (model_driven.loc[bpmac1, plane_mu] - model_driven.loc[exciter, plane_mu]) )
                     ) / TWOPI
        ) % .5 - .5
    except:
        psi = acc.get_elements_tfs().loc[bpmac1, plane_mu] - acc.get_elements_tfs().loc[exciter, plane_mu]
    return bpmac1, psi, k, exciter


def get_kick_from_bpm_list_w_ACdipole(MADTwiss_ac, bpm_list, measurements, plane):
    '''
    @author: F Carlier
    Function calculates kick from measurements with AC dipole using the amplitude of the main line. The main line
    amplitude is obtained from Drive/SUSSIX and is normalized with the model beta-function.

    Input:
        bpm_list:     Can be any list of bpms. Preferably only arc bpms for kick calculations, but other bad bpms may be
                      included as well.
        measurements: List of measurements when analyzing multiple measurements at once
        plane:        Either H or V
    Output:
        actions_sqrt:       is a list containing the actions of each measurement. Notice this is the square root of the
                            action, so sqrt(2JX) or sqrt(2JY) depending on the plane
        actions_sqrt_err:   is the list containing the errors for sqrt(2Jx/y) for each measurement.
    '''
    if plane == 'H':
        betmdl = MADTwiss_ac.loc[bpm_list.index, "BETX"].values
    if plane == 'V':
        betmdl = MADTwiss_ac.loc[bpm_list.index, "BETY"].values

    actions_sqrt = []
    actions_sqrt_err = []

    for meas in measurements:
        if plane == 'H':
            amp = 2 * meas.loc[bpm_list.index, "AMPX"].values
        if plane == 'V':
            amp = 2 * meas.loc[bpm_list.index, "AMPY"].values
        actions_sqrt.append(np.average(amp / np.sqrt(betmdl)))
        actions_sqrt_err.append(np.std(amp / np.sqrt(betmdl)))

    return actions_sqrt, actions_sqrt_err
