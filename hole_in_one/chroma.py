from __future__ import print_function
import logging
import numpy as np
import pandas as pd
from model import manager

LOGGER = logging.getLogger(__name__)

PI2I = 2 * np.pi * complex(0, 1)


def get_chroma(harpy_input, plane, bpm_names, bpm_samples, panda):
    arc_bpm_names, arc_bpm_samples = _get_lhc_arc_bpms(bpm_names, bpm_samples)
    tunez = np.mean(panda.set_index("NAME").loc[arc_bpm_names, "TUNEZ"])
    tunez_phase = _phase_mean(panda.set_index("NAME").loc[arc_bpm_names, "MUZ"])
    dpoverp_amp, pos, neg = get_dpoverp_amp(tunez, tunez_phase, bpm_samples)
    for bin in pos:
        mini, maxi = bin
        coefficients = panda[:, "AMP" + plane] * np.exp(PI2I * panda[:, "MU" + plane])
        frequencies = panda[:, "TUNE" + plane]
        new_signal = coefficients * np.exp(PI2I * np.outer(frequencies, np.arange(bpm_samples.shape[1])))
        # Remove synchro-betatron line?
        bpm_int = (bpm_samples - new_signal)[:, mini:maxi]
        # or subtract the synchrotron line? should be about the same
        bpm_int = (bpm_int.T-np.mean(bpm_int, axis=1)).T
        length = length + maxi - mini
    for bin in neg:
        mini, maxi = bin
        co = co - np.sum(bpm_samples[:, mini:maxi], axis=1)
        length = length + maxi - mini


def _get_positive_dpoverp_intervals(tunez, tunez_phase, turns):
    halfwidth = 0.1  # at 90 % of the maximum
    intervals = []
    start = (-halfwidth - tunez_phase) / tunez
    end = (halfwidth - tunez_phase) / tunez
    for i in range(int(turns * tunez) + 1):
        start_i = start + float(i) / tunez
        end_i = end + float(i) / tunez
        if start_i > 0 and end_i < turns:
            intervals.append((int(start_i), int(end_i)))
    return intervals


def _get_negative_dpoverp_intervals(tunez, tunez_phase, turns):
    return _get_positive_dpoverp_intervals(tunez, tunez_phase + 0.5, turns)


# # full period is one, phase is of the arc bpms in horizontal plane
def get_dpoverp_amp(model_tfs, panda, bpm_names, bpm_samples):
    # Interval is around integer phases for positive dpoverp
    # and around halfes for the negative dpoverp
    arc_bpm_names, arc_bpm_samples = _get_lhc_arc_bpms(bpm_names, bpm_samples)
    tunez = np.mean(panda.set_index("NAME").loc[arc_bpm_names, "TUNEZ"])
    tunez_phase = _phase_mean(panda.set_index("NAME").loc[arc_bpm_names, "MUZ"])
    turns = bpm_samples.shape[1]
    pos = _get_positive_dpoverp_intervals(tunez, tunez_phase, turns)
    neg = _get_negative_dpoverp_intervals(tunez, tunez_phase, turns)
    # We assume 3D kicks only happen in LHC for now.
    co = np.zeros_like(arc_bpm_names)
    model_dx = model_tfs.set_index("NAME").loc[arc_bpm_names, "DX"]
    length = 0
    for bin in pos:
        mini, maxi = bin
        co = co + np.sum(arc_bpm_samples[:, mini:maxi], axis=1)
        length = length + maxi - mini
    for bin in neg:
        mini, maxi = bin
        co = co - np.sum(arc_bpm_samples[:, mini:maxi], axis=1)
        length = length + maxi - mini
    co = co / 1e3  # Going from mm to m.
    codx = co * model_dx / length
    dx2 = model_dx ** 2
    return np.sum(codx) / np.sum(dx2), pos, neg


def _get_lhc_arc_bpms(bpm_names, bpm_samples):
    accel_cls = manager.get_accel_class("lhc")
    arc_bpms_mask = accel_cls.get_arc_bpms_mask(bpm_names)
    arc_bpm_samples = bpm_samples[arc_bpms_mask]
    arc_bpm_names = bpm_names[arc_bpms_mask]
    return arc_bpm_names, arc_bpm_samples


def _phase_mean(phases):
    return np.angle(np.sum(np.exp(PI2I * phases))) / (2 * np.pi)

#if plane == "X":
            #dpoverp_amp, pos, neg = get_dpoverp_amp(
             #   model_tfs, panda, bpm_names, bpm_matrix,
            #)

