"""
.. module: interaction_point

Created on 13/06/18

:author: Jaime Coello de Portugal


GetLLM.algorithms.interaction_point.py stores helper functions for phase calculations for GetLLM.
This module is not intended to be executed. It stores only functions.
"""
from numpy import sqrt, sin, cos, tan, pi
import pandas as pd
from utils import tfs_pandas

PI2 = 2 * pi
COLUMNS = ("IP", "BETASTAR", "EBETASTAR", "PHASEADV", "EPHASEADV",
           "MDLPHADV", "LSTAR")
PLANES = ("X", "Y")
MODE_TO_SUFFIX = {"D": "get_IP{}.out",
                  "F": "get_IP{}_free.out",
                  "F2": "get_IP{}_free2.out"}


def betastar_from_phase(accel, phase_d, model):
    """Writes the getIP files with the betastar computed using phase advance.

    Arguments:
        accel: The accelerator class to be used.
        phase_d: The GetLLM phase_d object, output of calculate_phase.
        mode: metaclass.twiss instance with a model of the machine.
    Returns:
        A nested dict with the same structure as the phase_d dict.
    """
    try:
        ips = accel.get_ips()
    except AttributeError:
        # TODO: Log no ips in the accelerator
        return None
    ip_dict = dict(zip(PLANES, ({}, {})))
    for plane in PLANES:
        for mode in phase_d[plane]:
            phases_df = phase_d[plane][mode]
            rows = []
            for ip_name, bpml, bpmr in ips:
                try:
                    phaseadv, ephaseadv, mdlphadv = _get_meas_phase(
                        bpml, bpmr, phases_df
                    )
                except KeyError:
                    # TODO: Log combination not in phase meas
                    continue  # Measurement on one of the BPMs not present.
                lstar = _get_lstar(bpml, bpmr, model)
                betastar, ebestar = phase_to_betastar(
                    lstar, PI2 * phaseadv, PI2 * ephaseadv
                )
                rows.append([ip_name, betastar, ebestar,
                            phaseadv, ephaseadv, mdlphadv,
                            lstar])
            ip_dict[plane][mode] = pd.DataFrame(columns=COLUMNS, data=rows)
    return ip_dict


def write_betastar_from_phase(ips_d, headers, output_dir):
    """ TODO
    """
    for plane in PLANES:
        for mode in ips_d[plane]:
            tfs_pandas.write_tfs(MODE_TO_SUFFIX[mode].format(plane.lower()),
                                 ips_d[plane][mode],
                                 headers_dict=headers,)


def phase_to_betastar(lstar, phase, errphase):
    """Return the betastar and its error given the phase advance across the IP.

    This function computes the betastar using the phase advance between the
    BPMs around the IP and their distance (lstar). The phase and error in the
    phase must be given in radians.

    Arguments:
        lstar: The distance between the BPMs and the IR.
        phase: The phase advance between the BPMs at each side of the IP. Must
            be given in radians.
        errphase: The error in phase. Must be given in radians.
    """
    return (_phase_to_betastar_value(lstar, phase),
            _phase_to_betastar_error(lstar, phase, errphase))


def _phase_to_betastar_value(l, ph):
    tph = tan(ph)
    return (l * (1 - sqrt(tph ** 2 + 1))) / tph


def _phase_to_betastar_error(l, ph, eph):
    return abs((eph * l * (abs(cos(ph)) - 1)) / (sin(ph) ** 2))


def _get_ip_tfs_writer(files_dict, filename, plane):
    tfs_writer = files_dict[filename]
    tfs_writer.add_column_names(COLUMNS)
    tfs_writer.add_column_datatypes(["%s"] +
                                    ["%le"] * (len(COLUMNS) - 1))
    tfs_writer.add_string_descriptor("PLANE", plane.upper())
    return tfs_writer


def _phase_d_combinations(phase_d, plane):
    for mode in ("D", "F", "F2"):
        try:
            yield (mode, phase_d[plane][mode]),
        except KeyError:
            continue


def _get_meas_phase(bpml, bpmr, phases_df):
    return (phases_df["MEAS"].loc[bpml, bpmr],
            phases_df["ERRMEAS"].loc[bpml, bpmr],
            phases_df["MODEL"].loc[bpml, bpmr])


def _get_lstar(bpml, bpmr, model):
    return abs(model.S[model.indx[bpml]] - model.S[model.indx[bpmr]]) / 2.
