'''
Created on 27 May 2013

@author: ?, vimaier

@version: 0.0.1

GetLLM.algorithms.interaction_point.py stores helper functions for phase calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''
from numpy import sqrt, sin, cos, tan, pi
from utils import logging_tools
PI2 = 2 * pi
COLUMNS = ("IP", "BETASTAR", "EBETASTAR", "PHASEADV", "EPHASEADV",
           "MDLPHADV", "LSTAR")
PLANES = ("x", "y")
HVPLANE = {"x": "H", "y": "V"}

LOGGER = logging_tools.get_logger(__name__)

def betastar_from_phase(accel, phase_d, model, files_dict):
    """Writes the getIP files with the betastar computed using phase advance.

    Arguments:
        accel: A string with the accelerator name, for now only LHCB1 or LHCB2.
        phase_d: The GetLLM phase_d object, output of calculate_phase.
        mode: metaclass.twiss instance with a model of the machine.
        files_dict: The GetLLM files_dict, storing the tfs_file_writers to use
            to write the results.
    """
    
    LOGGER.debug("Accelerator is: %s",accel)
    if "LHC" not in accel.upper():
        LOGGER.warn("Machine %s is not of LHC type, there is no IP's to have beta star.", accel)
        #print("Machine %s is not of LHC type, there is no IP's to have beta star"%accel)
        return
    
    for plane in PLANES:
        for filename, phases in _phase_d_combinations(phase_d, plane):
            tfs_writer = _get_ip_tfs_writer(files_dict, filename, plane)
            for ip_name, bpml, bpmr in _get_ips(accel):
                try:
                    phaseadv, ephaseadv, mdlphadv = _get_meas_phase(
                        bpml, bpmr, plane, phases
                    )
                except KeyError:
                    continue  # Measurement on one of the BPMs not present.
                lstar = _get_lstar(bpml, bpmr, model)
                betastar, ebestar = phase_to_betastar(
                    lstar, PI2 * phaseadv, PI2 * ephaseadv
                )
                tfs_writer.add_table_row([ip_name, betastar, ebestar,
                                          phaseadv, ephaseadv, mdlphadv,
                                          lstar])


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
    return (("getIP{}.out".format(plane),
             getattr(phase_d, "ph_{}".format(plane))),
            ("getIP{}_free.out".format(plane),
             getattr(phase_d, "{}_f".format(plane))),
            ("getIP{}_free2.out".format(plane),
             getattr(phase_d, "{}_f2".format(plane))))


LHC_IPS = ("1", "2", "5", "8")
NORMAL_IP_BPMS = "BPMSW.1{side}{ip}.B{beam}"
DOROS_IP_BPMS = "LHC.BPM.1{side}{ip}.B{beam}_DOROS"


# TODO: This should go in the accelerator class when available
def _get_ips(accel):
    if "LHC" not in accel:
        raise StopIteration
    beam = 1 if "B1" in accel else 2
    for ip in LHC_IPS:
        yield ("IP{}".format(ip),
               NORMAL_IP_BPMS.format(side="L", ip=ip, beam=beam),
               NORMAL_IP_BPMS.format(side="R", ip=ip, beam=beam))
        yield ("IP{}_DOROS".format(ip),
               DOROS_IP_BPMS.format(side="L", ip=ip, beam=beam),
               DOROS_IP_BPMS.format(side="R", ip=ip, beam=beam))


def _get_meas_phase(bpml, bpmr, plane, phases):
    combination_keys = (HVPLANE[plane] + bpml + bpmr,
                        HVPLANE[plane] + bpmr + bpml)
    for key in combination_keys:
        try:
            phadv, ephadvm, mdlphadv = phases[key]
            return abs(phadv), abs(ephadvm), abs(mdlphadv)
        except KeyError:
            continue
    raise KeyError("{0} {1} combination not in phases dict".format(bpml, bpmr))


def _get_lstar(bpml, bpmr, model):
    return abs(model.S[model.indx[bpml]] - model.S[model.indx[bpmr]]) / 2.
