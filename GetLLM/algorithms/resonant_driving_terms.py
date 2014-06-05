'''
Created on May 15, 2014

@author: rwestenb

@version: 0.0.1

GetLLM.algorithms.resonant_driving_terms.py stores helper functions for RDT calculations for GetLLM.
This module is not intended to be executed. It stores only functions.
'''

import sys

import Utilities.bpm
import helper

DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

def calculate_RDTs(mad_twiss, getllm_d, twiss_d, phase_d, tune_d, files_dict, pseudo_list_x, pseudo_list_y):
    '''
    Calculates line RDT amplitudes and phases and fills the following TfsFiles:
        f3000_line.out ...

    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            lhc_phase, accel, beam_direction and num_beams_for_coupling are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances.
    '''
    print "Calculating RDTs"

    """
    The rdt_set holds all RDTs which should be investigated with the parameters to call GetRDT()
    syntax is: rdt_set = [(plane, out_file, line), ...]
    with:
        plane in ["H", "V"]
        out_file in files_dict is the out file to write the data to (must be added to GetLLM.py)
        line in (int, int) is the corresponding line to the driving term
    """
    rdt_set = [
        ("H", files_dict["f3000_line.out"], (-2, 0)), # sextupolar
        ("H", files_dict["f4000_line.out"], (-3, 0))  # sextupolar
    ]
    for rdt in range(len(rdt_set)):
        _process_RDT(mad_twiss, phase_d, twiss_d, rdt_set[rdt])

def _process_RDT(mad_twiss, phase_d, twiss_d, (plane, out_file, line)):
    assert plane in ["H", "V"] # check user input plane

    # get plane corresponding phase and twiss data
    if plane == "H":
        phase_data = phase_d.ph_x
        list_zero_dpp = twiss_d.zero_dpp_x
    else:
        phase_data = phase_d.ph_y
        list_zero_dpp = twiss_d.zero_dpp_y

    dbpms = Utilities.bpm.intersect(list_zero_dpp)
    dbpms = Utilities.bpm.model_intersect(dbpms, mad_twiss)

    # init out file
    out_file.add_column_names(["NAME", "S", "COUNT", "AMP", "EAMP", "PHASE", "EPHASE"])
    out_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le"])

    for i in range(len(dbpms)-4):
        bpm1 = dbpms[i][1].upper()
        try:
            bpm_pair_data = _get_best_fitting_bpm(phase_data, bpm1, plane)
        except KeyError:
            print >> sys.stderr, "Could not find a BPM pair (%s, %s)!\n\t" % (plane, bpm1)
            continue
        bpm2 = bpm_pair_data[0]
        for j in range(0,len(list_zero_dpp)):
            amp_line, phase_line = _line_to_amp_and_phase_attr(line, list_zero_dpp[j])
            delta, edelta = bpm_pair_data[1][:2]
            amp1 = amp_line[list_zero_dpp[j].indx[bpm1]]
            amp2 = amp_line[list_zero_dpp[j].indx[bpm2]]
            phase1 = phase_line[list_zero_dpp[j].indx[bpm1]]
            phase2 = phase_line[list_zero_dpp[j].indx[bpm2]]
            line_amp, line_phase, line_amp_e, line_phase_e = helper.ComplexSecondaryLineExtended(delta,edelta, amp1,amp2, phase1,phase2)
            out_file.add_table_row([bpm1, dbpms[i][0], len(list_zero_dpp), line_amp, line_amp_e, line_phase, line_phase_e])

def _get_best_fitting_bpm(phase_d, bpm1, plane):
    '''Returns best fitting bpm_pair+phase_d of the given phase_d in the given plane for the given bpm1.'''
    # get all possible_pairs
    possible_pairs = {}
    for bpm_pair in phase_d:
        if bpm_pair[0] == plane and bpm_pair[1:].startswith(bpm1):
            possible_pairs[bpm_pair] = phase_d[bpm_pair]
    
    if not possible_pairs: raise KeyError
    # find best_fitting bpm_pair. We want the second bpm to be as close as possible to .25 offset
    #bpm_pair = min(possible_pairs, key=lambda bpm_pair: abs(float(possible_pairs[bpm_pair][0] - .25)))
    # for now we just pick the next bpm for pairing -> 0.0 phase offset
    bpm_pair = min(possible_pairs, key=lambda bpm_pair: abs(float(possible_pairs[bpm_pair][0])))
    return (bpm_pair[1:].replace(bpm1, ""), possible_pairs[bpm_pair])

def _line_to_amp_and_phase_attr(line, zero_dpp):
    '''To turn input line (-1,2) to (zero_dpp.AMP_12, zero_dpp.PHASE_12).'''
    line = (str(line[0])+str(line[1])).replace("-", "_")
    return (getattr(zero_dpp, "AMP"+line), getattr(zero_dpp, "PHASE"+line))