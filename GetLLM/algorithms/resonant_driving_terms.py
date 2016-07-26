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


RDT_LIST = ['f1001H', 'f3000H', 'f1002H', 
            'f4000H', 'f2010H', 'f1220H', 
            'f1020H', 'f1120H', 'f2003H', 
            'f2020H', 'f3001H', 'f1101H', 
            'f1200H', 
            'f0110V', 'f1020V', 'f1011V', 
            'f0111V', 'f2010V', 'f0120V', 
            'f0030V', 'f0040V', 'f0211V', 
            'f1012V', 'f1013V'
            ]


def determine_lines(rdt):
    r = list(rdt)
    j, k, l, m, plane = int(r[1]), int(r[2]), int(r[3]), int(r[4]), r[5]
    if plane == 'H':
        line = (1-j+k, m-l)
    elif plane == 'V':
        line = (k-j, 1-l+m)
    return line, plane


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

    for rdt in RDT_LIST:
        line, plane = determine_lines(rdt)
        _process_RDT(mad_twiss, phase_d, twiss_d, (plane, files_dict[rdt+'_line.out'], line))


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
            delta, edelta = bpm_pair_data[1:]
            amp1 = amp_line[list_zero_dpp[j].indx[bpm1]]
            amp2 = amp_line[list_zero_dpp[j].indx[bpm2]]
            phase1 = phase_line[list_zero_dpp[j].indx[bpm1]]
            phase2 = phase_line[list_zero_dpp[j].indx[bpm2]]
            line_amp, line_phase, line_amp_e, line_phase_e = helper.ComplexSecondaryLineExtended(delta,edelta, amp1,amp2, phase1,phase2)
            out_file.add_table_row([bpm1, dbpms[i][0], len(list_zero_dpp), line_amp, line_amp_e, line_phase, line_phase_e])


def _get_best_fitting_bpm(phase_d, bpm1, plane):
    ''' 
    @author: F Carlier
    phase_d is a dictionary containing two different sets of keys. 
    First there are the keys for the bpm pairs in the H & V planes, for example: 
    HBPM.33R3.B1BPM.32R3.B1 for which the values are [phase difference measurement, error phase difference, phase difference model]

    Second keys are the individual bpms ["BPM.33R3.B1","BPM.32R3.B1", etc..], with their values given by:
    [phase diff X, std phase diff X, phase diff Y, std phase diff Y, phase diff Model X, phase diff model Y, NEXT BPM]   

    Using the second key one iterates over the list of bpms to find the pair for which the phase difference is closest to 90 degrees (.25)
    The phase advances between a pair of bpms is always compared to the advance of the previous set. As soon as the absolute difference to .25 
    rises the last bpm is used for the pair. This iterations makes sure that the integer phase advance is taken into account. 
    
    plane_idx       : makes sure the right plane is used for the phase
    ph_advance      : total phase advance between first bpm and target bpm
    ph_advance_next : total phase advance between first bpm and next bpm
    value           : absolute difference between phase advance and .25 of first bpm set
    value_next      : absolute difference between phase advance and .25 of the next bpm set
    
    '''

    if plane == 'H':
        plane_idx = 0
    elif plane == 'V':
        plane_idx = 2
    else: 
        raise KeyError("No valid plane was found!")
    
    ph_advance = float(phase_d[bpm1][plane_idx])  
    ph_adv_err = float(phase_d[bpm1][plane_idx+1])  
    next_bpm = phase_d[bpm1][6]
    target_bpm = next_bpm
    ph_advance_next = ph_advance + float(phase_d[next_bpm][plane_idx])
    ph_advance_next_err = float(phase_d[next_bpm][plane_idx+1])
    value = abs(ph_advance - .25) 
    value_next = abs(ph_advance_next - .25)

    while value_next < value:
        target_bpm = next_bpm
        ph_adv_err = ph_advance**2*ph_adv_err + (ph_advance_next*ph_advance_next_err)**2 
        ph_advance = ph_advance_next
        next_bpm = phase_d[target_bpm][6]
        ph_advance_next += float(phase_d[next_bpm][plane_idx])
        ph_advance_next_err = float(phase_d[next_bpm][plane_idx+1])
        value = value_next
        value_next = abs(ph_advance_next - .25)
        
    return next_bpm , ph_advance, ph_adv_err


def _line_to_amp_and_phase_attr(line, zero_dpp):
    '''To turn input line (-1,2) to (zero_dpp.AMP_12, zero_dpp.PHASE_12).'''
    line = (str(line[0])+str(line[1])).replace("-", "_")
    return (getattr(zero_dpp, "AMP"+line), getattr(zero_dpp, "PHASE"+line))
