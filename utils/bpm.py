'''
.. module: utils.bpm

Created on 3 Jun 2013

This module contains helper functions concerning bpms. It contains functions to filter BPMs or to
intersect BPMs in multiple files or with a given model file.

.. moduleauthor:: gvanbavi, vimaier
'''

import sys
import pandas as pd

def filterbpm(list_of_bpms):
    '''Filter non-arc BPM.
       :returns: list -- a list with those bpms which start with name "BPM."
    '''
    return list_of_bpms.loc[list_of_bpms.index.str.match("^BPM\\.")]

def model_intersect(exp_bpms, model_twiss):
    '''
    Intersects BPMs from

    :param list exp_bpms: list with tuples (<S_value_i>,<bpm_i>)
    :param metaclass.Twiss model_twiss: Should be a Twiss object from model

    :returns: list with tuples: (<S_value_i>,<bpm_i>) -- A list with BPMs which are both in exp_bpms and model_twiss.
    '''

    common_index = exp_bpms.index.intersection(model_twiss.index)

    return exp_bpms.loc[common_index]

def intersect(list_of_twiss_files):
    '''
    Pure intersection of all bpm names in all files.

    :param list list_of_twiss_files: List of metaclass.Twiss objects with columns NAME and S.

    :returns: list with tuples: (<S_value_i>,<bpm_i>) -- bpm_i is in every twiss of list_of_twiss_files.
    '''
    if len(list_of_twiss_files) == 0:
        print >> sys.stderr, "Nothing to intersect!!!!"
        return []

    common_index = list_of_twiss_files[0].index
    for i in range(1, len(list_of_twiss_files)):
        common_index = common_index.intersection(list_of_twiss_files[i].index)

    commbpms = pd.DataFrame(list_of_twiss_files[0].loc[common_index, "S"])
    commbpms["NFILES"] = len(list_of_twiss_files)
    return commbpms


def intersect_with_bpm_list(exp_bpms, bpm_list):
    '''
    Intersects BPMs from

    :param list exp_bpms' list with tuples: (<S_value_i>,<bpm_i>)
    :param list bpm_list: List of bpm names

    :returns: list with tuples: (<S_value_i>,<bpm_i>) -- A list with BPMs which are both in exp_bpms and bpm_list.
    '''
    result = []

    for s_bpm_tupel in exp_bpms:
        bpm_name = s_bpm_tupel[1]
        if bpm_name in bpm_list:
            result.append(s_bpm_tupel)
    return result
