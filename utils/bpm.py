'''
.. module: utils.bpm

Created on 3 Jun 2013

This module contains helper functions concerning bpms. It contains functions to filter BPMs or to
intersect BPMs in multiple files or with a given model file.

.. moduleauthor:: gvanbavi, vimaier
'''

import sys

def filterbpm(list_of_bpms):
    '''Filter non-arc BPM.
       :returns: list -- a list with those bpms which start with name "BPM."
    '''
    result = []
    if len(list_of_bpms) == 0:
        print >> sys.stderr, "Nothing to filter!!!!"
        return result
    for b in list_of_bpms:
        if ('BPM.' in b[1].upper()):
            result.append(b)
    return result


def model_intersect(exp_bpms, model_twiss):
    '''
    Intersects BPMs from

    :param list exp_bpms: list with tuples (<S_value_i>,<bpm_i>)
    :param metaclass.Twiss model_twiss: Should be a Twiss object from model

    :returns: list with tuples: (<S_value_i>,<bpm_i>) -- A list with BPMs which are both in exp_bpms and model_twiss.
    '''
    return model_twiss.loc[exp_bpms.index, "S"]


def intersect(list_of_twiss_files):
    '''
    Pure intersection of all bpm names in all files.

    :param list list_of_twiss_files: List of metaclass.Twiss objects with columns NAME and S.

    :returns: list with tuples: (<S_value_i>,<bpm_i>) -- bpm_i is in every twiss of list_of_twiss_files.
    '''
    if len(list_of_twiss_files) == 0:
        print >> sys.stderr, "Nothing to intersect!!!!"
        return []

    names_list = list_of_twiss_files[0].index
    if len(names_list) == 0:
        print >> sys.stderr, "No exp BPMs..."
        sys.exit(1)
    for twiss_file in list_of_twiss_files:
        names_list = twiss_file.index.intersection(names_list)

    result = list_of_twiss_files[0].loc[names_list, "S"]
    return result

def get_list_of_tuples(bpms):
    """transforms the DataFrame bpms to a list of tuples to fake the old usage.
    """
    return [(bpms.loc[name], name) for name in bpms.index]


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
