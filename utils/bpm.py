'''
.. module: Utilities.bpm

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
    bpmsin = []
    #print "start Intersect, exp_bpms #:", len(exp_bpms)
    if len(exp_bpms) == 0:
        print >> sys.stderr, "Zero exp BPMs sent to model_intersect"
        return bpmsin

    for bpm in exp_bpms:
        try:
            model_twiss.indx[bpm[1].upper()]  # Check if bpm is in the model
            bpmsin.append(bpm)
        except KeyError:
            print >> sys.stderr, bpm, "Not in Model"

    if len(bpmsin) == 0:
        print >> sys.stderr, "Zero intersection of Exp and Model"
        print >> sys.stderr, "Please, provide a good Dictionary or correct data"
        print >> sys.stderr, "Now we better leave!"
        sys.exit(1)

    return bpmsin


def intersect(list_of_twiss_files):
    '''
    Pure intersection of all bpm names in all files.

    :param list list_of_twiss_files: List of metaclass.Twiss objects with columns NAME and S.

    :returns: list with tuples: (<S_value_i>,<bpm_i>) -- bpm_i is in every twiss of list_of_twiss_files.
    '''
    if len(list_of_twiss_files) == 0:
        print >> sys.stderr, "Nothing to intersect!!!!"
        return []

    names_list = list_of_twiss_files[0].NAME
    if len(names_list) == 0:
        print >> sys.stderr, "No exp BPMs..."
        sys.exit(1)
    for twiss_file in list_of_twiss_files:
        #TODO: have to use a set probably, does not detect duplicates! (vimaier)
        names_list = [b for b in twiss_file.NAME if b in names_list]

    result = [] # list of tupels (S, bpm_name)
    twiss_0 = list_of_twiss_files[0]
    for bpm in names_list:
        result.append((twiss_0.S[twiss_0.indx[bpm]], bpm))

    #SORT by S
    result.sort()
    return result


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
