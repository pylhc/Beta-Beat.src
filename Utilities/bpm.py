'''
Created on 3 Jun 2013

@author: vimaier

@version: 0.0.1

This module contains helper functions concerning bpms.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys


def modelIntersect(exp_bpms, model_twiss):
    '''
    Intersects BPMs from
    
    :Parameters:
        'exp_bpms': list with tuples: (<S_value_i>,<bpm_i>)
            list with tuples: (<S_value_i>,<bpm_i>).
        'model_twiss': metaclass.Twiss
            Should be a Twiss object from model
                
    :Return: list with tuples: (<S_value_i>,<bpm_i>)
        A list with BPMs which are both in exp_bpms and model_twiss.
        
    :Exception: SystemExit
        If the length of exp_bpms is 0 or the length of the resulting list of strings is 0 then 
        SystemExit will be raised (sys.exit(1) ),
    '''
    bpmsin = []
    #print "start Intersect, exp_bpms #:", len(exp_bpms)
    if len(exp_bpms) == 0:
        print >> sys.stderr, "Zero exp BPMs sent to modelIntersect"
        sys.exit(1)
    for bpm in exp_bpms:
        try:
            check_if_bpm_in_model = model_twiss.indx[bpm[1].upper()]  # @UnusedVariable
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
    
    :Parameters:
        'list_of_twiss_files': list
            List of metaclass.Twiss objects with columns NAME and S.
                
    :Return: list with tuples: (<S_value_i>,<bpm_i>)
        A list with following tuples: (<S_value_i>,<bpm_i>).
        bpm_i is in every twiss of list_of_twiss_files.
        
    :Exception: SystemExit
        If the length of list_of_twiss_files is 0 or there are no entries in column 'NAME' then 
        SystemExit will be raised (sys.exit(1) ),
    '''
    if len(list_of_twiss_files) == 0:
        print >> sys.stderr, "Nothing to intersect!!!!"
        sys.exit(1)
        
    z = list_of_twiss_files[0].NAME
    if len(z) == 0:
        print >> sys.stderr, "No exp BPMs..."
        sys.exit(1)
    for b in list_of_twiss_files:
        z=filter(lambda x: x in z   , b.NAME)
     
        
    result = []
    x0 = list_of_twiss_files[0]
    for bpm in z:
        result.append((x0.S[x0.indx[bpm]], bpm))

    #SORT by S
    result.sort()
    return result



