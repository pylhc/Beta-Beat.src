'''
.. module: Utilities.tfs_file

Created on 20 Aug 2013

This module contains some functions for TfsFiles

.. moduleauthor:: vimaier
'''
import numpy

import Python_Classes4MAD.metaclass
import Utilities.math
import Utilities.compare

def check_tunes_for_tfs_file(path_to_tfs_file):
    """
    Recalculates tunes(Q1, Q2, Q1RMS, Q2RMS, NATQ1, NATQ2, NATQ1RMS, NATQ2RMS) from columns(TUNEX,
    TUNEY, NATTUNEX, NATTUNEY) and compares with saved tunes.
    """
    print "Checking tune values for "+path_to_tfs_file
    tw = Python_Classes4MAD.metaclass.twiss(path_to_tfs_file)

    column = getattr(tw, "TUNEX", None)
    if not column is None:
        _compute_values(tw, column, "Q1")

    column = getattr(tw, "TUNEY", None)
    if not column is None:
        _compute_values(tw, column, "Q2")

    column = getattr(tw, "NATTUNEX", None)
    if not column is None:
        column = _clean_from_default_values(column)
        _compute_values(tw, column, "NATQ1")

    column = getattr(tw, "NATTUNEY", None)
    if not column is None:
        column = _clean_from_default_values(column)
        _compute_values(tw, column, "NATQ2")



def _compare_and_print(stored, computed, kind_of_value):
    if stored is None:
        print "\t"+kind_of_value+", attribute is not in twiss file. Calculated value: "+computed

    if not Utilities.compare.almost_equal_double(stored, computed):
        print "\t"+kind_of_value+" seems to be wrong! stored != computed: "+ str(stored) +" != "+ str(computed)

def _compute_values(tw, column, tune_str):
    stored_tune = getattr(tw, tune_str, None)
    stored_tune_rms = getattr(tw, tune_str+"RMS", None)
    computed_tune = Utilities.math.arithmetic_mean(column)
    # The values with 'RMS' are actually the standard deviation of the column and not the root
    # square mean
#     computed_tune_rms = Utilities.math.root_mean_square(column)
    computed_tune_rms = Utilities.math.standard_deviation(column)
    _compare_and_print(stored_tune, computed_tune, tune_str)
    _compare_and_print(stored_tune_rms, computed_tune_rms, tune_str+"RMS")

def _clean_from_default_values(column):
    """ Nattune column has sometimes default values(-100) which should not be included in calculation """
    default_values_indices = []
    for i in xrange(len(column)):
        if not column[i] > -100: # This comparison is enough
            default_values_indices.append(i)
    return numpy.delete(column, default_values_indices)


if __name__=="__main__":
    column = [-100, 1,2,3,-100,4,5,-100,6,-100]
    column = _clean_from_default_values(column)
    print column



