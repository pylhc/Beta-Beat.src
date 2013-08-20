'''
Created on 20 Aug 2013

@author: vimaier
'''
import metaclass
import Utilities.math
import Utilities.compare

def check_tunes_for_tfs_file(path_to_tfs_file):
    """
    Recalculates tunes(Q1, Q2, Q1RMS, Q2RMS, NATQ1, NATQ2, NATQ1RMS, NATQ2RMS) from columns(TUNEX,
    TUNEY, NATTUNEX, NATTUNEY) and compares with saved tunes. 
    """
    print "Checking tune values for "+path_to_tfs_file
    tw = metaclass.twiss(path_to_tfs_file)
    
    column = getattr(tw, "TUNEX", None)
    if not column is None:
        stored_q1 = getattr(tw, "Q1", None)
        stored_q1rms = getattr(tw, "Q1RMS", None)
        computed_q1 = Utilities.math.arithmetic_mean(column)
        computed_q1rms = Utilities.math.root_mean_square(column)
        _compare_and_print(stored_q1, computed_q1, "Q1")
        _compare_and_print(stored_q1rms, computed_q1rms, "Q1RMS")
    
    column = getattr(tw, "TUNEY", None)
    if not column is None:
        stored_q2 = getattr(tw, "Q2", None)
        stored_q2rms = getattr(tw, "Q2RMS", None)
        computed_q2 = Utilities.math.arithmetic_mean(column)
        computed_q2rms = Utilities.math.root_mean_square(column)
        _compare_and_print(stored_q2, computed_q2, "Q2")
        _compare_and_print(stored_q2rms, computed_q2rms, "Q2RMS")
    
    column = getattr(tw, "NATTUNEX", None)
    if not column is None:
        stored_natq1 = getattr(tw, "NATQ1", None)
        stored_natq1rms = getattr(tw, "NATQ1RMS", None)
        computed_natq1 = Utilities.math.arithmetic_mean(column)
        computed_natq1rms = Utilities.math.root_mean_square(column)
        _compare_and_print(stored_natq1, computed_natq1, "NATQ1")
        _compare_and_print(stored_natq1rms, computed_natq1rms, "NATQ1RMS")
    
    column = getattr(tw, "NATTUNEY", None)
    if not column is None:
        stored_natq2 = getattr(tw, "NATQ2", None)
        stored_natq2rms = getattr(tw, "NATQ2RMS", None)
        computed_natq2 = Utilities.math.arithmetic_mean(column)
        computed_natq2rms = Utilities.math.root_mean_square(column)
        _compare_and_print(stored_natq2, computed_natq2, "NATQ2")
        _compare_and_print(stored_natq2rms, computed_natq2rms, "NATQ2RMS")
    
      
        
def _compare_and_print(stored, computed, kind_of_value):
    if stored is None:
        print "\t"+kind_of_value+", attribute is not in twiss file. Calculated value: "+computed
        
    if not Utilities.compare.almost_equal_double(stored, computed):
        print "\t"+kind_of_value+" seems to be wrong! stored != computed: "+ str(stored) +" != "+ str(computed)
        
        
        