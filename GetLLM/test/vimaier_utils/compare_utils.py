'''
Created on 18 Jul 2013

@author: vimaier

GetLLM.test.vimaier_utils.compare_utils holds functions for comparing TfsFiles.

'''
import Python_Classes4MAD.metaclass as metaclass
import Utilities.ndiff

#===================================================================================================
# helper-functions
#===================================================================================================
def compare_tfs_files(name_valid, name_to_check):
    """ Compares both files. Whitespace does not matter.
        Returns True in success.
    """
    if name_valid.endswith(".gitignore"):
        return True

    if __file_not_valid(name_valid):
        return True

    return Utilities.ndiff.compare_tfs_files_and_ignore_header(name_valid, name_to_check )


def __file_not_valid(name_twiss):
    ''' Checks if the given twiss file is empty, thus not valid. '''
    try:
        tw_f = metaclass.twiss(name_twiss)
        return tw_f.has_no_bpm_data()
    except ValueError:
        # Probably empty file
        return True
