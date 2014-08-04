"""This module is capable of converting data from the old correction into data with the new 2015 correction polynomial.

This implementation follows the procedure described in the note:
"Geometrical non-linearity correction procedure of LHC beam position monitors", A. Nosych (BE/BI) (CERN EDMS Id:1342295)
"""
import copy

from SDDS_file import SDDS_file
import polynomial_correction

def _read_in_sdds_file(sdds_file_path):
    sdds_file = SDDS_file(sdds_file_path)
    sdds_file.read_in()
    return sdds_file

def _revert_old_correction(data, year):
    """Reverts the old correction of the given year to end up with uncorrected data."""
    for plane in data:
        for bpm in data[plane]:
            for i in range(len(data[plane][bpm])):
                bpm_type = bpm.split('.')[0]
                if bpm_type == 'BPMS': continue # we will skip skewed BPMS type BPMs as they need special treatment and are only a few
                
                coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.old_coefficients[year])
                if coefficients:
                    # the normalization by kf = coefficients[0] is important
                    data[plane][bpm][i] = polynomial_correction.revert_correction_polynomial(data[plane][bpm][i], coefficients) / coefficients[0]
                else:
                    print('Warning: No coefficients for bpm type: "%s" found! Skipping...' % bpm_type)
                    continue
    return data

def _apply_new_correction(data):
    """Applies the new 2015 correction to the data.
    
    Arguments:
    data -- Raw data without any correction applied (old correction should be already reverted at this point).
    
    Returns:
    return_data -- Corrected data using 2015 2D correction polynomial where possible. Where only data of one plane is available 1D fallback is used.
    """
    return_data = copy.deepcopy(data) # need a deepcopy of data here because we want to modify values but still need the old ones
    for plane in data:
        for bpm in data[plane]:
            bpm_type = bpm.split('.')[0]
            if bpm_type == 'BPMS': continue # we will skip skewed BPMS type BPMs as they need special treatment and are only a few
            
            use_1d_fallback = True if bpm not in data['X'] or bpm not in data['Y'] else False
            if use_1d_fallback:
                print('Warning: No corresponding BPM for %s found, using 1D fallback correction!' % bpm)

            for i in range(len(data[plane][bpm])):
                if not use_1d_fallback:
                    if plane == 'X':
                        u_raw = data['X'][bpm][i]
                        v_raw = data['Y'][bpm][i]
                    else: # bit confusing but we just want to replace x by y and vice versa for the correction or the other plane
                        u_raw = data['Y'][bpm][i]
                        v_raw = data['X'][bpm][i]
                    coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.new_coefficients)
                else:
                    u_raw = data[plane][bpm][i]
                    v_raw = 0
                    coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.new_coefficients_1D_fallback)

                return_data[plane][bpm][i] = polynomial_correction.new_P5(u_raw, v_raw, coefficients)
    return return_data

def _apply_old_correction(data, year):
    """Applies old correction to raw data. Normally not needed was implemented for testing the _revert_old_correction() function."""
    for plane in data:
        for bpm in data[plane]:
            for i in range(len(data[plane][bpm])):
                bpm_type = bpm.split('.')[0]
                u_raw = data[plane][bpm][i]
                coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.old_coefficients[year])

                data[plane][bpm][i] = polynomial_correction._old_P5(u_raw, coefficients)
    return data

def _output_new_data(corrected_bpm_data, sdds_file, new_sdds_file_path):
    """Writes the corrected_bpm_data to new_sdds_file_path.
    Old sdds_file is used to get the header and meta information for the new_sdds_file. Note that sdds_file is an SDDS_file class instance.
    """
    new_sdds_file = SDDS_file(new_sdds_file_path)
    new_sdds_file.file_header = sdds_file.file_header
    new_sdds_file.bpm_positions = sdds_file.bpm_positions
    new_sdds_file.bpm_data = corrected_bpm_data
    new_sdds_file.write_out()

def correct_non_linearities(in_sdds_file_path, out_sdds_file_path, year=2012, undo_old_correction=True, verbose=False):
    """Main correction algorithm reverting the old applied correction and applies the new one, then outputs the data.
    
    Arguments:
    in_sdds_file_path -- Complete path to the input sdds file.
    out_sdds_file_path -- Complete path to the ouput sdds file.
    year -- Year of the input data (2011, 2012 or 2013) (default 2012)
    
    Returns:
    None
    
    Note:
    It is assumed that the old data was already corrected by 2011, 2012 or 2013 correction polynomial. BPMS type BPMs are skipped (see _revert_old_correction() or _apply_new_correction()).
    """
    if verbose: print('Reading in: %s' % in_sdds_file_path)
    sdds_file = _read_in_sdds_file(in_sdds_file_path)
    data = sdds_file.get_bpm_data()
    
    if verbose: print('Got data. Reverting old correction...')
    
    if undo_old_correction:
        data = _revert_old_correction(data, year)

    if verbose: print('Reverted. Applying new correction...')
    
    data = _apply_new_correction(data)
    
    if verbose: print('Done. Writing out new data to: %s' % out_sdds_file_path)

    _output_new_data(data, sdds_file, out_sdds_file_path)
    print('Done!')

def parse_args():
    import sys
    from optparse import OptionParser

    usage = 'Usage: %prog -o OUT_FILE_PATH [options] INPUT_SDDS_FILE_PATH'
    parser = OptionParser(usage=usage)
    parser.add_option('-o', '--out-file-path', help='path of the resulting sdds out file', action='store', type='string', dest='out_sdds_file_path')
    parser.add_option('-y', '--year', help='measurement year of the input sdds file (for determining the old correction parameters) [default: %default]', action='store', type='int', dest='year', default=2012)
    parser.add_option('-n', '--not-undo-old-correction', help='do not undo old correction before applying new one', action='store_false', dest='undo_old_correction', default=True)
    parser.add_option('-v', '--verbose', help='extends output', action='store_true', dest='verbose')
    
    (options, args) = parser.parse_args()
    if not args:
        parser.print_help()
        print('Error: No input file given!')
        sys.exit(1)
    
    if len(args) > 1:
        parser.print_help()
        print('Error: Multiple input files given, expected only one!')
        sys.exit(1)
    
    if not options.out_sdds_file_path:
        options.out_sdds_file_path = args[0].replace('.sdds', '.nl_corrected.sdds')
        print('No output file path with -o given, using: %s' % options.out_sdds_file_path)
    
    return (options, args)

if __name__ == '__main__':
    options, args = parse_args()
    '''
    options.in_sdds_file_path = '/afs/cern.ch/user/r/rwestenb/workarea/11_12_injectionMD/measurementsB2inj_diagonal/analysis/manual_data_preparation/experimental/Corrected_8sig8sig_Beam2@Turn@2012_06_25@04_00_22_346_0.sdds.new'
    options.out_sdds_file_path = in_sdds_file_path.replace('.sdds', '.nl_corrected.sdds')
    options.year = 2012
    options.undo_old_correction = True
    options.verbose = True
    '''
    correct_non_linearities(args[0], options.out_sdds_file_path, year=options.year, undo_old_correction=options.undo_old_correction, verbose=options.verbose)