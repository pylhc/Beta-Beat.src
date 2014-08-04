from SDDS_file import SDDS_file
import polynomial_correction

DEBUG = True

def _read_in_sdds_file(sdds_file_path):
    sdds_file = SDDS_file(sdds_file_path)
    sdds_file.read_in()
    return sdds_file

def _output_new_data(corrected_bpm_data, sdds_file, new_sdds_file_path):
    '''Old sdds_file is used to get the header and meta information for the new_sdds_file.'''
    new_sdds_file = SDDS_file(new_sdds_file_path)
    new_sdds_file.file_header = sdds_file.file_header
    new_sdds_file.bpm_positions = sdds_file.bpm_positions
    new_sdds_file.bpm_data = corrected_bpm_data
    new_sdds_file.write_out()

def _test_print(data):
    bpm = 'BPM.33R3.B2'
    start_turn = 100
    if DEBUG:
        for i in range(3):
            print('%s (turn %d): (%.5f, %.5f)' % (bpm, start_turn+i, data['X'][bpm][start_turn+i], data['Y'][bpm][start_turn+i]))

def _to_bpm_data(data, year):
    for plane in data:
        for bpm in data[plane]:
            for i in range(len(data[plane][bpm])):
                bpm_type = bpm.split('.')[0]
                
                coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.old_coefficients[year])

                if coefficients:
                    data[plane][bpm][i] = polynomial_correction.revert_correction_polynomial(data[plane][bpm][i], coefficients)
                else:
                    print('Warning: No coefficients for bpm type: "%s" found! Skipping...' % bpm_type)
    return data

def from_ideal_to_bpm_data(sdds_file_path, new_sdds_file_path, year=2012):
    '''Year is needed in order to find the right coefficients for the non linear correction.'''
    if DEBUG: print('Reading in: %s' % sdds_file_path)
    sdds_file = _read_in_sdds_file(sdds_file_path)
    data = sdds_file.get_bpm_data()
    
    if DEBUG: print('Got data. Applying reverse correction...')
    _test_print(data)

    data = _to_bpm_data(data, year)
    print('Now we have raw BPM data.')

    _test_print(data)
    if DEBUG: print('Writing out new data to: %s' % new_sdds_file_path)

    _output_new_data(data, sdds_file, new_sdds_file_path)
    print('Done!')

if __name__ == '__main__':
    '''We get model data and transform it into raw bpm data by undoing the correction.'''
    sdds_file_path = '/afs/cern.ch/user/r/rwestenb/workarea/11_12_injectionMD/modelB2inj_correctedoptics/analysis/simulated_data/model.sdds'
    new_sdds_file_path = '/afs/cern.ch/user/r/rwestenb/workarea/11_12_injectionMD/modelB2inj_correctedoptics/analysis/bpm_data/model.bpm_data.sdds'

    from_ideal_to_bpm_data(sdds_file_path, new_sdds_file_path, year=2012)