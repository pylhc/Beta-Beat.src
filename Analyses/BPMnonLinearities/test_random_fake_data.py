import unittest
import os

class Test(unittest.TestCase):
    '''Test class to generate and validate random BPM data with all old corrections.
    
    TODO: Implement BPMS random data generation by using polynomial_correction.rotate_skewed_coordinates(x, y).
    '''

    def _create_fake_test_data(self):
        import textwrap, time, random
        import polynomial_correction

        self.sdds_header = textwrap.dedent('''\
            #SDDSASCIIFORMAT v1
            #Beam: %s
            #Created: %s By: Python unittest (NL correction)
            #bunchid :0
            #number of turns :%d
            #number of monitors :%d\n''' % (
                self.beam_name,
                time.strftime('%Y-%m-%d#%H-%M-%S', time.localtime()),
                self.number_of_turns,
                self.number_of_bpms
            )
        )

        # first create raw data and write to sdds file
        with open(self.valid_raw_path, 'w') as sdds_raw:
            # write header
            sdds_raw.write(self.sdds_header)

            # write data
            for bpm_num in range(self.number_of_bpms):
                random_bpm_position = random.random() * 27e3 # about 27km
                random_bpm_type = random.choice(self.bpm_types)
                bpm_name = '%s.%d.%s' % (random_bpm_type, bpm_num, self.beam_name)
                self.bpm_positions[bpm_name] = random_bpm_position
                for plane in ('X', 'Y'):
                    if plane == 'X': plane_num = 0
                    else: plane_num = 1
                    if bpm_name not in self.raw_data[plane]: self.raw_data[plane][bpm_name] = []

                    sdds_raw.write('%d %s      %.5f  ' % (plane_num, bpm_name, random_bpm_position))
                    for turn_num in range(self.number_of_turns):
                        random_amp = (random.random()*2 - 1)*self.max_raw_amp
                        self.raw_data[plane][bpm_name].append(random_amp)

                        format_string = '%.5f  '
                        if turn_num == self.number_of_turns - 1: format_string = '%.5f' # to avoid '  ' at the end of each line
                        sdds_raw.write(format_string % random_amp) # random number in (-self.max_raw_amp, +self.max_raw_amp)
                    sdds_raw.write('\n')

        # now we create old correction file with this data
        with open(self.valid_old_correction_path, 'w') as sdds_old:
            # write header
            sdds_old.write(self.sdds_header)

            # write data
            for plane in self.raw_data:
                if plane == 'X': plane_num = 0
                else: plane_num = 1

                for bpm_name in self.raw_data[plane]:
                    bpm_type = bpm_name.split('.')[0]
                    bpm_position = self.bpm_positions[bpm_name]
                    sdds_old.write('%d %s      %.5f  ' % (plane_num, bpm_name, bpm_position))
                    for turn_num in range(len(self.raw_data[plane][bpm_name])):
                        coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.old_coefficients[self.year_for_old_correction])
                        turn_amp = self.raw_data[plane][bpm_name][turn_num]
                        turn_amp *= coefficients[0]
                        turn_amp = polynomial_correction._old_P5(turn_amp, coefficients)

                        format_string = '%.5f  '
                        if turn_num == len(self.raw_data[plane][bpm_name]) - 1: format_string = '%.5f' # to avoid '  ' at the end of each line
                        sdds_old.write(format_string % turn_amp) # u_old = P5_old(kf * u_raw)
                    sdds_old.write('\n')

        # new correction file with this data
        with open(self.valid_new_correction_path, 'w') as sdds_new:
            # write header
            sdds_new.write(self.sdds_header)

            # write data
            for plane in self.raw_data:
                if plane == 'X': plane_num = 0
                else: plane_num = 1

                for bpm_name in self.raw_data[plane]:
                    bpm_type = bpm_name.split('.')[0]
                    bpm_position = self.bpm_positions[bpm_name]
                    sdds_new.write('%d %s      %.5f  ' % (plane_num, bpm_name, bpm_position))
                    for turn_num in range(len(self.raw_data[plane][bpm_name])):
                        coefficients = polynomial_correction.get_coefficients_of_bpm_type(bpm_type, polynomial_correction.new_coefficients)
                        turn_amp = self.raw_data[plane][bpm_name][turn_num]
                        if plane == 'X': v_raw = self.raw_data['Y'][bpm_name][turn_num]
                        else: v_raw = self.raw_data['X'][bpm_name][turn_num]
                        turn_amp = polynomial_correction.new_P5(self.raw_data[plane][bpm_name][turn_num], v_raw, coefficients)

                        format_string = '%.5f  '
                        if turn_num == len(self.raw_data[plane][bpm_name]) - 1: format_string = '%.5f' # to avoid '  ' at the end of each line
                        sdds_new.write(format_string % turn_amp) # u_new = P5_new(u_raw)
                    sdds_new.write('\n')
    
    def setUp(self):
        self.test_dir = './random_generated_test_data/'
        self.valid_dir = os.path.join(self.test_dir, 'valid/')
        self.to_test_dir = os.path.join(self.test_dir, 'to_test/')
        self.valid_raw_path = os.path.join(self.valid_dir, 'raw.sdds')
        self.valid_old_correction_path = os.path.join(self.valid_dir, 'old_correction.sdds')
        self.valid_new_correction_path = os.path.join(self.valid_dir, 'new_correction.sdds')
        self.raw_data = {'X':{}, 'Y':{}}
        self.bpm_positions = {}

        self.year_for_old_correction = 2012
        self.beam_name = 'LHCBtest'
        self.number_of_turns = 10
        self.number_of_bpms = len(self.bpm_types)*2
        self.max_raw_amp = 0.03 # maximum raw amplitude (.03 from 8sig8sig test)
        self.eps = 1.1e-5 # last digit is 1e-5
        if not os.path.exists(self.valid_dir): os.makedirs(self.valid_dir)
        if not os.path.exists(self.to_test_dir): os.makedirs(self.to_test_dir)

        self._create_fake_test_data()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.test_dir)

    def test_non_linearty_correction(self):
        import correct_non_linearities, SDDS_file

        self.to_test_new_correction_path = os.path.join(self.to_test_dir, 'new_correction.sdds')
        correct_non_linearities.correct_non_linearities(self.valid_old_correction_path, self.to_test_new_correction_path, self.year_for_old_correction)

        to_test_sdds = SDDS_file.SDDS_file(self.to_test_new_correction_path)
        to_test_sdds.read_in()
        valid_sdds = SDDS_file.SDDS_file(self.valid_new_correction_path)
        valid_sdds.read_in()

        to_test_bpm_data = to_test_sdds.get_bpm_data()
        valid_bpm_data = valid_sdds.get_bpm_data()
        for plane in to_test_bpm_data:
            for bpm in to_test_bpm_data[plane]:
                for turn_num in range(len(to_test_bpm_data[plane][bpm])):
                    try:
                        assert abs(to_test_bpm_data[plane][bpm][turn_num] - valid_bpm_data[plane][bpm][turn_num]) <= self.eps
                    except AssertionError:
                        raise AssertionError('Test failed for %s (%s): %.5f <-- valid | to test --> %.5f; diff: %.5f (eps: %.5f)' % (
                            bpm,
                            plane,
                            valid_bpm_data[plane][bpm][turn_num],
                            to_test_bpm_data[plane][bpm][turn_num],
                            abs(to_test_bpm_data[plane][bpm][turn_num] - valid_bpm_data[plane][bpm][turn_num]),
                            self.eps
                        ))

if __name__ == "__main__":
    unittest.main()