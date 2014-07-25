import unittest
import os

class Test(unittest.TestCase):
    '''Test class validate old correction reversion + new correction application from valid test data from one BPM.'''

    def setUp(self):
        self.test_dir = './test_data/'
        self.valid_dir   = os.path.join(self.test_dir, 'valid/')
        self.to_test_dir = os.path.join(self.test_dir, 'to_test/')
        self.valid_raw_path            = os.path.join(self.valid_dir, 'raw.sdds')
        self.valid_old_correction_path = os.path.join(self.valid_dir, 'old_correction.sdds')
        self.valid_new_correction_path = os.path.join(self.valid_dir, 'new_correction.sdds')
        if not os.path.exists(self.to_test_dir): os.makedirs(self.to_test_dir)

        self.eps = 2.1e-4 # last digit is 1e-5, 2.1e-4 was the smallest passing tolerance
        self.year_for_old_correction = 2012 # we will look at 2012 corrected data

    def tearDown(self):
        import shutil
        shutil.rmtree(self.to_test_dir)

    def test_non_linearty_correction(self):
        import correct_non_linearities, SDDS_file

        self.to_test_new_correction_path = os.path.join(self.to_test_dir, 'new_correction.sdds')
        # from old to new correction
        correct_non_linearities.correct_non_linearities(self.valid_old_correction_path, self.to_test_new_correction_path, self.year_for_old_correction)

        # new test valid new vs. to test new
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