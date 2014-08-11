'''
Created on Jul 18, 2014

This test class runs a tracking simulation for each file in "data/tracking_jobs", then it runs the Drive + GetLLM
algorithms and compares the beta outputs.

The maximum allowed relative error can be set in the constant MAX_BETA_REL_ERR.

@author: jcoellod
'''
import unittest
import os
from Python_Classes4MAD import madxrunner, metaclass
from Utilities import iotools
import sys
import subprocess
from numpy.ma.core import abs
import re

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
MAX_BETA_REL_ERR = 0.0004


class Test(unittest.TestCase):

    def setUp(self):
        self._set_up_paths()
        self._prepare_tracking_scripts()

    def tearDown(self):
        for output_path in self._output_paths:
            iotools.delete_item(output_path)

    def testName(self):
        for madx_file_path, output_path in zip(self._madx_file_paths, self._output_paths):
            print "Testing madx job " + os.path.basename(madx_file_path)
            self._run_tracking_script(madx_file_path)
            self._run_drive(output_path, madx_file_path)
            self._run_getLLM(output_path)
            self._compare_output(output_path)
            print "Done\n"

    def _set_up_paths(self):
        self._data_path = os.path.join(CURRENT_PATH, "data")
        self._tracking_jobs_path = os.path.join(self._data_path, "tracking_jobs")
        self._output_paths = []
        self._madx_template_paths = []
        self._madx_file_paths = []
        for file_name in iotools.get_all_filenames_in_dir(self._tracking_jobs_path):
            output_path = os.path.join(self._data_path, "output_" + file_name)
            iotools.create_dirs(output_path)
            self._output_paths.append(output_path)
            self._madx_template_paths.append(os.path.join(self._tracking_jobs_path, file_name))
            self._madx_file_paths.append(os.path.join(output_path, file_name))

        if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "linux3":
            self._path_to_drive = os.path.join(CURRENT_PATH, "..", "..", "drive", "Drive_God_lin")
        elif sys.platform == "win32":
            self._path_to_drive = os.path.join(CURRENT_PATH, "..", "..", "drive", "Drive_God_lin_win.exe")
        else:
            raise OSError("No drive version for given os: " + sys.platform)

        self._path_to_getllm = os.path.join(CURRENT_PATH, "..", "..", "GetLLM", "GetLLM.py")

    def _prepare_tracking_scripts(self):
        for template, output_job, output_dir in zip(self._madx_template_paths, self._madx_file_paths, self._output_paths):
            dict_to_replace = {"%DATA_PATH": self._data_path, "%OUTPUT_PATH": output_dir}
            self._replace_keywords_in_file(os.path.join(self._data_path, "parse_track_file.sh"),
                                           os.path.join(output_dir, "parse_track_file.sh"),
                                           dict_to_replace)
            self._replace_keywords_in_file(template, output_job, dict_to_replace)

    def _run_tracking_script(self, madx_file_path):
        print "Running tracking code..."
        errcode = madxrunner.runForInputFile(madx_file_path, stdout=subprocess.PIPE)
        self.assertEqual(errcode, 0, "Error running MADX tracking code.")

    def _run_drive(self, output_path, madx_file_path):
        print "Running drive..."
        self._prepare_drive_input(output_path, madx_file_path)
        self._replace_keywords_in_file(os.path.join(self._data_path, "DrivingTerms"),
                                       os.path.join(output_path, "DrivingTerms"),
                                       {"%OUTPUT_PATH": output_path})
        call_command = os.path.abspath(self._path_to_drive) + " " + os.path.abspath(output_path)

        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)

        (out_stream, err_stream) = process.communicate()
        errcode = process.returncode
        if 0 != errcode:
            print "Error running drive:", call_command
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream
        self.assertEqual(errcode, 0, "Error running drive.")

    def _prepare_drive_input(self, output_path, madx_file_path):
        tune_x, tune_y = 0.0, 0.0
        with open(madx_file_path, "r") as madx_file:
            for line in madx_file:
                if "mux=" in line and "muy=" in line:
                    tunes = re.findall("\d+.\d+", line)
                    tune_x = float(tunes[0]) % 1
                    tune_y = float(tunes[1]) % 1
        self.assertTrue(tune_x != 0.0 and tune_y != 0.0, "Cannot find tunes in " + madx_file_path)
        print "Tune x: " + str(tune_x) + ", Tune y: " + str(tune_y)
        self._replace_keywords_in_file(os.path.join(self._data_path, "Drive.inp"),
                                       os.path.join(output_path, "Drive.inp"),
                                       {"%TUNE_X": str(tune_x), "%TUNE_Y": str(tune_y)})

    def _run_getLLM(self, output_path):
        print "Running GetLLM..."
        call_command = sys.executable + " " + os.path.abspath(self._path_to_getllm) + \
        " --accel=LHCB1 --tbtana=SUSSIX --bpmu=mm " + \
        " -m " + os.path.join(output_path, "twiss.dat") + \
        " -f " + os.path.join(output_path, "ALLBPMs") + \
        " -o " + output_path

        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)

        (out_stream, err_stream) = process.communicate()
        errcode = process.returncode
        if 0 != errcode:
            print "Error running getLLM:", call_command
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream
        self.assertEqual(errcode, 0, "Error running getLLM.")

    def _compare_output(self, output_path):
        print "Comparing output..."
        beta_x_twiss = metaclass.twiss(os.path.join(output_path, "getbetax.out"))
        beta_y_twiss = metaclass.twiss(os.path.join(output_path, "getbetay.out"))

        for index in range(len(beta_x_twiss.NAME)):
            rel_error = abs((beta_x_twiss.BETX[index] - beta_x_twiss.BETXMDL[index]) / beta_x_twiss.BETXMDL[index])
            self.assertTrue(rel_error < MAX_BETA_REL_ERR,
                            "Relative error too big found in: " + beta_x_twiss.NAME[index] + " (" + str(rel_error) + ")")
        for index in range(len(beta_y_twiss.NAME)):
            rel_error = abs((beta_y_twiss.BETY[index] - beta_y_twiss.BETYMDL[index]) / beta_y_twiss.BETYMDL[index])
            self.assertTrue(rel_error < MAX_BETA_REL_ERR,
                            "Relative error too big found in: " + beta_y_twiss.NAME[index] + " (" + str(rel_error) + ")")

    def _replace_keywords_in_file(self, input_file, output_file, dict_to_replace):
        with open(output_file, "w") as output_data:
            with open(input_file, "r") as input_data:
                for line in input_data:
                    new_line = line
                    for key, value in dict_to_replace.iteritems():
                        new_line = new_line.replace(key, value)
                    output_data.write(new_line)

if __name__ == "__main__":
    unittest.main()
