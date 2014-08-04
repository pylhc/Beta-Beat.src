'''
Created on Jul 18, 2014

@author: jcoellod
'''
import unittest
import os
from Python_Classes4MAD import madxrunner, metaclass
from Utilities import iotools
import sys
import subprocess
from numpy.ma.core import abs

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
MAX_BETA_REL_ERR = 0.006


class Test(unittest.TestCase):

    def setUp(self):
        self._set_up_paths()
        self._prepare_tracking_script()

    def tearDown(self):
        iotools.delete_content_of_dir(self._output_path)
        pass

    def testName(self):
        self._run_tracking_script()
        self._run_drive()
        self._run_getLLM()
        self._compare_output()

    def _set_up_paths(self):
        self._data_path = os.path.join(CURRENT_PATH, "data")
        self._output_path = os.path.join(self._data_path, "output")
        self._madx_template_path = os.path.join(self._data_path, "job.tracking.madx")
        self._madx_file_path = os.path.join(self._output_path, "job.tracking.done.madx")

        if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "linux3":
            self._path_to_drive = os.path.join(CURRENT_PATH, "..", "..", "drive", "Drive_God_lin")
        elif sys.platform == "win32":
            self._path_to_drive = os.path.join(CURRENT_PATH, "..", "..", "drive", "Drive_God_lin_win.exe")
        else:
            raise OSError("No drive version for given os: " + sys.platform)

        self._path_to_getllm = os.path.join(CURRENT_PATH, "..", "..", "GetLLM", "GetLLM.py")

    def _prepare_tracking_script(self):
        dict_to_replace = {"%DATA_PATH": self._data_path, "%OUTPUT_PATH": self._output_path}
        self._replace_keywords_in_file(os.path.join(self._data_path, "parse_track_file.sh"),
                                       os.path.join(self._output_path, "parse_track_file.sh"),
                                       dict_to_replace)
        self._replace_keywords_in_file(self._madx_template_path, self._madx_file_path, dict_to_replace)

    def _run_tracking_script(self):
        print "Running tracking code..."
        errcode = madxrunner.runForInputFile(self._madx_file_path, stdout=subprocess.PIPE)
        self.assertEqual(errcode, 0, "Error running MADX tracking code.")

    def _run_drive(self):
        print "Running drive..."
        iotools.copy_item(os.path.join(self._data_path, "Drive.inp"), self._output_path)
        self._replace_keywords_in_file(os.path.join(self._data_path, "DrivingTerms"),
                                       os.path.join(self._output_path, "DrivingTerms"),
                                       {"%OUTPUT_PATH": self._output_path})
        call_command = os.path.abspath(self._path_to_drive) + " " + os.path.abspath(self._output_path)

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

    def _run_getLLM(self):
        print "Running GetLLM..."
        call_command = sys.executable + " " + os.path.abspath(self._path_to_getllm) + \
        " --accel=LHCB1 --tbtana=SUSSIX --bpmu=mm " + \
        " -m " + os.path.join(self._output_path, "twiss.dat") + \
        " -f " + os.path.join(self._output_path, "ALLBPMs") + \
        " -o " + self._output_path

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

    def _compare_output(self):
        print "Comparing output..."
        beta_x_twiss = metaclass.twiss(os.path.join(self._output_path, "getbetax.out"))
        beta_y_twiss = metaclass.twiss(os.path.join(self._output_path, "getbetay.out"))

        for index in range(len(beta_x_twiss.NAME)):
            rel_error = abs((beta_x_twiss.BETX[index] - beta_x_twiss.BETXMDL[index]) / beta_x_twiss.BETX[index])
            self.assertTrue(rel_error < MAX_BETA_REL_ERR,
                            "Relative error too big found in: " + beta_x_twiss.NAME[index] + " (" + str(rel_error) + ")")
        for index in range(len(beta_y_twiss.NAME)):
            rel_error = abs((beta_y_twiss.BETY[index] - beta_y_twiss.BETYMDL[index]) / beta_y_twiss.BETY[index])
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
