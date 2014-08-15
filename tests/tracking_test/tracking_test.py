'''
Created on Jul 18, 2014

This test class runs a tracking simulation for each file in "data/tracking_jobs", then it runs the Drive + GetLLM
algorithms and compares the beta outputs.

The maximum allowed relative error can be set in the constant MAX_BETA_REL_ERR.

@author: jcoellod
'''
import __init__  # @UnusedImport
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
SHORT_RUN = False

IP_SEGMENTS_B1 = "BPM.15L2.B1,BPM.15R2.B1,IP2"
IP_SEGMENTS_B2 = "BPM.15L2.B2,BPM.15R2.B2,IP2"


class Test(unittest.TestCase):

    def setUp(self):
        self._set_up_paths()
        self._prepare_tracking_scripts()

    def tearDown(self):
        try:
            os.unlink("ats")
            os.unlink("db")
            os.unlink("db5")
            os.unlink("ds")
        except:
            pass
        for output_path in self._output_paths:
            # iotools.delete_item(output_path)
            pass

    def testName(self):
        for madx_file_path, output_path in zip(self._madx_file_paths, self._output_paths):
            print "Testing madx job " + os.path.basename(madx_file_path)
            self._run_tracking_script(madx_file_path)
            self._run_drive(output_path, madx_file_path)
            self._run_getLLM(output_path)
            self._run_Segment_by_Segment(output_path)
            self._run_Segment_by_Segment_Match(output_path)
            self._compare_output(output_path)
            if SHORT_RUN:
                break
            print "Done\n"

    def _set_up_paths(self):
        self._data_path = os.path.join(CURRENT_PATH, "data")
        self._tracking_jobs_path = os.path.join(self._data_path, "tracking_jobs")
        self._output_paths = []
        self._madx_template_paths = []
        self._madx_file_paths = []
        for file_name in iotools.get_all_filenames_in_dir(self._tracking_jobs_path):
            if not file_name.startswith("~"):
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
        self._path_to_sbs = os.path.join(CURRENT_PATH, "..", "..", "SegmentBySegment", "SegmentBySegment.py")
        self._path_to_sbs_match = os.path.join(CURRENT_PATH, "..", "..", "SegmentBySegmentMatch", "SegmentBySegmentMatch.py")

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
        print "Running Drive..."
        self._prepare_drive_input(output_path, madx_file_path)
        self._replace_keywords_in_file(os.path.join(self._data_path, "DrivingTerms"),
                                       os.path.join(output_path, "DrivingTerms"),
                                       {"%OUTPUT_PATH": output_path})
        call_command = os.path.abspath(self._path_to_drive) + " " + os.path.abspath(output_path)

        self._run_outer_process(call_command, "Drive")

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

        self._run_outer_process(call_command, "GetLLM")

    def _run_Segment_by_Segment(self, output_path):
        print "Running Segment by Segment..."
        madx_bin_path = madxrunner.get_sys_dependent_path_to_mad_x()
        bb_source_path = iotools.get_absolute_path_to_betabeat_root()
        sbs_output_path = os.path.join(output_path, "sbs")
        twiss_model_path = os.path.join(output_path, "twiss.dat")

        ip_segments = None
        sequence = self._get_sequence(twiss_model_path)
        if sequence == "LHCB1":
            ip_segments = IP_SEGMENTS_B1
        elif sequence == "LHCB2":
            ip_segments = IP_SEGMENTS_B2
        self.assertTrue(ip_segments, "Unknown accelerator sequence.")

        call_command = sys.executable + " " + os.path.abspath(self._path_to_sbs) + \
        " --path " + output_path + \
        " --save " + sbs_output_path + \
        " --twiss " + twiss_model_path + \
        " --mad " + madx_bin_path + \
        " --bbsource " + bb_source_path + \
        " --accel " + sequence + " --cuts 10 --start " + ip_segments

        self._run_outer_process(call_command, "Segment by Segment")

    def _get_sequence(self, twiss_path):
        with open(twiss_path, "r") as twiss_file:
            for line in twiss_file:
                if "SEQUENCE" in line:
                    return line.split()[3].strip('"')

    def _run_Segment_by_Segment_Match(self, output_path):
        print "Running Segment by Segment Match..."
        beam1path = ""
        beam2path = ""
        if "beam1" in output_path:
            beam1path = output_path
            beam2path = beam1path.replace("beam1", "beam2")
        if "beam2" in output_path:
            beam2path = output_path
            beam1path = beam2path.replace("beam2", "beam1")
        if os.path.isdir(os.path.join(beam1path, "sbs")) and os.path.isdir(os.path.join(beam2path, "sbs")):
            temp_dir = beam1path.replace("beam1", "match")
            iotools.create_dirs(temp_dir)
            ip = 2
            call_command = sys.executable + " " + os.path.abspath(self._path_to_sbs_match) + \
            " --ip " + str(ip) + \
            " --beam1 " + beam1path + \
            " --beam2 " + beam2path + \
            " --temp " + temp_dir

            self._run_outer_process(call_command, "Segment by Segment Match")
        else:
            print "No data for both beams yet."

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

    def _run_outer_process(self, command, name):
        process = subprocess.Popen(command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)

        (out_stream, err_stream) = process.communicate()
        errcode = process.returncode
        if 0 != errcode:
            print "Error running", name + ":", command
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream
        self.assertEqual(errcode, 0, "Error running " + name)

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
