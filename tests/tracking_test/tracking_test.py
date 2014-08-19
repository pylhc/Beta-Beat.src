'''
Created on Jul 18, 2014

This test class runs a tracking simulation for each file in "data/tracking_jobs" that doesn't start with "~".
In this tracking jobs the "twiss_elements.dat" file must be generated.
The region where the optics are defined must be surrounded by the %OPTICS_START and %OPTICS_END tags (view examples),
because this text will be copied to the modifiers file needed by Segment by Segment.

Segment by segment match will only run if it finds information for both beams. The names of the jobs must be the same
only different in the words 'beam1' or 'beam2' depending on the case, for example:

- job.tracking.injection.beam1.madx
- job.tracking.injection.beam2.madx

If this happens, manual errors can be input in the madx jobs, marking them with the %ERROR tag:

- %ERROR kq8.l1b1 = kq8.l1b1 + 2e-5;

Then the test will check if segment by segment match has found the error.


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
MAX_CORRECTION_DEVIATION = 1e-5
SHORT_RUN = False

IP_SEGMENTS_B1 = "BPM.15L1.B1,BPM.15R1.B1,IP1"
IP_SEGMENTS_B2 = "BPM.15L1.B2,BPM.15R1.B2,IP1"


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
            iotools.delete_item(output_path)
        for output_path in self._sbs_match_output_paths:
            iotools.delete_item(output_path)

    def testName(self):
        for madx_file_path, output_path in zip(self._madx_file_paths, self._output_paths):
            print "Testing madx job " + os.path.basename(madx_file_path)
            self._run_tracking_script(madx_file_path)
            sequence = self._get_sequence(os.path.join(output_path, "twiss.dat"))
            self._run_drive(output_path, madx_file_path, sequence)
            self._run_getLLM(output_path, sequence)
            self._run_Segment_by_Segment(output_path, sequence)
            self._run_Segment_by_Segment_Match(output_path)
            self._compare_output_betas(output_path)
            if SHORT_RUN:
                break
            print "Done\n"
        print "Checking match errors..."
        self._compare_match_errors()
        print "All done"

    def _set_up_paths(self):
        self._data_path = os.path.join(CURRENT_PATH, "data")
        self._tracking_jobs_path = os.path.join(self._data_path, "tracking_jobs")
        self._output_paths = []
        self._madx_template_paths = []
        self._madx_file_paths = []
        self._sbs_match_output_paths = []
        self._errors = {}
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
            self._errors[output_dir] = {}
            modifiers = self._find_optics_and_errors_in_job(output_job, output_dir)
            iotools.write_string_into_new_file(os.path.join(output_dir, "modifiers.madx"), modifiers)

    def _find_optics_and_errors_in_job(self, madx_job, output_dir):
        new_lines = ""
        modifiers = ""
        in_modifiers = False
        with open(madx_job, "r") as madx_job_data:
            for line in madx_job_data:
                if "%OPTICS_END" in line and in_modifiers:
                    in_modifiers = False
                if in_modifiers:
                    modifiers += line
                if "%OPTICS_START" in line and not in_modifiers:
                    in_modifiers = True
                if not "%OPTICS_START" in line and not "%OPTICS_END" in line:
                    if "%ERROR" in line:
                        clean_line = line.replace("%ERROR", "")  # Example: %ERROR kq8.l1b1 = kq8.l1b1 + 2e-5;
                        var_value = clean_line.replace(";", "").replace("\n", "").split("=")[1].split("+")
                        var_value = var_value if len(var_value) == 2 else var_value[0].split("-")
                        self._errors[output_dir][var_value[0].strip()] = float(var_value[1])
                        new_lines += clean_line
                    else:
                        new_lines += line
        iotools.write_string_into_new_file(madx_job, new_lines)
        return modifiers

    def _run_tracking_script(self, madx_file_path):
        print "Running tracking code..."
        errcode = madxrunner.runForInputFile(madx_file_path, stdout=subprocess.PIPE)
        self.assertEqual(errcode, 0, "Error running MADX tracking code.")

    def _run_drive(self, output_path, madx_file_path, sequence):
        print "Running Drive..."
        self._prepare_drive_input(output_path, madx_file_path, sequence)
        self._replace_keywords_in_file(os.path.join(self._data_path, "DrivingTerms"),
                                       os.path.join(output_path, "DrivingTerms"),
                                       {"%OUTPUT_PATH": output_path})
        call_command = os.path.abspath(self._path_to_drive) + " " + os.path.abspath(output_path)

        self._run_outer_process(call_command, "Drive")

    def _prepare_drive_input(self, output_path, madx_file_path, sequence):
        tune_x, tune_y = 0.0, 0.0
        with open(madx_file_path, "r") as madx_file:
            for line in madx_file:
                if "mux=" in line and "muy=" in line:
                    tunes = re.findall("\d+.\d+", line)
                    tune_x = float(tunes[0]) % 1
                    tune_y = float(tunes[1]) % 1
        self.assertTrue(tune_x != 0.0 and tune_y != 0.0, "Cannot find tunes in " + madx_file_path)
        print "Tune x: " + str(tune_x) + ", Tune y: " + str(tune_y)
        if(sequence == "LHCB2"):  # TODO: Temporary sign change to make drive work properly with the tracking of beam 2
            tune_x = -tune_x
            tune_y = -tune_y
        self._replace_keywords_in_file(os.path.join(self._data_path, "Drive.inp"),
                                       os.path.join(output_path, "Drive.inp"),
                                       {"%TUNE_X": str(tune_x), "%TUNE_Y": str(tune_y)})

    def _run_getLLM(self, output_path, sequence):
        print "Running GetLLM..."
        call_command = sys.executable + " " + os.path.abspath(self._path_to_getllm) + \
        " --accel=" + sequence + " --tbtana=SUSSIX --bpmu=mm " + \
        " -m " + os.path.join(output_path, "twiss.dat") + \
        " -f " + os.path.join(output_path, "ALLBPMs") + \
        " -o " + output_path

        self._run_outer_process(call_command, "GetLLM")

    def _run_Segment_by_Segment(self, output_path, sequence):
        print "Running Segment by Segment..."
        madx_bin_path = madxrunner.get_sys_dependent_path_to_mad_x()
        bb_source_path = iotools.get_absolute_path_to_betabeat_root()
        sbs_output_path = os.path.join(output_path, "sbs")
        twiss_model_path = os.path.join(output_path, "twiss_elements.dat")

        ip_segments = None
        if sequence == "LHCB1":
            ip_segments = IP_SEGMENTS_B1
        elif sequence == "LHCB2":
            ip_segments = IP_SEGMENTS_B2
        self.assertTrue(ip_segments, "Unknown accelerator sequence.")

        call_command = sys.executable + " " + os.path.abspath(self._path_to_sbs) + \
        " --path=" + output_path + \
        " --save=" + sbs_output_path + \
        " --twiss=" + twiss_model_path + \
        " --mad=" + madx_bin_path + \
        " --bbsource=" + bb_source_path + \
        " --accel=" + sequence + " --start=" + ip_segments

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
            self._sbs_match_output_paths.append(temp_dir)
            ip = 1  # TODO: This has to fit the segments
            call_command = sys.executable + " " + os.path.abspath(self._path_to_sbs_match) + \
            " --ip " + str(ip) + \
            " --beam1 " + beam1path + \
            " --beam2 " + beam2path + \
            " --temp " + temp_dir

            self._run_outer_process(call_command, "Segment by Segment Match")
        else:
            print "No data for both beams yet."

    def _compare_output_betas(self, output_path):
        print "Comparing output..."
        if len(self._errors[output_path]) == 0:
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
        else:
            print "Manual errors in the madx job: skipping betas check"

    def _compare_match_errors(self):
        for sbs_match_output_path in self._sbs_match_output_paths:
            print "Checking path", sbs_match_output_path, "..."
            match_dir = os.path.join(sbs_match_output_path, "match")
            changeparameters_path = os.path.join(match_dir, "changeparameters.madx")
            beam1_path = sbs_match_output_path.replace("match", "beam1")
            beam2_path = sbs_match_output_path.replace("match", "beam2")
            errors_beam1 = self._errors[beam1_path]
            errors_beam2 = self._errors[beam2_path]
            with open(changeparameters_path, "r") as correction_lines:
                for correction_line in correction_lines:
                    correction_split = correction_line.replace(";", "").replace("d", "").replace("\n", "").split("=")
                    variable, value = correction_split[0].strip(), float(correction_split[1])
                    if variable in errors_beam1:
                        difference = abs(value - errors_beam1[variable])
                        self.assertTrue(difference < MAX_CORRECTION_DEVIATION,
                                               "Wrong correction for variable: " + variable)
                        print "Error for variable", variable, "corrected."
                    elif variable in errors_beam2:
                        difference = abs(value - errors_beam2[variable])
                        self.assertTrue(difference < MAX_CORRECTION_DEVIATION,
                                               "Wrong correction for variable: " + variable)
                        print "Error for variable", variable, "corrected."
                    else:
                        self.assertTrue(value < MAX_CORRECTION_DEVIATION, "Wrong correction for variable: " + variable)
        print "Done\n"

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
