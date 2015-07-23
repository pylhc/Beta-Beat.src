import __init__  # @UnusedImport
import unittest
import os
from Utilities import iotools
import sys
import subprocess


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

TOLERANCE = 1e-8


class TestOutput(unittest.TestCase):

    def setUp(self):
        self.path_to_valid = os.path.join(CURRENT_PATH, "data", "valid")
        self.path_to_check = os.path.join(CURRENT_PATH, "data", "to_check")
        self.path_to_input = os.path.join(CURRENT_PATH, "data", "input")
        self.path_to_modified_sbs_match = os.path.abspath(os.path.join(CURRENT_PATH, "..", "SegmentBySegmentMatch.py"))
        self._successful_test = False

    def tearDown(self):
        if self._successful_test:
            to_check_dirs = iotools.get_all_dir_names_in_dir(self.path_to_check)
            for to_check_dir in to_check_dirs:
                directory_to_check = os.path.join(self.path_to_check, to_check_dir)
                input_arguments = " clean --temp=" + directory_to_check
                call_command = os.path.abspath(self.path_to_modified_sbs_match) + " " + input_arguments
                call_command = sys.executable + " " + call_command
                process = subprocess.Popen(call_command,
                               shell=True)
                process.communicate()
                iotools.delete_item(directory_to_check)
            iotools.delete_item(os.path.join(self.path_to_check, "test_log.log"))
        else:
            print "Test wasn't successful, will not delete output"

    def testOutput(self):
        input_dirs = iotools.get_all_dir_names_in_dir(self.path_to_input)
        for input_dir_name in input_dirs:
            input_dir_path = os.path.join(self.path_to_input, input_dir_name)
            directory_to_check = os.path.join(self.path_to_check, input_dir_name)
            directory_valid = os.path.join(self.path_to_valid, input_dir_name)
            ips_to_test = self._get_ips_to_test(os.path.join(input_dir_path, "test.cfg"))
            print "Starting run for " + input_dir_path
            for ip in ips_to_test:
                print "    Running ip", ip
                self._run_sbsmatch_for_input(ip, input_dir_path, directory_to_check)
                print "    Checking results..."
                self._check_results_for_ip(ip, directory_to_check, directory_valid)
                print "    Done"
        self._successful_test = True

    def _get_ips_to_test(self, cfg_file_path):
        ip_list = []
        cfg_file = open(cfg_file_path)
        for line in cfg_file.readlines():
            if line.strip() != "":
                ip_list.append(int(line.strip()))
        return ip_list

    def _run_sbsmatch_for_input(self, ip, input_dir_path, directory_to_check):
        beam1_directory_path = os.path.join(input_dir_path, "beam1")
        beam2_directory_path = os.path.join(input_dir_path, "beam2")
        log_stream = open(os.path.join(self.path_to_check, "test_log.log"), "w")

        self._run_sbsmatch_variables_for_input(ip, directory_to_check, log_stream)
        self._run_sbsmatch_constraints_for_input(ip, directory_to_check, beam1_directory_path, beam2_directory_path, log_stream)
        self._run_sbsmatch_match_for_input(ip, directory_to_check, beam1_directory_path, beam2_directory_path, log_stream)

    def _run_sbsmatch_variables_for_input(self, ip, directory_to_check, log_stream):
        input_arguments = " variables"
        input_arguments += " --temp=" + directory_to_check
        input_arguments += " --ip=" + str(ip)

        var_command = self.path_to_modified_sbs_match + input_arguments
        self._run_subprocess(var_command, log_stream)

    def _run_sbsmatch_constraints_for_input(self, ip, directory_to_check, directory_beam1, directory_beam2, log_stream):
        input_arguments = " constraints --useerrors"
        input_arguments += " --temp=" + directory_to_check
        input_arguments += " --ip=" + str(ip)
        input_arguments += " --beam1=" + directory_beam1
        input_arguments += " --beam2=" + directory_beam2
        input_arguments += " --run=2"

        constraints_command = self.path_to_modified_sbs_match + input_arguments
        self._run_subprocess(constraints_command, log_stream)

    def _run_sbsmatch_match_for_input(self, ip, directory_to_check, directory_beam1, directory_beam2, log_stream):
        input_arguments = " --temp=" + directory_to_check
        input_arguments += " --ip=" + str(ip)
        input_arguments += " --beam1=" + directory_beam1
        input_arguments += " --beam2=" + directory_beam2
        input_arguments += " --run=2"

        match_command = self.path_to_modified_sbs_match + input_arguments
        self._run_subprocess(match_command, log_stream)

    def _run_subprocess(self, command, log_stream):
        command = "python " + command
        process = subprocess.Popen(command,
                           stdout=log_stream,
                           stderr=log_stream,
                           shell=True)
        process.communicate()

    def _check_results_for_ip(self, ip, directory_to_check, directory_valid):
        to_check_changeparameters_file = os.path.join(directory_to_check,
                                                      "match", "changeparameters.madx")
        valid_changeparameters_file = os.path.join(directory_valid, "match",
                                                   "changeparameters_ip" + str(ip) + ".madx")
        self.assertTrue(os.path.isfile(to_check_changeparameters_file),
                        "The changeparameters file hasn't been generated")
        self.assertTrue(os.path.isfile(valid_changeparameters_file),
                                       "There is not valid file for this IP")
        self.assertTrue(self._compare_changeparameter_files(
                            valid_changeparameters_file,
                            to_check_changeparameters_file),
                        "Differences found"
        )

    def _compare_changeparameter_files(self, valid_file_path, to_check_file_path):
        are_equal = True
        valid_delta_to_value = self._get_deltas_from_changeparameters(valid_file_path)
        to_check_delta_to_value = self._get_deltas_from_changeparameters(to_check_file_path)
        for name, valid_value in valid_delta_to_value.iteritems():
            if name not in to_check_delta_to_value:
                are_equal = False
                print >> sys.stderr, "Name", name, "not in file", to_check_file_path
            else:
                to_check_value = to_check_delta_to_value[name]
                difference = abs(to_check_value - valid_value)
                if  difference > TOLERANCE:
                    are_equal = False
                    print >> sys.stderr, "Too big difference found in", name, difference
        return are_equal

    def _get_deltas_from_changeparameters(self, changeparamenters_path):
        delta_to_value = {}
        with open(changeparamenters_path, "r") as changeparamenters_file:
            for line in changeparamenters_file:
                line = line.strip()
                if line == "":
                    continue
                parts = line.split()
                if len(parts) != 3:
                    continue
                delta_name = parts[0]
                value_str = parts[2].replace(";", "")
                try:
                    value = float(value_str)
                except ValueError:
                    continue
                delta_to_value[delta_name] = value
        return delta_to_value


if __name__ == "__main__":
    unittest.main()
