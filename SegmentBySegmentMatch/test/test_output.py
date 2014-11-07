import __init__  # @UnusedImport
import unittest
import os
from Utilities import iotools, ndiff
import sys
import subprocess


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class TestOutput(unittest.TestCase):

    path_to_valid = os.path.join(CURRENT_PATH, "data", "valid")
    path_to_check = os.path.join(CURRENT_PATH, "data", "to_check")
    path_to_input = os.path.join(CURRENT_PATH, "data", "input")
    path_to_modified_sbs_match = os.path.join(CURRENT_PATH, "..", "SegmentBySegmentMatch.py")

    def tearDown(self):
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

    def testOutput(self):
        input_dirs = iotools.get_all_dir_names_in_dir(self.path_to_input)
        for input_dir_name in input_dirs:
            input_dir_path = os.path.join(self.path_to_input, input_dir_name) + os.sep
            self._preprocess_input(input_dir_path)
            ips_to_test = self._get_ips_to_test(os.path.join(input_dir_path, "test.cfg"))
            for ip in ips_to_test:
                valid_dir_name = input_dir_name + "_IP" + ip
                directory_valid = os.path.join(self.path_to_valid, valid_dir_name)
                self.assertTrue(iotools.exists_directory(directory_valid),
                                "Valid directory for input file " + valid_dir_name + " doesn't exists")
                print "Starting run for " + valid_dir_name
                self._run_sbsmatch_for_input(ip, valid_dir_name, input_dir_path)
                print "Checking results..."
                self._check_results_for_input(valid_dir_name)
                print "Done"
            self._restore_input(input_dir_path)

    def _preprocess_input(self, input_dir_path):
        sbs_beam1_path = os.path.join(input_dir_path, "Beam1", "sbs")
        sbs_beam2_path = os.path.join(input_dir_path, "Beam2", "sbs")
        replace = [("__PATH__", input_dir_path)]
        self._replace_in_files(sbs_beam1_path, replace)
        self._replace_in_files(sbs_beam2_path, replace)

    def _restore_input(self, input_dir_path):
        sbs_beam1_path = os.path.join(input_dir_path, "Beam1", "sbs")
        sbs_beam2_path = os.path.join(input_dir_path, "Beam2", "sbs")
        replace = [(input_dir_path, "__PATH__")]
        self._replace_in_files(sbs_beam1_path, replace)
        self._replace_in_files(sbs_beam2_path, replace)

    def _replace_in_files(self, src, replace_pairs):
        src_files = iotools.get_all_filenames_in_dir(src)
        for file_name in src_files:
            full_file_name = os.path.join(src, file_name)
            self._replace_in_file(full_file_name, replace_pairs)

    def _replace_in_file(self, full_file_name, replace_pairs):  # TODO: Use a python function instead of sed
        sed_command = "sed -i "
        for pattern, replace in replace_pairs:
            full_command = sed_command + "'s" + "#" + pattern + "#" + replace + "#" + "g' " + full_file_name
            subprocess.call(full_command, shell=True)

    def _get_ips_to_test(self, cfg_file_path):
        ip_list = []
        cfg_file = open(cfg_file_path)
        for line in cfg_file.readlines():
            if line.strip() != "":
                ip_list.append(line.strip())
        return ip_list

    def _run_sbsmatch_for_input(self, ip, input_file_name, input_dir_path):
        directory_to_check = os.path.join(self.path_to_check, input_file_name)

        input_arguments = " --temp=" + directory_to_check
        input_arguments += " --ip=" + ip
        input_arguments += " --beam1=" + os.path.join(input_dir_path, "Beam1")
        input_arguments += " --beam2=" + os.path.join(input_dir_path, "Beam2")

        call_command = os.path.abspath(self.path_to_modified_sbs_match) + " " + input_arguments
        call_command = sys.executable + " " + call_command

        process = subprocess.Popen(call_command,
                           stdout=open("test_log.log", "w"),
                           stderr=subprocess.PIPE,
                           shell=True)
        process.communicate()

    def _check_results_for_input(self, input_file_name):
        directory_valid = os.path.join(self.path_to_valid, input_file_name)
        directory_to_check = os.path.join(self.path_to_check, input_file_name, "match")
        self._compare_directories_with_ndiff(directory_valid, directory_to_check)

    def _compare_directories_with_ndiff(self, directory_valid, directory_to_check):
        valid_item_list = os.listdir(directory_valid)
        for valid_item in valid_item_list:
            valid_item_path = os.path.join(directory_valid, valid_item)
            to_check_item_path = os.path.join(directory_to_check, valid_item)
            self.assertTrue(os.path.exists(to_check_item_path), "Unable to find item " + to_check_item_path)
            if(os.path.isdir(valid_item_path)):
                self._compare_directories_with_ndiff(valid_item_path, to_check_item_path)
            else:
                self._compare_files_with_ndiff(valid_item_path, to_check_item_path)

    def _compare_files_with_ndiff(self, valid_file_path, to_check_file_path):
        self.assertTrue(ndiff.compare_files(valid_file_path, to_check_file_path, os.path.join(CURRENT_PATH, "sbs_match_ndiff.cfg")),
                        "Differences found in file: " + to_check_file_path)


if __name__ == "__main__":
    unittest.main()
