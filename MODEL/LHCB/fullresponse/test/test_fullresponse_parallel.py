'''
Created on 20 Oct 2013

:author: vimaier

This module tests the output files of generateFullResponse_parallel.py . First it runs the 'valid' version in the dir test
and then it runs the modified version in the dir fullresponse. After that it compares the content of the dirs
test.data.output.valid and to_check.
'''
import unittest
import os
import sys
import subprocess
import filecmp
import shutil


import __init__ # @UnusedImport init will include paths
import Utilities.iotools

CURRENT_PATH = os.path.dirname(__file__)
# Path 'x/Beta-Beat.src/MODEL/LHCB/fullresponse/test'

_SHORT_RUN = False # If True, genFullResp will only run on first dir
_ARGUMENTS_FILE_NAME = "arguments.txt" # Optional file which is located inside the run dir to manipulate the arguments for genFullResp
_DELETE_VALID_OUTPUT = True

class TestGenFullRespParallel(unittest.TestCase):

    path_to_modified_script = os.path.join(CURRENT_PATH, "..", "generateFullResponse_parallel.py")
    path_to_valid_script = os.path.join(CURRENT_PATH, "valid_generateFullResponse_parallel.py")

    path_to_valid = os.path.join(CURRENT_PATH, "data", "output", "valid")
    path_to_to_check = os.path.join(CURRENT_PATH, "data", "output", "to_check")
    path_to_input = os.path.join(CURRENT_PATH, "data", "input")

    betabeat_root = Utilities.iotools.get_absolute_path_to_betabeat_root()

    successful = False
    __have_to_run_valid_file = False

    def setUp(self):
        self._run_dirs = Utilities.iotools.get_all_dir_names_in_dir(TestGenFullRespParallel.path_to_input)
        self._check_if_valid_output_exist_and_set_attribute()
        self._delete_modified_output()
        self._replace_keyword_in_masks()

    def tearDown(self):
        self._place_keywords_and_thus_produce_masks()
        if TestGenFullRespParallel.successful:
            self._delete_output()

    def testOutput(self):
        print "Start TestGenFullRespParallel"
        for index, run_dir in enumerate(self._run_dirs):
            if self._break_after_first_run(index):
                break
            self._run_valid_file_if_necesaary(run_dir)
            self._run_modified_file(run_dir)
            self._compare_output_dirs(run_dir)
        TestGenFullRespParallel.successful = True
        print "End TestGenFullRespParallel"

    #===============================================================================================
    # helper
    #===============================================================================================
    def _check_if_valid_output_exist_and_set_attribute(self):
        if self._no_valid_output_exists():
            print "No valid output. Have to run valid file."
            TestGenFullRespParallel.__have_to_run_valid_file = True


    def _no_valid_output_exists(self):
        """
        Returns true if cannot find valid output. Assuming that every subdir of path_to_valid has
        the files FullResponse, FullResponse_couple and FullResponse_chromcouple.
        """
        return not self._valid_output_exists()
    def _valid_output_exists(self):
        all_dirs = Utilities.iotools.get_all_dir_names_in_dir(TestGenFullRespParallel.path_to_valid)
        if 0 == len(all_dirs) or len(self._run_dirs) != len(all_dirs):
            return False
        is_valid = False
        for dir_in_valid_output in all_dirs:
            is_valid = self._dir_has_fullresponse_files(
                            os.path.join(TestGenFullRespParallel.path_to_valid, dir_in_valid_output)
                            )
            if not is_valid:
                return False
        return True
    def _dir_has_fullresponse_files(self, path_to_dir):
        """ Returns True if dir contains FullResponse, FullResponse_couple and FullResponse_chromcouple. """
        has_all_files = True
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "FullResponse"))
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "FullResponse_couple"))
        has_all_files = has_all_files and os.path.exists(os.path.join(path_to_dir, "FullResponse_chromcouple"))
        return has_all_files

    def _delete_output(self):
        if _DELETE_VALID_OUTPUT:
            Utilities.iotools.delete_content_of_dir(TestGenFullRespParallel.path_to_valid)
        self._delete_modified_output()

    def _delete_modified_output(self):
        ''' Deletes content in path_to_valid and path_to_to_check. '''
        Utilities.iotools.delete_content_of_dir(TestGenFullRespParallel.path_to_to_check)

    def _replace_keyword_in_masks(self):
        """
        job.iterate.madx and modifiers.madx files are masks with the keyword BB_SRC. It have to be replaced
        """
        for run_dirname in self._run_dirs:
            mask_file_path = os.path.join(TestGenFullRespParallel.path_to_input, run_dirname, "job.iterate.madx")
            replace_dict = {"BB_SRC" : TestGenFullRespParallel.betabeat_root}
            Utilities.iotools.replace_keywords_in_textfile(mask_file_path, replace_dict)
            mask_file_path = os.path.join(TestGenFullRespParallel.path_to_input, run_dirname, "modifiers.madx")
            Utilities.iotools.replace_keywords_in_textfile(mask_file_path, replace_dict)
    def _place_keywords_and_thus_produce_masks(self):
        """ Reverses changes of _replace_keyword_in_maks() """
        for run_dirname in self._run_dirs:
            mask_file_path = os.path.join(TestGenFullRespParallel.path_to_input, run_dirname, "job.iterate.madx")
            self._place_keyword_in_file(mask_file_path)
            mask_file_path = os.path.join(TestGenFullRespParallel.path_to_input, run_dirname, "modifiers.madx")
            self._place_keyword_in_file(mask_file_path)
    def _place_keyword_in_file(self, file_path):
        with open(file_path) as mask_file:
            all_lines = mask_file.read()
        replaced_lines = all_lines.replace(TestGenFullRespParallel.betabeat_root, "%(BB_SRC)s")
        with open(file_path, "w") as mask_file:
            mask_file.write(replaced_lines)


    def _break_after_first_run(self, index):
        return 0 < index and _SHORT_RUN


    def _run_valid_file_if_necesaary(self, run_dir):
        ''' Runs valid script for given dir '''
        if self.__have_to_run_valid_file:
            print "  Run valid script for "+run_dir
            self._run_in_src_dir(TestGenFullRespParallel.path_to_valid_script, run_dir, os.path.join(TestGenFullRespParallel.path_to_valid, run_dir))

    def _run_in_src_dir(self, path_to_script, run_dir_name, output_dst_path):
        """
        genFullResp produces files which have the output dir as path in it. Thus two produced files in different
        output dirs produce different results due to the saved path.
        Solution: Both scripts will be run in the input folder and the output files will be moved to the output_dst_path
        for comparing.
        """
        output_path = os.path.join(TestGenFullRespParallel.path_to_input, run_dir_name)
        self._run_script(path_to_script, output_path)
        self._move_produced_files(output_path, output_dst_path)


    def _run_modified_file(self, run_dir):
        print "  Run modified script for "+run_dir
        self._run_in_src_dir(TestGenFullRespParallel.path_to_modified_script, run_dir, os.path.join(TestGenFullRespParallel.path_to_to_check, run_dir))


    def _run_script(self, path_to_script, path_to_out_run_dir_root):
        ''' Runs given script with path_to_dir.
        path_to_dir has to contain 'job.iterate.madx' and 'modifiers.madx' and if needed the file ___ARGUMENTS_FILE_NAME.
        '''
        arguments_list = self._prepare_and_get_arguments_list(path_to_out_run_dir_root)
        call_command =sys.executable + " " + os.path.abspath(path_to_script) + " " + " ".join(arguments_list)

        process = subprocess.Popen(call_command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           shell=True)

        # wait for the process to terminate
        (out_stream, err_stream) = process.communicate()

        errcode = process.returncode

        if 0 != errcode:
            print "Error running command:", call_command
            print "Printing output:-------------------------"
            print out_stream
            print >> sys.stderr, "Printing error output:-------------------"
            print >> sys.stderr, err_stream
        self.assertEqual(0, errcode, "Execution failed")

    def _move_produced_files(self, src_path, dst_path):
        Utilities.iotools.create_dirs(dst_path)
        Utilities.iotools.delete_content_of_dir(dst_path) # Otherwise shutil.move() would raise an error
        for abs_filename in Utilities.iotools.get_all_absolute_filenames_in_dir_and_subdirs(src_path):
            if self._is_produced_output_file(abs_filename):
                shutil.move(abs_filename, dst_path)
    def _is_produced_output_file(self, file_path):
        return "job.iterate.madx" not in file_path and "modifiers.madx" not in file_path and "arguments.txt" not in file_path

    def _prepare_and_get_arguments_list(self, path_to_out_run_dir_root):
        args_dict = {}
        args_dict["--accel"] = "LHCB1"
        args_dict["--path"] = path_to_out_run_dir_root
        args_dict["--core"] = os.path.join(CURRENT_PATH, "..")
        args_dict["--deltak"] = "0.00002"

        self._set_args_from_file_in_run_dir(path_to_out_run_dir_root, args_dict)

        args_list = []
        for key in args_dict:
            args_list.append(key+"="+args_dict[key])

        return args_list


    def _set_args_from_file_in_run_dir(self, path_to_run, args_dict):
        """ Reads the arguments from file _ARGUMENTS_FILE_NAME
        File should have following syntax:
        --accel=LHCB2
        --deltak=0.00002
        # core and path should not be set.
        """
        arguments_file_path = os.path.join(path_to_run, _ARGUMENTS_FILE_NAME)
        if os.path.exists(arguments_file_path):
            with open(arguments_file_path) as arg_file:
                for arg_line in arg_file:
                    if "--accel" in arg_line:
                        args_dict["--accel"] = arg_line.split("=", 1)[1].strip()
                    if "--deltak" in arg_line:
                        args_dict["--deltak"] = arg_line.split("=", 1)[1].strip()


    def _compare_output_dirs(self, run_dir):
        ''' Compares output by using filecmp '''
        print "  Comparing output files"
        valid_dir = os.path.join(self.path_to_valid, run_dir)
        to_check_dir = os.path.join(self.path_to_to_check, run_dir)
        self._compare_binaries_in_dir_directly(valid_dir, to_check_dir)


    def _compare_dir_with_filecmp(self, valid_dir, to_check_dir):
        dir_compare = filecmp.dircmp(valid_dir, to_check_dir)
        dir_compare.report()
        self.assertEqual(0, len(dir_compare.diff_files), "Files are not equal.")
        self.assertEqual(0, len(dir_compare.left_only), "Files existing in only one dir")
        self.assertEqual(0, len(dir_compare.right_only), "Files existing in only one dir")

    def _compare_binaries_in_dir_directly(self, valid_dir_path, to_check_dir_path):
        for file_name in Utilities.iotools.get_all_filenames_in_dir_and_subdirs(valid_dir_path):
            valid_file_path = os.path.join(valid_dir_path, file_name)
            to_check_file_path = os.path.join(to_check_dir_path, file_name)
            self._compare_binary_files(valid_file_path, to_check_file_path)
    def _compare_binary_files(self, path1, path2):
        with open(path1,"rb") as file1:
            with open(path2, "rb") as file2:
                buffer1 = file1.read(1024)
                buffer2 = file2.read(1024)
                self.assertEqual(buffer1, buffer2, "Files are not equal. \n######################### Buffer1:\n"+
                                 buffer1+ "\n#########################\n != \n######################### Buffer2:\n"+buffer2+
                                 "\n#########################\nFile1: "+path1+"\nFile2: "+path2)



# END TestOutput -----------------------------------------------------------------------------------


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestOutput.testOutput']
    unittest.main()
