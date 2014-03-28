"""
Created on 28 Mar 2014

@author: jcoellod

Temporary extension of the test_output file, to test the drive results with different
amount of columns.

"""
import os
import unittest
from Utilities import ndiff
from drive.test import test_output


class TestOutputForExtension(test_output.TestOutput):

    MAXIMUM_ABS_ERR = 9e-7  # to compare the floating point values of the drive output
    MAXIMUM_REL_ERR = 9e-7  # to compare the floating point values of the drive output

    def _compare_output_dirs(self):
        print "  Comparing output files"
        self._delete_copied_files()
        self._search_and_compare_recursive()

    def _search_and_compare_recursive(self, relative_path=""):
        full_valid_path = os.path.join(self.path_to_valid, relative_path)
        full_to_check_path = os.path.join(self.path_to_to_check, relative_path)
        for index, directory in enumerate(os.listdir(full_valid_path)):
            if relative_path == "" and self._break_after_first_run(index):
                break
            new_valid_path = os.path.join(full_valid_path, directory)
            new_to_check_path = os.path.join(full_to_check_path, directory)
            if os.path.isdir(new_valid_path):
                self._search_and_compare_recursive(os.path.join(relative_path, directory))
            else:
                self.assertTrue(os.path.exists(new_to_check_path),
                                "File " + new_to_check_path + " don't exists in output folder")
                if directory.endswith("linx") or directory.endswith("liny"):
                    self._compare_files_by_value_table(new_valid_path, new_to_check_path)
                else:
                    ndiff.compare_files(new_valid_path, new_to_check_path)
        test_output.TestOutput.successful = True

    def _compare_files_by_value_table(self, valid_dir, to_check_dir):
        difference_count = 0
        to_check_data = self._gather_table_data_from_file(to_check_dir)
        valid_data = self._gather_table_data_from_file(valid_dir)
        for col_index in valid_data.keys():
            for row_index in range(len(valid_data[col_index])):
                fail = False
                if type(valid_data[col_index][row_index]) is float:
                    if (self._get_absolute_error(valid_data[col_index][row_index],
                                                 to_check_data[col_index][row_index]) > self.MAXIMUM_ABS_ERR and
                        self._get_relative_error(valid_data[col_index][row_index],
                                                 to_check_data[col_index][row_index]) > self.MAXIMUM_REL_ERR):
                        fail = True
                elif type(valid_data[col_index][row_index]) is str:
                    if valid_data[col_index][row_index] != to_check_data[col_index][row_index]:
                        fail = True
                else:
                    print('Unknown column data type, cant compare: %s' % col_index)
                    break
                if fail:
                    print('Diff: %s <valid---new> %s (line %s, column %s)' % (str(valid_data[col_index][row_index]),
                                                                              str(to_check_data[col_index][row_index]),
                                                                              row_index + 7,
                                                                              col_index))
                    difference_count += 1
        self.assertEqual(difference_count, 0, str(difference_count) + " differences found in file: " + to_check_dir)

    def _gather_table_data_from_file(self, outFileName):
        columnList = []
        dataFormat = {}
        data = {}
        with open(outFileName) as outFile:
            while True:  # i know this is ugly but there is a good reason to do this here because we sometimes want to read in two lines at once
                line = outFile.readline().strip()
                if line.startswith('@'):  # header comment
                    if line.startswith('@ FILES'):  # data file path
                        data['FILENAME'] = line.split('"')[1]
                elif line.startswith('*'):  # column names and formats
                    columnLine = line.split('*')[1].split()
                    formatLine = outFile.readline()
                    formatLine = formatLine.split('$')[1].split()
                    assert len(columnLine) == len(formatLine)

                    columnList = columnLine
                    for i in range(len(columnLine)):
                        dataFormat[columnLine[i]] = formatLine[i]
                        data[columnLine[i]] = []
                else:  # data
                    line = line.split()
                    if len(line) == 0:
                        break  # last line is empty

                    self.assertEqual(len(line),
                                     len(dataFormat),
                                     outFileName +
                                      ':\nError in the data! Wrong number of cols: ' +
                                      str(len(line)) + ' (expected: ' + str(len(dataFormat)) +
                                      ').')

                    for i in range(len(line)):
                        if dataFormat[columnList[i]] == '%s':
                            data[columnList[i]].append(line[i].replace('"', ''))
                        elif dataFormat[columnList[i]] == '%le':
                            data[columnList[i]].append(float(line[i]))
                        else:
                            raise ValueError('Unknown data format: %s' % dataFormat[columnList[i]])
        return data

    def _get_absolute_error(self, value1, value2):
        return abs(value1 - value2)

    def _get_relative_error(self, value1, value2):
        min_val = min(abs(value1), abs(value2))
        if min_val == 0:
            min_val = 1
        return abs(value1 - value2) / min_val


if __name__ == "__main__":
    unittest.main()
