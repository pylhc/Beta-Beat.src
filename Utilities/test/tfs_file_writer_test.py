'''
Created on 9 Oct 2013

@author: vimaier
'''
import unittest

import Utilities.iotools
import Utilities.tfs_file_writer
import Utilities.ndiff


class TestTfsFileWriter(unittest.TestCase):


    def setUp(self):
        # Saves created files
        self.__created_files = []


    def tearDown(self):
        # Delete all created files
        for created_file in self.__created_files:
            Utilities.iotools.delete_item(created_file)


    def testTfsFileWriter(self):
        file_from_tfs_writer = "test.out"
        tfs_file_writer = Utilities.tfs_file_writer.TfsFileWriter.open(file_from_tfs_writer)
        tfs_file_writer.add_string_descriptor("NAME", "TWISS")
        tfs_file_writer.add_float_descriptor("MASS", 0.938272013)
        tfs_file_writer.add_comment("I am a comment")
        tfs_file_writer.add_column_names("NAME S BETX ALFX BETY ALFY".split())
        tfs_file_writer.add_column_datatypes("%s %le %le %le %le %le".split())
        tfs_file_writer.add_table_row("BTVSS.6L2.B1  1.125  131.5873094  -1.899115044  67.61780908 1.699566347".split())
        tfs_file_writer.write_to_file()
        self.__created_files.append(tfs_file_writer.get_absolute_file_name_path())

        file_assumed = "test2.out"
        f_out = open(file_assumed, "w")
        f_out.write('@ NAME %s "TWISS"\n')
        f_out.write('@ MASS %le 0.938272013\n')
        f_out.write('# I am a comment\n')
        f_out.write('* NAME S BETX ALFX BETY ALFY\n')
        f_out.write('$ %s %le %le %le %le %le\n')
        f_out.write("BTVSS.6L2.B1 1.125 131.5873094 -1.899115044 67.61780908 1.699566347")
        f_out.close()
        self.__created_files.append(file_assumed)

        self.assertTrue(
                    Utilities.ndiff.compare_files_and_ignore_whitespace(tfs_file_writer.get_absolute_file_name_path(),
                                                                            file_assumed),
                        "Created tfs file is not equal to assumed file"
                        )


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()