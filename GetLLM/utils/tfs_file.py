"""
.. module:: tfs_file_writer

Created on 25 Apr 2013

This module contains the class GetllmTfsFile which handles the output files of GetLLM.

.. moduleauthor:: Viktor Maier <viktor.maier@cern.ch>
"""

import sys
import datetime

import Utilities.tfs_file_writer


class GetllmTfsFile(Utilities.tfs_file_writer.TfsFileWriter):
    '''
    It stores additionally to TfsFileWriter the following descriptors for :module:'GetLLM.py':
     - GetLLMVersion
     - Command
     - Date
     - FILES
    '''

    s_output_path = ""
    s_getllm_version = ""
    s_mad_filename = ""
    __s_current_date = datetime.datetime.today().strftime("%d. %B %Y, %H:%M:%S")#e.g.: 17. July 2013, 12:28:56
    __s_getllm_invocation_command = " ".join(sys.argv)

    def __init__(self, file_name, column_width=Utilities.tfs_file_writer.TfsFileWriter.DEFAULT_COLUMN_WIDTH):
        """
        Constructor

        :param string file_name: The file name without path where the file will be written.
        :param int column_width: Indicates the width of each column in the file.
        """
        super(GetllmTfsFile, self).__init__(file_name, column_width)
        self.set_outputpath(GetllmTfsFile.s_output_path)
        self.__getllm_srcfiles = []

        self.add_string_descriptor("GetLLMVersion", GetllmTfsFile.s_getllm_version)
        self.add_string_descriptor("Command", GetllmTfsFile.__s_getllm_invocation_command)
        self.add_string_descriptor("Date", GetllmTfsFile.__s_current_date)


    def add_filename_to_getllm_header(self, file_name):
        """ Adds a file to '@ FILES %s "<files-list>"\n' """
        self.__getllm_srcfiles.append(file_name)


    def write_to_file(self, formatted=True):
        self.add_string_descriptor("FILES", ",".join(self.__getllm_srcfiles))
        Utilities.tfs_file_writer.TfsFileWriter.write_to_file(self, formatted=formatted)




