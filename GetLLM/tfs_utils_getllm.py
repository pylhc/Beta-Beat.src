"""
.. module:: tfs_file_writer

Created on 25 Apr 2013

This module contains the class GetllmTfsFile which handles the output files of GetLLM.

.. moduleauthor:: Viktor Maier <viktor.maier@cern.ch>
"""

import os
import sys
import datetime

import utils.tfs_file_writer
from utils import tfs_pandas


class GetllmTfsFile(utils.tfs_file_writer.TfsFileWriter):
    '''
    It stores additionally to TfsFileWriter the following descriptors for :module:'GetLLM.py':
     - GetLLMVersion
     - Command
     - Current working directory
     - Date
     - FILES
    '''

    s_output_path = ""
    s_getllm_version = ""
    s_mad_filename = ""
    s_current_date = datetime.datetime.today().strftime("%d. %B %Y, %H:%M:%S")#e.g.: 17. July 2013, 12:28:56
    s_getllm_invocation_command = sys.executable+" '"+"' '".join([]+sys.argv)+"'"
    s_getllm_current_working_dir = os.getcwd()

    def __init__(self, file_name, column_width=utils.tfs_file_writer.TfsFileWriter.DEFAULT_COLUMN_WIDTH):
        """
        Constructor

        :param string file_name: The file name without path where the file will be written.
        :param int column_width: Indicates the width of each column in the file.
        """
        super(GetllmTfsFile, self).__init__(file_name, column_width)
        self.set_outputpath(GetllmTfsFile.s_output_path)
        self.__getllm_srcfiles = []

        self.add_string_descriptor("GetLLMVersion", GetllmTfsFile.s_getllm_version)
        self.add_string_descriptor("Command", GetllmTfsFile.s_getllm_invocation_command)
        self.add_string_descriptor("CWD", GetllmTfsFile.s_getllm_current_working_dir)
        self.add_string_descriptor("Date", GetllmTfsFile.s_current_date)


    def add_filename_to_getllm_header(self, file_name):
        """ Adds a file to '@ FILES %s "<files-list>"\n' """
        self.__getllm_srcfiles.append(file_name)


    def write_to_file(self, formatted=True):
        self.add_string_descriptor("FILES", ",".join(self.__getllm_srcfiles))
        utils.tfs_file_writer.TfsFileWriter.write_to_file(self, formatted=formatted)


class GetllmPandasTfs:
    """Class that mocks the behaviour of GetllmTfsFile.
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self.__getllm_srcfiles = []
        self.__headers = {}
        self.__df = None

        self.add_string_descriptor("GetLLMVersion", GetllmTfsFile.s_getllm_version)
        self.add_string_descriptor("Command", GetllmTfsFile.s_getllm_invocation_command)
        self.add_string_descriptor("CWD", GetllmTfsFile.s_getllm_current_working_dir)
        self.add_string_descriptor("Date", GetllmTfsFile.s_current_date)

    def add_header(self, key, value):
        self.__headers[key] = value

    def add_string_descriptor(self, key, value):
        self.add_header(key, value)

    def add_float_descriptor(self, key, value):
        self.add_header(key, value)

    def add_filename_to_getllm_header(self, file_name):
        """ Adds a file to '@ FILES %s "<files-list>"\n' """
        self.__getllm_srcfiles.append(file_name)

    def set_dataframe(self, dataframe):
        self.__df = dataframe

    def write_to_file(self, formatted=True):
        self.add_string_descriptor("FILES", ",".join(self.__getllm_srcfiles))
        if self.__df is None:
            raise AttributeError("DataFrame not set")
        tfs_pandas.write_tfs(
            os.path.join(GetllmTfsFile.s_output_path, self.file_name),
            self.__df, self.__headers)

