"""
Created on 25 Apr 2013

@author: vimaier

@version: 1.0.0

This module contains the class TfsFile which handles the output files of GetLLM.

Change history:
 -
 
"""

import os

DEFAULT_COLUMN_WIDTH = 17
""" Indicates width of columns in output file. """
INVALID_COLUMN_NUMBER = -1
""" Initial value for number of columns. """

# Indexes for the self.__table:
I_COLUMN_TYPES = 0
I_COLUMN_NAMES = 1
I_TABLE_LINES = 2

    
class TfsFile(object):
    '''
    This class represents a TFS file. It stores all header lines and the table in lists and write
    all the content at once by calling the write function.
    '''
    
    s_output_path = ""

    def __init__(self, file_name, column_width=DEFAULT_COLUMN_WIDTH):
        """
        Constructor
        
        :Parameters:
            'file_name': string
                The file name without path where the file will be written.
            'column_width': int
                Indicates the width of each column in the file.
        """
        
        self.__file_name = file_name
        self.__column_width = column_width
        
        self.__getllm_version = ""
        self.__getllm_madfile = ""
        self.__getllm_srcfiles = []
        
        self.__header = [];
        """ The header contains descriptor lines and comment lines. """
        
        self.__table = ([], [], [])
        """ A tupel whereby the first entry is the list of column names, the second. 
            entry is the list of the column data types and the third list contains lists of table
            row data.
        """
        
        self.__is_column_names_set = False
        self.__is_column_types_set = False
        self.__num_of_columns = INVALID_COLUMN_NUMBER
        
        
    def add_getllm_header(self, version=None, mad_file=None):
        """ Adds '@ GetLLMVersion %s "<version>"\n' and/or '@ MAD_FILE %s "<mad_file>"\n' 
            Returns self for concatenation reasons
        """
        if version is not None:
            self.__getllm_version = version
        if mad_file is not None:
            self.__getllm_madfile = mad_file
            
        return self

    
    def add_filename_to_getllm_header(self,file_name):
        """ Adds a file to '@ FILES %s "<files-list>"\n' """
        self.__getllm_srcfiles.append(file_name)
    
        
    def add_descriptor(self, name, data_type, data):
        """
        Adds the string "@ <name> <data_type> <data>" to the tfs header.
        """
        descriptor_line = "@ " + name + " " + data_type + " " + str(data)+"\n"
        
        self.__header.append(descriptor_line)


    def add_comment(self, comment):
        """
        Adds the string "# <comment>" to the tfs header.
        """
        self.__header.append("# "+comment+"\n")

        
    def add_column_names(self, list_names):
        """
        Adds the list of column names to the table header.
        If the number of columns is determined already(e.g. by adding column data types) and the 
        length of list_names does not match the number of columns a TypeError will be raised.
        
        :Parameters:
            'list_names': list
                Containing the names of the columns. Without prefix '$'
        """
        if INVALID_COLUMN_NUMBER != self.__num_of_columns:
            if self.__num_of_columns != len(list_names):
                raise AttributeError("Column number is set already but the length of the given list"+
                                        " does not match.("+self.__file_name+")")
        else:
            self.__num_of_columns = len(list_names)
            
        self.__table[I_COLUMN_NAMES].extend(list_names)        
        self.__is_column_names_set = True
        
        
    def add_column_datatypes(self, list_datatypes): 
        """
        Adds the list of column data types to the table header.
        If the number of columns is determined already(e.g. by adding column names) and the 
        length of list_datatypes does not match the number of columns a TypeError will be raised.
        
        :Parameters:
            'list_datatypes': list
                Containing the data type(%s, %le) of the columns. Without prefix '$'
        """
        if INVALID_COLUMN_NUMBER != self.__num_of_columns:
            if self.__num_of_columns != len(list_datatypes):
                raise AttributeError("Column number is set already but the length of the given list"+
                                        " does not match.("+self.__file_name+")")
        else:
            self.__num_of_columns = len(list_datatypes)
        
        self.__table[I_COLUMN_TYPES].extend(list_datatypes)
        self.__is_column_types_set = True
      
        
    def add_table_row(self, list_row_entries):
        """
        Adds the entries of one row to the table data.
        """
        if INVALID_COLUMN_NUMBER == self.__num_of_columns:
            raise TypeError("Before filling the table, set the names and datatypes.("+
                                self.__file_name+")")
        else:
            if self.__num_of_columns != len(list_row_entries):
                raise TypeError("Number of entries does not match the column number of the table.("+
                                self.__file_name+")")
                
        self.__table[I_TABLE_LINES].append(list_row_entries)
    
    
    def __write_getllm_header(self, tfs_file):
        """ Writes the getllm header to the file. """  
        tfs_file.write('@ GetLLMVersion %s "'+self.__getllm_version+'"\n')
        tfs_file.write('@ MAD_FILE %s "'+self.__getllm_madfile+'"\n')
        
        tfs_file.write('@ FILES %s "')
        tfs_file.write(" ".join(self.__getllm_srcfiles))
        tfs_file.write('"\n')
    
    
    def __write_formatted_table(self, tfs_file):
        """ Writes the table of this object formatted to file. """
        format_for_string = "{0:>"+str(self.__column_width)+"s}"
        
        # Write column names
        first_element = self.__table[I_COLUMN_NAMES][0]
        tfs_file.write( ("* {0:>"+str(self.__column_width-2)+"} ").format(first_element))
        remain =self.__table[I_COLUMN_NAMES][1:]
        tfs_file.write(" ".join(format_for_string.format(entry) for entry in remain))
        tfs_file.write("\n")
        
        # Write column types
        first_element = self.__table[I_COLUMN_TYPES][0]
        tfs_file.write( ("$ {0:>"+str(self.__column_width-2)+"} ").format(first_element))
        remain = self.__table[I_COLUMN_TYPES][1:]
        tfs_file.write(" ".join(format_for_string.format(entry) for entry in remain))
        tfs_file.write("\n")
        
        # Write table lines
        for table_line in self.__table[I_TABLE_LINES]:
            tfs_file.write(" ".join(format_for_string.format(str(entry)) for entry in table_line))
            tfs_file.write("\n")
                    
        
    def __write_unformatted_table(self, tfs_file):
        tfs_file.write("* ")  
        tfs_file.write(" ".join(self.__table[I_COLUMN_NAMES]))
                    
        tfs_file.write("\n$ ")
        
        tfs_file.write(" ".join(self.__table[I_COLUMN_TYPES]))
            
        tfs_file.write("\n")
        
        for row in self.__table[I_TABLE_LINES]:
            tfs_file.write(" ".join(str(entry) for entry in row))
            tfs_file.write("\n")
            
            
    def write_to_file(self, formatted = False):
        """ Writes the stored data to the file with the given filename. """
        #TODO: fix the thing with get[D](x|y) and uncomment        
        if not self.__is_column_names_set or not self.__is_column_types_set:
            formatted = False #Quickfix
#             raise AssertionError("Cannot write before column names and column types are set.("+
#                                 self.__file_name+")")
        
        path = os.path.join(TfsFile.s_output_path, self.__file_name)
        tfs_file = open(path,'w')  
        
        # Header
        self.__write_getllm_header(tfs_file)
        
        tfs_file.writelines(self.__header)
        
        # Table
        if formatted:
            self.__write_formatted_table(tfs_file)
        else:
            self.__write_unformatted_table(tfs_file)
            
        
        tfs_file.close()



        