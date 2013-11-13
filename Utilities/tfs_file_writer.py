"""
.. module:: tfs_file_writer


This module contains the class TfsFileWriter which is used to create easily TfsFiles.

.. moduleauthor:: Viktor Maier <viktor.maier@cern.ch>

Usage1::

    import Utilities.tfs_file_writer as tfs_writer
    ...
    tfs_file_writer = tfs_writer.TfsFileWriter.open("my_file.out")
    tfs_file_writer.set_column_width(15) # Specify desired column width
    tfs_file_writer.set_outputpath("/x/y/z/") # Specify a directory if desired

    tfs_file_writer.add_string_descriptor("NAME", "TWISS")
    tfs_file_writer.add_float_descriptor("MASS", 0.938272013)
    tfs_file_writer.add_comment("I am a comment")
    tfs_file_writer.add_column_names("NAME S BETX ALFX BETY ALFY".split())
    tfs_file_writer.add_column_datatypes("%s %le %le %le %le %le".split())
    tfs_file_writer.add_table_row("BTVSS.6L2.B1  1.125  131.5873094  -1.899115044  67.61780908 1.699566347".split())

    tfs_file_writer.write_to_file()

Usage2::

    tfs_file_writer = tfs_writer.TfsFileWriter.open("/x/y/z/my_file.out")
    #tfs_file_writer.set_outputpath("/x/y/z/") # Don't do this! Outputpath is already in file_name.

    ... # Add at least column_names, dolumn_datatypes and one table_row
    tfs_file_writer.write_to_file()

"""

import os

import Utilities.iotools
import Utilities.math



class TfsFileWriter(object):
    '''
    This class represents a TFS file. It stores all header lines and the table and write
    all the content formatted at once by calling the write function.
    '''

    DEFAULT_COLUMN_WIDTH = 17
    # Indicates width of columns in output file.
    MIN_COLUMN_WIDTH = 4

    @staticmethod
    def open(file_name):
        """ This function will create and return a TfsFileWriter object with the given filename.
            No file will be opened on the file system yet. An actual file will first be created
            after adding at least column_names, dolumn_datatypes and one table_row and the method
            call write_to_file().
         """
        return TfsFileWriter(file_name)

    def __init__(self, file_name, outputpath=None, column_width=DEFAULT_COLUMN_WIDTH):
        """
        Constructor

        :param string file_name: The file name without path where the file will be written (outputpath can be set via set_outputpath).
        :param int column_width: Indicates the width of each column in the file.
        """

        self.__file_name = ""
        self.__outputpath = ""
        self.__column_width = 0
        self.__tfs_header_lines = [] # Holds instances of subclasses of _TfsHeaderLine
        self.__tfs_table = _TfsTable(self)

        self.set_file_name(file_name)
        self.set_outputpath(outputpath)
        self.set_column_width(column_width)


    def set_file_name(self, file_name):
        if not isinstance(file_name, str) or 0 == len(file_name):
            raise ValueError("File name is not valid: "+ file_name)
        self.__file_name = file_name

    def set_outputpath(self, outputpath):
        if outputpath is None or not isinstance(outputpath, str):
            outputpath = os.path.abspath("./")
        if Utilities.iotools.not_exists_directory(outputpath):
            Utilities.iotools.create_dirs(outputpath)
        self.__outputpath = outputpath

    def set_column_width(self, column_width):
        if not isinstance(column_width, (int, long)) or column_width < TfsFileWriter.MIN_COLUMN_WIDTH:
            column_width = TfsFileWriter.DEFAULT_COLUMN_WIDTH
        self.__column_width = column_width


    def get_file_name(self):
        return self.__file_name


    def add_string_descriptor(self, name, str_value):
        """ Adds the string "@ <name> %s <data>" to the tfs header. """
        tfs_descriptor = _TfsDescriptor(name, str_value, _TfsDataType.get_new_string_instance())
        self.__tfs_header_lines.append(tfs_descriptor)


    def add_float_descriptor(self, name, float_value):
        """ Adds the string "@ <name> %le <data>" to the tfs header. """
        tfs_descriptor = _TfsDescriptor(name, float_value, _TfsDataType.get_new_float_instance())
        self.__tfs_header_lines.append(tfs_descriptor)


    def add_comment(self, comment):
        """
        Adds the string "# <comment>" to the tfs header.
        """
        self.__tfs_header_lines.append(_TfsComment(comment))


    def add_column_names(self, list_names):
        """
        Adds the list of column names to the table header.
        If the number of columns is determined already(e.g. by adding column data types) and the
        length of list_names does not match the number of columns a TypeError will be raised.

        :param list list_names: Containing the names of the columns. Without prefix '*'
        """
        self.__tfs_table.add_column_names(list_names)


    def add_column_datatypes(self, list_datatypes):
        """
        Adds the list of column data types to the table header.
        If the number of columns is determined already(e.g. by adding column names) and the
        length of list_datatypes does not match the number of columns a TypeError will be raised.

        :param list list_datatypes: Containing the data type(%s, %le) of the columns. Without prefix '$'
        """
        self.__tfs_table.add_column_datatypes(list_datatypes)


    def add_table_row(self, list_row_entries):
        """
        Adds the entries of one row to the table data.
        :param list list_row_entries: Values of one row. Datatypes will not be checked. Only length
                                    will be checked with the length of column names.
        """
        self.__tfs_table.add_table_row(list_row_entries)


    def get_absolute_file_name_path(self):
        return os.path.join(self.__outputpath, self.__file_name)

    def write_to_file(self, formatted = True):
        """ Writes the stored data to the file with the given filename. """
        if not self.__tfs_table.are_column_names_and_types_are_set():
            print self.__file_name+":", "Abort writing file. Cannot write file until column names and types are set."
            return

        if self.__tfs_table.is_empty():
            print self.__file_name+":", "Abort writing file. No rows in table."
            return

        path = self.get_absolute_file_name_path()
        tfs_file = open(path,'w')

        # Header
        tfs_file.writelines(x.get_line_as_string_with_newline() for x in self.__tfs_header_lines)

        # Table
        if formatted:
            self.__write_formatted_table(tfs_file)
        else:
            self.__write_unformatted_table(tfs_file)


        tfs_file.close()


    def __write_formatted_table(self, tfs_file):
        """ Writes the table of this object formatted to file. """
        format_for_string = "{0:>"+str(self.__column_width)+"s}"

        # Write column names
        list_column_names = self.__tfs_table.get_column_names()
        first_element = list_column_names[0]
        tfs_file.write( ("* {0:>"+str(self.__column_width-2)+"} ").format(first_element))
        remain = list_column_names[1:]
        tfs_file.write(" ".join(format_for_string.format(entry) for entry in remain))
        tfs_file.write("\n")

        # Write column types
        list_column_types = self.__tfs_table.get_column_data_types()
        first_element = list_column_types[0]
        tfs_file.write( ("$ {0:>"+str(self.__column_width-2)+"} ").format(first_element))
        remain = list_column_types[1:]
        tfs_file.write(" ".join(format_for_string.format(entry) for entry in remain))
        tfs_file.write("\n")

        # Write table lines
        for table_line in self.__tfs_table.get_data_rows():
            tfs_file.write(" ".join(format_for_string.format(str(entry)) for entry in table_line))
            tfs_file.write("\n")


    def __write_unformatted_table(self, tfs_file):
        tfs_file.write("* ")
        tfs_file.write(" ".join(self.__tfs_table.get_column_names()))

        tfs_file.write("\n$ ")

        tfs_file.write(" ".join(self.__tfs_table.get_column_data_types()))

        tfs_file.write("\n")

        for row in self.__tfs_table.get_data_rows():
            tfs_file.write(" ".join(str(entry) for entry in row))
            tfs_file.write("\n")


class _TfsHeaderLine(object):
    """ Abstract class which represents a header line.

        Subclasses: _TfsDescriptor, _TfsComment
    """
    def __init__(self):
        pass
    def get_line_as_string_with_newline(self):
        raise NotImplementedError()


class _TfsDescriptor(_TfsHeaderLine):
    """ Represents a descriptor in a TFS file. E.g.: '@ SomeText %s "Some Text"' """

    def __init__(self, name, value, tfs_data_type = None):
        super(_TfsDescriptor, self).__init__()
        if tfs_data_type is None:
            tfs_data_type = _TfsDataType.get_new_string_instance()

        self.__name = ""
        self.__tfs_data_type = None
        self.__value = ""

        self.set_name(name)
        self.set_tfs_data_type(tfs_data_type)
        self.set_value(value)

    def set_name(self, name_of_descriptor_as_string):
        if isinstance(name_of_descriptor_as_string, str):
            self.__name = name_of_descriptor_as_string
        else:
            raise ValueError(name_of_descriptor_as_string + "is not a string")

    def set_value(self, value_of_descriptor):
        if self.get_tfs_data_type().is_value_valid(value_of_descriptor):
            if _TfsDataType.TYPE_STRING == self.get_tfs_data_type_as_string():
                value_of_descriptor = '"'+value_of_descriptor+'"'
            self.__value = str(value_of_descriptor)
        else:
            raise ValueError(str(value_of_descriptor) + " does not correspond to tfs_data_type: "+
                             self.get_tfs_data_type_as_string())

    def set_tfs_data_type(self, tfs_data_type):
        if isinstance(tfs_data_type, _TfsDataType):
            self.__tfs_data_type = tfs_data_type
        else:
            raise ValueError(str(tfs_data_type) + " is not an instance of _TfsDataType")

    def get_name(self):
        return self.__name

    def get_tfs_data_type(self):
        return self.__tfs_data_type

    def get_tfs_data_type_as_string(self):
        return self.get_tfs_data_type().get_type_as_string()

    def get_value(self):
        return self.__value

    def get_line_as_string_with_newline(self):
        return '@ {0} {1} {2}\n'.format(self.get_name(),
                                          self.get_tfs_data_type_as_string(),
                                          self.get_value())


class _TfsComment(_TfsHeaderLine):

    def __init__(self, comment):
        super(_TfsComment, self).__init__()
        self.__comment = ""
        self.__set_comment(comment)

    def __set_comment(self, comment_as_string):
        if isinstance(comment_as_string, str):
            self.__comment = comment_as_string
        else:
            raise ValueError(str(comment_as_string) + " is not a string")

    def get_line_as_string_with_newline(self):
        return "# " + self.__comment + "\n"


class _TfsDataType:
    """ This class represents a data type for the descritpors and columns of the table of a TFS file.
        The following types are implemented: String, Float
        Usage::
            corresponding_type = _TfsDataType.get_type_from_string(type_as_string)
            string_type = _TfsDataType.get_new_string_instance()
            float_type = _TfsDataType.get_new_float_instance()

    """
    TYPE_STRING = "%s"
    TYPE_FLOAT = "%le"
    TYPE_INVALID = None

    def __init__(self):
        self.__type = _TfsDataType.TYPE_INVALID

    def set_type(self, tfs_type):
        all_types = [_TfsDataType.TYPE_STRING, _TfsDataType.TYPE_FLOAT]
        if tfs_type not in all_types:
            raise ValueError("Invalid type" + str(tfs_type))
        else:
            self.__type = tfs_type

    @staticmethod
    def get_type_from_string(type_as_string):
        if _TfsDataType.TYPE_STRING == type_as_string:
            return _TfsDataType.get_new_string_instance()
        elif _TfsDataType.TYPE_FLOAT == type_as_string:
            return _TfsDataType.get_new_float_instance()
        else:
            raise ValueError("Type in string not recognized: "+type_as_string)

    @staticmethod
    def get_new_string_instance():
        tfs_type = _TfsDataType()
        tfs_type.set_type(_TfsDataType.TYPE_STRING)
        return tfs_type

    @staticmethod
    def get_new_float_instance():
        tfs_type = _TfsDataType()
        tfs_type.set_type(_TfsDataType.TYPE_FLOAT)
        return tfs_type

    def get_type_as_string(self):
        return self.__type

    def is_value_valid(self, value):
        if _TfsDataType.TYPE_STRING == self.__type:
            return isinstance(value, str)
        elif _TfsDataType.TYPE_FLOAT == self.__type:
            return Utilities.math.can_str_be_parsed_to_number(value)
        else:
            raise TypeError("Type of _TfsDataType is unknown: "+ str(self.__type))

    def __str__(self):
        return self.__type

    def __get_type_as_string(self):
        return self.__str__()


class _TfsTable(object):
    """ Represents the table in a TFS file. """

    def __init__(self, tfs_file_writer):
        self.__tfs_file_writer = tfs_file_writer
        self.__num_of_columns = 0
        self.__list_of_column_data_types = []
        self.__list_of_column_names = []
        self.__list_of_table_rows = []

    def add_column_names(self, list_names):
        """
        :param list list_names: List of strings representing the names of the TFS columns. Without '*'.
        """
        if self.__column_data_types_are_set():
            if self.__length_is_not_equal_to_already_set_length(list_names):
                raise AttributeError("Column number is set already but the length of the given list"+
                                        " does not match.("+self.__tfs_file_writer.get_file_name()+")")
        else:
            self.__num_of_columns = len(list_names)

        self.__list_of_column_names.extend(list_names)

    def __length_is_not_equal_to_already_set_length(self, a_list):
        return self.__num_of_columns != len(a_list)

    def __column_data_types_are_set(self):
        return 0 != len(self.__list_of_column_data_types)

    def are_column_names_and_types_are_set(self):
        return self.__column_names_are_set() and self.__column_data_types_are_set()

    def is_empty(self):
        return 0 == len(self.__list_of_table_rows)

    def add_column_datatypes(self, list_datatypes):
        """
        Adds the list of column data types to the table header.
        If the number of columns is determined already(e.g. by adding column names) and the
        length of list_datatypes does not match the number of columns a TypeError will be raised.

        :param list list_datatypes: Containing the data type(%s, %le) of the columns. Without prefix '$'
        """
        if self.__column_names_are_set():
            if self.__length_is_not_equal_to_already_set_length(list_datatypes):
                raise AttributeError("Column number is set already but the length of the given list"+
                                        " does not match.("+self.__tfs_file_writer.get_file_name()+")")
        else:
            self.__num_of_columns = len(list_datatypes)

        self.__list_of_column_data_types.extend(list_datatypes)

    def __column_names_are_set(self):
        return 0 != len(self.__list_of_column_names)


    def add_table_row(self, list_row_entries):
        """
        Adds the entries of one row to the table data.
        """
        if not (self.__column_names_are_set() and self.__column_data_types_are_set()):
            raise TypeError("Before filling the table, set the names and datatypes("+
                                self.__tfs_file_writer.get_file_name()+").")
        else:
            if self.__num_of_columns != len(list_row_entries):
                raise TypeError("Number of entries does not match the column number of the table.("+
                                self.__tfs_file_writer.get_file_name()+")")

        self.__list_of_table_rows.append(list_row_entries)

    def get_column_names(self):
        return self.__list_of_column_names

    def get_column_data_types(self):
        return self.__list_of_column_data_types

    def get_data_rows(self):
        return self.__list_of_table_rows

