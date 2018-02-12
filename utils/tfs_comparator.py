import __init__  # @UnusedImport
import sys
from Python_Classes4MAD.metaclass import twiss
from utils import logging_tools

LOG = logging_tools.get_logger(__name__)


def compare_tfs(tfs_file_1, tfs_file_2, max_relative_error=1e-9, check_globals=False,
                check_same_columns=False, ignore_names=[]):
    """ Module to compare tfs files.

    It will return false and output the errors it found to
    the error output and true if the files are equal.

    Args:
        tfs_file_1: The first TFS file, either a path o a metaclass.twiss object
        tfs_file_2: The second TFS file.
        max_relative_error: Maximum relative error allowed in the numeric values.
        check_globals: If true, it will also compare the value of the global parameters.
        check_same_columns: If true, both files must have the same number and names of columns
        ignore_names: List of columns or global parameters names that will be
                      ignored while comparing.
    """
    ignore_names.append("indx")
    ignore_names.append("keys")
    ignore_names.append("filename")
    result = True
    if not isinstance(tfs_file_1, twiss):
        try:
            tfs_file_1 = twiss(tfs_file_1)
        except Exception:
            LOG.error("File '" + str(tfs_file_1) + "' is unreadable")
            return False
    if not isinstance(tfs_file_2, twiss):
        try:
            tfs_file_2 = twiss(tfs_file_2)
        except Exception:
            LOG.error("File '", str(tfs_file_2), "' is unreadable")
            return False

    attr_list = ([], [])
    tfs_files = (tfs_file_1, tfs_file_2)
    [attr_list[0].append(col_name) for col_name in dir(tfs_file_1) if not col_name.startswith("_")]
    [attr_list[1].append(col_name) for col_name in dir(tfs_file_2) if not col_name.startswith("_")]

    if check_same_columns and len(attr_list[0]) != len(attr_list[1]):
        LOG.error("The files don't have the same number of columns or global parameters.")
        return False

    longest_list = 0
    shortest_list = 1
    if len(attr_list[1]) > len(attr_list[0]):
        longest_list = 1
        shortest_list = 0

    for item_name in attr_list[longest_list]:
        if item_name not in ignore_names:
            if item_name not in attr_list[shortest_list] and check_same_columns:
                LOG.error("The column '" + item_name + "' is not present in both files.")
                return False
            elif item_name in attr_list[shortest_list]:
                next_value = getattr(tfs_files[shortest_list], item_name)
                try:
                    iter(next_value)
                    if type(next_value) is str:
                        raise TypeError
                    result = _check_column_values(item_name, next_value,
                                                  getattr(tfs_files[longest_list], item_name),
                                                  max_relative_error)
                except TypeError:
                    if check_globals:
                        if not callable(next_value):
                            result = _check_single_value(item_name, next_value,
                                                         getattr(tfs_files[longest_list],
                                                                 item_name),
                                                         max_relative_error)
            if not result:
                break
    return result


def _check_column_values(column_name, column1, column2, max_relative_error):
    if len(column1) != len(column2):
        LOG.error("The columns named '" + column_name + "'  do not have the same size.")
        return False
    line = 1
    result = True
    for value1, value2 in zip(column1, column2):
        new_result = _check_single_value(column_name, value1, value2, max_relative_error, line)
        result = result and new_result
        line += 1
    if not result:
        LOG.error("Differences found in column '" + column_name + "'. Aborting.")
    return result


def _check_single_value(value_name, value1, value2, max_relative_error, line=-1):
    if type(value1) != type(value2):
        line_string = "" if line == -1 else "(row " + str(line) + ")"
        LOG.error("Values of different type found on " + value_name + line_string)
        return False
    try:
        result = _compare_numeric_values(float(value1), float(value2), max_relative_error)
        if not result:
            line_string = "" if line == -1 else "(row " + str(line) + ")"
            LOG.error("Difference found on " +  value_name + line_string +
                      ":" + str(value1) + "!=" + str(value2))
            return False
    except (ValueError, TypeError):
        result = str(value1) == str(value2)
        if not result:
            line_string = "" if line == -1 else "(row " + str(line) + ")"
            LOG.error("Difference found on" + value_name + line_string +
                      ":" + str(value1) + "!=" + str(value2))
            return False
    return result


def _compare_numeric_values(value1, value2, max_relative_error):
    bigest = value1
    if abs(value2) > abs(value1):
        bigest = value2
    if bigest == 0:
        return True
    return abs((value1 - value2) / bigest) < max_relative_error


if __name__ == "__main__":
    compare_tfs("twiss_coll.dat", "twiss_coll2.dat")
