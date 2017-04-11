from __future__ import print_function
from collections import OrderedDict
import sys
import os
import logging
import pandas
import numpy as np
import tfs_file_writer

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())

HEADER = "@"
NAMES = "*"
TYPES = "$"

ID_TO_TYPE = {
    "%s": np.str,
    "%bpm_s": np.str,
    "%le": np.float64,
    "%hd": np.int,
    "%d": np.int,
}

TYPE_TO_ID = {
    np.str: "%s",
    np.float64: "%le",
    float: "%le",
    np.int: "%d",
}


def read_tfs(tfs_path):
    """
    Parses the TFS table present in tfs_path and returns a Pandas DataFrame,
    adding a "headers" attribute:
    data_frame.headers[header_name] = header_value
    """
    LOGGER.debug("Reading path: " + tfs_path)
    headers = OrderedDict()
    column_names = column_types = None
    rows_list = []
    with open(tfs_path, "r") as tfs_data:
        for line in tfs_data:
            parts = line.split()
            if len(parts) == 0:
                continue
            if parts[0] == HEADER:
                headers[parts[1]] = _parse_header(
                    parts[2], " ".join(parts[3:]))
            elif parts[0] == NAMES:
                LOGGER.debug("Setting column names.")
                column_names = np.array(parts[1:])
            elif parts[0] == TYPES:
                LOGGER.debug("Setting column types.")
                column_types = _compute_types(parts[1:])
            else:
                if column_names is None:
                    raise TfsFormatError("Column names have not been set.")
                if column_types is None:
                    raise TfsFormatError("Column types have not been set.")
                parts = [part.strip('"') for part in parts]
                rows_list.append(parts)
    return _create_data_frame(column_names, column_types, rows_list, headers)


def write_tfs(data_frame, headers_dict, tfs_path):
    """
    Writes the Pandas DataFrame data_frame into tfs_path with the headers_dict
    as headers dictionary. If you want to keep the order of the headers, use
    collections.OrderedDict.
    """
    tfs_name = os.path.basename(tfs_path)
    tfs_dir = os.path.dirname(tfs_path)
    LOGGER.debug("Attempting to write file: " + tfs_name + " in " + tfs_dir)
    tfs_writer = tfs_file_writer.TfsFileWriter(tfs_name, outputpath=tfs_dir)
    column_names = _get_column_names(data_frame)
    column_types = _get_column_types(data_frame)
    for head_name, head_value in headers_dict.iteritems():
        if type(head_value) is str:
            tfs_writer.add_string_descriptor(head_name, head_value)
        else:
            tfs_writer.add_float_descriptor(head_name, head_value)
    tfs_writer.add_column_names(column_names)
    tfs_writer.add_column_datatypes(column_types)
    for _, row in data_frame.iterrows():
        tfs_writer.add_table_row(row)
    tfs_writer.write_to_file()


def add_coupling(data_frame):
    """
    Computes the coupling for data_frame adding 3 columns to it:
    - f1001
    - f1010
    - gamma
    """
    df = data_frame
    j = np.array([[0., 1.],
                  [-1., 0.]])
    rs = np.reshape(df.as_matrix(columns=["R11", "R12",
                                          "R21", "R22"]),
                    (len(df), 2, 2))
    cs = np.einsum("ij,kjn,no->kio",
                   -j, np.transpose(rs, axes=(0, 2, 1)), j)
    cs = np.einsum("k,kij->kij", (1 / np.sqrt(1 + np.linalg.det(rs))), cs)

    g11a = 1 / np.sqrt(df.loc[:, "BETX"])
    g12a = np.zeros(len(df))
    g21a = df.loc[:, "ALFX"] / np.sqrt(df.loc[:, "BETX"])
    g22a = np.sqrt(df.loc[:, "BETX"])
    gas = np.reshape(np.array([g11a, g12a,
                               g21a, g22a]).T,
                     (len(df), 2, 2))

    ig11b = np.sqrt(df.loc[:, "BETY"])
    ig12b = np.zeros(len(df))
    ig21b = -df.loc[:, "ALFY"] / np.sqrt(df.loc[:, "BETY"])
    ig22b = 1. / np.sqrt(df.loc[:, "BETY"])
    igbs = np.reshape(np.array([ig11b, ig12b,
                                ig21b, ig22b]).T,
                      (len(df), 2, 2))
    cs = np.einsum("kij,kjl,kln->kin", gas, cs, igbs)
    gammas = np.sqrt(1 - np.linalg.det(cs))
    data_frame.loc[:, "gamma"] = gammas
    data_frame.loc[:, "f1001"] = ((cs[:, 0, 0] + cs[:, 1, 1]) * 1j +
                                  (cs[:, 0, 1] - cs[:, 1, 0])) / 4 / gammas
    data_frame.loc[:, "f1010"] = ((cs[:, 0, 0] - cs[:, 1, 1]) * 1j +
                                  (-cs[:, 0, 1]) - cs[:, 1, 0]) / 4 / gammas


class TfsFormatError(Exception):
    """
    Raised when wrong format is detected in the TFS file.
    """
    pass


def _create_data_frame(column_names, column_types, rows_list, headers):
    data_frame = pandas.DataFrame(data=np.array(rows_list),
                                  columns=column_names)
    _assign_column_types(data_frame, column_names, column_types)
    data_frame.headers = headers
    return data_frame


def _assign_column_types(data_frame, column_names, column_types):
    names_to_types = dict(zip(column_names, column_types))
    for name, type_f in names_to_types.iteritems():
        data_frame[name] = data_frame[name].astype(type_f)


def _compute_types(str_list):
    return [_id_to_type(string) for string in str_list]


def _parse_header(type_str, value_str):
    return _id_to_type(type_str)(value_str.strip('"'))


def _id_to_type(type_str):
    try:
        return ID_TO_TYPE[type_str]
    except KeyError:
        if type_str.startswith("%") and type_str.endswith("s"):
            return str
        _raise_unknown_type(type_str)


def _type_to_id(type_f):
    try:
        return TYPE_TO_ID[type_f]
    except KeyError:
        return "%s"


def _get_column_names(data_frame):
    return data_frame.columns.values


def _get_column_types(data_frame):
    types = []
    for column in data_frame.columns:
        type_f = data_frame[column].dtype
        types.append(_type_to_id(type_f.type))
    return types


def _raise_unknown_type(name):
    raise TfsFormatError("Unknown data type: " + name)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    LOGGER.debug(read_tfs(sys.argv[1]))
