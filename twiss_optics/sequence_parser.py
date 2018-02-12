"""
    Functions to create mapping from variables to magnets based on sequence_file.
     Assumptions:

     - magnets named "mb.***" are defined by their bending radius and are hence immutable!
     - Sequence is saved in a way, that magnets definitions contain 'knl:={}', 'ksl:={}' or 'K#:='
     - There is only one value in each knl/ksl array.
     - zero-offset, i.e. no fixed number summation (see hint below)
     - linearity, i.e. variables do not multiply with each other (will result in zeros)
     - the variable-name is final (i.e. it is not reassigned by ':=' somewhere else)
     - the variable-name does not start with "l." (reserved for length multiplier)
     - the length multiplier is given by l.<magnettype>
     - apart from the length, there are no other named multipliers (numbers are fine)

     If a magnet is redefined, the last definition is used. (not really an assumption)

    HINT: 'magnet := #.### + var' WILL NOT RESULT IN A VALID MAPPING !!!
    (e.g. '1 + var' is indistinguishable from '2 * var', which is not what you want!)
"""

import cPickle as pickle
import os
import re
from collections import OrderedDict

import pandas as pd

from utils import logging_tools as logtool
from utils import tfs_pandas as tfs
from utils.contexts import timeit
from utils.iotools import json_dumps_readable

LOG = logtool.get_logger(__name__)

EXT = "varmap"  # Extension Standard

"""
=============================   Main   =============================
"""
DEFAULT = {
    'return': "dictionary",
    'ret_choices': ["dataframe", "dictionary"],
}


def parse_variable_mapping(seqfile_path, ret=DEFAULT['return']):
    """ Main Function for creating the variable mapping

    Args:
        seqfile_path: Saved Sequence file
        ret: return format, either "tfs" or "dictionary"

    Returns:
        Dictionary of orders containing either a
        DataFrame of magnets(index) vs. variables(columns)
        containing the variable-coefficients
        or a
        Dictionary of all variables containing magnet-coefficient Series
    """
    LOG.info("Parsing file '{:s}'".format(seqfile_path))

    _check_ret(ret)

    with timeit(lambda t:
                LOG.debug("  File parsed in {:f}s".format(t))):
        magnet_strings, length_constants = _read_file_for_magnets(seqfile_path)
        if ret.lower() == "dataframe":
            return _build_variable_mapping_df(magnet_strings, length_constants)

        elif ret.lower() == "dictionary":
            return _build_variable_mapping_dict(magnet_strings, length_constants)


def load_variable_mapping(seqfile_path, ret=DEFAULT['return']):
    """ Load mapping from file(s).

    Args:
        seqfile_path: Saved Sequence file
        ret: return format, either "tfs" or "dictionary"

    Returns:
        Dictionary of orders containing either a
        DataFrame of magnets(index) vs. variables(columns)
        containing the variable-coefficients
        or a
        Dictionary of all variables containing magnet-coefficient Series
    """
    _check_ret(ret)
    varmapfile_path = seqfile_path.replace(".seq", "").replace("." + EXT, "")
    if ret == "dictionary":
        full_file_path = "{f:s}.{e:s}".format(f=varmapfile_path, e=EXT)
        with open(full_file_path, "rb") as varmapfile:
            mapping = pickle.load(varmapfile)
        LOG.debug("Loaded mapping from varmap file '{:s}'.".format(full_file_path))

    elif ret == "dataframe":
        varmapfile_name = os.path.basename(varmapfile_path)
        order = []
        for f in os.listdir(os.path.dirname(varmapfile_path)):
            if f.startswith(varmapfile_name) and f.endswith(EXT):
                order += [f.replace(varmapfile_name + ".", "").replace("." + EXT, "")]

        if len(order) == 0:
            raise IOError("Could not find varmap files of scheme: '{f:s}.{o:s}.{e:s}'".format(
                f=varmapfile_path, o="(order)", e=EXT))

        mapping = dict.fromkeys(order)
        for o in order:
            full_file_path = "{f:s}.{o:s}.{e:s}".format(f=varmapfile_path, o=o, e=EXT)
            mapping[o] = tfs.read_tfs(full_file_path)
            LOG.debug("Loaded mapping from varmap file '{:s}'.".format(full_file_path))
    return mapping


def save_variable_mapping(mapping, outfile_path, format=DEFAULT['return']):
    """ Save mapping to file(s).

    Args:
        mapping: The mapping to save
        outfile_path: Output File, either extension "",".seq",".varmap" will be changed to ".varmap"
        format: mapping format, either "tfs" or "dictionary"
    """
    _check_ret(format)
    varmapfile_path = outfile_path.replace(".seq", "").replace("." + EXT, "")

    if format == "dictionary":
        full_file_path = "{f:s}.{e:s}".format(f=varmapfile_path, e=EXT)
        with open(full_file_path, "wb") as varmapfile:
            pickle.dump(mapping, varmapfile, -1)
        LOG.debug("Saved Variable mapping into file '{:s}'".format(full_file_path))

    elif format == "dataframe":
        for order in mapping:
            full_file_path = "{f:s}.{o:s}.{e:s}".format(f=varmapfile_path, o=order, e=EXT)
            tfs.write_tfs(full_file_path,
                          mapping[order],
                          save_index=True)
            LOG.debug("Saved Variable mapping into file '{:s}'".format(full_file_path))


def load_or_parse_variable_mapping(seqfile_path, ret=DEFAULT['return']):
    """ Load mapping, or parse if not found. Convenience wrapper for parse and load functions.

    Loads variable mapping from a file. If not found it will do the parsing instead and saves
    the results for later use.

    Args:
        seqfile_path: Saved Sequence file
        ret: return format, either "tfs" or "dictionary"

    Returns:
        Dictionary of orders containing either a
        DataFrame of magnets(index) vs. variables(columns)
        containing the variable-coefficients
        or a
        Dictionary of all variables containing magnet-coefficient Series

    """
    _check_ret(ret)
    try:
        mapping = load_variable_mapping(seqfile_path, ret=ret)
    except IOError:
        mapping = parse_variable_mapping(seqfile_path, ret=ret)
        try:
            save_variable_mapping(mapping, seqfile_path, format=ret)
        except IOError as e:
            LOG.warn("  IOError: {:s}.".format(e.message))
    return mapping


def varmap_variables_to_json(varmap_or_file, outfile_path=None, format=DEFAULT['return']):
    """ Saves all variable names from mapping to json file.

        The variables will be saved by their order in the file.

        Args:
            varmap_or_file: varmap or saved varmap-file (sequence file is okay as well)
            format: varmap format, either "tfs" or "dictionary"
    """
    _check_ret(format)
    LOG.debug("Converting varmap to json-file.")
    if isinstance(varmap_or_file, basestring):
        mapping = load_or_parse_variable_mapping(varmap_or_file, ret=format)
        if outfile_path is None:
            outfile_path = varmap_or_file.replace(".seq", "").replace("." + EXT, "") + "_all_list.json"
    else:
        mapping = varmap_or_file
        if outfile_path is None:
            IOError("Outputfile not given!")

    json_dict = OrderedDict.fromkeys(sorted(mapping.keys()))
    for order in mapping:
        if format == "dictionary":
            json_dict[order] = sorted(mapping[order].keys())
        elif format == "dataframe":
            json_dict[order] = sorted(mapping[order].columns.tolist())

    json_dict["all"] = sorted(list(set([var for order in json_dict for var in json_dict[order]])))

    json_dumps_readable(outfile_path, json_dict)
    LOG.debug("Variables saved to '{:s}'.".format(outfile_path))


"""
=============================   Read Sequence   =============================
"""


def _read_file_for_magnets(sequence_file):
    """ Read the file and return the magnet definitions and magnet lengths """
    LOG.debug("  Reading File")
    length_constants = {}
    magnet_strings = {}
    with open(sequence_file, 'r') as f_seq:
        for line in f_seq:
            var_and_value = _find_element_length(line)
            if var_and_value is not None:
                length_constants[var_and_value[0]] = var_and_value[1]
            else:
                var_and_value = _find_magnet_strength(line)
                if var_and_value is not None:
                    magnet_strings[var_and_value[0]] = var_and_value[1]
    return magnet_strings, length_constants


def _find_element_length(line):
    """ Search for length variable in line """
    match = re.match(r"const\s(l\.[^;]+)", line)
    if match is not None:
        eq = match.group(1).replace(" ", "")
        return eq.split("=")
    else:
        return None


def _find_magnet_strength(line):
    """ Search for magnet strengths in line
    """
    line = line.lower()
    matches = list(re.finditer(
        r",\s*k((?P<s>[ns])l:=\{(?P<knl>[^\}]+)\}|(?P<n>\d+s?):=(?P<k>[^,]+))", line))

    if len(matches) > 0:
        magnet = re.match(r"[\w.]*", line).group(0)

        knl_dict = {}
        for match in matches:
            if match.group("knl") is not None:
                skew = "S" if match.group('s') == "s" else ""
                knls = match.group('knl').split(',')
                for n, knl in enumerate(knls):
                    try:
                        float(knl)  # check could also be "len(knl) > 1"
                    except ValueError:
                        ########## HACK TO AVOID DIPOLES AS THEY ARE DEFINED BY LRAD!
                        # TODO: Find a way to change dipoles in MADX!?
                        if n == 0 and not re.search(r":\s*multipole\s*,", line):
                            return None
                        ##############################################################
                        order = "K{n:d}{skew:s}L".format(n=n, skew=skew)
                        knl_dict[order] = knl.replace(" ", "")
            else:
                if match.group("n") in ['0', '0s']:
                    # dipole strength are defined by their angles
                    knl = re.search(r"angle\s*:=\s*([^,]+)", line).group(1)
                else:
                    length = "l." + re.search(r":(?!=)\s*([^,]+)", line).group(1)
                    knl = "({kn:s}) * {l:s}".format(kn=match.group("k"), l=length)

                order = "K{:s}L".format(match.group("n").upper())
                knl_dict[order] = knl.replace(" ", "")

        return magnet, knl_dict
    else:
        return None


"""
=============================   Build Mapping   =============================
"""


def _build_variable_mapping_df(magnet_strings, length_constants):
    """ Build the data frame representing the mapping variables to magnets and return as df
        SLOW!!
    """
    LOG.debug("  Building Dataframe Mapping")
    var_to_mag = {}
    for magnet in magnet_strings:
        for order, value_string in magnet_strings[magnet].iteritems():
            if order not in var_to_mag:
                var_to_mag[order] = tfs.TfsDataFrame()

            k_dict = _eval_magnet_strength(value_string, length_constants)
            var_to_mag[order] = var_to_mag[order].append(
                tfs.TfsDataFrame([k_dict.values()],
                                 index=[magnet],
                                 columns=k_dict.keys()
                                 )).fillna(0)
    return var_to_mag


def _build_variable_mapping_dict(magnet_strings, length_constants):
    """ Build the data frame representing the mapping variables to magnets and return as dict
        Faster!
    """
    LOG.debug("  Building Dictionary Mapping")
    var_to_mag = {}
    for magnet in magnet_strings:
        for order, value_string in magnet_strings[magnet].iteritems():
            if order not in var_to_mag:
                var_to_mag[order] = {}

            k_dict = _eval_magnet_strength(value_string, length_constants)
            for var in k_dict:
                try:
                    var_to_mag[order][var].loc[magnet] = k_dict[var]
                except KeyError:
                    var_to_mag[order][var] = pd.Series(k_dict[var], index=[magnet])
    return var_to_mag


def _eval_magnet_strength(k_string, length_constants):
    """ Evaluate the magnet-strength string and return dictionary with values """
    # replace element-length multiplier by their value
    el_length = re.search(r"l\.[\w.]+", k_string)
    if el_length is not None:
        l_var = el_length.group(0)
        k_string = k_string.replace(l_var, length_constants[l_var])

    # get variable names
    variables = re.findall(r"[a-zA-Z][\w.]*", k_string)
    d = {}
    for variable in variables:
        # replace current variable by 1 and all others by 0
        knl_temp = k_string
        for other in variables:
            knl_temp = knl_temp.replace(other, "1" if variable == other else "0")

        # evaluate coefficient (let python sort out parenthesis, order etc.)
        d[variable] = eval(knl_temp)

    return d


"""
=============================   Helper   =============================
"""


def _check_ret(ret):
    """ Checks if the given 'ret' value is valid """
    if ret not in DEFAULT["ret_choices"]:
            ValueError("Return format '{ret:s}' unknown, choices: {choices:s}".format(
                ret=ret, choices=", ".join(DEFAULT["ret_choices"])
            ))


if __name__ == '__main__':
    # # Testing:
    # import os
    # seq_name = "lhcb1_tmp.seq"
    # df = load_or_parse_variable_mapping(os.path.join("tests", seq_name))
    # varmap_variables_to_json(df, seq_name)
    raise EnvironmentError("{:s} is not meant to be run as main.".format(__file__))

