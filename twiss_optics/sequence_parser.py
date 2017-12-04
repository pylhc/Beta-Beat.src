import re
from Utilities import logging_tools as logtool
from Utilities.contexts import timeit
from Utilities import tfs_pandas as tfs
import cPickle as pickle
import pandas as pd

LOG = logtool.get_logger(__name__)


"""
=============================   Main   =============================
"""


def get_variable_mapping(sequence_file, ret="dictionary"):
    """
    Create Mapping from variables to magnets based on sequence_file.
     assumptions:
     - Sequence is saved in a way, that magnets definitions contain 'knl:={}', 'ksl:={}' or 'K#:='
     - There is only one value in each knl/ksl array.
     - Dipole strengths (if not given as knl/ksl) are defined by their angle instead of k0
     - zero-offset, i.e. no fixed number summation (see hint below)
     - linearity, i.e. variables do not multiply with each other (will result in zeros)
     - the variable-name is final (i.e. it is not reassigned by ':=' somewhere else)
     - the variable-name does not start with "l." (reserved for length multiplier)
     - the length multiplier is given by l.<magnettype>
     - apart from the length, there are no other named multipliers (numbers are fine)

     HINT: 'magnet := #.### + var' WILL NOT RESULT IN A VALID MAPPING !!!
     (e.g. '1 + var' is indistinguishable from '2 * var', which is not what you want!)

    :param sequence_file: Saved Sequence file
    :param ret: return format, either "tfs" or "dictionary"
    :return: TFS of magnets(index) vs. variables(columns)
             containing the variable-coefficients
             or
             Dictionary of all variables containing magnet-coefficient Series
    """
    LOG.info("Parsing file '{:s}'".format(sequence_file))

    _check_ret(ret)

    with timeit(lambda t:
                LOG.debug("  File parsed in {:f}s".format(t))):
        magnet_strings, length_constants = _read_file_for_magnets(sequence_file)
        if ret.lower() == "dataframe":
            return _build_variable_mapping_df(magnet_strings, length_constants)
        elif ret.lower() == "dictionary":
            return _build_variable_mapping_dict(magnet_strings, length_constants)


def load_or_calc_variable_mapping(seqfile_path, ret="dictionary"):
    """
    Convenience function for loading variable mapping from a file or calculating and saving
    it if not found
    :param seqfile_path: Saved Sequence file
    :return:
    """
    _check_ret(ret)
    varmapfile_path = seqfile_path.replace(".seq", "").replace(".varmap", "") + ".varmap"

    try:
        if ret == "dictionary":
            with open(varmapfile_path, "rb") as varmapfile:
                mapping = pickle.load(varmapfile)
        elif ret == "dataframe":
            mapping = tfs.read_tfs(varmapfile_path)
        LOG.debug("Using mapping from varmap file '{:s}'.".format(varmapfile_path))

    except IOError:
        mapping = get_variable_mapping(seqfile_path, ret=ret)

        try:
            if ret == "dictionary":
                with open(varmapfile_path, "wb") as varmapfile:
                    pickle.dump(mapping, varmapfile, -1)
            elif ret == "dataframe":
                tfs.write_tfs(varmapfile_path, mapping, save_index=True)
            LOG.debug("Saved Variable mapping into file '{:s}'".format(varmapfile_path))

        except IOError as e:
            LOG.debug("Could not save mapping to '{:s}',".format(varmapfile_path))
            LOG.debug("  IOError: {:s}.".format(e.message))

    return mapping


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
        name_and_value = eq.split("=")
        return name_and_value
    else:
        return None


def _find_magnet_strength(line):
    """ Search for magnet strengths in line
        [If one ever needs n exchange '\d' for '(?P<N>\d) and use enumerate(knls)]
    """
    match = re.search(r",\s*k([ns]l:=\{(?P<knl>[^\}]+)\}|(?P<n>\d):=(?P<k>[^,]+))", line)
    if match is not None:
        magnet = re.match(r"[\w.]*", line).group(0)

        if match.group("knl") is not None:
            knls = match.group('knl').split(',')
            for knl in knls:
                try:
                    float(knl)  # check could also be "len(knl) > 1"
                except ValueError:
                    return magnet, knl.replace(" ", "")
        else:
            if match.group("n") == '0':
                # dipole strength are defined by their angles
                knl = re.search(r"angle\s*:=\s*([^,]+)", line).group(1)
            else:
                length = "l." + re.search(r":(?!=)\s*([^,]+)", line).group(1)
                knl = "({kn:s}) * {l:s}".format(kn=match.group("k"), l=length)
            return magnet, knl.replace(" ", "")
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
    var_to_mag = tfs.TfsDataFrame()
    for magnet, value_string in magnet_strings.iteritems():
        k_dict = _eval_magnet_strength(value_string, length_constants)
        var_to_mag = var_to_mag.append(
            tfs.TfsDataFrame([k_dict.values()],
                             index=[magnet],
                             columns=k_dict.keys()
                             ))
    return var_to_mag.fillna(0)


def _build_variable_mapping_dict(magnet_strings, length_constants):
    """ Build the data frame representing the mapping variables to magnets and return as dict
        Faster!
    """
    LOG.debug("  Building Dictionary Mapping")
    var_to_mag = {}
    for magnet, value_string in magnet_strings.iteritems():
        k_dict = _eval_magnet_strength(value_string, length_constants)
        for var in k_dict:
            try:
                var_to_mag[var].loc[magnet] = k_dict[var]
            except KeyError:
                var_to_mag[var] = pd.Series(k_dict[var], index=[magnet])
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
    ret_choices =["dataframe", "dictionary"]
    if ret not in ret_choices:
            ValueError("Return format '{ret:s}' unknown, choices: {choices:s}".format(
                ret=ret, choices=", ".join(ret_choices)
            ))


if __name__ == '__main__':
    pass

    # # Testing:
    # import os
    # df = get_variable_mapping(os.path.join("tests", "lhcb1_tmp.seq"))
    # print df.loc['mq.26r5.b1', df.loc['mq.26r5.b1', :] != 0]
