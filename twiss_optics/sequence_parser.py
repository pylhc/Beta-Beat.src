import re
from Utilities import logging_tools as logtool
from Utilities.contexts import timeit
from Utilities import tfs_pandas as tfs


LOG = logtool.get_logger(__name__)


"""
=============================   Main   =============================
"""


def parse_sequence_file(sequence_file):
    """
    Create Mapping from variables to magnets based on sequence_file.
     assumptions:
     - Sequence is saved in a way, that magnets definitions contain 'knl:={}' or 'K#:='
     - we don't care about the N in KNL (could be found quickly, though)
     - linearity: variables do not multiply with each other (will result in zeros)
     - the variable-name is final (i.e. it is not reassigned by ':=' somewhere else)
     - the variable-name does not start with "l." (reserved for length multiplier)
     - the length multiplier is given by l.<magnettype>
     - apart from the length, there are no other named multipliers (numbers are fine)

    :param sequence_file: Saved Sequence file
    :return: TFS_DF containing a matrix -  magnets(index) vs. variables(columns)
             containing the variable-coefficients
    """
    LOG.info("Parsing file '{:s}'".format(sequence_file))

    with timeit(lambda t:
                LOG.debug("  File parsed in {:f}s".format(t))):
        magnet_strings, length_constants = _read_file(sequence_file)
        var_to_mag = _build_dataframe(magnet_strings, length_constants)
    return var_to_mag


"""
=============================   Read Sequence   =============================
"""


def _read_file(sequence_file):
    """ Read the file and return the magnet definitions and magnet lengths """
    length_constants = {}
    magnet_strings = {}
    with open(sequence_file, 'r') as f_seq:
        for line in f_seq:
            var_and_value = _find_length(line)
            if var_and_value is not None:
                length_constants[var_and_value[0]] = var_and_value[1]
            else:
                var_and_value = _find_strength(line)
                if var_and_value is not None:
                    magnet_strings[var_and_value[0]] = var_and_value[1]
    return magnet_strings, length_constants


def _find_length(line):
    """ Search for length variable in line """
    match = re.match(r"const\s(l\.[^;]+)", line)
    if match is not None:
        eq = match.group(1).replace(" ", "")
        name_and_value = eq.split("=")
        return name_and_value
    else:
        return None


def _find_strength(line):
    """ Search for magnet strengths in line
        [If one ever needs n exchange '\d' for '(?P<N>\d) and use enumerate(knls)]
    """
    match = re.search(r",\s*k(nl:=\{(?P<knl>[^\}]+)\}|\d:=(?P<k>[^,]+))", line)
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
            length = "l." + re.search(r":(?!=)\s*([^,]+)", line).group(1)
            knl = "({kn:s}) * {l:s}".format(kn=match.group("k"), l=length)
            return magnet, knl.replace(" ", "")
    else:
        return None


"""
=============================   Build Mapping   =============================
"""


def _build_dataframe(magnet_strings, length_constants):
    """ Build the data frame representing the mapping variables to magnets """
    var_to_mag = tfs.TfsDataFrame()
    for magnet, value_string in magnet_strings.iteritems():
        k_dict = _eval_strength(value_string, length_constants)
        var_to_mag = var_to_mag.append(
            tfs.TfsDataFrame([k_dict.values()],
                             index=[magnet],
                             columns=k_dict.keys()
                             ))
    return var_to_mag.fillna(0)


def _eval_strength(k_string, length_constants):
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


if __name__ == '__main__':
    import os
    df = parse_sequence_file(os.path.join("tests", "lhcb1_tmp.seq"))
    print df.loc['mq.26r5.b1', df.loc['mq.26r5.b1', :] != 0]
