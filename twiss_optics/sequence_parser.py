import re
from Utilities import logging_tools as logtool
from Utilities.contexts import timeit
from Utilities import tfs_pandas as tfs


LOG = logtool.get_logger(__name__)


def parse_sequence_file(sequence_file):
    """
    Create Mapping from knobs to magnets based on sequence_file.
     assumptions:
     - we don't care about the N in KNL
     - linearity: knobs do not multiply with each other (will result in zeros)
     - the knob-name is final (i.e. it is not reassigned by ':=' somewhere else)
     - the knob-name does not start with "l." (reserved for length multiplier)
     - if there is a '(0.####)' at the end of the KNL it comes from splitting the element

    :param sequence_file: Saved Sequence file
    :return: TFS_DF containing a matrix magnets vs. knobs containing the knob-coefficients
    """
    LOG.info("Parsing file '{:s}'".format(sequence_file))

    with timeit(lambda t:
                LOG.debug("  File parsed in {:f}s".format(t))):
        m_t_k = tfs.TfsDataFrame()
        with open(sequence_file, 'r') as f_seq:
            for line in f_seq:
                match = re.search(r"knl:=\{(?P<knl>[^\}]+)\}", line)
                if match is not None:
                    magnet = re.match(r"([\w]|\.(?!\.))*", line).group(0)
                    if magnet not in m_t_k.index:
                        knls = match.group('knl').split(',')
                        for n, knl in enumerate(knls):
                            try:
                                float(knl)  # check could also be "len(knl) > 1"
                            except ValueError:
                                # if one ever needs it: n is available as well!
                                k_dict = _parse_knl_element(knl)
                                m_t_k = m_t_k.append(
                                    tfs.TfsDataFrame([k_dict.values()],
                                                     index=[magnet],
                                                     columns=k_dict.keys()
                                                     ))
    return m_t_k.fillna(0)


def _parse_knl_element(knl):
    """ Finds knob-names and coefficients in string 'knl' """
    knl = knl.replace(" ", "")  # to simplify regex (no \s*)

    # replace those splitting values by 1
    num_split = re.search(r"\((0\.\d*)\)+$", knl)
    if num_split is not None:
        knl = _replace_chars_at(knl, num_split.start(1), num_split.end(1), "1")

    # replace element-length multiplier by 1
    el_length = re.search(r"l\.[\w.]+", knl)
    if el_length is not None:
        knl = knl.replace(el_length.group(0), "1")

    # get knob names
    knobs = re.findall(r"[a-zA-Z][\w.]*", knl)
    d = {}
    for knob in knobs:
        # replace current knob by 1 and all others by 0
        knl_temp = knl
        for other in knobs:
            knl_temp = knl_temp.replace(other, "1" if knob == other else "0")

        # evaluate coefficient (let python sort out parenthesis, order etc.)
        d[knob] = eval(knl_temp)

    return d


def _replace_chars_at(s, start, end, replace=""):
    """ kind of s[start:end] = replace, where replace can have any length """
    return s[:start] + replace + s[end:]


if __name__ == '__main__':
    import os
    df = parse_sequence_file(os.path.join("tests", "lhcb1.seq"))
    print df