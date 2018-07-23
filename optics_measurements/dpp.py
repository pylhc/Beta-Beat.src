"""
.. module: dpp

Created on 2017

:author: Elena Fol

Arranges \frac{\Delta p}{p} of given files
"""

import logging
import numpy as np

DPP_TOLERANCE = 0.0001
LOGGER = logging.getLogger(__name__)


def arrange_dpp(list_of_tfs):
    """
    Grouping of dpp-values in the given linx,liny-files and computing new values
    """
    list_of_tfs_arranged = []
    for tfs_file in list_of_tfs:
        if "DPP" not in tfs_file.headers:
            tfs_file.headers["DPP"] = 0.0
    if len(list_of_tfs) == 1:
        only_dpp = list_of_tfs[0].headers["DPP"]
        if np.abs(only_dpp) > DPP_TOLERANCE:
            LOGGER.warn(
                'It looks like the file you are analyzing has too '
                'high momentum deviation {}. Optics parameters might '
                'be wrong.'.format(only_dpp)
            )
        list_of_tfs[0].headers["DPP"] = 0.0
        return list_of_tfs
    dpp_values = [tfs_file.DPP for tfs_file in list_of_tfs]
    if 0.0 in dpp_values:
        LOGGER.warn('Exact 0.0 found, the dp/p values are probably already grouped.')
        return list_of_tfs
    closest_to_zero = np.argmin(np.absolute(dpp_values))
    ordered_indices = np.argsort(dpp_values)
    ranges = _compute_ranges(list_of_tfs, ordered_indices)
    offset_range = _find_range_with_element(ranges, closest_to_zero)
    offset_dpps = _values_in_range(offset_range, dpp_values)
    LOGGER.debug("dp/p closest to zero is {}".format(dpp_values[closest_to_zero]))
    zero_offset = np.mean(offset_dpps)
    LOGGER.debug("Detected dpp differences, aranging as: {0}, zero offset: {1}."
                 .format(ranges, zero_offset))
    for idx in range(len(dpp_values)):
        range_to_use = _find_range_with_element(ranges, idx)
        dpps_from_range = _values_in_range(range_to_use, dpp_values)
        range_mean = np.mean(dpps_from_range)
        list_of_tfs[idx].headers["DPP"] = range_mean - zero_offset
        list_of_tfs_arranged.append(list_of_tfs[idx])
    return list_of_tfs_arranged


def _values_in_range(range_to_use, dpp_values):
    dpps_from_range = []
    for dpp_idx in range_to_use:
        dpps_from_range.append(dpp_values[dpp_idx])
    return dpps_from_range


def _find_range_with_element(ranges, element):
    range_with_element = None
    for dpp_range in ranges:
        if element in dpp_range:
            range_with_element = dpp_range
    return range_with_element


def _compute_ranges(list_of_tfs, ordered_indices):
    list_of_ranges = []
    last_range = None
    for idx in ordered_indices:
        if (list_of_ranges and
                _is_in_same_range(list_of_tfs[last_range[0]].DPP,
                                  list_of_tfs[idx].DPP)):
            last_range.append(idx)
        else:
            new_range = [idx]
            list_of_ranges.append(new_range)
            last_range = new_range
    return list_of_ranges


def _is_in_same_range(a, b):
    return a + DPP_TOLERANCE >= b >= a - DPP_TOLERANCE
