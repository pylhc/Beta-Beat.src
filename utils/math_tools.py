"""
Module utils.math_tools
------------------------

Mathematical helper functions for everyone.
"""
import numpy as np


def get_next_scientific_exponent(number):
    """ Returns the next exponent of 10 which is a multiple of 3, so that
        the coefficient is <100.

        Args:
            number: number to convert

        Returns:
            Tupel of coefficient and exponent.
    """
    exp = np.log10(abs(number))
    exp = int(np.ceil(exp - 2 * (exp < 0))/3.)*3
    coeff = number * 10**(- exp)
    return coeff, exp


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
        Source: https://stackoverflow.com/a/23535934/5609590

        Args:
            arr: numpy array
    """
    arr = np.ma.array(arr).compressed()  # in case of MaskedArray.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

