'''
.. module: Utilities.math

Created on 20 Aug 2013

math holds some mathematical functions.
This functions aims at understandable and readable code.
Use this functions to make you code more readable.

Usage::

    import Python_Classes4MAD.metaclass
    import Utilities.math

    twiss_file = Python_Classes4MAD.metaclass.twiss("twiss.C.dat")
    column = getattr(twiss_file, "BETX", None)
    if not column is None:
        print Utilities.math.arithmetic_mean(column)
        print Utilities.math.root_mean_square(column)
        print Utilities.math.standard_deviation(column)


.. moduleauthor:: vimaier
'''
from __future__ import division # To get float division instead of integer division

import numpy

def arithmetic_mean(numbers):
    """ Calculates the arithmetic mean """
    return numpy.mean(numbers)

def root_mean_square(numbers):
    """ Calculates the root mean square"""
    return numpy.linalg.norm(numbers) / numpy.sqrt(len(numbers))

def standard_deviation(numbers):
    """ Calculates the standard deviation """
    return numpy.std(numbers)

def can_str_be_parsed_to_number(str_to_test):
    """ Checks if the given string is a float.
    "str" --> False,
    "23.45" --> True,
    "2E3" --> True,
    "2e-3" --> True
    """
    try:
        float(str_to_test)
        return True
    except ValueError:
        return False

