'''
.. module: Utilities.math

Created on 20 Aug 2013

math holds some mathematical functions

.. moduleauthor:: vimaier
'''
from __future__ import division # To get float division instead of integer division

import numpy

def arithmetic_mean(numbers):
    return numpy.mean(numbers)

def root_mean_square(numbers):
    return numpy.linalg.norm(numbers) / numpy.sqrt(len(numbers))

def standard_deviation(numbers):
    return numpy.std(numbers)

def can_str_be_parsed_to_number(str_to_test):
    try:
        float(str_to_test)
        return True
    except ValueError:
        return False