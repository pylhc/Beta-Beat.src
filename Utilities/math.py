'''
Created on 20 Aug 2013

@author: vimaier
'''
from __future__ import division # To get float division instead of integer division

import numpy

def arithmetic_mean(numbers):
    return sum(numbers)/len(numbers)
#     sum = 0.0
#     for i in xrange(0, len(numbers)):
#         sum += numbers[i]
#     return sum / len(numbers)

def root_mean_square(numbers):
    return numpy.sqrt(sum(n*n for n in numbers)/len(numbers))
#     sum = 0.0
#     for i in xrange(0, len(numbers)):
#         sum += numbers[i] * numbers[i]
#     return numpy.sqrt(sum / len(numbers))


def standard_deviation(numbers):
    return numpy.std(numbers)