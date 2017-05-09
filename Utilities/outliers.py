r'''
.. module: Utilities.outliers

Created on 08/05/17

:author: Lukas Malina  

It filter the array of values which are meant to be constant 
or a linear function of the other array if that is provided
Returns a filter mask for the original array
'''

import numpy as np
from scipy.stats import t


def get_filter_mask(data, x_data=None, limit=0.0, niter=20, nsig=None): # limit for not being cleaned
    mask = np.ones_like(data, dtype=bool)
    if x_data is not None:
        if not len(data)==len(x_data):
            print "datasets are not equally long, will not use the second dataset"
            x_data=None
      
    if nsig is None:
        ns=_get_significance_cut_from_length(len(data))
    else:
        ns = nsig
    prevlen = np.sum(mask) + 1 # to fulfill the condition for the first iteration  
    iteration = 0    

    while (np.sum(mask) < prevlen) and (np.sum(mask) > 2) and (iteration < niter):
        prevlen = np.sum(mask)
        iteration+=1      
        y1, y1_orig = _get_data(mask,data) if x_data is None else _get_data_without_slope(mask, x_data, data)
        av, st = _get_moments(y1)
        mask = np.abs(y1_orig - av) < np.max([limit, ns*st])
        if nsig is None:
            ns=_get_significance_cut_from_length(len(data))
        
    return mask


def _get_moments(data):
    return np.mean(data), np.std(data)


def _get_data_without_slope(mask, x, y):
    m,b=np.polyfit(x[mask],y[mask],1)
    return y[mask]-b-m*x[mask], y-b-m*x


def _get_data(mask,data):
    return data[mask], data

# set the sigma cut, that expects 1 value to be cut if it is sample of normal distribution
def _get_significance_cut_from_length(length): 
    return t.ppf([1 - 0.5/float(length)],length)

