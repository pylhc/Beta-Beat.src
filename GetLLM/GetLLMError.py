#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:25:51 2017

@author: awegsche
"""

import __init__
    
class GetLLMError(Exception):
    """
    Raised when an error occured in GetLLM or one of its algorithms, like
    missing input files.
    """
    pass