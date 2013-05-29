'''
Created on 27 May 2013

@author: vimaier

@version: 1.0.1

This module serves as an example for good code style.

Every module should have a short/abstract description of what it does, how to use it, where it is 
used and any other relevant information.
Comments in this document which start with '##' are explanations

The module name should contain only lower letters and words should be divided by an underscore '_'.
Good examples: 'segment_by_segment.py','io_utils.py', 'get_phase.py'
Bad examples: 'SegmentBySegment.py', 'IoUtils.py', 'getPhase.py'

There should also be a change history with an entry for every important change.

Change history:
 - 1.0.1, vimaier, 27th May 2013: 
    Added totally new feature x in method y.
 - <version>, <author>, <date>:
    <description>
'''

## Divide the imports into three sections
## System imports(Python Standard Library)
import sys
import os

## Third party libraries
import numpy as np

## Own modules or libraries
import metaclass


## 'import * from x' should never be used!


class ExampleClass(object):
    ''' Document the purpose and usage of the class with docstrings. '''
    
    def __init__(self, real_part=0.0, imag_part=0.0):
        '''Constructor
        
        :Parameters:
            'real_part': float
                Real part of a complex number.
            'imag_part': float
                Imaginary part of a complex number.
        '''
        ## Define all attributes(data which belongs to the object) only in the constructor
        ## Use meaningful names. 'a' is not an identifier!
        self.real_part = real_part
        self.imag_part = imag_part 
        
    def get_modulus(self):
        '''  Document the method if it is not obvious. '''
        return np.sqrt(self.real_part**2 + self.imag_part**2)
    
    def do_awesome_stuff(self, awesome_num, awesome_vector, num_of_rounds):
        '''
        Document complex methods and functions. Since Python is not a strong typed language
        describe which parameters a function assume.
        
        :Parameters:
            'awesome_num': float
                An awesome number which will be added to awesome_vector.
            'awesome_vector': array
                A one dimensional numpy array.
            'num_of_rounds': int
                Defines the number of rounds which adds an awesome number to awesome_vector.
                
        :Return: array
            An awesome one dimensional array.     
        '''
        # Since awesome_vector is a reference we will copy the array
        awesome_vector = np.array(awesome_vector)        
        
        ## Use short identifier if it improves readability of math functions but describe them 
        ## before!
        re = self.real_part
        im = self.imag_part
        z = self.get_modulus()
        
        # another awesome number
        x = 0.0

        # Calculate awesome stuff        
        for i in range(num_of_rounds):
            x = (z+re-im)*i + awesome_num
            awesome_vector += x
            
        return awesome_vector 
        
# END class ExampleClass ---------------------------------------------------------------------------
        
    
