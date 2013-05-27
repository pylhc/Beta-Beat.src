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
import optparse
import traceback

## Third party libraries
import numpy as np

## Own modules or libraries
import metaclass


## 'import * from x' should never be used!


## Try to avoid global variables. If you really need one use the prefix 'g_' for the identifier.
g_my_global_variable = 0.0
''' Needed for ... '''

## If you need constants indicate it by upper letters and underscores
SPECIAL_TIMEOUT = 1000
''' Timeout in ms for... '''

## Never mess up the global area for calculations. Instead of this use a main function.
## You can divide a module in the following sections:
## - parse_args()-function
## - main()-function
## - helper-functions
## - main invocation


#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
    ''' Parses the arguments, checks for valid input and returns tupel '''
    parser = optparse.OptionParser()
    parser.add_option("-a", "--accel",
                    help="Which accelerator: LHCB1 LHCB2 LHCB4? SPS RHIC TEVATRON",
                    metavar="ACCEL", default="LHCB1",dest="accel")

    (options, args) = parser.parse_args()
    
    ## Check arguments to be sure to have valid input
    if options.accel not in ("LHCB1", "LHCB2", "LHCB4", "SPS", "RHIC", "TEVATRON"):
        raise ValueError("Invalid accelerator: "+options.accel)
    
    return options,args

#===================================================================================================
# main()-function
#===================================================================================================
def main(accel):
    '''
    :Parameters:
        'accel': string
            Indicates the used accelerator
        'names_dict': dict: string --> float
            <description>
    :Return: int
        0 if execution was successful otherwise !=0
    '''
    
    ## Never use exec() or execFile()
    
    ## If you have to use try/catch then catch specific exceptions but never all
    try:
        twiss_file = metaclass.twiss("I_do_not_exist.dat")
    except IOError:
        traceback.print_exc()
        return 1
      
    print twiss_file
    print accel
    
    return 0


#===================================================================================================
# helper-functions
#===================================================================================================
def do_x(values_list, names_dict,):
    '''
    Does x...
    
    :Parameters:
        'values_list': list of float
            <description>
        'names_dict': dict: string --> float
            <description>
    :Return: float
        <description>     
    '''
    pass

def do_y():
    '''
    Used by do_x() to do y...
    
    :Return: float
        <description>     
    '''
    pass

#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":
    (options, args) = parse_args()
    
    return_value = main(accel = options.accel)
    
    sys.exit(return_value)
    
