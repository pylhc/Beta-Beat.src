'''
.. module: test.test_template

Created on 27 May 2013

This module serves as an example for good code style.

Every module should have a short/abstract description of what it does, how to use it, where it is
used and any other relevant information.
Comments in this document which start with '##' are explanations

The module name should contain only lower letters and words should be divided by an underscore '_'.
Good examples: 'segment_by_segment.py','io_utils.py', 'get_phase.py'
Bad examples: 'SegmentBySegment.py', 'IoUtils.py', 'getPhase.py'

Since we use Git we need no change history. Every change should be stated in the commit message.

.. moduleauthor:: vimaier
'''

## Divide the imports into three sections
## System imports(Python Standard Library)
import sys
import argparse
import traceback

## Third party libraries
import numpy as np

## Own modules or libraries
import Python_Classes4MAD.metaclass


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
    ''' Parses the arguments and returns args '''
    ## Argparse is first available with 2.7
    ## However, it is possible to install the argparse package for Python >=2.3
    ## All important server have argparse
    ## How-to install argparse: https://pypi.python.org/pypi/argparse
    ## Argparse tutorial: http://docs.python.org/2/howto/argparse.html
    parser = argparse.ArgumentParser(description='A description for the script.')

    parser.add_argument("x", type=int, help="the base")
    parser.add_argument("y", type=int, help="the exponent")
    parser.add_argument("-a", "--accel", choice=["LHCB1", "LHCB2", "SPS", "RHIC"])
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    return parser.parse_args()

#===================================================================================================
# main()-function
#===================================================================================================
def main(accel):
    '''
    :param string accel: Indicates the used accelerator
    :param dict names_dict: string --> float -- <description>
    :returns: int -- 0 if execution was successful otherwise !=0
    '''

    ## Never use exec() or execFile()

    ## If you have to use try/catch then catch specific exceptions but never all(except: ....)
    try:
        twiss_file = Python_Classes4MAD.metaclass.twiss("I_do_not_exist.dat")
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

    :param list values_list: float list -- <description>
    :param dict names_dict: string --> float -- <description>
    :returns: float -- <description>
    '''
    pass

def do_y():
    '''
    Used by do_x() to do y...

    :returns: float -- <description>
    '''
    pass

#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":
    args = parse_args()

    return_value = main(accel = args.accel)

    sys.exit(return_value)

