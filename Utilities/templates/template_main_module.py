'''
Created on 27 May 2013

@author: vimaier

@version: 0.0.1

templates.template_main_module <description>

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import optparse


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
    
    # Check arguments
    if options.accel not in ("LHCB1", "LHCB2", "LHCB4", "SPS", "RHIC", "TEVATRON"):
        raise ValueError("Invalid accelerator: "+options.accel)
    
    return options,args


#===================================================================================================
# main()-function
#===================================================================================================
def main():
    '''
    :Return: int
        0 if execution was successful otherwise !=0
    '''
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
    

#===================================================================================================
# main invocation
#===================================================================================================
if __name__ == "__main__":
    (options, args) = parse_args()
    
    return_value = main(accel = options.accel)
    
    sys.exit(return_value)