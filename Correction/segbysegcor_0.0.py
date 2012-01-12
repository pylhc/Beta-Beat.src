################################################################
#                                                              #
#  @ Glenn Vanbavinckhove  => Date: 14/09/09                   #
#                                                              #
################################################################
#
#  !=> segbysegcor_0.0.py : Construction of the program
#
#
#

###### imports
from optparse import OptionParser
from metaclass import twiss
import os,sys




###### optionparser
parser = OptionParser()
parser.add_option("-a", "--accel",
                help="Which accelerator: LHCB1 LHCB2 SPS RHIC SOLEIL",
                metavar="ACCEL", default="LHCB1",dest="accel")


(options, args) = parser.parse_args()
