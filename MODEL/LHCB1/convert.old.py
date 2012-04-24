#from Numeric import *
from string import split, replace
import sys, os
#f11='/afs/cern.ch/user/r/rcalaga/scratch1/bblr.rdt/'
#f11='/afs/cern.ch/eng/sl/lintrack/ForRam/BetaBeat2006-07-15/multit/'
#f22='37gev000amp1'
f22=sys.argv[1]

filename=f22

output=f22+'.sdds'

os.system('cp '+filename+' '+ output)

