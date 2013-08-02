import sys
import os

if "/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/" not in sys.path: # add internal path for python scripts to current environment (tbach)
    sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')

if sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..")) not in sys.path: # Root of Beta-Beat.src
	sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..",".."))
