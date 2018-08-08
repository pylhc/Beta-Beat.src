import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("/afs/cern.ch/work/f/fcarlier/public/for/Jaime/sps/spstools/")
from MULTIT import sdds_to_dict 

def read_file(path):
  data_dict = sdds_to_dict(path) 
  namesex, namesey = [], []
  datax, datay = [], []
  
  for key, data in data_dict.iteritems():
    if key.startswith('SPS.'):
      if key.endswith('H-SUM'):
        new_key = key.replace('SPS.', '').replace('.H/H-SUM','')
        namesex.append(new_key)
        datax.append(data)
      elif key.endswith('V-SUM'):
        new_key = key.replace('SPS.', '').replace('.V/V-SUM','')
        namesey.append(new_key)
        datay.append(data)

  return namesex, np.array(datax), namesey, np.array(datay), None 
