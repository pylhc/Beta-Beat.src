
from metaclass import *


model=twiss('twiss.dat')

names=model.NAME


filee=open('./observe_bpms_new','w')

for i in range(len(names)):

    filee.write('ptc_observe, place= '+names[i]+' ;\n')


filee.close
