import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import math

fin=open("ALLBPMs", "r")
lines=fin.readlines()
fin.close()

position_x=[]
position_y=[]
for line in lines:
    if line.startswith("#"):
        continue
    elif line.startswith("0") and "bpmw.4r3.b1" in line:
        i=3
        while i<8003:
            position_x.append(line.split(" ")[i])
            i=i+1
    elif line.startswith("1") and "bpmw.4r3.b1" in line:
        i=3
        while i<8003:
            position_y.append(line.split(" ")[i])
            i=i+1
    else:
        continue
 
fout=open("singleBPMsignal.dat", "w")
for i in range(len(position_x)):
    dataout=[position_x[i], position_y[i]]
    csvfout=csv.writer(fout,delimiter='\t')
    csvfout.writerow(dataout)
fout.close()
