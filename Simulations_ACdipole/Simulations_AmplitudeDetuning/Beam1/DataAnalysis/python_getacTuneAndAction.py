import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import math

amplitude = [0.5, 0.75, 1.0, 1.25, 1.5]

action_x=[]
action_y=[]
erraction_x=[]
erraction_y=[]
tune_x=[]
tune_y=[]
errtune_x=[]
errtune_y=[]
for AMP in amplitude:
    fin=open("Amplitude_"+str(AMP)+"/getkickac.out","r")
    lines=fin.readlines()
    fin.close()
    for line in lines:
        if line.startswith("@") or line.startswith("#") or line.startswith("$") or line.startswith("*"):
            continue
        else:
            tune_x.append(float((line.lstrip()).split()[5]))
            errtune_x.append(float((line.lstrip()).split()[6]))
            tune_y.append(float((line.lstrip()).split()[7]))
            errtune_y.append(float((line.lstrip()).split()[8]))
            action_x.append(float((line.lstrip()).split()[13]))
            erraction_x.append(float((line.lstrip()).split()[14]))
            action_y.append(float((line.lstrip()).split()[15]))
            erraction_y.append(float((line.lstrip()).split()[16]))


fout=open("det_with_amp_ac.dat","w")
fout.write("#Qx\t\tQxSTD\t\tQy\t\tQySTD\t\t2Jx\t\t\t2JxSTD\t\t\t2Jy\t\t\t2JySTD\n ")
for i in range(len(amplitude)):
    dataout = [tune_x[i], errtune_x[i], tune_y[i], errtune_y[i], action_x[i], erraction_x[i], action_y[i], erraction_y[i]]
    csvfout=csv.writer(fout,delimiter='\t')
    csvfout.writerow(dataout)
fout.close()
