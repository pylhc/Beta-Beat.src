import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import math

amplitudes = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6]
ampx=[]
ampy=[]
p2p_x_flattop=[]
p2p_y_flattop=[]
p2p_x_end=[]
p2p_y_end=[]

for amp_x in amplitudes:
    for amp_y in amplitudes:

        fin=open("Amplitude_"+str(amp_x)+"-"+str(amp_y)+"/acdipole.dynappert.1.dat.obs0002.p0001", "r")
        lines=fin.readlines()
        fin.close()
        x_flattop=[]
        y_flattop=[]
        x_end=[]
        y_end=[]
        for line in lines:
            if line.startswith("@") or line.startswith("*") or line.startswith("$"):
                continue
            elif float((line.lstrip()).split()[1])>2000 and float((line.lstrip()).split()[1])<=2500:
                x_flattop.append(float((line.lstrip()).split()[2]))
                y_flattop.append(float((line.lstrip()).split()[4]))
            elif float((line.lstrip()).split()[1])>10000:
                x_end.append(float((line.lstrip()).split()[2]))
                y_end.append(float((line.lstrip()).split()[4]))
            else:
                continue
        
        if len(x_flattop)==0:
            p2p_x_flattop.append(0.0)
        else:
            p2p_x_flattop.append(max(x_flattop)-min(x_flattop))


        if len(y_flattop)==0:
            p2p_y_flattop.append(0.0)
        else:
            p2p_y_flattop.append(max(y_flattop)-min(y_flattop))

        if len(x_end)==0:
            p2p_x_end.append(0.0)
        else:
            p2p_x_end.append(max(x_end)-min(x_end))

        if len(y_end)==0:
            p2p_y_end.append(0.0)
        else:
            p2p_y_end.append(max(y_end)-min(y_end))

        ampx.append(amp_x)
        ampy.append(amp_y)

fout=open("dynamic_apperture.dat", "w")
fout.write("x[mm]\ty[mm]\tampx_meas[mm]\tampy_meas[mm]\tp2p_x[m]\tp2p_y[m]\n")
for k in range(len(p2p_x_end)):
    if k==0 or ampx[k]==ampx[k-1]:
        dataout = [ampx[k], ampy[k], 0.5*1000*p2p_x_flattop[k], 0.5*1000*p2p_y_flattop[k], p2p_x_end[k], p2p_y_end[k]]
        csvfout=csv.writer(fout,delimiter='\t')
        csvfout.writerow(dataout)
    else:
        fout.write("\n")
        dataout = [ampx[k], ampy[k], 0.5*1000*p2p_x_flattop[k], 0.5*1000*p2p_y_flattop[k], p2p_x_end[k], p2p_y_end[k]]
        csvfout=csv.writer(fout,delimiter='\t')
        csvfout.writerow(dataout)
fout.close()
