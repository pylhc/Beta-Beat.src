import __init__
import numpy as np
import matplotlib.pyplot as plt
import metaclass
import os.path
from pylab import *
import matplotlib
import os



def plot_fitting(fitx_L,fitx_R,fity_L,fity_R,left_data,right_data,path):
    karray = np.linspace(-0.02,0.02,100)
    lx = fitx_L[0]*karray+fitx_L[1]
    ly = fity_L[0]*karray+fity_L[1]
    rx = fitx_R[0]*karray+fitx_R[1]
    ry = fity_R[0]*karray+fity_R[1]

    fig1 = plt.figure(figsize=(11,8))
    fig1.patch.set_facecolor('white')
    t1 = fig1.add_subplot(111)
    t1.set_ylim([min(left_data.TUNEX)-0.001,max(left_data.TUNEX)+0.001])
    t1.set_xlim([min(left_data.K)*1e3-0.001,max(left_data.K)*1e3+0.001])
    t1.set_xlabel(r'K [$10^{-3}$]')
    t1.set_ylabel(r'Qx [tune units]')
    t1.scatter(left_data.K*1e3, left_data.TUNEX)
    t1.plot(karray*1e3, lx,'r', label='Qx fit: Left Data')
    fig1.tight_layout()
    savefig(os.path.join(path, 'left.x.pdf'))

    fig2 = plt.figure(figsize=(11,8))
    fig2.patch.set_facecolor('white')
    t2 = fig2.add_subplot(111)
    t2.set_ylim([min(left_data.TUNEY)-0.001,max(left_data.TUNEY)+0.001])
    t2.set_xlim([min(left_data.K)*1e3-0.001,max(left_data.K*1e3)+0.001])
    t2.set_xlabel(r'K [$10^{-3}$]')
    t2.set_ylabel(r'Qy [tune units]')
    t2.scatter(left_data.K*1e3, left_data.TUNEY)
    t2.plot(karray*1e3, ly,'r', label='Qy fit: Left Data')
    fig2.tight_layout()
    savefig(os.path.join(path, 'left.y.pdf'))

    fig3 = plt.figure(figsize=(11,8))
    fig3.patch.set_facecolor('white')
    t3 = fig3.add_subplot(111)
    t3.set_ylim([min(right_data.TUNEX)-0.001,max(right_data.TUNEX)+0.001])
    t3.set_xlim([min(right_data.K)*1e3-0.001,max(right_data.K)*1e3+0.001])
    t3.set_xlabel('K [$10^{-3}$]')
    t3.set_ylabel('Qx [tune units]')
    t3.scatter(right_data.K*1e3, right_data.TUNEX)
    t3.plot(karray*1e3, rx,'r', label='Qx fit: Right Data')
    fig3.tight_layout()
    savefig(os.path.join(path, 'right.x.pdf'))

    fig4 = plt.figure(figsize=(11,8))
    fig4.patch.set_facecolor('white')
    t4 = fig4.add_subplot(111)
    t4.set_ylim([min(right_data.TUNEY)-0.001,max(right_data.TUNEY)+0.001])
    t4.set_xlim([min(right_data.K)*1e3-0.001,max(right_data.K)*1e3+0.001])
    t4.set_xlabel(r'K [$10^{-3}$]')
    t4.set_ylabel(r'Qy [tune units]')
    t4.scatter(right_data.K*1e3, right_data.TUNEY)
    t4.plot(karray*1e3, ry,'r', label='Qy fit: Right Data')
    fig4.tight_layout()
    savefig(os.path.join(path, 'right.y.pdf'))
    plt.show()



if __name__=='__main__':
	''' insert data '''
	make_fit_plots(fitx_L,fitx_R,fity_L,fity_R,left_data,right_data,path)
