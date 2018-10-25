import os.path
from pylab import *
import os


def plot_fitting(fitx_L,fitx_R,fity_L,fity_R,left_data,right_data,path):
    karray = np.linspace(-0.02,0.02,100)
    lx = fitx_L[0]*karray+fitx_L[1]
    ly = fity_L[0]*karray+fity_L[1]
    rx = fitx_R[0]*karray+fitx_R[1]
    ry = fity_R[0]*karray+fity_R[1]

    fig1 = plt.figure(figsize=(11,8))
    fig1.patch.set_facecolor('white')
    t1 = fig1.add_subplot(221)
    t2 = fig1.add_subplot(222)
    t3 = fig1.add_subplot(223)
    t4 = fig1.add_subplot(224)

    t1.set_ylim([min(left_data.TUNEX)-0.001,max(left_data.TUNEX)+0.001])
    t1.set_xlim([min(left_data.K)*1e3-0.001,max(left_data.K)*1e3+0.001])
    t1.set_title('Left x fit')
    t1.set_xlabel(r'K [$10^{-3}$]')
    t1.set_ylabel(r'Qx [tune units]')
    t1.errorbar(left_data.K*1e3, left_data.TUNEX, left_data.TUNEX_ERR, fmt='o')
    t1.plot(karray*1e3, lx,'r', label='Qx fit: Left Data')

    t2.set_ylim([min(right_data.TUNEX)-0.001,max(right_data.TUNEX)+0.001])
    t2.set_xlim([min(right_data.K)*1e3-0.001,max(right_data.K)*1e3+0.001])
    t2.set_title('Right x fit')
    t2.set_xlabel('K [$10^{-3}$]')
    t2.set_ylabel('Qx [tune units]')
    t2.errorbar(right_data.K*1e3, right_data.TUNEX, right_data.TUNEX_ERR, fmt='o')
    t2.plot(karray*1e3, rx,'r', label='Qx fit: Right Data')

    t3.set_ylim([min(left_data.TUNEY)-0.001,max(left_data.TUNEY)+0.001])
    t3.set_xlim([min(left_data.K)*1e3-0.001,max(left_data.K*1e3)+0.001])
    t3.set_title('Left y fit')
    t3.set_xlabel(r'K [$10^{-3}$]')
    t3.set_ylabel(r'Qy [tune units]')
    t3.errorbar(left_data.K*1e3, left_data.TUNEY, left_data.TUNEY_ERR, fmt='o')
    t3.plot(karray*1e3, ly,'r', label='Qy fit: Left Data')

    t4.set_ylim([min(right_data.TUNEY)-0.001,max(right_data.TUNEY)+0.001])
    t4.set_xlim([min(right_data.K)*1e3-0.001,max(right_data.K)*1e3+0.001])
    t4.set_title('Right y fit')
    t4.set_xlabel(r'K [$10^{-3}$]')
    t4.set_ylabel(r'Qy [tune units]')
    t4.errorbar(right_data.K*1e3, right_data.TUNEY, right_data.TUNEY_ERR, fmt='o')
    t4.plot(karray*1e3, ry,'r', label='Qy fit: Right Data')

    fig1.tight_layout()
    savefig(os.path.join(path, 'fit_plots.pdf'))
    savefig(os.path.join(path, 'fit_plots.png'))
    plt.show()


if __name__=='__main__':
    # insert data
    make_fit_plots(fitx_L,fitx_R,fity_L,fity_R,left_data,right_data,path)
