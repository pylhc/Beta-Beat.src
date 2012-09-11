from Timber import *
import sys, pylab
import numpy as np
from optparse import OptionParser

def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default="None",dest="FILE")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    return opt, args

def SGfilter(y, win, order, deriv=0):
    try:
        window_size = np.abs(np.int(win))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] \
                for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m, y, mode='valid')

def kalman(sz,R=0.01,Q=2.0e-4):
    n_iter=len(sz);xhat=np.zeros_like(sz);P=np.copy(xhat)         
    xhatminus=np.copy(xhat); Pminus=np.copy(xhat)    
    K=np.copy(xhat); xhat[0]=0.0; P[0]=1.0;
    
    for k in range(1,n_iter):
        # time update
        xhatminus[k] = xhat[k-1]
        Pminus[k] = P[k-1]+Q

        # measurement update
        K[k] = Pminus[k]/( Pminus[k]+R )
        xhat[k] = xhatminus[k]+K[k]*(sz[k]-xhatminus[k])
        P[k] = (1-K[k])*Pminus[k]
    return xhat

if __name__ == "__main__":
    opt,args=parsse()
    data=parseout(opt.FILE);dat=[]
    B2H=np.array([float(j[0]) for j in data['LHC.BQBBQ.UA43.FFT1_B2:TUNE_H'][1]])
    B1H=np.array([float(j[0]) for j in data['LHC.BQBBQ.UA47.FFT1_B1:TUNE_H'][1]])
    B2V=np.array([float(j[0]) for j in data['LHC.BQBBQ.UA43.FFT1_B2:TUNE_V'][1]])
    B1V=np.array([float(j[0]) for j in data['LHC.BQBBQ.UA47.FFT1_B1:TUNE_V'][1]])

    #-- beam 1 SG filter
    a=SGfilter(B1H,win=51,order=4)
    b=SGfilter(B1V,win=51,order=4)
    pylab.subplot(2,2,1);pylab.title('SG Filter (4th order)')
    pylab.plot(B1H); pylab.plot(a,linewidth=2,label='Beam 1');
    pylab.plot(B1V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune")

    #-- beam 1 kalman filter
    a=kalman(B1H,R=0.01,Q=2.0e-4)
    b=kalman(B1V,R=0.01,Q=2.0e-4)
    pylab.subplot(2,2,2);pylab.ylim(0.3,0.35)
    pylab.title('Kalman Filter');
    pylab.plot(B1H); pylab.plot(a,linewidth=2,label='Beam 1');
    pylab.plot(B1V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune")

    #-- beam 2 SG filter
    a=SGfilter(B2H,win=51,order=4)
    b=SGfilter(B2V,win=51,order=4)
    pylab.subplot(2,2,3);
    pylab.plot(B2H); pylab.plot(a,linewidth=2,label='Beam 2');
    pylab.plot(B2V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune");

    #-- beam 2 kalman filter
    a=kalman(B2H,R=0.01,Q=2.0e-4)
    b=kalman(B2V,R=0.01,Q=2.0e-4)
    pylab.subplot(2,2,4);
    pylab.plot(B2H); pylab.plot(a,linewidth=2,label='Beam 2');
    pylab.plot(B2V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune");
    
    pylab.show()

    

