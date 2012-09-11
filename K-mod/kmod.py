#--- read timber data to calc fitted tune values (including)
#--- smoothed data correponding to the Qtrims
#--- R. Calaga, Mar 14, 2011

from Timber import parseout
import sys, pylab,re, time, dateutil
import numpy as np
from string import split, replace
from optparse import OptionParser
from scipy.optimize import leastsq

def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default="None",dest="FILE")
    parser.add_option("-p", "--ip",
                      help="calculate tunes",
                      metavar="IP", default=None,dest="IP")
    parser.add_option("-s", "--smooth",
                      help="smooth data [0 or 1, default=1]",
                      metavar="SMOOTH", default="1",dest="SMOOTH")
    parser.add_option("-t", "--threshold",
                      help="time thres to match Qtrim-BBQ [def=1sec]",
                      metavar="THRES", default="1.0",dest="THRES")
    parser.add_option("-i", "--shift",
                      help="subset between start-end points [def=0]",
                      metavar="SHIFT", default="0",dest="SHIFT")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    return opt, args

def MAVG(y,window=None):
    #--- moving average  but specify window region
    #--- also available 'hanning', 'hamming', 'bartlett', 'blackman' windows
    y=y[st:ed]
    win=len(y);w=np.ones(win,'d')
    s=np.r_[2*y[0]-y[win-1::-1],y,2*y[-1]-y[-1:-win:-1]]
    if window: w=eval('np.'+window+'(win)')
    x=np.convolve(w/w.sum(),s,mode='same')
    return x[win:-win+1]


def ST(y,order,dropVals=False):
    #moving triangle smoothing, variable order."""
    triangle=np.array(range(order)+[order]+range(order)[::-1])+1
    smoothed=[]
    for i in range(order,len(y)-order*2):
        point=y[i:i+len(triangle)]*triangle
        smoothed.append(sum(point)/sum(triangle))
    if dropVals: return smoothed
    smoothed=[smoothed[0]]*(order+order/2)+smoothed
    while len(smoothed)<len(y):smoothed.append(smoothed[-1])
    return smoothed

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

fit=lambda p,x: p[0]+p[1]*x 
err=lambda p,x,y:fit(p, x)-y 

def tt(f1,f2):
    return abs(time.mktime(time.strptime(f1,"%Y-%m-%d %H:%M:%S"))-\
               time.mktime(time.strptime(f2,"%Y-%m-%d %H:%M:%S")))

def tunes(x,st,ed,diff=1.0):
    for ky in x.keys():
        if 'TUNE' in ky:
            pl=replace(ky[-9:],':TUNE_','')
            aa=np.array([float(j[0]) for j in x[ky][1]])
            if opt.SMOOTH=="1":aa=SGfilter(aa,win=51,order=4);
            st1=None; ed1=None
            for dx in range(len(x[ky][0])):
                if tt(st,x[ky][0][dx][:-4])<diff: st1=dx;
                if tt(ed,x[ky][0][dx][:-4])<diff: ed1=dx 
            if not (st1 and ed1): 
                print pl, "--> No Matching time within",diff,"sec"
            else:
                xnew=aa[st1:ed1]; x0=np.arange(len(xnew))+st1
                p2,cov,info,mesg,suc=leastsq(err,[0.0,0.0],\
                                             args=(x0,xnew),\
                                             full_output=True)
                chisq=sum(info["fvec"]**2);dof=len(x0)-2
                print "%3s %12s  %12s %10f +/- %10f" \
                      %(pl,split(x[ky][0][st1])[1],split(x[ky][0][ed1])[1],\
                        p2[0],np.sqrt(cov[0,0])*np.sqrt(chisq/dof))
    
def IMEAS(y,key2,thres=1.0e-3):
    indi=[];st=[];ed=[];fac=1.5
    for key in y.keys():
        if 'RPMBC' in key and key2 in key:
            print '#--',replace(key[-9:],':I_MEAS','')
            cur=[float(j[0]) for j in y[key][1]]
            dt=[dateutil.parser.parse(j) for j in y[key][0]]
            dcur=np.diff(cur);acur=np.average(cur)
            for k in range(1,len(cur)):
                if (cur[k]-acur>fac or cur[k]-acur<-fac) and \
                        (dcur[k-1]<thres and dcur[k-1]>-thres):
                    indi.append(k)
            st.append(indi[0])
            for k in range(1,len(indi)):
                if indi[k]-indi[k-1]>10:
                    st.append(indi[k-1]);st.append(indi[k])
            st.append(indi[-1]);st=np.array(st); 
            st=np.reshape(st,(len(st)/2,2))
            for sta,end in st:
                #-- remove the last three digits
                sta1=dt[sta+sh];end1=dt[end-sh]
                sta=y[key][0][sta+sh][:-4];end=y[key][0][end-sh][:-4]
                tunes(y,sta,end,diff=thr)
                #pylab.plot_date(pylab.date2num([sta1,end1]),[0.0,0.0]);
            #pylab.plot_date(pylab.date2num(dt),cur,fmt='-',marker='');
            #pylab.plot_date(pylab.date2num(dt)[1:],dcur,ls='-',marker='');
            #pylab.show()
            
    
if __name__ == "__main__":
    opt,args=parsse(); data=parseout(opt.FILE);
    thr=float(opt.THRES); sh=int(opt.SHIFT)
    #--- for each IP find Q1 trims + calculate smoothed tunes
    if opt.IP:
        IMEAS(data,"L"+opt.IP);
        IMEAS(data,"R"+opt.IP)
        sys.exit()
    
    #-- beam 1 SG filter
    for key in data.keys():
        if 'B1:TUNE_H' in key:
            B1H=np.array([float(j[0]) for j in data[key][1]])
        if 'B1:TUNE_V' in key:
            B1V=np.array([float(j[0]) for j in data[key][1]])
        if 'B2:TUNE_H' in key:
            B2H=np.array([float(j[0]) for j in data[key][1]])
        if 'B2:TUNE_V' in key:
            B2V=np.array([float(j[0]) for j in data[key][1]])
            
    a=SGfilter(B1H,win=51,order=4)
    b=SGfilter(B1V,win=51,order=4)
    pylab.subplot(2,1,1);pylab.title('SG Filter (4th order)')
    pylab.plot(B1H); pylab.plot(a,linewidth=2,label='Beam 1');
    pylab.plot(B1V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune")

    #-- beam 2 SG filter
    a=SGfilter(B2H,win=51,order=7)
    b=SGfilter(B2V,win=51,order=4)
    pylab.subplot(2,1,2);
    pylab.plot(B2H); pylab.plot(a,linewidth=2,label='Beam 2');
    pylab.plot(B2V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune");

    pylab.show()
    sys.exit()

    #-- beam 1 kalman filter
    a=kalman(B1H,R=0.01,Q=2.0e-4)
    b=kalman(B1V,R=0.01,Q=2.0e-4)
    pylab.subplot(2,2,2);pylab.ylim(0.3,0.35)
    pylab.title('Kalman Filter');
    pylab.plot(B1H); pylab.plot(a,linewidth=2,label='Beam 1');
    pylab.plot(B1V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune")

    #-- beam 2 kalman filter
    a=kalman(B2H,R=0.01,Q=2.0e-4)
    b=kalman(B2V,R=0.01,Q=2.0e-4)
    pylab.subplot(2,2,4);
    pylab.plot(B2H); pylab.plot(a,linewidth=2,label='Beam 2');
    pylab.plot(B2V); pylab.plot(b,linewidth=2);pylab.legend(loc=3)
    pylab.ylim(0.3,0.325);pylab.ylabel("Tune");
    
    pylab.show()
    #-- intervals for 2nd scan
    #sted2=[[50,330],[430,500],[536,580],[634,771],[898,1042],\
    #          [1160,1300],[1432,1562],[1740,1878],[2052,2178],\
    #          [2320,2467],[2563,2716]]
    #sted1=[[132,340],[402,464],[510,548],[595,741],[870,1022],\
    #          [1148,1269],[1396,1523],[1704,1825],[2006,2140],\
    #          [2297,2430],[2534,2680]]

    

