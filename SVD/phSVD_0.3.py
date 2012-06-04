#--- new ver of phSVD.py with multi-threading for SVD
#--- Slightly better peak detection included
#--- R. Calaga, Feb 24, 2011
import sys,os,time

if "/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/" not in sys.path: # add internal path for python scripts to current environment (tbach, 2012/05)
  sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')

from numpy import *; svd=linalg.svd; fftt=fft.fft
from rhicdata25 import rhicdata;
from metaclass25 import twiss
from string import split, replace
from operator import mod
from optparse import OptionParser
from handythread import foreach

tp=transpose;ctn=concatenate


def parsse():
    mode={"0":"uncoupled", "1":"coupled", "2":"Concatenated"} 
    parser = OptionParser()
    parser.add_option("-a", "--accel", 
                      help="Accelerator: LHCB1 LHCB2 SPS RHIC",
                      metavar="ACCEL", default="LHCB1",dest="ACCEL")
    parser.add_option("-p", "--path",
                      help="Path to bpm file",
                      metavar="PATH", default="./",dest="PATH")
    parser.add_option("-m", "--mode",
                      help="Uncoupled/Coupled/Concatenated 0/1/2",
                      metavar="MODE", default="0",dest="MODE")
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default="None",dest="FILE")
    parser.add_option("-c", "--calculate",
                      help="optics/modes-only 0/1",
                      metavar="CALC", default="0",dest="CALC")
    parser.add_option("-u", "--usemode",
                      help="specify 4 modes to use, or auto",
                      metavar="USE", default="0,1,0,1",dest="USE")
    parser.add_option("-n", "--noise",
                      help="flag to calculate noise, 0/1",
                      metavar="FLAG", default="1",dest="FLAG")
    (opt, args) = parser.parse_args()
    print "Selected accelerator:", opt.ACCEL
    print "Selected path:", opt.PATH
    print "Selected mode:", mode[opt.MODE]
    return opt, args

def checkTune(q1,q2):
    if q1>0.5: q1=1-q1
    if q2>0.5: q2=1-q2
    if q1<0.0: q1=q1%1.0
    if q2<0.0: q2=q2%1.0 
    return q1, q2

def rinput(drr='./'):
    global qx0, qy0, win
    if os.path.isfile(drr+'/DrivingTerms') and \
           os.path.isfile(drr+'/Drive.inp'):
        ff=open(drr+'/DrivingTerms').readline()
        FLE=ff.split()[0];NTR=float(ff.split()[2])
        for j in open(drr+'/Drive.inp'):
            if 'ISTUN' in j: win=float(j.split('=')[1])
            if 'TUNE' and 'X' in j: qx0=float(j.split('=')[1])
            if 'TUNE' and 'Y' in j: qy0=float(j.split('=')[1])
        print "Input {Qx, Qy}:", qx0, qy0
        qx0,qy0=checkTune(qx0,qy0)
    else:
        print "DrivingTerms/Drive.inp not in",drr, '\n exiting...'
        sys.exit()
    return FLE, NTR

def bMat(aa):
    xD=tp(array([j.data for j in aa.H]))
    yD=tp(array([j.data for j in aa.V]))
    sx=array([j.location for j in aD.H])
    sy=array([j.location for j in aD.V])
    return xD,yD,sx,sy

def findpk(v,delta,x=None):
    #-- min/max peak finder, but only real part
    maxtab=[]; mintab=[]; v=asarray(v)
    if x is None: x=arange(len(v))
    if len(v) != len(x):sys.exit('Input vx!= x')
    if not isscalar(delta):sys.exit('Input delta not scalar')
    if delta <= 0:sys.exit('Input delta not positive')
    mn,mx=Inf,-Inf;mnpos,mxpos=NaN,NaN; lookformax=True
    for i in arange(len(v)):
        this = v[i]
        if this>mx: mx=this; mxpos=x[i]
        if this<mn: mn=this; mnpos=x[i]
        if lookformax: #-- only real part
            if this<mx-delta and mxpos<500:
                maxtab.append((mxpos, mx))
                mn=this; mnpos=x[i]; lookformax=False
        else: #--only real part
            if this > mn+delta and mxpos<500:
                mintab.append((mnpos,mn))
                mx=this; mxpos=x[i]; lookformax=True
    return array(maxtab)

def peakInterpolate(dt,indx):
    lg=len(dt);freq=fft.fftfreq(lg);
    #--- peak and 2 points (before & after)
    x0=freq[indx]; xm1=freq[indx-1]; x1=freq[indx+1]
    ym1=dt[indx-1];y0=dt[indx];y1=dt[indx+1]
    if ym1!=0.0 and y0!=0.0 and y1!=0.0:
        #--- quad interpolation
        p = x0+(y1-ym1)/(2*(2*y0-y1-ym1))/lg
        y = y0-0.25*(ym1-y1)*p
        a = 0.5*(ym1-2*y0+y1)
    else: p = 0.0; y=0.0; a=0.0
    return p,y

def findTunes(data,q0,tol=0.1):#-- pk finder + quad interolate
    row,col=shape(data); txx=[];l=range(col)
    amp=abs(fftt(data-average(data,0),axis=0))
    #fr=argsort(amp,axis=0);txx=fr[-1,:]
    for j in range(col):
        fr=array(findpk(amp[:,j],tol))
        indx=argsort(fr[:,1]);fr=take(fr,indx,0)
        pk,am=peakInterpolate(amp[:,j],fr[-1,0])
        txx.append(pk)
        if q0-win>pk>q0+win:print "Tunes outside window", pk
    return txx

def findmodes(FREQ,q0,tol=0.1):
    mm=[];row,col=shape(FREQ); nm=10
    for j in range(nm):
        fr=array(findpk(FREQ[:,j],tol))
        indx=argsort(fr[:,1]);fr=take(fr,indx,0)[-1,0]
        if q0-win<fr/row<q0+win:
            mm.append(j)
        else:  print "Tunes outside window", fr/row, "mode",j
        if len(mm)==2: break
    
    return mm
    
def ptp(data):
    co=[]; ptop=[]; row,col=shape(data)
    for j in range(col):
        co.append(average(data[:,j]))
        ptop.append(max(data[:,j])-min(data[:,j]))
    return array(co), array(ptop)

def computeSVD(data,level=20):
    dat=(data-average(data,0))/sqrt(len(data))
    U,S,V=svd(dat,full_matrices=0)
    if opt.FLAG!="0": NOISE=noise(U,S,V)
    else: NOISE=data[0,:]*0.0
    return U,S,tp(V),NOISE

def phBeta(L1,L2): 
    #L1=S1*V1;L2=S2*V2
    phase=arctan2(L1,L2)/2./pi
    amp=(L1**2+L2**2)
    return array(phase), array(amp)

def checkphase(phase,location,OA):
    indx=argsort(location);lg=len(phase);p=0;m=0
    ph=take(phase,indx);ss=take(location,indx)
    for i in range(lg-1):
        if ph[i+1]>ph[i]: p+=1
        else: m+=1
    if   OA=="LHCB1" and p<m: phase=-phase
    elif OA=="LHCB2" and p>m: phase=-phase
    return phase

def noise(U,S,V,level=20):
    Sn=copy(S); Sn[:level]=0.0
    BN=dot(U,dot(diag(Sn),V))
    nrms=sqrt(sum(BN**2,0))
    return nrms

def nHist(BBX,BBY,filename,nbinx=70,nbiny=70):
    fx=open(filename+'_nhistx','w')
    fx.write("*"+' %5s %5s' % (" RMSNX", " FREQ")+"\n")
    fx.write("$"+' %5s  %5s' % (" %le"," %le")+"\n")

    fy=open(filename+'_nhisty','w')
    fy.write("*"+' %5s %5s' % (" RMSNY", " FREQ")+"\n")
    fy.write("$"+' %5s  %5s' % (" %le"," %le")+"\n")
    ax=histogram(BBX,bins=nbinx)
    ay=histogram(BBY,bins=nbiny)
    for j in range(nbinx): print >> fx, ax[1][j], ax[0][j]
    for j in range(nbinx): print >> fy, ay[1][j], ay[0][j]
    fx.close(); fy.close()

def pUSV(U,V,ss,st=0,nm=2):
    import pylab
    F=abs(fftt(U,axis=0));F0=fft.fftfreq(len(F))
    k=0;jj=arange(nm*3)+1;
    for j in range(st,st+nm):        
        pylab.subplot(nm,3,jj[k]);pylab.plot(U[:,j]);k+=1
        pylab.subplot(nm,3,jj[k]);pylab.plot(F0,F[:,j]);k+=1
        pylab.subplot(nm,3,jj[k]);pylab.plot(V[j,:]);k+=1    
    pylab.show()

def wrtUSV(filename,U,S,V,ss,pl=None,nm=20):
    fV=open(filename+'_V'+pl,'w');fF=open(filename+'_F'+pl,'w')
    fS=open(filename+'_S'+pl,'w');fU=open(filename+'_U'+pl,'w')
    fS.write('* INDEX   S\n$ %le %le\n')    
    if nm>len(V): nm=len(V); print '# of modes =',nm
    F=abs(fftt(U,axis=0));F0=fft.fftfreq(len(F))
    fV.write('*  S  ');fF.write('*  INDEX  ');
    for j in range(1,nm+1): fV.write('MOD_'+str(j)+'  ')
    for j in range(1,nm+1): fF.write('MOD_'+str(j)+'  ')
    fV.write('\n$   %le   ');fF.write('\n$   %le   ')
    for j in range(1,nm+1): fV.write('%le   ')
    for j in range(1,nm+1): fF.write('%le   ')
    fV.write('\n'); fF.write('\n')
    for k in range(len(ss)):
        fV.write(str(ss[k])+' ')
        for j in range(nm):fV.write(str(V[j,k])+' ')
        fV.write('\n')
    for k in range(len(F)):
        fF.write(str(F0[k])+' ')
        for j in range(nm):fF.write(str(F[k,j])+' ')
        fF.write('\n')
    for k in range(len(U)):
        fU.write(str(k)+' ')
        for j in range(nm):fU.write(str(U[k,j])+' ')
        fU.write('\n')
    for j in range(len(S)):
        fS.write(str(j)+'   '+str(S[j])+'\n')
    fV.close();fF.close();fS.close();fU.close()
    return F

def writexy(filename, aDD):
    print filename
    f=open(filename+'_svdx','w')    
    f.write('@ Q1 %le '+str(Q1)+'\n'+'@ Q1RMS %le '+str(Q1RMS)+'\n')
    f.write('* NAME  S   BINDEX SLABEL  TUNEX   MUX  AMPX'+\
            '  NOISE   PK2PK  AMP01 PHASE01 CO C11 C12 C21 C22\n')
    f.write('$ %s  %le   %le    %le  %le  %le  %le  %le  %le'+\
            '  %le  %le  %le  %le  %le  %le  %le  \n')
    for j in range(len(aDD.H)):
        print >> f, aDD.H[j].name, aDD.H[j].location, BLABEL,\
              SLABEL,TUNEX[j],PHX[j],sqrt(AMPX[j]),NOISEX[j],PK2PKX[j],\
              AMP01,PHASE01,COX[j],C11[j],C12[j],C21[j],C22[j]
    f.close()

    f=open(filename+'_svdy','w') 
    f.write('@ Q2 %le '+str(Q2)+'\n'+'@ Q2RMS %le '+str(Q2RMS)+'\n')
    f.write('* NAME  S   BINDEX SLABEL  TUNEY   MUY  AMPY'+\
            '  NOISE   PK2PK  AMP10 PHASE10 CO C11 C12 C21 C22\n')
    f.write('$ %s  %le   %le    %le  %le  %le  %le  %le  %le'+\
            '  %le  %le  %le  %le  %le  %le  %le\n')
    for j in range(len(aDD.V)):
        print >> f, aDD.V[j].name, aDD.V[j].location, BLABEL,\
              SLABEL,TUNEY[j],PHY[j],sqrt(AMPY[j]),NOISEY[j],PK2PKY[j],\
              AMP01,PHASE01,COY[j],C11[j],C12[j],C21[j],C22[j]
    f.close()    

def plotOptics(sx,px,sy,py,modelFile):
    import pylab; phm=[]; phx=[]; phy=[]
    a=twiss(modelFile);Q1=a.Q1%1.0; Q2=a.Q2%1.0
    phm.append([a.MUX[0],a.MUY[0]]);
    phx.append((Q1-(px[-1]-px[0]))/2)
    phy.append((-Q2+(py[-1]-py[0]))/2)
    for j in range(1,len(a.MUX)):
        phm.append([a.MUX[j]-a.MUX[j-1], a.MUY[j]-a.MUY[j-1]])
    for j in range(1,len(px)):
        phx.append(mod(px[j]-px[j-1],1))
    for j in range(1,len(py)):
        phy.append(mod(py[j]-py[j-1],1))
    phm=array(phm)
    pylab.subplot(2,1,1);
    pylab.scatter(sx,phx);pylab.plot(a.S,phm[:,0]);
    pylab.subplot(2,1,2);
    pylab.scatter(sy,phy);pylab.plot(a.S,phm[:,1]);
    pylab.show()
    
    
if __name__ == "__main__":
    t0=time.time()
    #-- Input options & read data
    opt,args=parsse()
    file, NTURNS=rinput(opt.PATH)
    print "Input File:", file    
    BLABEL=0;SLABEL=0;AMP01=0;PHASE01=0
    aD=rhicdata(file);bx,by,sx,sy=bMat(aD)
    
    #-- Tunes, CO, PTOP & Noise & write modes
    TUNEX=findTunes(bx,qx0,tol=0.1);
    TUNEY=findTunes(by,qy0,tol=0.1);    
    Q1=average(TUNEX);Q2=average(TUNEY)
    print "Av. Measured {Qx, Qy}:", round(Q1,4), round(Q2,4)
    Q1RMS=sqrt(average(TUNEX)**2); Q2RMS=sqrt(average(TUNEY)**2)
    COX, PK2PKX=ptp(bx); COY, PK2PKY=ptp(by)
    #Ux,Sx,Vx,NOISEX=computeSVD(bx);Uy,Sy,Vy,NOISEY=computeSVD(by);
    [Ux,Sx,Vx,NOISEX],[Uy,Sy,Vy,NOISEY]=foreach(computeSVD,[bx,by],threads=2,return_=True)
    Fx=wrtUSV(file,Ux,Sx,Vx,sx,pl='x');
    Fy=wrtUSV(file,Uy,Sy,Vy,sy,pl='y');

    #--- find modes & compute phadv(x,y)
    if opt.MODE=="0": #--- ph/amp from uncoupled case
        x1,x2=findmodes(Fx,qx0); y1,y2=findmodes(Fy,qy0)
        PHX,AMPX=phBeta(Sx[x1]*Vx[:,x1],Sx[x2]*Vx[:,y2])
        PHY,AMPY=phBeta(Sy[y1]*Vy[:,y1],Sy[y2]*Vy[:,y2])
        C11=zeros(len(COX)+len(COY));C12=C11;C21=C11;C22=C11
    else:
        print "Coupled SVD not available yet"
        print "Only SVD modes written out, exiting.."
        sys.exit()

    #--- check phase for mode inversion
    if opt.ACCEL[:3]=="LHC":#-- ryoichi's phase check
        PHX=checkphase(PHX,sx,opt.ACCEL);
        PHY=checkphase(PHY,sy,opt.ACCEL)

    #-- write svdx & svdy
    writexy(file, aD)
    nHist(NOISEX, NOISEY,file,nbinx=45,nbiny=140)
    
    print "Total time", round(time.time()-t0,2),'s'
    
    sys.exit()
