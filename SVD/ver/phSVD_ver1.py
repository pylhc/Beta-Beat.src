from numpy import *
import os,pylab,sys,re
from rhicdata25 import rhicdata
from operator import mod
from metaclass25 import twiss

def writexy(filename, aDD):
    f=open(filename+'_svdx','w')    
    f.write('@ Q1 %le '+str(Q1)+'\n'+'@ Q1RMS %le '+str(Q1RMS)+'\n')
    f.write('* NAME  S   BINDEX SLABEL  TUNEX   MUX  AMPX'+\
            '  NOISE   PK2PK  AMP01 PHASE01 CO\n')
    f.write('$ %s  %le   %le    %le  %le  %le  %le  %le  %le'+\
            '  %le  %le  %le  %le\n')
    for j in range(len(aDD.H)):
        print >> f, aDD.H[j].name, aDD.H[j].location, BLABEL,\
              SLABEL, TUNEX[j], PHX[j], AMPX[j], NOISEX[j], \
              PK2PKX[j], AMP01, PHASE01, COX[j]
    f.close()

    f=open(filename+'_svdy','w') 
    f.write('@ Q2 %le '+str(Q2)+'\n'+'@ Q2RMS %le '+str(Q2RMS)+'\n')
    f.write('* NAME  S   BINDEX SLABEL  TUNEY   MUY  AMPY'+\
            '  NOISE   PK2PK  AMP10 PHASE10 CO\n')
    f.write('$ %s  %le   %le    %le  %le  %le  %le  %le  %le'+\
            '  %le  %le  %le  %le\n')
    for j in range(len(aDD.V)):
        print >> f, aDD.V[j].name, aDD.V[j].location, BLABEL,\
              SLABEL, TUNEY[j], PHY[j], AMPY[j], NOISEY[j], \
              PK2PKY[j], AMP01, PHASE01, COY[j]
    f.close()

def nHist(BBX,BBY,filename,nbinx=70,nbiny=70):
    fx=open(filename+'_nhistx','w')
    fx.write("*"+'%5s %5s' % ("RMSNX", "FREQ")+"\n")
    fx.write("$"+'%5s  %5s' % ("%le","%le")+"\n")

    fy=open(filename+'_nhisty','w')
    fy.write("*"+'%5s %5s' % ("RMSNY", "FREQ")+"\n")
    fy.write("$"+'%5s  %5s' % ("%le","%le")+"\n")
    ax=histogram(BBX,bins=nbinx)
    ay=histogram(BBY,bins=nbiny)
    for j in range(nbinx): print >> fx, ax[1][j], ax[0][j]
    for j in range(nbinx): print >> fy, ay[1][j], ay[0][j]
    fx.close(); fy.close()
     
def bMat(aa):
    xD=[];yD=[]
    for j in aa.H: xD.append(j.data)
    for j in aa.V: yD.append(j.data)
    return transpose(array(xD)),transpose(array(yD))

def phBeta(S1,S2,V1,V2): 
    L1=S1*V1;L2=S2*V2
    phase=arctan2(L1,L2)/2./pi
    amp=(L1**2+L2**2)
    return array(phase), array(amp)

def phAdv(phase,bm=1):
    adv=[]
    for j in range(1,len(phase)):
        if bm==2: phadv=phase[j-1]-phase[j]
        else: phadv=phase[j]-phase[j-1]
        while phadv<0: phadv=phadv+1.0
        adv.append(phadv)
    return adv
    
def ptp(BB):
    co=[]; ptop=[]
    for j in range(shape(BB)[1]):
        co.append(average(BB[:,j]))
        ptop.append(max(BB[:,j])-min(BB[:,j]))
    return array(co), array(ptop)

def noise(BB,level=5):
    U,S,V = linalg.svd(BB);S[:level]=0.0
    Sigma=zeros_like(BB);n=min(BB.shape);Sigma[:n,:n]=diag(S)
    BN=dot(U,dot(Sigma,V))
    nrms=sqrt(sum(BN**2,0)/len(BN))
    return nrms

def fTune(BB):
    txx=[]; tyy=[]
    tx=transpose(array([abs(fft.rfft(j.data-average(j.data)))\
                        for j in BB.H]))
    ty=transpose(array([abs(fft.rfft(j.data-average(j.data)))\
                        for j in BB.V]))
    q1=int((qx0-win)*len(tx)*2);q2=int((qx0+win)*len(tx)*2)
    for j in range(shape(tx)[1]):
        aa=argmax(tx[q1:q2,j])+q1; aa1=argmax(tx[:,j])
        if aa!=aa1:
            print "warning, Qx outside window:",aD.H[j].name
        txx.append(aa/2.0/len(tx))
    q1=int((qy0-win)*len(ty)*2);q2=int((qy0+win)*len(ty)*2)
    for j in range(shape(ty)[1]):
        aa=argmax(ty[q1:q2,j])+q1; aa1=argmax(ty[:,j])
        if aa!=aa1:
            print "warning, Qy outside window:", aD.V[j].name
        tyy.append(aa/2.0/len(ty))
    return txx, tyy

def rinput(drr='./'):
    global qx0, qy0, win
    try:
        ff=open(drr+'DrivingTerms').readline()
        FLE=ff.split()[0];NTR=float(ff.split()[2])
        for j in open(drr+'Drive.inp'):
            if 'ISTUN' in j: win=float(j.split('=')[1])
            if 'TUNE' and 'X' in j: qx0=float(j.split('=')[1])
            if 'TUNE' and 'Y' in j: qy0=float(j.split('=')[1])
        print "Input {Qx, Qy}:", qx0, qy0
    except IOError:
        print "Drive.inp -or- DrivingTerms not in",dir
        print 'exiting...'
        sys.exit()        
    return FLE, NTR

def fmod(BB,md1=1,md2=2, tol=10.0):
    BB=(BB-average(BB,0))/sqrt(len(BB))    
    U,S,V = linalg.svd(BB,full_matrices=0)
    if S[md1]/S[md2]>tol:
        print "warning: singular vals are very different"
    return S[md1],S[md2],V[md1,:],V[md2,:]

def pmod(BB,ss,st=0):
    BB=(BB-average(BB,0))/sqrt(len(BB))    
    U,S,V=linalg.svd(BB,full_matrices=0); print S[:5]
    F=transpose(array([abs(fft.rfft(U[:,j])) for j in range(len(S))]))
    #pylab.scatter(arange(len(S)),S);pylab.show()
    F0=arange(len(F))/2.0/len(F)
    k=0; jj=arange(12)+1
    for j in range(st,st+4):        
        pylab.subplot(4,3,jj[k]);pylab.plot(U[:,j]);k+=1
        pylab.subplot(4,3,jj[k]);pylab.plot(F0,F[:,j]);k+=1
        pylab.subplot(4,3,jj[k]);pylab.plot(V[j,:]);k+=1    
    pylab.show()
   
if __name__ == "__main__":
    #-- specify directory of drive.inp, default ./
    file, NTURNS=rinput() 
    BLABEL=0;SLABEL=0;AMP01=0;PHASE01=0
    
    #--- read tbt input
    aD=rhicdata(file);bx,by=bMat(aD)
    sx=array([j.location for j in aD.H])/1.0e3
    sy=array([j.location for j in aD.V])/1.0e3
    
    #-- Tunes, CO, PTOP & Noise
    TUNEX, TUNEY=fTune(aD)
    Q1=average(TUNEX); Q2=average(TUNEY)
    print "Average {Qx, Qy}:", round(Q1,4), round(Q2,4)
    Q1RMS=std(TUNEX); Q2RMS=std(TUNEY)
    COX, PK2PKX=ptp(bx); COY, PK2PKY=ptp(by)
    NOISEX=noise(bx); NOISEY=noise(by)

    #--- find modes & compute phadv(x,y)
    s1,s2,v1,v2=fmod(bx,md1=0,md2=1)
    print "Sing vals, X:", round(s1,3),round(s2,3)
    PHX,AMPX=phBeta(s1,s2,v1,v2)
    s1,s2,v1,v2=fmod(by,md1=0,md2=1)
    print "Sing vals, Y:", round(s1,3),round(s2,3)
    PHY,AMPY=phBeta(s1,s2,v1,v2)
    
    #-- plot stuff
    #pmod(by,sy,st=0);sys.exit()
    #adx=phAdv(PHX,bm=2)
    #pylab.scatter(sx[1:],adx);pylab.show()
    #pylab.plot(sx,AMPX);pylab.show()
    #sys.exit()
    
    #-- write svdx & svdy
    writexy(file, aD)
    nHist(NOISEX, NOISEY,file,nbinx=45,nbiny=140)
    sys.exit()

