from numpy import *
import os,pylab,sys,re
from rhicdata25 import rhicdata
from operator import mod
from metaclass25 import twiss
from optparse import OptionParser
from SussixBS import *
tp=transpose

def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--accel", 
                      help="Accelerator: LHCB1 LHCB2 SPS RHIC",
                      metavar="ACCEL", default="LHCB1",dest="ACCEL")
    parser.add_option("-p", "--path",
                      help="Path to bpm file",
                      metavar="PATH", default="./",dest="PATH")
    parser.add_option("-m", "--mode",
                      help="Uncoupled/Coupled 0/1",
                      metavar="MODE", default="1",dest="MODE")
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default="None",dest="FILE")
    (opt, args) = parser.parse_args()
    print "Selected accelerator:", opt.ACCEL
    print "Selected path:", opt.PATH
    print "Selected mode:", opt.MODE
    return opt, args

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
    return tp(array(xD)),tp(array(yD))
    
def phBeta(L1,L2): 
    #L1=S1*V1;L2=S2*V2
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

def fTune(BX,BY):
    txx=[]; tyy=[]
    tx=abs(fft.rfft(BX-average(BX,0), axis=0))
    ty=abs(fft.rfft(BY-average(BY,0), axis=0))
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

def fmod(BB,m1=0,m2=1,m3=2,m4=3,tol=10.0):
    BB=(BB-average(BB,0))/sqrt(len(BB))    
    U,S,V = linalg.svd(BB,full_matrices=0)
    if S[m1]/S[m2]>tol:
        print "warning: singular vals are very different"
    S=take(S,[m1,m2,m3,m4]);U=take(U,[m1,m2,m3,m4],1)
    V=tp(take(V,[m1,m2,m3,m4],0))
    return U,S,V

def pUSV(U,V,ss,st=0,nm=4):
    F=abs(fft.rfft(U,axis=0));F0=arange(len(F))/2.0/len(F)
    k=0;jj=arange(nm*3)+1
    for j in range(st,st+nm):        
        pylab.subplot(nm,3,jj[k]);pylab.plot(U[:,j]);k+=1
        pylab.subplot(nm,3,jj[k]);pylab.plot(F0,F[:,j]);k+=1
        pylab.subplot(nm,3,jj[k]);pylab.plot(V[j,:]);k+=1    
    pylab.show()

def bMatXY(aa):
    xD={};yD={};cKeys=[];BBX=[];BBY=[]
    for j in aa.H: xD[j.name]=j.data
    for j in aa.V: yD[j.name]=j.data
    for j in xD.keys():
        if j in yD.keys():cKeys.append(j)
    for j in cKeys: BBX.append(xD[j])
    for j in cKeys: BBY.append(yD[j])
    BBX=tp(array(BBX));BBY=tp(array(BBY))
    return concatenate((BBX,BBY),1)

def fmodXY(BB,m1=0,m2=1,m3=2,m4=4,tol=10.0):
    BB=(BB-average(BB,0))/sqrt(len(BB))    
    U,S,V = linalg.svd(BB,full_matrices=0)
    if S[m1]/S[m2]>tol or S[m3]/S[m4]>tol:
        print "warning: singular vals are very different"
    S=take(S,[m1,m2,m3,m4]);U=take(U,[m1,m2,m3,m4],1)
    V=tp(take(V,[m1,m2,m3,m4],0))
    return U,S,V

def fPeak(UU, win=0.01):
    sussix_inp(ir=1,turns=len(UU),tunex=qx0,tuney=qy0,\
               istun=0.02,idam=2,narm=100)
    Sus=sussix(UU,UU,zeros(len(UU)))
    n1=Sus.tunexy[0];n2=Sus.tunexy[1]
    print "Mode Tunes:", round(n1,4), round(n2,4)
    return n1,n2

def projA(BB,U,S,V):
    n1,n2=fPeak(U[:,0])
    NTR=len(U); NBP=len(V);CS=zeros_like(U);BBR=[]
    for k in range(NTR):
        CS[k,0]=cos(2*pi*n1*k);CS[k,1]=sin(2*pi*n1*k)
        CS[k,2]=cos(2*pi*n2*k);CS[k,3]=-sin(2*pi*n2*k)
    for k in range(4):
        BBR.append([dot(BB[:,j],CS[:,k]) for j in range(NBP)])
    BBR=array(BBR)
    
    try: BBR=dot(linalg.inv(dot(BBR,tp(BBR))), BBR)
    except: raise ValueError, "Singular Matrix"
    SVT=dot(S*identity(4), tp(V))
    O=dot(SVT, tp(BBR))
    DET = abs(linalg.det(O))
    if DET != 0.0: O = O/math.pow(DET,.25)
    print 'det rotation |O|:', linalg.det(O)
    U1=dot(U, O);V1=dot(transpose(O),SVT)
    return U1,V1

if __name__ == "__main__":
    #-- Input options
    opt,args=parsse()
    if opt.FILE=="None": file, NTURNS=rinput()
    else:file=opt.PATH+opt.FILE
    print "Input File:", file    
    BLABEL=0;SLABEL=0;AMP01=0;PHASE01=0
    
    #--- read tbt input
    aD=rhicdata(file);bx,by=bMat(aD)
    sx=array([j.location for j in aD.H])/1.0e3
    sy=array([j.location for j in aD.V])/1.0e3
    
    #-- Tunes, CO, PTOP & Noise
    TUNEX, TUNEY=fTune(bx,by)
    Q1=average(TUNEX); Q2=average(TUNEY)
    print "Average {Qx, Qy}:", round(Q1,4), round(Q2,4)
    Q1RMS=std(TUNEX); Q2RMS=std(TUNEY)
    COX, PK2PKX=ptp(bx); COY, PK2PKY=ptp(by)
    NOISEX=noise(bx); NOISEY=noise(by)
    
    #--- find modes & compute phadv(x,y)
    if opt.MODE=="0": #--- ph/amp from uncoupled case
        Ux,Sx,Vx=fmod(bx,m1=0,m2=1,m3=2,m4=3)
        print "Sing vals, X:", round(Sx[0],3),round(Sx[1],3)
        PHX,AMPX=phBeta(Sx[0]*Vx[:,0],Sx[1]*Vx[:,1])
        Uy,Sy,Vy=fmod(by,m1=0,m2=1,m3=2,m4=3)
        print "Sing vals, Y:", round(Sy[0],3),round(Sy[1],3)
        PHY,AMPY=phBeta(Sy[0]*Vy[:,0],Sy[1]*Vy[:,1])
        #pUSV(u1,v1,sx,nm=4)
    elif opt.MODE=="1": #-- ph/amp from rotated vectors
        Ux,Sx,Vx=fmod(bx,m1=0,m2=1,m3=2,m4=3)
        print "Sing vals, X:", round(Sx[0],3),round(Sx[1],3)
        Uy,Sy,Vy=fmod(by,m1=0,m2=1,m3=2,m4=3)
        print "Sing vals, Y:", round(Sy[0],3),round(Sy[1],3)
        u1,v1=projA(bx,Ux,Sx,Vx);#pUSV(u1,v1,sx,nm=4)
        PHX,AMPX=phBeta(v1[0,:],v1[1,:])
        u2,v2=projA(by,Uy,Sy,Vy);#pUSV(u1,v1,sy,nm=4)
        PHY,AMPY=phBeta(v2[0,:],v2[1,:])
    elif opt.MODE=="2": 
        bxy=bMatXY(aD);sxy=concatenate((sx,sy))
        u,s,v=fmodXY(bxy,m1=0,m2=1,m3=2,m4=3)
        print "Sing vals, X:", round(s[0],3),round(s[1],3)
        print "Sing vals, Y:", round(s[2],3),round(s[3],3)
        u1,v1=projA(bxy,u,s,v); xlen=shape(v1)[1]/2
        PHX,AMPX=phBeta(v1[0,:xlen],v1[1,:xlen])
        PHY,AMPY=phBeta(v1[2,xlen:],v1[3,xlen:])
        #pUSV(u1,v1,sxy,nm=4)
    else: print "no mode selected, exiting..."; sys.exit()

    #-- plot stuff
    #adx=phAdv(PHX,bm=2)
    #pylab.scatter(sx[1:],adx);pylab.show()
    #pylab.plot(sx,AMPX);pylab.show()
    #sys.exit()
    
    #-- write svdx & svdy
    writexy(file, aD)
    nHist(NOISEX, NOISEY,file,nbinx=45,nbiny=140)
    sys.exit()

