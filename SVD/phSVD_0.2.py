#!/afs/cern.ch/eng/sl/lintrack/Python-2.5_32bit/Python-2.5_32bit/bin/python

import sys,os,re
try:
    from numpy import *; svd=linalg.svd; fftt=fft.fft
    from rhicdata25 import rhicdata;
    from metaclass import twiss
except:
    from Numeric import *; from FFT import fft as fftt 
    from rhicdata import rhicdata
    from metaclass import twiss
    from LinearAlgebra import singular_value_decomposition as svd
    
    
from string import split, replace
from operator import mod
from optparse import OptionParser
tp=transpose;ctn=concatenate

def parsse():
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
    (opt, args) = parser.parse_args()
    print "Selected accelerator:", opt.ACCEL
    print "Selected path:", opt.PATH
    print "Selected mode:", opt.MODE
    return opt, args

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
     
def bMat(aa):
    xD=[];yD=[]
    for j in aa.H: 
      xD.append(j.data)
    for j in aa.V: 
      yD.append(j.data)
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


def noise(U,S,V,level=20):
    #U,S,V = svd(BB);S[:level]=0.0
    #Sigma=zeros_like(BB);n=min(BB.shape);Sigma[:n,:n]=diag(S)
    Sn=1.0*S; Sn[:level]=0.0
    BN=dot(U,dot(diag(Sn),V))
    #nrms=sqrt(sum(BN**2,0)/len(BN))
    nrms=sqrt(sum(BN**2,0))
    return nrms

def noiseold(BB,level=20): #--older way w/o subtracting C.orbit
    U,S,V = svd(BB);S[:level]=0.0
    Sigma=zeros_like(BB);n=min(BB.shape);Sigma[:n,:n]=diag(S)
    BN=dot(U,dot(Sigma,V))
    nrms=sqrt(sum(BN**2,0)/len(BN))
    return nrms

def window(seq, n=2): #--- not finished
    "Returns a sliding window"
    import itertools
    it=iter(seq)
    result=tuple(itertools.islice(it,n))
    if len(result)==n: yield result    
    for elem in it:
        result=result[1:]+(elem,)
        yield result    

def fTuneOld(BX,BY):#--dumb pk finder in window
    txx=[];tx=abs(fftt(BX-average(BX,0), axis=0))
    tyy=[];ty=abs(fftt(BY-average(BY,0), axis=0))
    q1=int((qx0-win)*len(tx));q2=int((qx0+win)*len(tx))
    for j in range(shape(tx)[1]):
        aa=argmax(tx[q1:q2,j])+q1; aa1=argmax(tx[:,j])
        if aa!=aa1:
            print "warning, Qx outside window:",aD.H[j].name
        txx.append(aa/1.0/len(tx))
    q1=int((qy0-win)*len(ty));q2=int((qy0+win)*len(ty))
    for j in range(shape(ty)[1]):
        aa=argmax(ty[q1:q2,j])+q1; aa1=argmax(ty[:,j])
        if aa!=aa1:
            print "warning, Qy outside window:", aD.V[j].name
        tyy.append(aa/1.0/len(ty))
    return txx, tyy

def fTune(BX,BY):#-- pk finder in window + quad interolate
    txx=[];tx=abs(fftt(BX-average(BX,0), axis=0))
    tyy=[];ty=abs(fftt(BY-average(BY,0), axis=0))
    qs=int((qx0-win)*len(tx));qe=int((qx0+win)*len(tx))
    for j in range(shape(tx)[1]):
        pk,am=peak(tx[:,j],qs,qe)
        if qx0-win>pk>qx0+win:
            print "Warning: Qx outside window:",aD.H[j].name
        txx.append(pk)
    qs=int((qy0-win)*len(ty));qe=int((qy0+win)*len(ty))
    for j in range(shape(ty)[1]):
        pk,am=peak(ty[:,j],qs,qe)
        if qy0-win>pk>qy0+win:
            print "Warning: Qy outside window:",aD.V[j].name
        tyy.append(pk)
    return txx,tyy

def peak(dt,q1,q2):
    indx=argmax(dt[q1:q2])+q1;lg=len(dt)*1.0
    #--- peak and 2 points (before & after)
    x0=indx/lg; xm1=(indx-1)/lg; x1=(indx+1)/lg
    ym1=dt[indx-1];y0=dt[indx];y1=dt[indx+1]
    if ym1!=0.0 and y0!=0.0 and y1!=0.0:
        #--- quad interpolation
        p = x0+(y1-ym1)/(2*(2*y0-y1-ym1))/lg 
        y = y0-0.25*(ym1-y1)*p
        a = 0.5*(ym1-2*y0+y1)
    else: p = 0.0; y=0.0; a=0.0
    return p,y
        
def checkTune(q1,q2):
    if q1> 0.5: q1=1-q1;
    if q2> 0.5: q2=1-q2
    if q1< 0.0: q1=q1 % 1.0;
    if q2< 0.0: q2=q2 % 1.0 
    return q1, q2

def rinput(drr='./'):
    global qx0, qy0, win
    try:
        ff=open(drr+'/DrivingTerms').readline()
        FLE=ff.split()[0];NTR=float(ff.split()[2])
        for j in open(drr+'/Drive.inp'):
            if 'ISTUN' in j: win=float(j.split('=')[1])
            if 'TUNE' and 'X' in j: qx0=float(j.split('=')[1])
            if 'TUNE' and 'Y' in j: qy0=float(j.split('=')[1])
        print "Input {Qx, Qy}:", qx0, qy0
        qx0,qy0=checkTune(qx0,qy0)
    except IOError:
        print "Drive.inp -or- DrivingTerms not in",drr
        print 'exiting...'
        sys.exit()        
    return FLE, NTR

def findModes(FF,qs,qe,nm=20):
    mm=[]
    for j in range(nm):
        pk,ap=peak(FF[:,j],qs,qe)
        pk0=argmax(FF[:,j])*1.0/len(FF)
        print win,pk,pk0,argmax(FF[:,j])*1.0/len(FF),pk
        #if pk0>0.5:pk0=1-pk0
        if pk-win<pk0<pk+win: mm.append(j)
        if len(mm)==2: break
    if len(mm)==0:
        print "Window is out of range, please change the window and check your tunes \n system exit"
        sys.exit()

    # by me (glenn)
    if len(mm)<2:
        mode1=0;mode2=2
    else:
        mode1=mm[0];mode2=mm[1]
    
    #print "Suggested Modes "+pl+":",mm[0],mm[1]
    return mode1,mode2

def fmod(U,S,V,F,md,tol=10.0):
    #BB=(BB-average(BB,0))/sqrt(len(BB))    
    #U,S,V=svd(BB,full_matrices=0);F=abs(fftt(U,axis=0))
    if S[m1]/S[m2]>tol:
        print "warning: singular vals are very different"
    S=take(S,md);U=take(U,md,1)
    V=tp(take(V,md,0));F=take(F,md,1)
    return U,S,V,F

def pUSV(U,V,ss,st=0,nm=4):
    F=abs(fftt(U,axis=0));F0=arange(len(F))/1.0/len(F)
    k=0;jj=arange(nm*3)+1
    for j in range(st,st+nm):        
        pylab.subplot(nm,3,jj[k]);pylab.plot(U[:,j]);k+=1
        pylab.subplot(nm,3,jj[k]);pylab.plot(F0,F[:,j]);k+=1
        pylab.subplot(nm,3,jj[k]);pylab.plot(V[j,:]);k+=1    
    pylab.show()

def wrtUSV(filename,BB,ss,nm=20):
    fV=open(filename+'_V'+pl,'w');fF=open(filename+'_F'+pl,'w')
    fS=open(filename+'_S'+pl,'w');fU=open(filename+'_U'+pl,'w')
    fS.write('* INDEX   S\n$ %le %le\n')
    BB=(BB-average(BB,0))/sqrt(len(BB))
    U,S,V = svd(BB,full_matrices=0)
    if nm>len(V): nm=len(V); print '# of modes =',nm
    F=abs(fftt(U,axis=0));F0=arange(len(F))/1.0/len(F)
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
    return U,S,V,F

def bMatXY(aa):
    xD={};yD={};cKeys=[];BBX=[];BBY=[]
    for j in aa.H: xD[j.name]=j.data
    for j in aa.V: yD[j.name]=j.data
    for j in xD.keys():
        if j in yD.keys():cKeys.append(j)
    for j in cKeys: BBX.append(xD[j])
    for j in cKeys: BBY.append(yD[j])
    BBX=tp(array(BBX));BBY=tp(array(BBY))
    return ctn((BBX,BBY),1)

def fmodXY(BB,m1=0,m2=1,m3=2,m4=4,tol=10.0):
    BB=(BB-average(BB,0))/sqrt(len(BB))    
    U,S,V = svd(BB,full_matrices=0);F=abs(fftt(U,axis=0))
    if S[m1]/S[m2]>tol or S[m3]/S[m4]>tol:
        print "warning: singular vals are very different"
    S=take(S,[m1,m2,m3,m4]);U=take(U,[m1,m2,m3,m4],1)
    V=tp(take(V,[m1,m2,m3,m4],0));F=take(F,[m1,m2,m3,m4],1)
    return U,S,V,F

def fPeak(UU):
    try: from SussixBS import *
    except: from SussixBS32 import * 
    print "Sussix input:",qx0,qy0,win
    sussix_inp(ir=1,turns=len(UU),tunex=qx0,tuney=qy0,\
               istun=win,idam=2,narm=100)
    Sus=sussix(UU,UU,zeros(len(UU)))
    n1=Sus.tunexy[0];n2=Sus.tunexy[1]
    print "Mode Tunes:", round(n1,4), round(n2,4)
    return n1,n2

def projA(BB,U,S,V,pln):
    if pln==0 or pln==2: n1,n2=fPeak(U[:,0])
    if pln==1: n2,n1=fPeak(U[:,0])
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

def calC21(k11,k12,k22,ph1,ph2):
    k21=[];PH1=phAdv(ph1,bm=2);PH2=phAdv(ph2,bm=2)
    for j in range(len(k12)-1):
        k21.append((-k11[j]*cos(PH1[j])*sin(PH2[j])\
                   +k22[j]*sin(PH1[j])*cos(PH2[j])\
                   +k12[j]*cos(PH1[j])*cos(PH2[j])\
                   -k12[j+1])/(sin(PH1[j])*sin(PH2[j])))
    j=len(k12)-2 #--- for first/last bpms
    k21.append((-k11[j]*cos(PH1[j])*sin(PH2[j])\
               + k22[j]*sin(PH1[j])*cos(PH2[j])\
               + k12[j]*cos(PH1[j])*cos(PH2[j])\
               - k12[0])/(sin(PH1[j])*sin(PH2[j])))
    return k21

def cMatXY(V1,XL):    
    print 'calculating C Matrix ...'
    A=sqrt(sum(V1[0:2,:]**2,0));B=sqrt(sum(V1[2:4,:]**2,0))
    Ra = V1[0:2,:]/reshape(ctn((A,A)),(2, 2*XL))
    Rb = V1[2:4,:]/reshape(ctn((B,B)),(2, 2*XL))
    Bt=A[XL:];A=A[:XL]; At=B[:XL];B=B[XL:]
    J = reshape(array([0,1,-1,0]),(2,2))
    sCb = tp(diagonal(dot(tp(Ra[:,:XL]),dot(J, Ra[:,XL:]))))
    sCa = tp(diagonal(dot(tp(Rb[:,XL:]),dot(J, Rb[:,:XL]))))
    cCb = tp(diagonal(dot(tp(Ra[:,:XL]), Ra[:,XL:])))
    cCa = tp(diagonal(dot(tp(Rb[:,XL:]),Rb[:,:XL])))

    c12 = sqrt(abs(At*Bt*sCa*sCb/(A*B)))*sign(sCa)
    c11=c12*cCa/sCa; c22=-c12*cCb/sCb
    c21=calC21(c11,c12,c22,PHX,PHY)
    #Jratio = abs((A/At/sCa)/(B/Bt/sCb))
    return c11,c12,c21,c22
    

def cMat(U1,V1,U2,V2):
    print 'calculating C matrix ...'    
    nbx=shape(V1)[1];nby=shape(V2)[1]
    A=sqrt(sum(V1[0:2,:]**2,0));At=sqrt(sum(V1[2:4,:]**2,0))
    B=sqrt(sum(V2[0:2,:]**2,0));Bt=sqrt(sum(V2[2:4,:]**2,0))
    Ra=V1[0:2,:]/reshape(ctn((A,A)),(2,nbx))
    Rca=V1[2:4,:]/reshape(ctn((At,At)),(2, nbx))
    Rb=V2[0:2,:]/reshape(ctn((B,B)),(2,nby))
    Rcb=V2[2:4,:]/reshape(ctn((Bt,Bt)),(2,nby))
        
    J = reshape(array([0,1,-1,0]),(2,2))
    cCa=tp(diagonal(dot(tp(Rb),Rca)))
    sCa=tp(diagonal(dot(tp(Rb), dot(J, Rca))))
    cCb=tp(diagonal(dot(tp(Ra),Rcb)))
    sCb=tp(diagonal(dot(tp(Ra), dot(J, Rcb))))        
    O=dot(tp(U1),ctn((U2[:,2:4], U2[:,0:2]),axis=1))
        
    O[0:2,0:2] = O[0:2,0:2]/sqrt(abs(linalg.det(O[0:2,0:2])))
    O[2:4,2:4] = O[2:4,2:4]/sqrt(abs(linalg.det(O[2:4,2:4])))
    ppp = shape(ctn((cCa, sCa),0))[0]
    
    R=dot(tp(O[2:4,2:4]),reshape(ctn((cCa, sCa)),(2,ppp/2)))
    cCa = R[0,:];sCa = R[1,:]
    R=dot(tp(O[0:2,0:2]),reshape(ctn((cCb, sCb)),(2,ppp/2)))
    cCb = R[0,:];sCb = R[1,:]
    
    c12=sqrt((At*Bt)*sCa*sCb/(A*B))*sign(sCa)
    c11=c12*cCa/sCa; c22 = -c12*cCb/sCb
    c21=calC21(c11,c12,c22,PHX,PHY)
    return c11,c12,c21,c22

def findPairs(xx,yy,sr=0):
    if len(xx)<len(yy):
        sr=1;pry=arange(len(yy));prx=[]
        for j in range(len(yy)):
            indx=argmin(abs(xx-yy[j]))
            prx.append(indx)
    else:
        prx=arange(len(xx));pry=[]
        for j in range(len(xx)):
            indx=argmin(abs(yy-xx[j]))
            pry.append(indx)
    return prx, pry

def mods():
    m4i=opt.USE.split(',')
    if len(m4i)==4:
        m1i= int(m4i[0]); m2i=int(m4i[1]);
        m3i=int(m4i[2]);m4i=int(m4i[3])
    else:
        print "WARNING: # of mode input should be 4!"
        print "will use default modes 0,1,2,3 instead"
        m1i=0;m2i=1;m3i=0;m4i=1
    return m1i,m2i,m3i,m4i

def checkphase(PH,Z,aD,OA):
    PH1=[]
    if Z=='X': PH1=[[aD.H[i].location,PH[i]] for i in range(len(aD.H))]
    else     : PH1=[[aD.V[i].location,PH[i]] for i in range(len(aD.V))]
    PH1.sort()
    p=0;m=0
    for i in range(len(PH1)-1):
        if PH1[i+1][1]>PH1[i][1]: p+=1
        else                    : m+=1
    if   OA=="LHCB1" and p<m: PH=-PH
    elif OA=="LHCB2" and p>m: PH=-PH
    return PH


if __name__ == "__main__":
    
    #-- Input options
    opt,args=parsse()
    file, NTURNS=rinput(opt.PATH)
    print "Input File:", file    
    BLABEL=0;SLABEL=0;AMP01=0;PHASE01=0

    #--- read tbt input
    aD=rhicdata(file);bx,by=bMat(aD)
    sx=array([j.location for j in aD.H])/1.0e3
    sy=array([j.location for j in aD.V])/1.0e3

    #--- define tune windows for x & y
    qsx=int((qx0-win)*len(bx));qex=int((qx0+win)*len(bx))
    qsy=int((qy0-win)*len(by));qey=int((qy0+win)*len(by))
    
    #---- write USV modes only
    pl='x'; Ux,Sx,Vx,Fx=wrtUSV(file,bx,sx);
    pl='y'; Uy,Sy,Vy,Fy=wrtUSV(file,by,sy);
    
    #--- find x & y modes in each plane
    m1,m2=findModes(Fx,qsx,qex,nm=20);
    m3,m4=findModes(Fy,qsy,qey,nm=20)
    
    if opt.CALC=="1": print "Writing SVD modes only";sys.exit()
    
    #-- Tunes, CO, PTOP & Noise
    TUNEX, TUNEY=fTune(bx,by);   
    Q1=average(TUNEX);Q2=average(TUNEY)
    print "Average {Qx, Qy}:", round(Q1,4), round(Q2,4)
    # rms or std?
    Q1RMS=sqrt(average(TUNEX)**2); Q2RMS=sqrt(average(TUNEY)**2)
    COX, PK2PKX=ptp(bx); COY, PK2PKY=ptp(by)
    #NOISEX=noiseold(bx); NOISEY=noiseold(by) #-- using old and slow way to calc
    NOISEX=noise(Ux,Sx,Vx); NOISEY=noise(Uy,Sy,Vy)

    #--- user input of modes, else default
    if opt.USE!='auto': m1,m2,m3,m4=mods()
    
    
    #--- find modes & compute phadv(x,y)
    if opt.MODE=="0": #--- ph/amp from uncoupled case
        #  NO FREEDOM FOR 3 and 4 modes
        mdx=[m1,m2,m1+2,m2+2]; mdy=[m3,m4,m3+2,m4+2]
        print "Modes Used X:",m1,m2; print "Modes Used Y:",m3,m4
        Ux,Sx,Vx,Fx=fmod(Ux,Sx,Vx,Fx,mdx)  
        print "Sing vals, X:", round(Sx[0],3),round(Sx[1],3)
        PHX,AMPX=phBeta(Sx[0]*Vx[:,0],Sx[1]*Vx[:,1])
        Uy,Sy,Vy,Fy=fmod(Uy,Sy,Vy,Fy,mdy)  
        print "Sing vals, Y:", round(Sy[0],3),round(Sy[1],3)
        PHY,AMPY=phBeta(Sy[0]*Vy[:,0],Sy[1]*Vy[:,1])
        #pUSV(Uy,tp(Vy),sx,nm=4)
        #--- just to avoid writing problems, extra zeros
        C11=zeros(len(COX)+len(COY));C12=C11;C21=C11;C22=C11
    elif opt.MODE=="1": #-- ph/amp from rotated vectors
        m11,m22=findModes(Fx,qsy,qey,nm=20);
        m33,m44=findModes(Fy,qsx,qex,nm=20)
        print "Modes Used X:",m1,m2,m11,m22;
        print "Modes Used Y:",m3,m4,m33,m44
        mdx=[m1,m2,m11,m22]; mdy=[m3,m4,m33,m44]
        Ux,Sx,Vx,Fx=fmod(Ux,Sx,Vx,Fx,mdx)  
        print "Sing vals, X:", round(Sx[0],3),round(Sx[1],3)
        Uy,Sy,Vy,Fy=fmod(Uy,Sy,Vy,Fy,mdy)  
        print "Sing vals, Y:", round(Sy[0],3),round(Sy[1],3)
        u1,v1=projA(bx,Ux,Sx,Vx,0);#pUSV(u1,v1,sx,nm=4)
        PHX,AMPX=phBeta(v1[0,:],v1[1,:])
        u2,v2=projA(by,Uy,Sy,Vy,1);#pUSV(u1,v1,sy,nm=4)
        PHY,AMPY=phBeta(v2[0,:],v2[1,:])
        try: C11,C12,C21,C22=cMat(u1,v1,u2,v2)
        except:
            C11=zeros(len(COX)+len(COY));C12=C11;C21=C11;C22=C11
            print "Coupling info not calculated, error in modes(?)"
    elif opt.MODE=="2":
        print "Mode 2 not working"; sys.exit()
        bxy=bMatXY(aD);sxy=ctn((sx,sy))
        u,s,v,f=fmodXY(bxy,m1,m2,m3,m4)
        print "Sing vals, X:", round(s[0],3),round(s[1],3)
        print "Sing vals, Y:", round(s[2],3),round(s[3],3)
        u1,v1=projA(bxy,u,s,v,2); xlen=shape(v1)[1]/2
        PHX,AMPX=phBeta(v1[0,:xlen],v1[1,:xlen])
        PHY,AMPY=phBeta(v1[2,xlen:],v1[3,xlen:])
        #pUSV(u1,v1,sxy,nm=4)
        try: C11,C12,C21,C22=cMatXY(v1,xlen)
        except:
            C11=zeros(len(COX)+len(COY));C12=C11;C21=C11;C22=C11
            print "Coupling info not available"
    else: print "no mode selected, exiting..."; sys.exit()

    #-- plot stuff
    #adx=phAdv(PHX,bm=2)
    #pylab.scatter(sx[1:],adx);pylab.show()
    #pylab.plot(sx,AMPX);pylab.show()
    #sys.exit()

    #-------------------------------------------------#
    # Temp fix of the mode issue (Ryoichi)            #
    # - Works only for a small phase advance per BPM. #
    # - Currenty only applied to the LHC.             #
    # - May have to apply to PHASE01 as well.         #
    #-------------------------------------------------#
    if opt.ACCEL[:3]=="LHC":
        PHX=checkphase(PHX,'X',aD,opt.ACCEL)
        PHY=checkphase(PHY,'Y',aD,opt.ACCEL)

    #-- write svdx & svdy
    writexy(file, aD)
    nHist(NOISEX, NOISEY,file,nbinx=45,nbiny=140)
    sys.exit()

