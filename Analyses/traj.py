from numpy import *
import os, sys,gzip,pylab, re, time
from string import split
from YASPmetaclassLHC import YASPtwiss
from operator import mod
from metaclass import twiss

def readYASP(filename):
    B1H={};B2H={};B1V={};B2V={};ctM=0
    if '.gz' in filename: f=gzip.open(filename, 'rb')
    else:f=open(filename, 'r')
    for line in f:
        cl=split(line)
        if '#' in line and "MONITOR" in line:ctM+=1
        if '*' in line and ctM<2:
            label=cl[1:]
            for j in cl[1:]:
                B1H[j]=[];B1V[j]=[];B2H[j]=[];B2V[j]=[]

        if ("BPM" in line and "B1" in line and "H" in line) :
            for j in range(2,len(cl)-1):
                B1H[label[j]].append(float(cl[j]))
        if ("BPM" in line and "B1" in line and "V" in line) :
            for j in range(2,len(cl)-1):
                B1V[label[j]].append(float(cl[j]))
        if ("BPM" in line and "B2" in line and "H" in line) :
            for j in range(2,len(cl)-1):
                B2H[label[j]].append(float(cl[j]))
        if ("BPM" in line and "B2" in line and "V" in line) :
            for j in range(2,len(cl)-1):
                B2V[label[j]].append(float(cl[j]))
    return B1H, B2H, B1V, B2V

def ttime(par):
    if "-" in par: ttm=time.strptime(par,"%H-%M-%S")[3:6]
    if ":" in par: ttm=time.strptime(par,"%H:%M:%S")[3:6]
    return ttm

def phAdv(A,b1=0,b2=1):
    U,S,V=linalg.svd(A); V=transpose(V)
    print S[:10]
    L1=S[b1]*V[:,b1]; L2=S[b2]*V[:,b2]
    #L1=V[:,b1]; L2=V[:,b2]
    ph=arctan2(L2,L1)
    ph = ph[1:]-ph[0:-1]
    for j in range(1,len(ph)): ph[j]=mod(ph[j],1)
    #ph=add.accumulate(ph[1:])
    beta = L1**2 + L2**2
    return ph/2/pi, beta


def pplot(aa):
    U,S,V=linalg.svd(aa); V=transpose(V)
    print S[:10]
    pylab.subplot(611);pylab.plot(V[:,0])
    pylab.subplot(612);pylab.plot(V[:,1])
    pylab.subplot(613);pylab.plot(V[:,2])
    pylab.subplot(614);pylab.plot(V[:,3])
    pylab.subplot(615);pylab.plot(V[:,4])
    pylab.subplot(616);pylab.plot(V[:,5])
    pylab.show()

def scoord(bb,sx,sy):
    ssx=[];ssy=[]
    for k in range(len(sx)):
        ssx.append(bb.S[bb.indx[sx[k]]])
    for k in range(len(sy)):
        ssy.append(bb.S[bb.indx[sy[k]]])
    return ssx, ssy

def sxbx(aa):
    SX=[]; BETX=[]
    for k in range(st,len(aa.HNAME)):
        SX.append(b.S[b.indx[aa.HNAME[k]]])
        BETX.append(b.BETX[b.indx[aa.HNAME[k]]])
    return array(SX), array(BETX)

def syby(aa):
    SY=[]; BETY=[]
    for k in range(st,len(aa.VNAME)):
        SY.append(b.S[b.indx[aa.VNAME[k]]])
        BETY.append(b.BETY[b.indx[aa.VNAME[k]]])
    return array(SY), array(BETY)

def phase(b):
    mux=b.MUX[1:]-b.MUX[0:-1]
    muy=b.MUY[1:]-b.MUY[0:-1]
    return mux, muy

if __name__ == "__main__":
    fileM='/afs/cern.ch/user/r/rcalaga/beta.beat/lhc.numpy/twiss.base.b1'
    fileM='twiss.ti2lhcb1.b1'
    b=twiss(fileM)
    aper=twiss('aper.ti2lhc.tfs.new')
    st=2; ex=1e-6; ey=0.5e-6
    app=[]
    for j in range(len(aper.NAME)):
        if 0<aper.APER_1[j]<1 and 0<aper.APER_2[j]<1:
            app.append([aper.S[j], aper.APER_1[j], aper.APER_2[j]])
    app=array(app)


    dir='/user/slops/data/LHC_DATA/OP_DATA/FILL_DATA/5/'
    x=[];y=[];sx=[];sy=[]
    for j in os.listdir(dir):
        if re.match('TRAJ_LHCRING_10-08-08_04',j)\
               or re.match('TRAJ_LHCRING_10-08-08_05-[0-4]',j)\
               or re.match('TRAJ_LHCRING_10-08-08_03-[3-5]',j):
            a=YASPtwiss(dir+j)
            if len(a.HCO)!=0:
                sx,betax=sxbx(a)
                orbx=a.HCO[st:]#/sqrt(2*betax*ex*1e6)
                pylab.subplot(211);pylab.plot(sx/1e3,orbx)
                #pylab.ylabel("Hor Orbit [sigx units]",fontsize=20)
                pylab.ylabel("Hor Orbit [mm]",fontsize=20)
                pylab.xlabel("Longitudinal Position [km]",fontsize=20)

                if len(a.HCO)==73:
                    x.append(a.HCO[6:-7])
        if re.match('TRAJ_LHCRING_10-08-08_0[2-3]',j) \
               or re.match('TRAJ_LHCRING_10-08-08_01-[2-5]',j):
            a=YASPtwiss(dir+j)
            if len(a.VCO)!=0:
                case=0; print len(a.VCO)
                try:
                    a.VCO[a.Vindx['BPM.33R2.B1']]=-a.VCO[a.Vindx['BPM.33R2.B1']]
                    a.VCO[a.Vindx['BPM.6L3.B1']]=-a.VCO[a.Vindx['BPM.6L3.B1']]
                except: print "BPM not found"; case=1
                sy,betay=syby(a)
                orby=a.VCO[st:]#/sqrt(2*betay*ey*1e6)
                pylab.subplot(212);pylab.plot(sy/1e3,orby);

                #pylab.ylabel("Ver Orbit [sigy units]",fontsize=20)
                pylab.ylabel("Ver Orbit [mm]",fontsize=20)
                pylab.xlabel("Longitudinal Position [km]",fontsize=20)
                if len(a.VCO)==73 and case==0:
                    y.append(a.VCO[26:-7])

    pylab.subplot(211);pylab.plot(app[:,0]/1e3,app[:,1]*1e3,'k')
    pylab.plot(app[:,0]/1e3,-app[:,1]*1e3,'k',ls='-')
    pylab.axis([3.5,7,-25,25])
    pylab.text(4,17.5,"Hor Aper Scan LHC Sector 2-3",fontsize=20)
    pylab.subplot(212);pylab.plot(app[:,0]/1e3,app[:,2]*1e3,'k')
    pylab.plot(app[:,0]/1e3,-app[:,2]*1e3,'k',ls='steps')
    pylab.axis([3.5,7,-25,25])
    pylab.text(4,12,"Ver Aper Scan LHC Sector 2-3",fontsize=20)
    #pylab.text(4,12,"Ver Bump Scan LHC IR3 region",fontsize=20)
    pylab.show();sys.exit()

    x=array(x);y=array(y);x=x[:41,:]
    print shape(x), shape(y)
    #ssx,ssy=scoord(b,sx,sy)
    #mux,muy=phase(b)
    #phix,betax=phAdv(x,1,0)
    #phiy,betay=phAdv(y,1,0)
    #pylab.plot(ssx[1:], phix)
    #pylab.plot(b.S[1:], mux);pylab.axis([1500,10000,0,.5])
    #pylab.subplot(212);pylab.plot(ssx, betax/22)
    #pylab.plot(b.S[140:230], b.BETX[140:230])
    #pylab.subplot(412);pylab.plot(phiy)

    #pylab.subplot(414);pylab.plot(betay)
    #pylab.show()
    pplot(concatenate((x,y),1))
