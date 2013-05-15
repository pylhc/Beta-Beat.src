#!/usr/bin/env pythonafs


import pickle
from Numeric import *
from LinearAlgebra import *
import sys
from os import system
from metaclass import twiss
from random import gauss

#########################
def writeparams(deltafamilie):
#########################
    global variables
    g = open ('changeparameters', 'w')
    i=0
    for var in variables:
        g.write(var+' ='+ var+'+ ('+str(deltafamilie[i])+');\n')
        i +=1
    g.close()
    return

#########################
def writeseed(j):
#########################
    g = open ('SEEDR', 'w')
    g.write('SEEDR ='+str(j)+';\n')
    g.close()
    return

#########################
def betabeat(a,b):
#########################
    global bpmjump
    rmsx=sqrt(sum(((a.BETX-b.BETX)/b.BETX)**2)/len(b.BETX))
    rmsy=sqrt(sum(((a.BETY-b.BETY)/b.BETY)**2)/len(b.BETY))
    peakx=max(abs((a.BETX-b.BETX)/b.BETX))
    peaky=max(abs((a.BETY-b.BETY)/b.BETY))
    aveby=sum(a.BETY)/len(a.BETY)
    avebx=sum(a.BETX)/len(a.BETX)
    avebyr=sum(b.BETY)/len(b.BETY)
    avebxr=sum(b.BETX)/len(b.BETX)
    avedxbx=sum(b.DX/sqrt(b.BETX))/len(b.DX)
    avedxbxr=sum(a.DX/sqrt(a.BETX))/len(a.DX)
    avedxr=sum(a.DX)/len(a.DX)
    avedx=sum(b.DX)/len(b.DX)
    ndxbeat=sqrt(sum((b.DX/sqrt(b.BETX)-a.DX/sqrt(a.BETX))**2)/len(b.BETX))
    
    dphixa=[]
    dphiya=[]
    dphixb=[]
    dphiyb=[]
    for i in range(bpmjump,len(a.MUX)):
        dphixa.append(a.MUX[i]-a.MUX[i-bpmjump])
        dphiya.append(a.MUY[i]-a.MUY[i-bpmjump])
        dphixb.append(b.MUX[i]-b.MUX[i-bpmjump])
        dphiyb.append(b.MUY[i]-b.MUY[i-bpmjump])
    dphixa=array(dphixa)
    dphiya=array(dphiya)
    dphixb=array(dphixb)
    dphiyb=array(dphiyb)
    rmsphix=sqrt(sum((dphixa-dphixb)**2)/len(dphixa))
    rmsphiy=sqrt(sum((dphiya-dphiyb)**2)/len(dphiya))
    peakphix=max(abs(dphixa-dphixb))
    peakphiy=max(abs(dphiya-dphiyb))
    return array([rmsx,rmsy, peakx, peaky, rmsphix, rmsphiy,peakphix,peakphiy, avebx, aveby, avebxr, avebyr, avedxbx,avedxbxr, avedx,avedxr , ndxbeat])




#R=transpose(pickle.load(open('sensitivity_matrix','r')))
#zerovector=array(pickle.load(open('zerovector','r')))
#print R
#print len(R), len(R[0]), len(variables)


y=twiss('twiss.dat')
beatbefore=[]
beatafter=[]
seedstart=0
seedend=1
sigmaphase=0.25/360.
for seed in range(seedstart,seedend):
    delta=zeros(210)*1.0
#    writeparams(delta)
#    writeseed(seed)
#    system('madx < job.err4example.madx > scum ')
    x=twiss('twiss.4G.dat')
    #dphix=[]
    #dphiy=[]
    #dphixb=[]
    #dphiyb=[]
    bpmjump=2
    #for i in range(bpmjump,len(x.MUX)):
    #    dphix.append(x.MUX[i]-x.MUX[i-bpmjump]+2*gauss(0,sigmaphase))
    #    dphiy.append(x.MUY[i]-x.MUY[i-bpmjump]+2*gauss(0,sigmaphase))
    #    dphixb.append(y.MUX[i]-y.MUX[i-bpmjump])
    #    dphiyb.append(y.MUY[i]-y.MUY[i-bpmjump])

    

    #vector=array(concatenate([dphix, dphiy]))
    #delta=-matrixmultiply(generalized_inverse(R,0.04), vector-zerovector)
    #writeparams(delta)


    #system('madx < job.errors.madx > scum  ')

    #z=twiss('twiss.dat')
    bb=betabeat(x,y)
    f=open("res.dat","a")
    print >>f, bb[0], bb[1], bb[2],bb[3], bb[4],bb[5], bb[6],bb[7], bb[8],bb[9], bb[10], bb[11], bb[12],bb[13], bb[14],bb[15],bb[16]
    f.close()
    #beatafter.append(betabeat(z,y))

#for k in range(0,len(dphix)):
#    print k, dphix[k]-dphixb[k], dphiy[k]-dphiyb[k], (x.BETX[k]-y.BETX[k])/y.BETX[k],  (x.BETY[k]-y.BETY[k])/y.BETY[k]

#print  beatbefore
