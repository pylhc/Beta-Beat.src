from metaclass import *
from Numeric import *
from math import *

prefix='cern8_6kV_4kV_3.txt'

flin=twiss(prefix+'_linx')
k={}
for i in range(0,len(flin.NAME)):
    bn1=flin.NAME[i]
    A30un=flin.AMPX[flin.indx[bn1]]*flin.AMP_30[flin.indx[bn1]] # unnormalized A30
    A=flin.AMPX[flin.indx[bn1]]-3.0*A30un
    k[bn1]=4.0*A30un/A**3 # temporal k, closed orbit distortion is ignored
    k1=k[bn1]
    A=flin.AMPX[flin.indx[bn1]]-3.0*A30un-3.0*k[bn1]*flin.CO[flin.indx[bn1]]**2.0
    k[bn1]=4.0*A30un/A**3 # approximated k
    k2=k[bn1]
    A=flin.AMPX[flin.indx[bn1]]-3.0*A30un-3.0*k[bn1]*flin.CO[flin.indx[bn1]]**2.0
    k[bn1]=4.0*A30un/A**3 # iteration
    k3=k[bn1]
    print k1,k2,k3,flin.CO[flin.indx[bn1]]
    
    


fbpm=open(prefix,'r')
fbpmo=open(prefix+'_filter','w')

#print k

for line in fbpm:
    sline=line.split()
    if ('#' not in line) and sline[0]=='0':
        wline=' '+sline[0]+' '+sline[1]+' '+sline[2]
        for i in range(3,len(sline)):
            xout=float(sline[i])
            if k[sline[1]]==0.0:
                xin=xout
            else:
                #Fontana
                p=1.0/3.0/k[sline[1]]
                q=-1.0/2.0/k[sline[1]]*xout
                a=(-q+sqrt(q**2.0+p**3.0))
                b=(-q-sqrt(q**2.0+p**3.0))
                c=abs(a)**0.3333333
                if a<0: c=-c
                d=abs(b)**0.3333333
                if b<0: d=-d
                xin=c+d
            wline=wline+' '+str(xin)
        wline=wline+'\n'
    else:
        wline=line
    fbpmo.write(wline)
