import sys
import os

import __init__  # @UnusedImport adds path to Beta-Beat.src
from Python_Classes4MAD.metaclass import *


path=sys.argv[1]


t2=twiss(path+'/twiss_cor.dat')
t3=twiss(path+'/twiss_no.dat')
t2.Cmatrix()
t3.Cmatrix()



#############################
#
# normal quad
#
#############################

f1=open(path+"/bbx.out","w")
f2=open(path+"/bby.out","w")

print >> f1,"NAME S MEA ERROR MODEL"
print >> f1,"%s %le %le %le %le"

print >> f2,"NAME S MEA ERROR MODEL"
print >> f2,"%s %le %le %le %le"

if os.path.exists(path+'/getbetax_free.out'):
    t1=twiss(path+'/getbetax_free.out')
else:
    t1=twiss(path+'/getbetax.out')

for i in range(len(t1.NAME)):
    el=t1.NAME[i]
    a=1
    try:
        el2=t2.NAME[t2.indx[el]]
    except:
        print "No ", el
        a=0
    if a==1:
        j=t2.indx[el]
        print>> f1,el, t1.S[i], (t1.BETX[i]-t1.BETXMDL[i])/t1.BETXMDL[i], t1.STDBETX[i]/t1.BETXMDL[i],(t2.BETX[j]-t3.BETX[j])/t3.BETX[j]

if os.path.exists(path+'/getbetay_free.out'):
    t1=twiss(path+'/getbetay_free.out')
else:
    t1=twiss(path+'/getbetay.out')

for i in range(len(t1.NAME)):
    el=t1.NAME[i]
    a=1
    try:
        el2=t2.NAME[t2.indx[el]]
    except:
        print "No ", el
        a=0
    if a==1:
        j=t2.indx[el]
        print>> f2,el, t1.S[i], (t1.BETY[i]-t1.BETYMDL[i])/t1.BETYMDL[i], t1.STDBETY[i]/t1.BETYMDL[i],(t2.BETY[j]-t3.BETY[j])/t3.BETY[j]

f1.close()
f2.close()

dx=1
try:
    t1=twiss(path+'/getDx.out')
except:
    print "NO dispersion"
    dx=0

if dx==1:
    f3=open(path+"/dx.out","w")

    print >> f3,"NAME S MEA ERROR MODEL"
    print >> f3,"%s %le %le %le %le"

    for i in range(len(t1.NAME)):
        el=t1.NAME[i]
        a=1
        try:
            el2=t2.NAME[t2.indx[el]]
        except:
            print "No ", el
            a=0
        if a==1:
            j=t2.indx[el]
            print>> f3,el, t1.S[i], (t1.DX[i]-t1.DXMDL[i]), t1.STDDX[i],(t2.DX[j]-t3.DX[j])

    f3.close()

#############################
#
# skew quad
#
#############################
f1=open(path+"/couple.out","w")


print >> f1,"NAME S F1001re F1001im F1001e F1001re_m F1001im_m"
print >> f1,"%s %le %le %le %le %le %le"



if os.path.exists(path+'/getcouple_free.out'):
    t1=twiss(path+'/getcouple_free.out')
else:
    t1=twiss(path+'/getcouple.out')


for i in range(len(t1.NAME)):
    el=t1.NAME[i]
    a=1
    try:
        el2=t2.NAME[t2.indx[el]]
    except:
        print "No ", el
        a=0
    if a==1:
        j=t2.indx[el]
        print>> f1,el, t1.S[i],t1.F1001R[i],t1.F1001I[i],t1.FWSTD1[i],t2.f1001[j].real,t2.f1001[j].imag




dy=1
try:
    t1=twiss(path+'/getDy.out')
except:
    print "NO dispersion"
    dy=0

if dy==1:
    f3=open(path+"/dy.out","w")

    print >> f3,"NAME S MEA ERROR MODEL"
    print >> f3,"%s %le %le %le %le"

    for i in range(len(t1.NAME)):
        el=t1.NAME[i]
        a=1
        try:
            el2=t2.NAME[t2.indx[el]]
        except:
            print "No ", el
            a=0
        if a==1:
            j=t2.indx[el]
            print>> f3,el, t1.S[i], (t1.DY[i]-t1.DYMDL[i]), t1.STDDY[i],(t2.DY[j]-t3.DY[j])

    f3.close()
