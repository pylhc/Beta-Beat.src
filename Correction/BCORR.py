r"""
.. module: Python_Classes4MAD.BCORR

Created on ??

set of python routines:
 1. best corrector for beta-beat (July 28, 2008)
 2. Best N Corrector ~ like micado
 3. changeparams contains all corrs, even if zero value

Both Numeric + Numpy should work --> Not true anymore. Changed everything to numpy (vimaier)


.. moduleauthor:: Unknown
"""

import datetime
import time
import os

import numpy as np
from numpy import dot as matrixmultiply
from numpy.linalg import pinv as  generalized_inverse

default = 1


def wrtparOLD(dfam, app=0, path="./"):
    if (app == 0):
        mode = 'w'
    if (app == 1):
        mode = 'a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'changeparameters', mode)
    f = open (path+'changeparameters.tfs', mode)
    print >> f, "@", "APP", "%le", app
    print >> f, "@", "PATH", "%s", path
    print >> f, "@", "DATE", "%s", a.ctime()
    print >> f, "*", "NAME    ", "    DELTA"
    print >> f, "$", "%s    ", "    %le"
    if default == 0:
        items = sorted(dfam.iteritems(), key=lambda (k, v):(np.abs(v), k), reverse=True)
    if default == 1:
        items = [(k, v) for (k, v) in dfam.items()]
    for k, v in items:
        g.write("%12s = %12s + ( %e );\n" % (k, k, v))
        f.write("%12s   %e\n" % (k, v))
    g.close()
    f.close()


def wrtpar(corr, dfam, app=0, path="./"):
    if (app == 0):
        mode = 'w'
    if (app == 1):
        mode = 'a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open( os.path.join(path, "changeparameters"), mode )
    f = open( os.path.join(path, "changeparameters.tfs"), mode )
    print >> f, "@", "APP", "%le", app
    print >> f, "@", "PATH", "%s", path
    print >> f, "@", "DATE", "%s", a.ctime()
    print >> f, "*", "NAME    ", "    DELTA"
    print >> f, "$", "%s    ", "    %le"
    for k in corr:
        if k in dfam.keys():
            g.write("%12s = %12s + ( %e );\n" % (k, k, dfam[k]))
            #--- note minus sign for change.tfs to do corr
            f.write("%12s   %e\n" % (k, -dfam[k]))
        else:
            f.write("%12s   %e\n" % (k, 0.0))
    g.close()
    f.close()

def calcRMS(R):
    rms = np.sqrt(  np.sum(R**2) / len(R)  )
    ptp = np.max(R) - np.min(R)
    return rms, ptp

def calcRMSNumeric(R):
    rms = np.sqrt(  sum(R**2) / len(R)  )
    ptp = max(R)-min(R)
    return rms, ptp

def sortDict(adict):
    sorted_list = sorted(adict.items(), key=lambda (k, v): (v, k))
    for item in sorted_list:
        print item[0], ":", item[1]

def bCorr(X, Y, DX, beat_input, cut=0.001, app=0, path="./"):
    R = np.transpose(beat_input.sensitivity_matrix)
    b = beat_input.computevectorEXP(X, Y, DX) - beat_input.zerovector
    corr = beat_input.varslist
    m, n = np.shape(R)
    if len(b) == m and len(corr) == n:
        rms, ptop = calcRMS(b)
        inva = np.linalg.pinv(R)
        print "initial {RMS, Peak}: { %e , %e } mm" % (rms, ptop)
        print "finding best over", n, "correctors"
        for i in range(n):
            dStren = np.dot(inva[i, :], b)
            bvec = b - matrixmultiply(R[:, i], dStren)
            rm, ptp = calcRMS(bvec)
            if rm < rms:
                rms = rm
                rbest = i
                rStren = dStren
            if ptp < ptop:
                ptop = ptp
                pbest = i
                pStren = dStren
        print "final {RMS, Peak}: { %e , %e }" % (rms, ptop)
        if rbest == pbest:
            print "best corr:", corr[rbest], '%e'% (rStren), "1/m"
        else:
            print "--- warning: best corr for rms & peak are not same"
            print "RMS best corr:", corr[rbest], '%e' % (rStren), "1/m"
            print "Peak best corr:", corr[pbest], '%e' % (pStren), "1/m"
    else:
        print "dimensional mismatch in input variables"

def bCorrNumeric(X, Y, DX, beat_input, cut=0.001, app=0, path="./"):
    R = np.transpose(beat_input.sensitivity_matrix)
    b = beat_input.computevectorEXP(X, Y, DX) - beat_input.zerovector
    corr = beat_input.varslist
    m, n = np.shape(R)
    if len(b) == m and len(corr) == n:
        rms, ptop = calcRMSNumeric(b)
        inva = generalized_inverse(R, cut)
        print "initial {RMS, Peak}: { %e , %e } mm" % (rms, ptop)
        print "finding best over", n, "correctors"
        for i in range(n):
            dStren = matrixmultiply(inva[i, :], b)
            bvec = b - matrixmultiply(R[:, i], dStren)
            rm, ptp = calcRMSNumeric(bvec)
            if rm < rms:
                rms = rm
                rbest = i
                rStren = dStren
            if ptp < ptop:
                ptop = ptp
                pbest = i
                pStren = dStren
        print "final {RMS, Peak}: { %e , %e }" % (rms, ptop)
        if rbest == pbest:
            print "best corr:", corr[rbest], '%e' % (rStren), "1/m"
        else:
            print "--- warning: best corr for rms & peak are not same"
            print "RMS best corr:", corr[rbest], '%e' % (rStren), "1/m"
            print "Peak best corr:", corr[pbest], '%e' % (pStren), "1/m"
    else:
        print "dimensional mismatch in input variables"


def bNCorr(X, Y, DX, beat_input, cut=0.001, ncorr=3, app=0, tol=1e-9, path="./"):
    R = np.transpose(beat_input.sensitivity_matrix)
    n = np.shape(R)[1]
    b = beat_input.computevectorEXP(X, Y, DX) - beat_input.zerovector
    corr = beat_input.varslist
    inva = np.linalg.pinv(R, cut)
    RHO2 = {}
    rmss, ptopp = calcRMS(b)
    for ITER in range(ncorr):
        for j in range(n):
            if j not in RHO2:
                RHO = [k for k in RHO2.keys()]
                RHO.append(j)
                RR = np.take(R, RHO, 1)
                invaa = np.take(inva, RHO, 0)
                dStren = matrixmultiply(invaa, b)
                bvec = b - matrixmultiply(RR, dStren)
                rm, ptp = calcRMS(bvec)
                if rm < rmss:
                    rmss = rm
                    ptopp = ptp
                    rbest = j
                    rStren = dStren
        print 'ITER:', ITER+1, '  RMS,PK2PK:', rmss, ptopp
        RHO2[rbest] = (rmss, ptopp)
        if (rm < tol):
            print "RMS converged with", ITER, "correctors"
            break
        if (rm < rmss):
            print "stopped after", ITER, "correctors"
            break
    itr = 0
    RHO3 = {} #-- make dict RHO3={corr:strength,...}
    for j in RHO2:
        RHO3[corr[j]] = rStren[itr]
        itr += 1

    print '\n', sortDict(RHO3)
    wrtpar(corr, RHO3, app, path)
    return RHO3


def itrSVD(X, Y, DX, beat_input, cut=0.001, num_iter=1, app=0, tol=1e-9, path="./"):
    R = np.transpose(beat_input.sensitivity_matrix)
    b = beat_input.computevectorEXP(X, Y, DX) - beat_input.zerovector
    corr = beat_input.varslist
    inva = generalized_inverse(R, cut)
    rmss, ptopp = calcRMSNumeric(b)
    RHO2 = {}
    dStren = np.zeros(len(corr))
    print 'Initial Phase-Beat:', '{RMS,PK2PK}', rmss, ptopp
    for ITER in range(num_iter):  # @UnusedVariable
        dStren = dStren + matrixmultiply(inva, b)
        bvec = b - matrixmultiply(R, dStren)
        rm, ptp = calcRMSNumeric(bvec)
        print 'ITER', num_iter, '{RMS,PK2PK}', rm, ptp
    for j in range(len(corr)):
        RHO2[corr[j]] = dStren[j]
    wrtpar(corr, RHO2, app, path)


def bNCorrNumeric(X, Y, DX, beat_input, cut=0.001, ncorr=3, app=0, tol=1e-9, path="./", beta_x=None, beta_y=None):
    R = np.transpose(beat_input.sensitivity_matrix)
    n = np.shape(R)[1]
    b = beat_input.computevectorEXP(X, Y, DX, beta_x, beta_y) - beat_input.zerovector
    corr = beat_input.varslist
    inva = generalized_inverse(R, cut)
    RHO2 = {}
    rmss, ptopp = calcRMSNumeric(b)
    print 'Initial Phase-beat {RMS,PK2PK}:', rmss, ptopp
    for ITER in range(ncorr):
        for j in range(n):
            if j not in RHO2:
                RHO = [k for k in RHO2.keys()]
                RHO.append(j)
                RR = np.take(R, RHO, 1)
                invaa = np.take(inva, RHO, 0)
                dStren = matrixmultiply(invaa, b)
                #--- calculate residual due to 1..nth corrector
                bvec = b - matrixmultiply(RR, dStren)
                rm, ptp = calcRMSNumeric(bvec)
                if rm < rmss:
                    rmss = rm
                    ptopp = ptp
                    rbest = j
                    rStren = dStren
        print 'ITER:', ITER+1, '  RMS,PK2PK:', rmss, ptopp
        RHO2[rbest] = (rmss, ptopp)
        if (rm < tol):
            print "RMS converged with", ITER, "correctors"
            break
        if (rm < rmss):
            print "stopped after", ITER, "correctors"
            break
    itr = 0
    RHO3 = {} #-- make dict RHO3={corr:strength,...}
    for j in RHO2:
        RHO3[corr[j]] = rStren[itr]
        itr += 1
    print '\n', sortDict(RHO3)
    wrtpar(corr, RHO3, app, path)
    return RHO3


def bNCorrNumericSim(a, beat_input, cut=0.1, ncorr=3, app=0, tol=1e-9, path="./"):
    R = np.transpose(beat_input.sensitivity_matrix)
    n = np.shape(R)[1]
    b = beat_input.computevector(a) - beat_input.zerovector
    corr = beat_input.varslist
    inva = generalized_inverse(R, cut)
    RHO2 = {}
    rmss, ptopp = calcRMSNumeric(b)
    for ITER in range(ncorr):
        for j in range(n):
            if j not in RHO2:
                RHO = [k for k in RHO2.keys()]
                RHO.append(j)
                RR = np.take(R, RHO, 1)
                invaa = np.take(inva, RHO, 0)
                dStren = matrixmultiply(invaa, b)
                bvec = b - matrixmultiply(RR, dStren)
                rm, ptp = calcRMSNumeric(bvec)
                if rm < rmss:
                    rmss = rm
                    ptopp = ptp
                    rbest = j
                    rStren = dStren
        print 'ITER:', ITER+1, '  RMS,PK2PK:', rmss, ptopp
        RHO2[rbest] = (rmss, ptopp)
        if (rm < tol):
            print "RMS converged with", ITER, "correctors"
            break
        if (rm < rmss):
            print "stopped after", ITER, "correctors"
            break
    itr = 0
    RHO3 = {} #-- make dict RHO3={corr:strength,...}
    for j in RHO2:
        RHO3[corr[j]] = rStren[itr]
        itr += 1
    print RHO3
    wrtpar(corr, RHO3, app, path)
    return RHO3

