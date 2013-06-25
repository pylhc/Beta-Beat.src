'''
Created on 27 May 2013

@author: vimaier

@version: 0.0.1

GetLLM.algorithms.coupling.py stores helper functions for coupling calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import math

import numpy as np
from numpy import sin, cos

import Utilities.bpm
import phase
import helper


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


#===================================================================================================
# helper-functions
#===================================================================================================

def GetCoupling1(MADTwiss, list_zero_dpp_x, list_zero_dpp_y, Q1, Q2):

    # not applicable to db=-1 for the time being...

    tp=2.0*np.pi

    # find operation point
    try:
        #TODO: There is no global outputpath. Will always crash. See Github issue #9 (vimaier)
        fdi=open(outputpath+'Drive.inp','r')  # Drive.inp file is normally in the outputpath directory in GUI operation
        for line in fdi:
            if "TUNE X" in line:
                fracxinp = line.split("=")
                fracx = fracxinp[1]
            if "TUNE Y" in line:
                fracyinp = line.split("=")
                fracy = fracyinp[1]
        fdi.close()
    except:
        fracx = Q1 # Otherwise, the fractional parts are assumed to be below 0.5
        fracy = Q2

    if fracx < 0.0 :
        fracx = 1.0 - Q1
    else:
        fracx = Q1
    if fracy < 0.0 :
        fracx = 1.0 - Q2
    else:
        fracy=Q2

    if fracx > fracy:
        sign_QxmQy = 1.0
    else:
        sign_QxmQy = -1.0

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement.
    if len(list_zero_dpp_x)!=len(list_zero_dpp_y):
        print >> sys.stderr, 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]


    XplusY=list_zero_dpp_x+list_zero_dpp_y
    dbpms=Utilities.bpm.intersect(XplusY)
    dbpms=Utilities.bpm.model_intersect(dbpms, MADTwiss)


    # caculate fw and qw, exclude bpms having wrong phases

    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)):
        bn1=str.upper(dbpms[i][1])

        fij=[]
        q1j=[]
        q2j=[]
        badbpm=0
        for j in range(0,len(list_zero_dpp_x)):
            jx=list_zero_dpp_x[j]
            jy=list_zero_dpp_y[j]
            C01ij=jx.AMP01[jx.indx[bn1]]
            C10ij=jy.AMP10[jy.indx[bn1]]
            fij.append(0.5*math.atan(math.sqrt(C01ij*C10ij)))

            #q1=(jx.MUX[jx.indx[bn1]]-jy.PHASE10[jy.indx[bn1]]+0.25)%1.0 # note that phases are in units of 2pi
            #q2=(jx.PHASE01[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]]-0.25)%1.0
            #q1=(0.5-q1)%1.0 # This sign change in the real part is to comply with MAD output
            #q2=(0.5-q2)%1.0
            q1j.append((jx.MUX[jx.indx[bn1]]-jy.PHASE10[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
            q2j.append((jx.PHASE01[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]]-0.25)%1.0)
            q1j[j]=(0.5-q1j[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2j[j]=(0.5-q2j[j])%1.0

            #if abs(q1-q2)<0.25:
            #       qij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            #       qij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
            #       badbpm=1
            #       countBadPhase += 1
            #       #print "Bad Phases in BPM ",bn1, "total so far", countBadPhase
        q1j=np.array(q1j)
        q2j=np.array(q2j)
        q1=np.average(q1j)
        q2=np.average(q2j)

        if abs(q1-q2)<0.25:  # Very rough cut !!!!!!!!!!!!!!!!!!!
            qi=(q1+q2)/2.0
        elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            qi=q1 # Note that q1 and q2 are confined 0. to 1.
        else:
            badbpm=1
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1, "total so far", countBadPhase



        if badbpm==0:
            fij=np.array(fij)
            fi=np.average(fij)
            fistd=math.sqrt(np.average(fij*fij)-(np.average(fij))**2.0+2.2e-16)
            #qij=np.array(qij)
            #qi=np.average(qij)
            #qistd=math.sqrt(np.average(qij*qij)-(np.average(qij))**2.0+2.2e-16)
            qistd=math.sqrt(np.average(q1j*q1j)-(np.average(q1j))**2.0+2.2e-16) # Not very exact...
            fi=fi*complex(cos(tp*qi),sin(tp*qi))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])
            # Trailing "0,0" in following lists because of compatibility. 
            # See issue on github pylhc/Beta-Beat.src#3
            # --vimaier
            fwqw[bn1]=[[fi,fistd,0,0],[qi,qistd,0,0]]


    dbpms=dbpmt


    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)):
        jx=list_zero_dpp_x[j]
        jy=list_zero_dpp_y[j]
        bn1=str.upper(dbpms[i][1])
        CG=CG+math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(jx.MUX[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]])


    CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
    QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]


    return [fwqw,dbpms]


def GetCoupling2(MADTwiss, list_zero_dpp_x, list_zero_dpp_y, Q1, Q2, phasex, phasey, bd, oa):
    # find operation point
    try:
        fdi = open(outputpath+'Drive.inp','r')  # Drive.inp file is normally in the outputpath directory in GUI operation
        for line in fdi:
            if "TUNE X" in line:
                fracxinp=line.split("=")
                fracx=fracxinp[1]
            if "TUNE Y" in line:
                fracyinp=line.split("=")
                fracy=fracyinp[1]
        fdi.close()
    except:
        fracx=Q1 # Otherwise, the fractional parts are assumed to be below 0.5
        fracy=Q2

    if fracx<0.0 :
        fracx=1.0-Q1
    else:
        fracx=Q1

    if fracy<0.0 :
        fracx=1.0-Q2
    else:
        fracy=Q2

    if fracx>fracy:
        sign_QxmQy=1.0
    else:
        sign_QxmQy=-1.0

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(list_zero_dpp_x)!=len(list_zero_dpp_y):
        print >> sys.stderr, 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]

    XplusY=list_zero_dpp_x+list_zero_dpp_y
    dbpms=Utilities.bpm.intersect(XplusY)
    dbpms=Utilities.bpm.model_intersect(dbpms, MADTwiss)

    # caculate fw and qw, exclude bpms having wrong phases

    tp=2.0*np.pi
    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)-1):
        bn1=str.upper(dbpms[i][1])
        bn2=str.upper(dbpms[i+1][1])

        delx= phasex[bn1][0] - 0.25  # Missprint in the coupling note
        dely= phasey[bn1][0] - 0.25

        f1001ij=[]
        #q1001ij=[]
        f1010ij=[]
        #q1010ij=[]
        q1js=[]
        q2js=[]
        q1jd=[]
        q2jd=[]
        badbpm=0
        for j in range(0,len(list_zero_dpp_x)):
            jx=list_zero_dpp_x[j]
            jy=list_zero_dpp_y[j]
            [SA0p1ij,phi0p1ij] = helper.ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]],
                    jx.PHASE01[jx.indx[bn1]], jx.PHASE01[jx.indx[bn2]])
            [SA0m1ij,phi0m1ij] = helper.ComplexSecondaryLine(delx, jx.AMP01[jx.indx[bn1]], jx.AMP01[jx.indx[bn2]],
                    -jx.PHASE01[jx.indx[bn1]], -jx.PHASE01[jx.indx[bn2]])
            [TBp10ij,phip10ij] = helper.ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]],
                    jy.PHASE10[jy.indx[bn1]], jy.PHASE10[jy.indx[bn2]])
            [TBm10ij,phim10ij] = helper.ComplexSecondaryLine(dely, jy.AMP10[jy.indx[bn1]], jy.AMP10[jy.indx[bn2]],
                    -jy.PHASE10[jy.indx[bn1]], -jy.PHASE10[jy.indx[bn2]])


            #print SA0p1ij,phi0p1ij,SA0m1ij,phi0m1ij,TBp10ij,phip10ij,TBm10ij,phim10ij
            f1001ij.append(0.5*math.sqrt(TBp10ij*SA0p1ij/2.0/2.0))
            f1010ij.append(0.5*math.sqrt(TBm10ij*SA0m1ij/2.0/2.0))

            if bd==1:
                q1jd.append((phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append((-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0)
            elif bd==-1:
                q1jd.append((phi0p1ij-jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append(-(-phip10ij+jx.MUX[jx.indx[bn1]]-0.25)%1.0)
            #print q1,q2
            q1jd[j]=(0.5-q1jd[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2jd[j]=(0.5-q2jd[j])%1.0


            #if abs(q1-q2)<0.25:
                #q1001ij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
                #q1001ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
                #badbpm=1
                #q1001ij.append(q1)
                #countBadPhase += 1
                #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

            if bd==1:
                q1js.append((phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append((phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0)
            if bd==-1:
                q1js.append((phi0m1ij+jy.MUY[jy.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append(-(phim10ij+jx.MUX[jx.indx[bn1]]+0.25)%1.0)
            #print q1,q2
            q1js[j]=(0.5-q1js[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2js[j]=(0.5-q2js[j])%1.0

            #if abs(q1-q2)<0.25:
                #q1010ij.append((q1+q2)/2.0)
            #elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
                #q1010ij.append(q1) # Note that q1 and q2 are confined 0. to 1.
            #else:
                #badbpm=1
                #if (oa=="SPS" or oa=="RHIC"):
                #       badbpm=0
                #q1010ij.append(q1)
                #countBadPhase += 1
                #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

        q1jd = np.array(q1jd)
        q2jd = np.array(q2jd)
        q1d = phase.phase_mean(q1jd,1.0)
        q2d = phase.phase_mean(q2jd,1.0)

        q1js = np.array(q1js)
        q2js = np.array(q2js)
        q1s = phase.phase_mean(q1js,1.0)
        q2s = phase.phase_mean(q2js,1.0)

        if min(abs(q1d-q2d),1.0-abs(q1d-q2d))>0.25 or min(abs(q1s-q2s),1.0-abs(q1s-q2s))>0.25:
            badbpm=1
            countBadPhase += 1

        if (oa=="SPS" or oa=="RHIC"):
            # No check for the SPS or RHIC
            badbpm=0
            q1010i=q1d
            q1010i=q1s
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1,bn2, "total so far", countBadPhase

        if badbpm==0:

            f1001ij=np.array(f1001ij)
            f1001i=np.average(f1001ij)
            f1001istd=math.sqrt(np.average(f1001ij*f1001ij)-(np.average(f1001ij))**2.0+2.2e-16)
            f1010ij=np.array(f1010ij)
            f1010i=np.average(f1010ij)
            f1010istd=math.sqrt(np.average(f1010ij*f1010ij)-(np.average(f1010ij))**2.0+2.2e-16)

            q1001i = phase.phase_mean(np.array([q1d,q2d]),1.0)
            q1010i = phase.phase_mean(np.array([q1s,q2s]),1.0)
            q1001istd = phase.calc_phase_std(np.append(q1jd,q2jd),1.0)
            q1010istd = phase.calc_phase_std(np.append(q1js,q2js),1.0)

            f1001i=f1001i*complex(cos(tp*q1001i),sin(tp*q1001i))
            f1010i=f1010i*complex(cos(tp*q1010i),sin(tp*q1010i))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])

            if bd==1:
                fwqw[bn1]=[[f1001i,f1001istd,f1010i,f1010istd],[q1001i,q1001istd,q1010i,q1010istd]]
            elif bd==-1:
                fwqw[bn1]=[[f1010i,f1010istd,f1001i,f1001istd],[q1010i,q1010istd,q1001i,q1001istd]]


    dbpms=dbpmt

    # possible correction ??
    #bn0=str.upper(dbpms[0][1])
    #up1=fwqw[bn0][0][0]
    #up2=fwqw[bn0][0][2]
    #for i in range(1,len(dbpms)):
        #bn0=str.upper(dbpms[i-1][1])
        #bn1=str.upper(dbpms[i][1])
        #df1001=math.sqrt(fwqw[bn0][0][0].real**2+fwqw[bn0][0][0].imag**2)/math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        #df1010=math.sqrt(fwqw[bn0][0][2].real**2+fwqw[bn0][0][2].imag**2)/math.sqrt(fwqw[bn1][0][2].real**2+fwqw[bn1][0][2].imag**2)
        #fwqw[bn0][0][0]=up1
        #fwqw[bn0][0][2]=up2
        #up1=complex(df1001*fwqw[bn1][0][0].real,fwqw[bn1][0][0].imag)
        #up2=complex(df1010*fwqw[bn1][0][2].real,fwqw[bn1][0][2].imag)

    #fwqw[bn1][0][0]=up1
    #fwqw[bn1][0][2]=up2
    # end of possible correction

    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)-1):
        jx=list_zero_dpp_x[0]
        jy=list_zero_dpp_y[0]
        bn1=str.upper(dbpms[i][1])
        CG=CG+math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(jx.MUX[jx.indx[bn1]]-jy.MUY[jy.indx[bn1]])

    if len(dbpms)==0:
        print >> sys.stderr, 'Warning: There is no BPM to output linear coupling properly... leaving Getcoupling.'
        fwqw['Global']=[CG,QG] #Quick fix Evian 2012
        return [fwqw,dbpms]
    else:
        CG=abs(4.0*(Q1-Q2)*CG/len(dbpms))
        QG=(QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]

    return [fwqw,dbpms]


def getCandGammaQmin(fqwq,bpms,tunex,tuney,twiss):
    # Cut the fractional part of Q1 and Q2
    QQ1 = float( int(twiss.Q1) )
    QQ2 = float( int(twiss.Q2) )

    tunex=float(tunex)+QQ1
    tuney=float(tuney)+QQ2

    tunefactor=(cos(2*np.pi*tunex)-cos(2*np.pi*tuney))/(np.pi*(sin(2*np.pi*tunex)+sin(2*np.pi*tuney)))

    coupleterms={}
    Qmin=[]

    if len(bpms)==0:
        print >> sys.stderr, "No bpms in getCandGammaQmin. Returning emtpy stuff"
        return coupleterms,0,0,bpms

    for bpm in bpms:
        bpmm=bpm[1].upper()
        detC=1-(1/(1+4*(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2)))
        check2=0.25+abs(fqwq[bpmm][0][0])**2

        if check2>abs(fqwq[bpmm][0][2])**2: # checking if sum or difference resonance is dominant!
            gamma=math.sqrt(1/(1/(1+4*(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2))))
            ffactor= 2*gamma*tunefactor*math.sqrt(abs(detC)) # cannot take abs
            C11=-(fqwq[bpmm][0][0].imag-fqwq[bpmm][0][2].imag)*2*gamma
            C12=-(fqwq[bpmm][0][0].real+fqwq[bpmm][0][2].real)*2*gamma
            C21=(fqwq[bpmm][0][0].real+fqwq[bpmm][0][2].real)*2*gamma
            C22=(fqwq[bpmm][0][0].imag-fqwq[bpmm][0][2].imag)*2*gamma
        else: # negative gamma
            gamma=-1
            ffactor=-1
            C11=C12=C21=C22=-1

        Qmin.append(ffactor)

        if (abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2)>0.0:
            err=(2*((abs(fqwq[bpmm][0][1])*abs(fqwq[bpmm][0][0]))+(abs(fqwq[bpmm][0][3])*abs(fqwq[bpmm][0][2]))))/(abs(fqwq[bpmm][0][0])**2-abs(fqwq[bpmm][0][2])**2)
        else:
            err=-1

        coupleterms[bpmm]=[detC,err,gamma,err,C11,C12,C21,C22]

    if gamma==-1:
        print "WARN: Sum resonance is dominant! "

    Qmin=np.array(Qmin)

    Qminerr=math.sqrt(np.average(Qmin*Qmin)-(np.average(Qmin))**2+2.2e-16)
    Qminav=np.average(Qmin)

    return coupleterms,Qminav,Qminerr,bpms


#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def getFreeCoupling(tunefreex,tunefreey,tunedrivenx,tunedriveny,fterm,twiss,bpms):
    if DEBUG:
        print "Calculating free fterms"
    couple={}
    couple['Global']=[fterm['Global'][0],fterm['Global'][1]]

    QQ1=float(str(twiss.Q1).split('.')[0])
    QQ2=float(str(twiss.Q2).split('.')[0])

    if(tunefreey>0.50):
        tunefreey=1-tunefreey
        tunefreey=abs(QQ2+tunefreey)
    else:
        tunefreey=abs(QQ2+abs(tunefreey))
    if(tunefreex>0.50):
        tunefreex=1-float(tunefreex)
        tunefreex=abs(QQ1+tunefreex)
    else:
        tunefreex=abs(QQ1+abs(tunefreex))

    if(tunedrivenx>0.50):
        tunedrivenx=1-tunedrivenx
    if(tunedriveny>0.50):
        tunedriveny=1-tunedriveny

    tunedrivenx=abs(QQ1+abs(tunedrivenx))
    tunedriveny=abs(QQ2+abs(tunedriveny))


    # diff f1001
    factor_top_diff=math.sqrt(sin(np.pi*(tunedrivenx-tunefreey))*sin(np.pi*(tunefreex-tunedriveny)))
    factor_bottom_diff=sin(np.pi*(tunefreex-tunefreey))

    factor_diff=abs((factor_top_diff/factor_bottom_diff))
    
    if DEBUG:
        print "Factor for coupling diff ",factor_diff

    # sum f1010
    factor_top_sum=math.sqrt(sin(np.pi*(tunedrivenx+tunefreey))*sin(np.pi*(tunefreex+tunedriveny)))
    factor_bottom_sum=sin(np.pi*(tunefreex+tunefreey))

    factor_sum=abs((factor_top_sum/factor_bottom_sum))

    if DEBUG:
        print "Factor for coupling sum ",factor_sum

    for bpm in bpms:

        bpmm=bpm[1].upper()
        [amp,phase]=fterm[bpmm]

        ampp=[amp[0]*factor_diff,amp[1],amp[2]*factor_sum,amp[3]]
        pphase=[phase[0]*factor_diff,phase[1],phase[2]*factor_sum,phase[3]]

        couple[bpmm]=[ampp,pphase]

    return couple,bpms





