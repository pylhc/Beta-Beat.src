'''
Created on 5 Jun 2013

@author: ?, vimaier

@version: 0.0.2

GetLLM.algorithms.helper contains some helper functions which are used by GetLLM.py.
This module is not intended to be executed. It stores only functions.

Change history:
 - 0.0.2, vimaier, 27/06/2013:
    Extracted reated functions to new modules(phase,beta...)
'''

import sys

import Python_Classes4MAD.metaclass as metaclass
import os
import math

from numpy import sin, cos, tan
import numpy as np

import phase
import beta
import Utilities.bpm


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)




#===================================================================================================
# helper-functions
#===================================================================================================



#-------------------------

def calculate_orbit(mad_twiss, list_of_files):

    commonbpms=Utilities.bpm.intersect(list_of_files)
    commonbpms=Utilities.bpm.model_intersect(commonbpms, mad_twiss)
    co={} # Disctionary for output
    for i in range(0,len(commonbpms)):
        bpm_name=str.upper(commonbpms[i][1])
        coi=0.0
        coi2=0.0

        for tw_file in list_of_files:
            coi=coi + tw_file.CO[tw_file.indx[bpm_name]]
            coi2=coi2 + tw_file.CO[tw_file.indx[bpm_name]]**2

        coi = coi/len(list_of_files)
        corms = math.sqrt(coi2/len(list_of_files)-coi**2+2.2e-16)
        co[bpm_name] = [coi,corms]

    return [co, commonbpms]




def ComplexSecondaryLine(delta, cw, cw1, pw, pw1):
    tp=2.0*np.pi
    a1=complex(1.0,-tan(tp*delta))
    a2=cw*complex(cos(tp*pw),sin(tp*pw))
    a3=-1.0/cos(tp*delta)*complex(0.0,1.0)
    a4=cw1*complex(cos(tp*pw1),sin(tp*pw1))
    SL=a1*a2+a3*a4
    sizeSL=math.sqrt(SL.real**2+SL.imag**2)
    phiSL=(np.arctan2(SL.imag , SL.real)/tp) %1.0
    #SL=complex(-SL.real,SL.imag)    # This sign change in the real part is to comply with MAD output
    return [sizeSL,phiSL]


def ComplexSecondaryLineExtended(delta,edelta, amp1,amp2, phase1,phase2):
    '''
     Input : - delta: phase advance between two BPMs
             - edelta: error on the phase advance between two BPMs
             - amp1: amplitude of secondary line at ith BPM
             - amp2: amplitude of secondary line at i+1th BPM
             - phase1: phase of secondary line at ith BPM
             - phase2: phase of secondary line at i+1th BPM
     Return: - amp: amplitude of the complex signal
             - phase: phase of the complex signal
             - eamp: error on amplitude of the complex signal
             - ephase: error on phase of the complex signal
    '''


    # functions
    tp=2.0*np.pi
    C=cos(delta*tp)
    T=tan(delta*tp)
    SC=sin(delta*tp)/((cos(delta*tp*2)+1)/2)

    # signal
    cs1=cos(tp*phase1)
    ss1=sin(tp*phase1)
    cs2=cos(tp*phase2)
    ss2=sin(tp*phase2)

    sig1=amp1*complex(cs1,ss1)
    sig2=amp2*complex(cs2,ss2)

    # computing complex secondary line (h-)
    sig=sig1*complex(1,-T)-sig2*complex(0,1)*(1/C)

    amp=abs(sig)
    phase=(np.arctan2(sig.imag,sig.real)/tp) %1.0

    # computing error secondary line (h-)
    esig=(sig1*complex(1,-(1/C))-sig2*complex(0,1)*(SC))*edelta

    eamp=abs(esig)
    ephase=(np.arctan2(esig.imag,esig.real)/tp) %1.0

    return [amp,phase,eamp,ephase]


def pseudo_double_plane_monitors(mad_twiss, list_of_zero_dpp_x, list_of_zero_dpp_y, bpm_dictionary):
    execfile(mad_twiss.filename.replace("twiss.dat","BPMpair.py"))

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(list_of_zero_dpp_x)!=len(list_of_zero_dpp_y):
        print 'Leaving pseudo_double_plane_monitors as linx and liny files seem not correctly paired...'
        dum0={}
        dum1=[]
        return [dum0,dum1]

    bpmh=Utilities.bpm.intersect(list_of_zero_dpp_x)
    bpmv=Utilities.bpm.intersect(list_of_zero_dpp_y)
    bpmh=Utilities.bpm.model_intersect(bpmh, mad_twiss)
    bpmv=Utilities.bpm.model_intersect(bpmv, mad_twiss)


    fbpmx=[]
    fbpmy=[]
    for i in range(0,len(list_of_zero_dpp_x)):
        filex='temp'+str(i)+'_linx'
        filey='temp'+str(i)+'_liny'
        fbpmxi=open(filex,'w')
        fbpmyi=open(filey,'w')
        fbpmx.append(fbpmxi)
        fbpmy.append(fbpmyi)
        fbpmx[i].write('* NAME   S      TUNEX  MUX    AMPX   AMP01  PHASE01\n')
        fbpmy[i].write('* NAME   S      TUNEY  MUY    AMPY   AMP10  PHASE10\n')
        fbpmx[i].write('$ %s     %le    %le    %le    %le    %le    %le\n')
        fbpmy[i].write('$ %s     %le    %le    %le    %le    %le    %le\n')

        # bpmhp will be used for storing in this section. Not used further. Replaced by
        # 'tentative solution' with bpmpair. See some lines below.
        # --vimaier
#     bpmhp=[]
#     for i in range(0,len(bpmh)):
#         smin=1.0e10
#         jsave=0
#         for j in range (0,len(bpmv),10):
#             sdiff=abs(bpmh[i][0]-bpmv[j][0])
#             if sdiff<smin:
#                 smin=sdiff
#                 jsave=j

#         jlower=jsave-9
#         jupper=jsave+9
#         if jupper > len(bpmv):
#             jupper=len(bpmv)
#         for j in range (jlower,jupper):
#             sdiff=abs(bpmh[i][0]-bpmv[j][0])
#             if sdiff<smin:
#                 smin=sdiff
#                 jsave=j
#
#         bpmhp.append([bpmh[i][0],bpmh[i][1],bpmv[jsave][1],0])


    #bpmvp=[]
    #for i in range(0,len(bpmv)):
        #smin=1.0e10
        #jsave=0
        #for j in range (0,len(bpmh),10):
            #sdiff=abs(bpmv[i][0]-bpmh[j][0])
            #if sdiff<smin:
                #smin=sdiff
                #jsave=j
        #jlower=jsave-9
        #jupper=jsave+9
        #if jupper > len(bpmh):
            #jupper=len(bpmh)
        #for j in range (jlower,jupper):
            #sdiff=abs(bpmv[i][0]-bpmh[j][0])
            #if sdiff<smin:
                #smin=sdiff
                #jsave=j

        #bpmvp.append([bpmv[i][0],bpmv[i][1],bpmh[jsave][1],1])


    #dbpms=combinebpms(bpmhp,bpmvp)
    # dbpms is replaced through 'tentative solution' (vimaier)
#     dbpms=bpmhp

    # tentative solution
    dbpms = bpmpair() # model BPM name
    count_of_missing_bpms = 0
    for i in xrange(0, len(dbpms)):
        wname = str.upper(dbpms[i][1]) # horizontal BPM basis of the pairing (model name)
        pname = str.upper(dbpms[i][2]) # vertical pair of the horizontal as in SPSBPMpairs (model name)
        ws = dbpms[i][0]  # Location
        #Check whether the inputs (linx/y) have BPM name of model or experiment
        try:
            exwname = bpm_dictionary[wname][0] #Experimental BPM name of horizontal To be paired
            expname = bpm_dictionary[pname][1] #Experimental BPM name of vertical  (one of them does not exist!) to be paired

        except KeyError:
            if len(bpm_dictionary)!=0:
                count_of_missing_bpms = count_of_missing_bpms + 1
                print wname, "or", pname, "not found in the BPMdictionary. Total so far = ", count_of_missing_bpms
        try:
            for j in xrange(0, len(list_of_zero_dpp_x)):
                twiss_x = list_of_zero_dpp_x[j]
                twiss_y = list_of_zero_dpp_y[j]
                #if dbpms[i][3]==0:
                # dphix is used only in commented out code beneath (vimaier)
#                 dphix=mad_twiss.MUX[mad_twiss.indx[str.upper(pname)]]-mad_twiss.MUX[mad_twiss.indx[str.upper(wname)]]
                dphiy = mad_twiss.MUY[mad_twiss.indx[str.upper(pname)]]-mad_twiss.MUY[mad_twiss.indx[str.upper(wname)]]
                # Going to try using model names, to be able to use simulation data
                try:
                    wampx = twiss_x.AMPX[twiss_x.indx[wname]]
                    wampy = twiss_y.AMPY[twiss_y.indx[pname]]
                    wamp01 = twiss_x.AMP01[twiss_x.indx[wname]]
                    wamp10 = twiss_y.AMP10[twiss_y.indx[pname]]
                    wtunex = twiss_x.TUNEX[twiss_x.indx[wname]]
                    wtuney = twiss_y.TUNEY[twiss_y.indx[pname]]
                    wmux = twiss_x.MUX[twiss_x.indx[wname]]
                    wmuy = (twiss_y.MUY[twiss_y.indx[pname]]-dphiy)%1.0
                    if (wmuy > 0.5):
                        wmuy = wmuy-1.0
                    wphase01 = twiss_x.PHASE01[twiss_x.indx[wname]]
                    wphase10 = (twiss_y.PHASE10[twiss_y.indx[pname]]-dphiy)%1.0
                    if (wphase10 > 0.5):
                        wphase10 = wphase10-1.0
                # This seems to be experiment data, going to try with experimental names
                except:
                    wampx = twiss_x.AMPX[twiss_x.indx[exwname]]
                    wampy = twiss_y.AMPY[twiss_y.indx[expname]]
                    wamp01 = twiss_x.AMP01[twiss_x.indx[exwname]]
                    wamp10 = twiss_y.AMP10[twiss_y.indx[expname]]
                    wtunex = twiss_x.TUNEX[twiss_x.indx[exwname]]
                    wtuney = twiss_y.TUNEY[twiss_y.indx[expname]]
                    wmux = twiss_x.MUX[twiss_x.indx[exwname]]
                    wmuy = (twiss_y.MUY[twiss_y.indx[expname]]-dphiy)%1.0
                    if (wmuy > 0.5):
                        wmuy = wmuy-1.0
                    wphase01 = twiss_x.PHASE01[twiss_x.indx[exwname]]
                    wphase10 = (twiss_y.PHASE10[twiss_y.indx[expname]]-dphiy)%1.0
                    if (wphase10 > 0.5):
                        wphase10 = wphase10-1.0
                #elif dbpms[i][3]==1:
                    #wampx=twiss_x.AMPX[twiss_x.indx[pname]]
                    #wampy=twiss_y.AMPY[twiss_y.indx[wname]]
                    #wamp01=twiss_x.AMP01[twiss_x.indx[pname]]
                    #wamp10=twiss_y.AMP10[twiss_y.indx[wname]]
                    #wtunex=twiss_x.TUNEX[twiss_x.indx[pname]]
                    #wtuney=twiss_y.TUNEY[twiss_y.indx[wname]]
                    #dphix=mad_twiss.MUX[mad_twiss.indx[str.upper(pname)]]-mad_twiss.MUX[mad_twiss.indx[str.upper(wname)]]
                    #dphiy=mad_twiss.MUY[mad_twiss.indx[str.upper(pname)]]-mad_twiss.MUY[mad_twiss.indx[str.upper(wname)]]
                    #wmux=(twiss_x.MUX[twiss_x.indx[pname]]-dphix)%1.0
                    #if (wmux > 0.5): wmux=wmux-1
                    #wmuy=twiss_y.MUY[twiss_y.indx[wname]]
                    #wphase01=(twiss_x.PHASE01[twiss_x.indx[pname]]-dphix)%1.0
                    #wphase10=twiss_y.PHASE10[twiss_y.indx[wname]]
                    #if (wphase01 > 0.5): wphase01=wphase01-1
                #elif dbpms[i][3]==2:
                    #wampx=twiss_x.AMPX[twiss_x.indx[wname]]
                    #wampy=twiss_y.AMPY[twiss_y.indx[wname]]
                    #wamp01=twiss_x.AMP01[twiss_x.indx[wname]]
                    #wamp10=twiss_y.AMP10[twiss_y.indx[wname]]
                    #wtunex=twiss_x.TUNEX[twiss_x.indx[wname]]
                    #wtuney=twiss_y.TUNEY[twiss_y.indx[wname]]
                    #wmux=twiss_x.MUX[twiss_x.indx[wname]]
                    #wmuy=twiss_y.MUY[twiss_y.indx[wname]]
                    #wphase01=twiss_x.PHASE01[twiss_x.indx[wname]]
                    #wphase10=twiss_y.PHASE10[twiss_y.indx[wname]]
                fbpmx[j].write('"'+wname+'" '+str(ws)+' '+str(wtunex)+' '+str(wmux)+' '+str(wampx)+' '+str(wamp01)+' '+str(wphase01)+'\n')
                fbpmy[j].write('"'+wname+'" '+str(ws)+' '+str(wtuney)+' '+str(wmuy)+' '+str(wampy)+' '+str(wamp10)+' '+str(wphase10)+'\n')
        except:
            if len(bpm_dictionary)!=0:
                count_of_missing_bpms = count_of_missing_bpms + 1
                print wname, "or", pname, "not found in the DATA. Total so far = ", count_of_missing_bpms


    pseudo_list_x = []
    pseudo_list_y = []
    for j in range(0, len(list_of_zero_dpp_x)):
        fbpmx[j].close()
        fbpmy[j].close()
        filex = 'temp' + str(j) + '_linx'
        filey = 'temp' + str(j) +'_liny'
        pseudo_list_x.append(metaclass.twiss(filex))
        pseudo_list_y.append(metaclass.twiss(filey))

        # Delete temp files again. (vimaier)
        os.remove(filex)
        os.remove(filey)


    return [pseudo_list_x, pseudo_list_y]

#----------------------- for finding the lines of the sextupoles (@ Glenn Vanbavinckhove)
def _convert_f_term_to_h_term(amp,ampphase,termj,factor,term,M2M):
    # conversion to include _convert_f_term_to_h_term
    tp=2.0*np.pi
    H=(amp/(termj*factor))*(1-np.e**complex(0,term*tp))

    Ampi=H.imag
    Ampr=H.real
    Amp=abs(H)/M2M
    phase=(math.atan2(Ampi,Ampr))/tp


    fh=[Ampi,Ampr,Amp,phase]

    return fh


def Getsextupole(MADTwiss,amp20list,phase,tune,j,k):
    '''
    function written to calculate resonance driving terms
    '''
    # constructing complex amplitude and phase using two BPM method

    bpms=Utilities.bpm.intersect(amp20list)
    bpms=Utilities.bpm.model_intersect(bpms,MADTwiss)

    # Since beta,rmsbb(return_value[0:2]) are not used, slice return value([2:4])(vimaier)
    [bpms,invariantJx] = ( beta.beta_from_amplitude(MADTwiss,amp20list,'H') )[2:4]
    sqrt2jx=invariantJx[0]

    Q=tune+float(str(MADTwiss.Q1).split(".")[0])

    afactor=(1-cos(2*(j-k)*np.pi*Q))#(2*sin(np.pi*(j-k)*Q))
    #print (2*sin(np.pi*(j-k)*Q)),(1-cos(6*np.pi*Q))
    #sys.exit()
    pfactor=(np.pi*(j-k)*Q)

    htot={}

    for i in range(len(bpms)):

        if i<(len(bpms)-1):
            bpm=bpms[i][1]
            bpm1=bpms[i+1][1]
            s=bpms[i][0]
        else:
            bpm=bpms[i][1]
            bpm1=bpms[0][1]
            s=bpms[i][0]

        amp_i_list=[]
        phase_i_list=[]

        hlist=[]
        hplist=[]

        flist=[]
        fplist=[]

        for fileamp in amp20list:
            amp_201=fileamp.AMP_20[fileamp.indx[bpm]]*fileamp.AMPX[fileamp.indx[bpm]]
            amp_202=fileamp.AMP_20[fileamp.indx[bpm1]]*fileamp.AMPX[fileamp.indx[bpm1]]

            phase_201=fileamp.PHASE_20[fileamp.indx[bpm]]
            phase_202=fileamp.PHASE_20[fileamp.indx[bpm1]]

            delta=phase[bpm.upper()][0]-0.25
            edelta=phase[bpm.upper()][1]

            #computing complex line
            # Since eampi,ephasei(return_value[2:4]) are not used, slice return value([0:1])(vimaier)
            ampi,phasei = ( ComplexSecondaryLineExtended(delta,edelta,amp_201,amp_202,phase_201,phase_202) )[0:2]

            if ampi!=0.0:

                amp_i_list.append(ampi)
                phase_i_list.append(phasei)

                if (j==3 and k==0):
                    factor=math.sqrt(2)### factor
                    fterm=ampi/(factor*2*j*sqrt2jx**2)
                    pterm=(phasei-phase[bpm.upper()][0]+0.25)%1

                    hterm=fterm/afactor

                    hpterm=(pterm-pfactor)%1

                elif (j==2 and k==1):
                    factor=math.sqrt(2)### factor
                    fterm=ampi/(factor*2*j*sqrt2jx**2)
                    pterm=(phasei-phase[bpm][0]+0.25)%1

                    hterm=fterm/afactor

                    hpterm=(pterm-pfactor)%1

                flist.append(fterm)
                fplist.append(pterm)
                hlist.append(hterm)
                hplist.append(hpterm)

        if len(amp_i_list)!=0.0:
            al=np.mean(amp_i_list)
            alstd=np.std(amp_i_list)

            pl=np.mean(phase_i_list)
            plstd=np.mean(phasei)

            fl=np.mean(flist)
            fstd=np.std(flist)

            fpl=np.mean(fplist)
            fpstd=np.std(fplist)

            hl=np.mean(hlist)
            hstd=np.std(hlist)

            hpl=np.mean(hplist)
            hpstd=np.std(hplist)


            htot[bpm]=[bpm,s,al,alstd,pl,plstd,fl,fstd,fpl,fpstd,hl,hstd,hpl,hpstd]

    return htot,afactor,pfactor


def Getoctopole(MADTwiss,plane,twiss_files,phaseI,Q,fname,fM,NAMES):
    '''
    for finding secondary lines of the octuple (@ Glenn Vanbavinckhove)
    '''
    # intersects BPMs
    dbpms=Utilities.bpm.intersect(twiss_files[0])
    dbpms=Utilities.bpm.model_intersect(dbpms,MADTwiss)

    # value definition
    hMODELT=[]
    hMODELTi=[]
    hMODELTr=[]
    h_phase_MODELT=[]

    AT=[]
    A_RMST=[]

    phaseT=[]
    phase_RMST=[]

    hT=[]
    hTi=[]
    hTr=[]
    h_RMST=[]

    h_phaseT=[]
    h_phase_RMST=[]

    invarianceJx=[]
    invarianceJy=[]

    # finding the invariances
    for j in range(0,len(twiss_files[0])):
        singleFilex=[twiss_files[0][j]]
        singleFiley=[twiss_files[1][j]]

        # Since beta,rmsbb,bpms(return_value[0:3]) are not used, slice return value([3])(vimaier)
        invariantJx = ( beta.beta_from_amplitude(MADTwiss,singleFilex,'H') )[3]

        # Since beta,rmsbb,bpms(return_value[0:3]) are not used, slice return value([3])(vimaier)
        invariantJy = ( beta.beta_from_amplitude(MADTwiss,singleFiley,'V') )[3]

        invarianceJx.append(invariantJx)
        invarianceJy.append(invariantJy)

    # for the model
    for i in range(0,len(dbpms)):

        bpm=str.upper(dbpms[i][1])

        #TODO: think about the thing with bpm_name_according_to_tw_file = tw.NAME[tw.indx[bpm_name]]
        #Maybe change name in twiss ctor to upper.
        bpmC=MADTwiss.NAME[MADTwiss.indx[bpm]]


        for j in range(0,len(NAMES)):
            try:
                name=NAMES[j]

                if name==bpmC:

                    amp=abs(fM[j])
                    ampr=fM[i].real
                    ampi=fM[j].imag
                    phase=np.arctan2(ampi,ampr)%1

                    hMODELT.append(amp)
                    hMODELTr.append(ampr)
                    hMODELTi.append(ampi)
                    h_phase_MODELT.append(phase)

            except:
                print 'name '+str(NAMES[j])+' is not found in dictionary'
            hMODEL=[hMODELT,hMODELTi,hMODELTr,h_phase_MODELT]

    #calculation of f,q,h,qh
    for i in range(0,len(dbpms)-1):
        bn1=str.upper(dbpms[i][1])
        bn2=str.upper(dbpms[i+1][1])

        dell= phaseI[bn1][0] - 0.25

        # internal value definition
        AS=[]
        A_SRMS=[]
        phaseS=[]
        phase_RMSS=[]

        hS=[]
        hSi=[]
        hSr=[]
        h_RMSS=[]
        h_phaseS=[]
        h_phase_RMSS=[]

        for j in range(0,len(twiss_files[0])):

            single_twiss = twiss_files[0][j]

            # for f4000
            if fname=='f4000':

                [A,phi]=ComplexSecondaryLine(dell, single_twiss.AMP_30[single_twiss.indx[bn1]], single_twiss.AMP_30[single_twiss.indx[bn2]], single_twiss.PHASE_30[single_twiss.indx[bn1]], single_twiss.PHASE_30[single_twiss.indx[bn2]])

                factor=float(8*invarianceJx[j][0]**1.5)   # 1 to fit with model
                term=float(4*Q[0])
                termj=4
                M2M=0.5

            #------ converting
            h = _convert_f_term_to_h_term(A,phi,termj,factor,term,M2M)

            #----- adding the terms
            AS.append(A)
            phaseS.append(phi)
            hSi.append(h[0])
            hSr.append(h[1])
            hS.append(h[2])
            h_phaseS.append(h[3])

        # array and taking average for all the input files for one BPM
        AS=np.array(AS)
        A_SRMS=math.sqrt(np.average(AS*AS)-(np.average(AS))**2+2.2e-16)

        phaseS=np.array(phaseS)
        try:
            phase_RMSS=math.sqrt(np.average(phaseS*phaseS)-(np.average(phaseS))**2+2.2e-16)
        except:
            phase_RMSS=0

        hS=np.array(hS)
        hSi=np.array(hSi)
        hSr=np.array(hSr)
        try:
            h_RMSS=math.sqrt(np.average(hS*hS)-(np.average(hS))**2+2.2e-16)
        except:
            h_RMSS=0

        h_phaseS=np.array(h_phaseS)
        try:
            phase_rms=np.average(h_phaseS*h_phaseS)-(np.average(h_phaseS))**2+2.2e-16
        except:
            phase_rms=0
        h_phase_RMSS=math.sqrt(phase_rms)

        # real output
        AT.append(np.average(AS))
        A_RMST.append(A_SRMS)

        phaseT.append(np.average(phaseS))
        phase_RMST.append(phase_RMSS)

        hT.append(np.average(hS))
        hTi.append(np.average(hSi))
        hTr.append(np.average(hSr))
        h_RMST.append(h_RMSS)

        h_phaseT.append(np.average(h_phaseS))
        h_phase_RMST.append(h_phase_RMSS)

        A=[AT,A_RMST,phaseT,phase_RMST]
        h=[hT,hTi,hTr,h_RMST,h_phaseT,h_phase_RMST]

    return [A, h, hMODEL, dbpms]


#---- construct offmomentum
def construct_off_momentum_model(mad_twiss, dpp, dictionary):

    twi = mad_twiss # abbreviation
    bpms = Utilities.bpm.intersect([mad_twiss])

    q_x = twi.Q1+dpp*twi.DQ1
    q_y = twi.Q2+dpp*twi.DQ2

    ftemp_name = "./TempOffMomentumTwiss.dat"
    ftemp = open(ftemp_name,"w")
    ftemp.write("@ Q1 %le "+str(q_x)+"\n")
    ftemp.write("@ Q2 %le "+str(q_y)+"\n")
    ftemp.write("@ DPP %le "+str(dpp)+"\n")
    ftemp.write("* NAME S BETX BETY ALFX ALFY MUX MUY\n")
    ftemp.write("$ %s %le %le  %le  %le  %le  %le %le\n")


    for bpm in bpms:
        b_n = str.upper(bpm[1]) # bpm name
        bpm_s = bpm[0]

        # dbeta and dalpha will be extract via metaclass. As it is for the time being.
        a_x = twi.WX[twi.indx[b_n]]*cos(2.0*np.pi*twi.PHIX[twi.indx[b_n]])
        b_x = twi.WX[twi.indx[b_n]]*sin(2.0*np.pi*twi.PHIX[twi.indx[b_n]])
        b_x1 = b_x+twi.ALFX[twi.indx[b_n]]*a_x
        nbetx = twi.BETX[twi.indx[b_n]]*(1.0+a_x*dpp)
        nalfx = twi.ALFX[twi.indx[b_n]]+b_x1*dpp
        nmux = twi.MUX[twi.indx[b_n]]+twi.DMUX[twi.indx[b_n]]*dpp

        a_y = twi.WY[twi.indx[b_n]]*cos(2.0*np.pi*twi.PHIY[twi.indx[b_n]])
        b_y = twi.WY[twi.indx[b_n]]*sin(2.0*np.pi*twi.PHIY[twi.indx[b_n]])
        b_y1 = b_y+twi.ALFY[twi.indx[b_n]]*a_y
        nbety = twi.BETY[twi.indx[b_n]]*(1.0+a_y*dpp)
        nalfy = twi.ALFY[twi.indx[b_n]]+b_y1*dpp
        nmuy = twi.MUY[twi.indx[b_n]]+twi.DMUY[twi.indx[b_n]]*dpp

        ftemp.write('"'+b_n+'" '+str(bpm_s)+" "+str(nbetx)+" "+str(nbety)+" "+str(nalfx)+" "+str(nalfy)+" "+str(nmux)+" "+str(nmuy)+"\n")

    ftemp.close()
    dpp_twiss = metaclass.twiss(ftemp_name, dictionary)

    # Delete temp file again(vimaier)
    os.remove(ftemp_name)

    return dpp_twiss


#---- finding kick
def getkick(files,mad_twiss,beta_d,bbthreshold,errthreshold):

    #One source for each beta-function source
    Sources = ['model','amp','phase']

    #To count number of BPM rej in each calculation
    bpmrejx={}
    bpmrejy={}

    #meansqrt_2j-> will be mean{sqrt{2Jx}} and mean_2J-> will be mean{2J} where mean is taken over all BPMs not rejected
    meansqrt_2jx={}
    meansqrt_2jy={}
    mean_2jx={}
    mean_2jy={}

    for source in Sources:
        bpmrejx[source] = []
        bpmrejy[source] = []
        meansqrt_2jx[source] = []
        meansqrt_2jy[source] = []
        mean_2jx[source] = []
        mean_2jy[source] = []

    tunex = []
    tuney = []
    tunexRMS = []
    tuneyRMS = []

    nat_tunex = []
    nat_tuney = []
    nat_tunexRMS = []
    nat_tuneyRMS = []

    dpp=[]
    #thresholds for rejecting, input by user
    bbthreshold = float(bbthreshold)
    errthreshold = float(errthreshold)

    for j in range(0,len(files[0])):
        tw_x = files[0][j]
        tw_y = files[1][j]
        #Loop uses gen_kick_calc to get action for each beta function source in each plane
        #Note that the result for source='amp' is not currently being used
        for source in Sources:
            meansqrt_2jx_temp, mean_2jx_temp, rejected_bpm_countx = gen_kick_calc([tw_x], mad_twiss, beta_d, source, 'H', bbthreshold, errthreshold)
            meansqrt_2jy_temp, mean_2jy_temp, rejected_bpm_county = gen_kick_calc([tw_y], mad_twiss, beta_d, source, 'V', bbthreshold, errthreshold)
            meansqrt_2jx[source].append(meansqrt_2jx_temp)
            meansqrt_2jy[source].append(meansqrt_2jy_temp)
            mean_2jx[source].append(mean_2jx_temp)
            mean_2jy[source].append(mean_2jy_temp)
            bpmrejx[source].append(rejected_bpm_countx)
            bpmrejy[source].append(rejected_bpm_county)


        dpp.append(getattr(tw_x, "DPP", 0.0))

        tunex.append(getattr(tw_x, "Q1", 0.0))
        tuney.append(getattr(tw_y, "Q2", 0.0))
        tunexRMS.append(getattr(tw_x, "Q1RMS", 0.0))
        tuneyRMS.append(getattr(tw_y, "Q2RMS", 0.0))

        nat_tunex.append(getattr(tw_x, "NATQ1", 0.0))
        nat_tuney.append(getattr(tw_y, "NATQ2", 0.0))
        nat_tunexRMS.append(getattr(tw_x, "NATQ1RMS", 0.0))
        nat_tuneyRMS.append(getattr(tw_y, "NATQ2RMS", 0.0))

    tune_values_list = [tunex, tunexRMS, tuney, tuneyRMS, nat_tunex, nat_tunexRMS, nat_tuney, nat_tuneyRMS]
    return [meansqrt_2jx, meansqrt_2jy, mean_2jx, mean_2jy, tune_values_list, dpp, bpmrejx, bpmrejy]


def gen_kick_calc(list_of_files,mad_twiss,beta_d,source, plane, bbthreshold,errthreshold):
    '''
    Gets action for a given beta function source and plane
     :Parameters:
        'list_of_files': list
            list of either linx or liny files input, produced by drive
        'mad_twiss': twiss
            The model file.
        'beta_d': BetaData
            Contains measured beta functions calculated in GetLLM
        'source': string
            From where the beta function was calcuated, either model, phase, or amp
        returns:  list, list, int
            The value and uncertainty for sqrt(2J) and 2J resp., along #BPM rej
    '''
    commonbpms = Utilities.bpm.intersect(list_of_files)
    if source == 'model':
        commonbpms = Utilities.bpm.model_intersect(commonbpms, mad_twiss)
    elif source == 'amp':
        if plane == 'H':
            commonbpms = Utilities.bpm.intersect_with_bpm_list(commonbpms, beta_d.x_amp.keys())
        elif plane == 'V':
            commonbpms = Utilities.bpm.intersect_with_bpm_list(commonbpms, beta_d.y_amp.keys())
    elif source == 'phase':
        if plane == 'H':
            commonbpms = Utilities.bpm.intersect_with_bpm_list(commonbpms, beta_d.x_phase.keys())
        elif plane == 'V':
            commonbpms = Utilities.bpm.intersect_with_bpm_list(commonbpms, beta_d.y_phase.keys())

    meansqrt2j = []
    mean2j = []
    rejbpmcount = 0

    #Two different action calcs follow, one is by finding mean(sqrt(2j))-->meansqrt2j and the other is by finding mean(2j)-->mean2j
    for i in range(0, len(commonbpms)):
        bn1 = str.upper(commonbpms[i][1])
        #Gets beta function for use in action calculation- performs checks on its quality, rejects if above threshold
        if plane == 'H':
            if source == 'model':
                tembeta = mad_twiss.BETX[mad_twiss.indx[bn1]]
            elif source == 'amp':
                if beta_d.x_amp[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif abs((beta_d.x_amp[bn1][0] - mad_twiss.BETX[mad_twiss.indx[bn1]])/mad_twiss.BETX[mad_twiss.indx[bn1]]) > bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.x_amp[bn1][1]/beta_d.x_amp[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.x_amp[bn1][0]
            elif source == 'phase':
                if 	beta_d.x_phase[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif abs((beta_d.x_phase[bn1][0] - mad_twiss.BETX[mad_twiss.indx[bn1]])/mad_twiss.BETX[mad_twiss.indx[bn1]]) > bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.x_phase[bn1][1]/beta_d.x_phase[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.x_phase[bn1][0]
        elif plane == 'V':
            if source == 'model':
                tembeta = mad_twiss.BETY[mad_twiss.indx[bn1]]
            elif source == 'amp':
                if 	beta_d.y_amp[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif abs((beta_d.y_amp[bn1][0] - mad_twiss.BETY[mad_twiss.indx[bn1]])/mad_twiss.BETY[mad_twiss.indx[bn1]]) > bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.y_amp[bn1][1]/beta_d.y_amp[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.y_amp[bn1][0]
            elif source == 'phase':
                if 	beta_d.y_phase[bn1][0] <= 0:
                    rejbpmcount += 1
                    continue
                elif (abs(beta_d.y_phase[bn1][0] - mad_twiss.BETY[mad_twiss.indx[bn1]])/mad_twiss.BETY[mad_twiss.indx[bn1]])> bbthreshold:
                    rejbpmcount += 1
                    continue
                elif (beta_d.y_phase[bn1][1]/beta_d.y_phase[bn1][0]) > errthreshold:
                    rejbpmcount += 1
                    continue
                tembeta = beta_d.y_phase[bn1][0]

        meansqrt2j_i = 0.0
        mean2j_i = 0.0

        for tw_file in list_of_files:
            meansqrt2j_i += tw_file.PK2PK[tw_file.indx[bn1]]/2.
            mean2j_i = (tw_file.PK2PK[tw_file.indx[bn1]]/2)**2

        meansqrt2j_i = meansqrt2j_i/len(list_of_files)
        meansqrt2j.append(meansqrt2j_i/math.sqrt(tembeta))
        mean2j_i = mean2j_i/len(list_of_files)
        mean2j.append(mean2j_i/tembeta)

    if len(commonbpms) == rejbpmcount:
        print "Beta function calculated from", source, "in plane", plane, "is no good, no kick or action will be calculated from it."
        j_old, mean_2j = [0,0], [0,0]
        return j_old, mean_2j, rejbpmcount

    meansqrt2j = np.array(meansqrt2j)
    meansqrt2j_ave = np.average(meansqrt2j)
    meansqrt2j_std = math.sqrt(np.average(meansqrt2j*meansqrt2j)-meansqrt2j_ave**2+2.2e-16) #this is std dev, not RMS
    mean2j = np.array(mean2j)
    mean2j_ave = np.average(mean2j)
    mean2j_std = math.sqrt(np.average(mean2j*mean2j)-mean2j_ave**2+2.2e-16) #this is std dev, not RMS

    meansqrt_2j = [meansqrt2j_ave, meansqrt2j_std]
    mean_2j = [ mean2j_ave, mean2j_std]

    return meansqrt_2j, mean_2j, rejbpmcount
