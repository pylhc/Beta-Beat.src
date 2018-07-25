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

import beta
import utils.bpm


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)


def ComplexSecondaryLine(delta, cw, cw1, pw, pw1):
    tp = 2.0*np.pi
    a1 = complex(1.0,-tan(tp*delta))
    a2 = cw*complex(cos(tp*pw),sin(tp*pw))
    a3 = -1.0/cos(tp*delta)*complex(0.0,1.0)
    a4 = cw1*complex(cos(tp*pw1),sin(tp*pw1))
    SL = a1*a2+a3*a4
    sizeSL = math.sqrt(SL.real**2+SL.imag**2)
    phiSL = (np.arctan2(SL.imag , SL.real)/tp) %1.0
    #SL=complex(-SL.real,SL.imag)    # This sign change in the real part is to comply with MAD output
    return [sizeSL,phiSL]


def ComplexSecondaryLineSTD(delta, cw, cw1, pw, pw1, std, std1):
    '''Calculates the propagated error on ComplexSecondaryLine'''
    tp = 2.0*np.pi
    # Aij corresponds to ai*aj from above without cw/cw1
    A12 = complex(1.0,-tan(tp*delta))*complex(cos(tp*pw),sin(tp*pw))
    A34 = -1.0/cos(tp*delta)*complex(0.0,1.0)*complex(cos(tp*pw1),sin(tp*pw1))
    # calculate the propagated STD on the output
    # the first np.abs in the sqrt suppresses an error but canNOT modify the result!
    sigma = 0.5/np.abs(A12*cw+A34*cw1)*math.sqrt(np.abs((std*(2*np.abs(A12)**2*cw+cw1*(A12*np.conj(A34)+np.conj(A12)*A34)))**2+((std1*(2*np.abs(A34)**2*cw1+cw*(np.conj(A12)*A34+A12*np.conj(A34)))))**2))
    
    return sigma


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
    S=sin(delta*tp)
    T=tan(delta*tp)
    S2_1 = 1/S**2
    # signal
    cs1=cos(tp*phase1)
    ss1=sin(tp*phase1)
    cs2=cos(tp*phase2)
    ss2=sin(tp*phase2)

    sig1=amp1*complex(cs1,ss1)
    sig2=amp2*complex(cs2,ss2)

    # computing complex secondary line (h-)
    sig=sig1*complex(1,1/T)-sig2*complex(0,1/S)

    amp=abs(sig)/2.
    phase=(np.arctan2(sig.imag,sig.real)/tp) %1.0

    # computing error secondary line (h-)
    esig = (sig1 * complex(0, S2_1 ) + sig2*complex(0, -S2_1*C))*edelta

    eamp=abs(esig)/2.
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

    bpmh=utils.bpm.intersect(list_of_zero_dpp_x)
    bpmv=utils.bpm.intersect(list_of_zero_dpp_y)
    bpmh=utils.bpm.model_intersect(bpmh, mad_twiss)
    bpmv=utils.bpm.model_intersect(bpmv, mad_twiss)


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


def _calculate_getsextupoles(twiss_d, phase_d, mad_twiss, files_dict, q1f):
    '''
    Fills the following TfsFiles:
     - getsex3000.out

    :returns: dict string --> GetllmTfsFile -- The same instace of files_dict to indicate that the dict was extended.
    '''
    print "Calculating getsextupoles"
    # For getsex1200.out andgetsex2100.out take a look at older revisions. (vimaier)

    htot, afactor, pfactor = Getsextupole(mad_twiss, twiss_d.zero_dpp_x, phase_d.ph_x, q1f, 3, 0)

    tfs_file = files_dict["getsex3000.out"]
    tfs_file.add_float_descriptor("f2h_factor", afactor)
    tfs_file.add_float_descriptor("p_f2h_factor", pfactor)
    tfs_file.add_column_names(["NAME", "S", "AMP_20", "AMP_20std", "PHASE_20", "PHASE_20std", "f3000", "f3000std", "phase_f_3000", "phase_f_3000std", "h3000", "h3000_std", "phase_h_3000", "phase_h_3000_std"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    for bpm_key in htot:
        li = htot[bpm_key]
        list_row_entries = [li[0], li[1], li[2], li[3], li[4], li[5], li[6], li[7], li[8], li[9], li[10], li[11], li[12], li[13]]
        tfs_file.add_table_row(list_row_entries)

    return files_dict

def Getsextupole(MADTwiss,amp20list,phase,tune,j,k):
    '''
    function written to calculate resonance driving terms
    '''
    # constructing complex amplitude and phase using two BPM method

    bpms=utils.bpm.intersect(amp20list)
    bpms=utils.bpm.model_intersect(bpms, MADTwiss)

    # Since beta,rmsbb(return_value[0:2]) are not used, slice return value([2:4])(vimaier)
    [bpms,invariantJx] = (beta.beta_from_amplitude(MADTwiss, amp20list, 'H'))[2:4]
    sqrt2jx=invariantJx[0]

    Q=tune+float(str(MADTwiss.Q1).split(".")[0])

    afactor=(1-cos(2*(j-k)*np.pi*Q))
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
    dbpms=utils.bpm.intersect(twiss_files[0])
    dbpms=utils.bpm.model_intersect(dbpms, MADTwiss)

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
        invariantJx = (beta.beta_from_amplitude(MADTwiss, singleFilex, 'H'))[3]

        # Since beta,rmsbb,bpms(return_value[0:3]) are not used, slice return value([3])(vimaier)
        invariantJy = (beta.beta_from_amplitude(MADTwiss, singleFiley, 'V'))[3]

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
    bpms = utils.bpm.intersect([mad_twiss])

    q_x = twi.Q1+dpp*twi.DQ1
    q_y = twi.Q2+dpp*twi.DQ2

    ftemp_name = "./TempOffMomentumTwiss.dat"
    ftemp = open(ftemp_name,"w")
    ftemp.write("@ Q1 %le "+str(q_x)+"\n")
    ftemp.write("@ Q2 %le "+str(q_y)+"\n")
    ftemp.write("@ DPP %le "+str(dpp)+"\n")
    ftemp.write("* NAME S BETX BETY ALFX ALFY MUX MUY K1L K2L\n")
    ftemp.write("$ %s %le %le  %le  %le  %le  %le %le %le %le\n")


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
        k1l = twi.K1L[twi.indx[b_n]]  # TODO: Check this formula
        k2l = twi.K2L[twi.indx[b_n]]  # TODO: Check this formula

        ftemp.write('"'+b_n+'" '+str(bpm_s)+" "+str(nbetx)+" "+str(nbety)+" "+str(nalfx)+" "+str(nalfy)+" "+str(nmux)+" "+str(nmuy)+" "+str(k1l)+" "+str(k2l)+"\n")

    ftemp.close()
    dpp_twiss = metaclass.twiss(ftemp_name, dictionary)

    # Delete temp file again(vimaier)
    os.remove(ftemp_name)

    return dpp_twiss

