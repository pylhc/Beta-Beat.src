'''
Created on 27 May 2013

@author: ?, vimaier

@version: 0.0.1

GetLLM.algorithms.coupling.py stores helper functions for coupling calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import traceback
import math

import numpy as np
from numpy import sin, cos

import Utilities.bpm
import phase
import helper
import compensate_ac_effect


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# main part
#===================================================================================================

def calculate_coupling(getllm_d, twiss_d, phase_d, tune_d, mad_twiss, mad_ac, files_dict, pseudo_list_x, pseudo_list_y):
    '''
    Calculates coupling and fills the following TfsFiles:
        getcouple.out        getcouple_free.out        getcouple_free2.out        getcoupleterms.out

    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            lhc_phase, accel, beam_direction and num_beams_for_coupling are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In/Out-param, values will be read and set)
            Holds tunes and phase advances. q1, mux, q2 and muy will be set if
            "num_beams_for_coupling == 2" and accel is 'SPS' or 'RHIC'.

    :Return: _TuneData
        the same instance as param tune_d to indicate that tunes will be set.
    '''
    print "Calculating coupling"

    if twiss_d.has_zero_dpp_x() and twiss_d.has_zero_dpp_y():
        #-- Coupling in the model
        try:
            mad_twiss.Cmatrix()
        except:
            traceback.print_exc()
        #-- Main part
        if getllm_d.num_beams_for_coupling == 1:
            # Avoids crashing the programm(vimaier)
            fwqwf = None
            fwqwf2 = None
            [fwqw, bpms] = GetCoupling1(mad_twiss, twiss_d.zero_dpp_x, twiss_d.zero_dpp_y, tune_d.q1, tune_d.q2, getllm_d.outputpath)
            tfs_file = files_dict['getcouple.out']
            tfs_file.add_float_descriptor("CG", fwqw['Global'][0])
            tfs_file.add_float_descriptor("QG", fwqw['Global'][1])
            tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD", "Q1001W", "QWSTD", "MDLF1001R", "MDLF1001I"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                try:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), (math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2)), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, mad_twiss.f1001[mad_twiss.indx[bn1]].real, mad_twiss.f1001[mad_twiss.indxy[bn1]].imag, mad_ac.f1010[mad_ac.indx[bn1]].real, mad_ac.f1010[mad_ac.indx[bn1]].imag]
                #-- Output zero if the model does not have couping parameters
                except AttributeError:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), (math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2)), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, 0.0, 0.0]
                tfs_file.add_table_row(list_row_entries)

        elif getllm_d.num_beams_for_coupling == 2:
            if getllm_d.accel == "SPS" or "RHIC" in getllm_d.accel:
                [phasexp, tune_d.q1, tune_d.mux, bpmsx] = phase.get_phases(getllm_d, mad_twiss, pseudo_list_x, None, 'H')
                [phaseyp, tune_d.q2, tune_d.muy, bpmsy] = phase.get_phases(getllm_d, mad_twiss, pseudo_list_y, None, 'V')
                [fwqw, bpms] = GetCoupling2(mad_twiss, pseudo_list_x, pseudo_list_y, tune_d.q1, tune_d.q2, phasexp, phaseyp, getllm_d.beam_direction, getllm_d.accel, getllm_d.outputpath)
            else:
                [fwqw, bpms] = GetCoupling2(mad_twiss, twiss_d.zero_dpp_x, twiss_d.zero_dpp_y, tune_d.q1, tune_d.q2, phase_d.ph_x, phase_d.ph_y, getllm_d.beam_direction, getllm_d.accel, getllm_d.outputpath)
            tfs_file = files_dict['getcouple.out']
            tfs_file.add_float_descriptor("CG", fwqw['Global'][0])
            tfs_file.add_float_descriptor("QG", fwqw['Global'][1])
            tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD1", "F1001R", "F1001I", "F1010W", "FWSTD2", "F1010R", "F1010I", "MDLF1001R", "MDLF1001I", "MDLF1010R", "MDLF1010I"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(len(bpms)):
                bn1 = str.upper(bpms[i][1])
                bns1 = bpms[i][0]
                try:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), (math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2)), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, math.sqrt(fwqw[bn1][0][2].real ** 2 + fwqw[bn1][0][2].imag ** 2), fwqw[bn1][0][3], fwqw[bn1][0][2].real, fwqw[bn1][0][2].imag, mad_ac.f1001[mad_ac.indx[bn1]].real, mad_ac.f1001[mad_ac.indx[bn1]].imag, mad_ac.f1010[mad_ac.indx[bn1]].real, mad_ac.f1010[mad_ac.indx[bn1]].imag]
                #-- Output zero if the model does not have couping parameters
                except AttributeError:
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), math.sqrt(fwqw[bn1][0][0].real ** 2 + fwqw[bn1][0][0].imag ** 2), fwqw[bn1][0][1], fwqw[bn1][0][0].real, fwqw[bn1][0][0].imag, math.sqrt(fwqw[bn1][0][2].real ** 2 + fwqw[bn1][0][2].imag ** 2), fwqw[bn1][0][3], fwqw[bn1][0][2].real, fwqw[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0]
                tfs_file.add_table_row(list_row_entries)

            #-- ac to free coupling
            if getllm_d.with_ac_calc:
                #-- analytic eqs
                try:
                    [fwqwf, bpmsf] = compensate_ac_effect.GetFreeCoupling_Eq(mad_twiss, twiss_d.zero_dpp_x, twiss_d.zero_dpp_y, tune_d.q1, tune_d.q2, tune_d.q1f, tune_d.q2f, phase_d.acphasex_ac2bpmac, phase_d.acphasey_ac2bpmac, getllm_d.beam_direction)
                    tfs_file = files_dict['getcouple_free.out']
                    tfs_file.add_float_descriptor("CG", fwqw['Global'][0])
                    tfs_file.add_float_descriptor("QG", fwqw['Global'][1])
                    tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD1", "F1001R", "F1001I", "F1010W", "FWSTD2", "F1010R", "F1010I", "MDLF1001R", "MDLF1001I", "MDLF1010R", "MDLF1010I"])
                    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                    for i in range(len(bpmsf)):
                        bn1 = str.upper(bpmsf[i][1])
                        bns1 = bpmsf[i][0]
                        try:
                            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), math.sqrt(fwqwf[bn1][0][0].real ** 2 + fwqwf[bn1][0][0].imag ** 2), fwqwf[bn1][0][1], fwqwf[bn1][0][0].real, fwqwf[bn1][0][0].imag, math.sqrt(fwqwf[bn1][0][2].real ** 2 + fwqwf[bn1][0][2].imag ** 2), fwqwf[bn1][0][3], fwqwf[bn1][0][2].real, fwqwf[bn1][0][2].imag, mad_twiss.f1001[mad_twiss.indx[bn1]].real, mad_twiss.f1001[mad_twiss.indx[bn1]].imag, mad_twiss.f1010[mad_twiss.indx[bn1]].real, mad_twiss.f1010[mad_twiss.indx[bn1]].imag] #-- Output zero if the model does not have couping parameters
                        except:
                            traceback.print_exc()
                            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), math.sqrt(fwqwf[bn1][0][0].real ** 2 + fwqwf[bn1][0][0].imag ** 2), fwqwf[bn1][0][1], fwqwf[bn1][0][0].real, fwqwf[bn1][0][0].imag, math.sqrt(fwqwf[bn1][0][2].real ** 2 + fwqwf[bn1][0][2].imag ** 2), fwqwf[bn1][0][3], fwqwf[bn1][0][2].real, fwqwf[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0]
                        tfs_file.add_table_row(list_row_entries)

                except:
                    traceback.print_exc()

                #-- global factor
                [fwqwf2, bpmsf2] = getFreeCoupling(tune_d.q1f, tune_d.q2f, tune_d.q1, tune_d.q2, fwqw, mad_twiss, bpms)
                tfs_file = files_dict['getcouple_free2.out']
                tfs_file.add_float_descriptor("CG",  fwqw['Global'][0])
                tfs_file.add_float_descriptor("QG",  fwqw['Global'][1])
                tfs_file.add_column_names(["NAME", "S", "COUNT", "F1001W", "FWSTD1", "F1001R", "F1001I", "F1010W", "FWSTD2", "F1010R", "F1010I", "MDLF1001R", "MDLF1001I", "MDLF1010R", "MDLF1010I"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(len(bpmsf2)):
                    bn1 = str.upper(bpmsf2[i][1])
                    bns1 = bpmsf2[i][0]
                    try:
                        list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), math.sqrt(fwqwf2[bn1][0][0].real ** 2 + fwqwf2[bn1][0][0].imag ** 2), fwqwf2[bn1][0][1], fwqwf2[bn1][0][0].real, fwqwf2[bn1][0][0].imag, math.sqrt(fwqwf2[bn1][0][2].real ** 2 + fwqwf2[bn1][0][2].imag ** 2), fwqwf2[bn1][0][3], fwqwf2[bn1][0][2].real, fwqwf2[bn1][0][2].imag, mad_twiss.f1001[mad_twiss.indx[bn1]].real, mad_twiss.f1001[mad_twiss.indx[bn1]].imag, mad_twiss.f1010[mad_twiss.indx[bn1]].real, mad_twiss.f1010[mad_twiss.indx[bn1]].imag] #-- Output zero if the model does not have couping parameters
                    except:
                        traceback.print_exc()
                        list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), math.sqrt(fwqwf2[bn1][0][0].real ** 2 + fwqwf2[bn1][0][0].imag ** 2), fwqwf2[bn1][0][1], fwqwf2[bn1][0][0].real, fwqwf2[bn1][0][0].imag, math.sqrt(fwqwf2[bn1][0][2].real ** 2 + fwqwf2[bn1][0][2].imag ** 2), fwqwf2[bn1][0][3], fwqwf2[bn1][0][2].real, fwqwf2[bn1][0][2].imag, 0.0, 0.0, 0.0, 0.0]
                    tfs_file.add_table_row(list_row_entries)

        else:
            raise ValueError('Number of monitors for coupling analysis should be 1 or 2 (option -n)')
        #-- Convert to C-matrix:
        if getllm_d.with_ac_calc and (fwqwf is not None or fwqwf2 is not None):
            try:
                [coupleterms, q_minav, q_minerr, bpms] = getCandGammaQmin(fwqwf, bpmsf, tune_d.q1f, tune_d.q2f, mad_twiss)
            except:
                traceback.print_exc()
                [coupleterms, q_minav, q_minerr, bpms] = getCandGammaQmin(fwqwf2, bpmsf2, tune_d.q1f, tune_d.q2f, mad_twiss)
        else:
            [coupleterms, q_minav, q_minerr, bpms] = getCandGammaQmin(fwqw, bpms, tune_d.q1f, tune_d.q2f, mad_twiss)
        tfs_file = files_dict['getcoupleterms.out']
        tfs_file.add_float_descriptor("DQMIN", q_minav)
        tfs_file.add_float_descriptor("DQMINE", q_minerr)
        tfs_file.add_column_names(["NAME", "S", "DETC", "DETCE", "GAMMA", "GAMMAE", "C11", "C12", "C21", "C22"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for bpm in bpms:
            bps = bpm[0]
            bpmm = bpm[1].upper()
            list_row_entries = [bpmm, bps, coupleterms[bpmm][0], coupleterms[bpmm][1], coupleterms[bpmm][2], coupleterms[bpmm][3], coupleterms[bpmm][4], coupleterms[bpmm][5], coupleterms[bpmm][6], coupleterms[bpmm][7]]
            tfs_file.add_table_row(list_row_entries)

    return tune_d
# END calculate_coupling ---------------------------------------------------------------------------

#===================================================================================================
# helper-functions
#===================================================================================================

def GetCoupling1(MADTwiss, list_zero_dpp_x, list_zero_dpp_y, tune_x, tune_y, outputpath):

    # not applicable to db=-1 for the time being...

    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement.
    if len(list_zero_dpp_x)!=len(list_zero_dpp_y):
        print >> sys.stderr, 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0 = {"Global":[0.0,0.0]}
        dum1 = []
        return [dum0,dum1]


    XplusY = list_zero_dpp_x+list_zero_dpp_y
    dbpms = Utilities.bpm.intersect(XplusY)
    dbpms = Utilities.bpm.model_intersect(dbpms, MADTwiss)


    # caculate fw and qw, exclude bpms having wrong phases

    fwqw = {}
    dbpmt = []
    countBadPhase = 0
    for i in range(0,len(dbpms)):
        bn1 = str.upper(dbpms[i][1])

        fij = []
        q1j = []
        q2j = []
        badbpm = 0
        for j in range(0,len(list_zero_dpp_x)):
            tw_x = list_zero_dpp_x[j]
            tw_y = list_zero_dpp_y[j]
            C01ij = tw_x.AMP01[tw_x.indx[bn1]]
            C10ij = tw_y.AMP10[tw_y.indx[bn1]]
            fij.append(0.5*math.atan(math.sqrt(C01ij*C10ij)))

            q1j.append((tw_x.MUX[tw_x.indx[bn1]]-tw_y.PHASE10[tw_y.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
            q2j.append((tw_x.PHASE01[tw_x.indx[bn1]]-tw_y.MUY[tw_y.indx[bn1]]-0.25)%1.0)
            q1j[j] = (0.5-q1j[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2j[j] = (0.5-q2j[j])%1.0

        q1j = np.array(q1j)
        q2j = np.array(q2j)
        q1 = np.average(q1j)
        q2 = np.average(q2j)

        if abs(q1-q2)<0.25:  # Very rough cut !!!!!!!!!!!!!!!!!!!
            qi = (q1+q2)/2.0
        elif abs(q1-q2)>0.75: # OK, for example q1=0.05, q2=0.95 due to measurement error
            qi = q1 # Note that q1 and q2 are confined 0. to 1.
        else:
            badbpm = 1
            countBadPhase += 1
            #print "Bad Phases in BPM ",bn1, "total so far", countBadPhase



        if badbpm == 0:
            fij = np.array(fij)
            fi = np.average(fij)
            fistd = math.sqrt(np.average(fij*fij)-(np.average(fij))**2.0+2.2e-16)
            qistd = math.sqrt(np.average(q1j*q1j)-(np.average(q1j))**2.0+2.2e-16) # Not very exact...
            fi = fi*complex(cos(2.0*np.pi*qi), sin(2.0*np.pi*qi))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])
            # Trailing "0,0" in following lists because of compatibility.
            # See issue on github pylhc/Beta-Beat.src#3
            # --vimaier
            fwqw[bn1] = [[fi,fistd,0,0],[qi,qistd,0,0]]


    dbpms = dbpmt


    # compute global values
    CG = 0.0
    QG = 0.0
    for i in range(0,len(dbpms)):
        tw_x=list_zero_dpp_x[j]
        tw_y=list_zero_dpp_y[j]
        bn1=str.upper(dbpms[i][1])
        CG=CG+math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(tw_x.MUX[tw_x.indx[bn1]]-tw_y.MUY[tw_y.indx[bn1]])

    # find operation point
    sign_QxmQy = _find_sign_QxmQy(outputpath, tune_x, tune_y)

    CG = abs(4.0*(tune_x-tune_y)*CG/len(dbpms))
    QG = (QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    fwqw['Global']=[CG,QG]


    return [fwqw,dbpms]


def GetCoupling2(MADTwiss, list_zero_dpp_x, list_zero_dpp_y, tune_x, tune_y, phasex, phasey, beam_direction, accel, outputpath):
    # check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(list_zero_dpp_x)!=len(list_zero_dpp_y):
        print >> sys.stderr, 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0={"Global":[0.0,0.0]}
        dum1=[]
        return [dum0,dum1]

    XplusY=list_zero_dpp_x+list_zero_dpp_y
    dbpms=Utilities.bpm.intersect(XplusY)
    dbpms=Utilities.bpm.model_intersect(dbpms, MADTwiss)

    # caculate fw and qw, exclude bpms having wrong phases
    fwqw={}
    dbpmt=[]
    countBadPhase=0
    for i in range(0,len(dbpms)-1):
        bn1=str.upper(dbpms[i][1])
        bn2=str.upper(dbpms[i+1][1])

        delx= phasex[bn1][0] - 0.25  # Missprint in the coupling note
        dely= phasey[bn1][0] - 0.25

        f1001ij=[]
        f1010ij=[]
        q1js=[]
        q2js=[]
        q1jd=[]
        q2jd=[]
        badbpm=0
        for j in range(0,len(list_zero_dpp_x)):
            tw_x = list_zero_dpp_x[j]
            tw_y = list_zero_dpp_y[j]
            [SA0p1ij,phi0p1ij] = helper.ComplexSecondaryLine(delx, tw_x.AMP01[tw_x.indx[bn1]], tw_x.AMP01[tw_x.indx[bn2]],
                    tw_x.PHASE01[tw_x.indx[bn1]], tw_x.PHASE01[tw_x.indx[bn2]])
            [SA0m1ij,phi0m1ij] = helper.ComplexSecondaryLine(delx, tw_x.AMP01[tw_x.indx[bn1]], tw_x.AMP01[tw_x.indx[bn2]],
                    -tw_x.PHASE01[tw_x.indx[bn1]], -tw_x.PHASE01[tw_x.indx[bn2]])
            [TBp10ij,phip10ij] = helper.ComplexSecondaryLine(dely, tw_y.AMP10[tw_y.indx[bn1]], tw_y.AMP10[tw_y.indx[bn2]],
                    tw_y.PHASE10[tw_y.indx[bn1]], tw_y.PHASE10[tw_y.indx[bn2]])
            [TBm10ij,phim10ij] = helper.ComplexSecondaryLine(dely, tw_y.AMP10[tw_y.indx[bn1]], tw_y.AMP10[tw_y.indx[bn2]],
                    -tw_y.PHASE10[tw_y.indx[bn1]], -tw_y.PHASE10[tw_y.indx[bn2]])


            #print SA0p1ij,phi0p1ij,SA0m1ij,phi0m1ij,TBp10ij,phip10ij,TBm10ij,phim10ij
            f1001ij.append(0.5*math.sqrt(TBp10ij*SA0p1ij/2.0/2.0))
            f1010ij.append(0.5*math.sqrt(TBm10ij*SA0m1ij/2.0/2.0))

            if beam_direction == 1:
                q1jd.append((phi0p1ij-tw_y.MUY[tw_y.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append((-phip10ij+tw_x.MUX[tw_x.indx[bn1]]-0.25)%1.0)
            elif beam_direction == -1:
                q1jd.append((phi0p1ij-tw_y.MUY[tw_y.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append(-(-phip10ij+tw_x.MUX[tw_x.indx[bn1]]-0.25)%1.0)
            #print q1,q2
            q1jd[j]=(0.5-q1jd[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2jd[j]=(0.5-q2jd[j])%1.0

            if beam_direction==1:
                q1js.append((phi0m1ij+tw_y.MUY[tw_y.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append((phim10ij+tw_x.MUX[tw_x.indx[bn1]]+0.25)%1.0)
            if beam_direction==-1:
                q1js.append((phi0m1ij+tw_y.MUY[tw_y.indx[bn1]]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append(-(phim10ij+tw_x.MUX[tw_x.indx[bn1]]+0.25)%1.0)
            #print q1,q2
            q1js[j]=(0.5-q1js[j])%1.0 # This sign change in the real part is to comply with MAD output
            q2js[j]=(0.5-q2js[j])%1.0

        q1jd = np.array(q1jd)
        q2jd = np.array(q2jd)
        q1d = phase.calc_phase_mean(q1jd,1.0)
        q2d = phase.calc_phase_mean(q2jd,1.0)

        q1js = np.array(q1js)
        q2js = np.array(q2js)
        q1s = phase.calc_phase_mean(q1js,1.0)
        q2s = phase.calc_phase_mean(q2js,1.0)

        if min(abs(q1d-q2d),1.0-abs(q1d-q2d))>0.25 or min(abs(q1s-q2s),1.0-abs(q1s-q2s))>0.25:
            badbpm=1
            countBadPhase += 1

        if (accel == "SPS" or accel == "RHIC"):
            # No check for the SPS or RHIC
            badbpm=0
            q1010i=q1d
            q1010i=q1s
            countBadPhase += 1

        if badbpm==0:

            f1001ij=np.array(f1001ij)
            f1001i=np.average(f1001ij)
            f1001istd=math.sqrt(np.average(f1001ij*f1001ij)-(np.average(f1001ij))**2.0+2.2e-16)
            f1010ij=np.array(f1010ij)
            f1010i=np.average(f1010ij)
            f1010istd=math.sqrt(np.average(f1010ij*f1010ij)-(np.average(f1010ij))**2.0+2.2e-16)

            q1001i = phase.calc_phase_mean(np.array([q1d,q2d]),1.0)
            q1010i = phase.calc_phase_mean(np.array([q1s,q2s]),1.0)
            q1001istd = phase.calc_phase_std(np.append(q1jd,q2jd),1.0)
            q1010istd = phase.calc_phase_std(np.append(q1js,q2js),1.0)

            f1001i=f1001i*complex(cos(2.0*np.pi*q1001i),sin(2.0*np.pi*q1001i))
            f1010i=f1010i*complex(cos(2.0*np.pi*q1010i),sin(2.0*np.pi*q1010i))
            dbpmt.append([dbpms[i][0],dbpms[i][1]])

            if beam_direction==1:
                fwqw[bn1]=[[f1001i,f1001istd,f1010i,f1010istd],[q1001i,q1001istd,q1010i,q1010istd]]
            elif beam_direction==-1:
                fwqw[bn1]=[[f1010i,f1010istd,f1001i,f1001istd],[q1010i,q1010istd,q1001i,q1001istd]]


    dbpms=dbpmt

    # compute global values
    CG=0.0
    QG=0.0
    for i in range(0,len(dbpms)-1):
        tw_x=list_zero_dpp_x[0]
        tw_y=list_zero_dpp_y[0]
        bn1=str.upper(dbpms[i][1])
        CG=CG+math.sqrt(fwqw[bn1][0][0].real**2+fwqw[bn1][0][0].imag**2)
        QG=QG+fwqw[bn1][1][0]-(tw_x.MUX[tw_x.indx[bn1]]-tw_y.MUY[tw_y.indx[bn1]])

    # find operation point
    sign_QxmQy = _find_sign_QxmQy(outputpath, tune_x, tune_y)

    if len(dbpms)==0:
        print >> sys.stderr, 'Warning: There is no BPM to output linear coupling properly... leaving Getcoupling.'
        fwqw['Global']=[CG,QG] #Quick fix Evian 2012
        return [fwqw,dbpms]
    else:
        CG=abs(4.0*(tune_x-tune_y)*CG/len(dbpms))
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
# helper
#===================================================================================================
def _find_sign_QxmQy(outputpath, tune_x, tune_y):
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
    except IOError:
        fracx=tune_x # Otherwise, the fractional parts are assumed to be below 0.5
        fracy=tune_y

    if fracx<0.0 :
        fracx=1.0-tune_x
    else:
        fracx=tune_x

    if fracy<0.0 :
        fracx=1.0-tune_y
    else:
        fracy=tune_y

    if fracx>fracy:
        return 1.0
    else:
        return -1.0

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





