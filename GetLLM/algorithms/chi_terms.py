'''
Created on 27 May 2013

@author: ?, vimaier

@version: 0.0.1

GetLLM.algorithms.chi_terms.py stores helper functions for coupling calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import math

import numpy as np
from numpy import cos, tan

import utils.bpm
import phase
import beta


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# main part
#===================================================================================================

def calculate_chiterms(getllm_d, twiss_d, mad_twiss, files_dict):
    '''
    Fills the following TfsFiles:
        getchi3000.out        getchi1010.out        getchi4000.out
        
    :Return: dict: string --> TfsFile
        The same instace of files_dict to indicate that the dict was extended

    #-> 1) chi3000
    #-> 2) chi1010
    #-> 2) chi4000
    '''
    # Designed to work only with zero_dpp_x AND zero_dpp_y files
    if not (twiss_d.has_zero_dpp_x() and twiss_d.has_zero_dpp_y()):
        return files_dict
    
    print "Calculating chiterms"
    # 1) chi3000
    tfs_file = files_dict['getchi3000.out']
    tfs_file.add_column_names(["NAME", "S", "S1", "S2", "X3000", "X3000i", "X3000r", "X3000RMS", "X3000PHASE", "X3000PHASERMS", "X3000M", "X3000Mi", "X3000Mr", "X3000MPHASE"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    files = [twiss_d.zero_dpp_x, twiss_d.zero_dpp_y]
    name = 'chi3000'
    plane = 'H'
    [dbpms, pos, xi_tot, xi_model] = get_chi_terms(mad_twiss, files, plane, name, twiss_d.zero_dpp_x, twiss_d.zero_dpp_y)
    for i in range(0, len(dbpms) - 2):
        bpm_name = str.upper(dbpms[i][1])
        list_row_entries = ['"' + bpm_name + '"', pos[0][i], pos[1][i], pos[2][i], xi_tot[0][i], xi_tot[1][i], xi_tot[2][i], xi_tot[3][i], xi_tot[4][i], xi_tot[5][i], xi_model[0][i], xi_model[1][i], xi_model[2][i], xi_model[3][i]]
        tfs_file.add_table_row(list_row_entries)
    
    # 2) chi1010
    if getllm_d.accel != 'SPS':
        tfs_file = files_dict['getchi1010.out']
        tfs_file.add_column_names(["NAME", "S", "X1010", "X1010RMS", "X1010PHASE", "X1010PHASERMS", "X1010M", "X1010MPHASE"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        name = 'chi1010'
        plane = 'H'
        [dbpms, xi_tot] = getchi1010(mad_twiss, plane, name, twiss_d.zero_dpp_x, twiss_d.zero_dpp_y)
        for i in range(len(dbpms) - 2):
            bpm_name = str.upper(dbpms[i][1])
            bns = dbpms[i][0]
            list_row_entries = ['"' + bpm_name + '"', bns, xi_tot[0][i], xi_tot[1][i], xi_tot[2][i], xi_tot[3][i], 0, 0]
            tfs_file.add_table_row(list_row_entries)
    
    # 1) chi4000
    # for getchi4000.out take a look at previous version of GetLLM. Was commented out so I deleted it (vimaier)
    
    return files_dict
# END calculate_chiterms ---------------------------------------------------------------------------


#===================================================================================================
# helper-functions
#===================================================================================================

def _compute_chi_terms(amp,phase_20,phase,terms,J,plane,ima,rea):
    ''' for finding the chi terms '''

    #computes the chiterms for different inputs
    twoPi=2*np.pi

    delta1=((phase[1]-phase[0]-0.25)*twoPi)
    delta2=((phase[2]-phase[1]-0.25)*twoPi)

    inp=0.13 # ????
    #term1=((amp[0]*np.e**complex(0,phase_20[0]*twoPi)))/cos(delta1)
    #term2=((amp[1]*np.e**complex(0,phase_20[1]*twoPi)))*(tan(delta1)+tan(delta2))
    #term3=((amp[2]*np.e**complex(0,phase_20[2]*twoPi)))/cos(delta2)
    term1=((amp[0]*np.e**complex(0,(phase_20[0]+inp)*twoPi)))/cos(delta1)
    term2=((amp[1]*np.e**complex(0,(phase_20[1]+inp)*twoPi)))*(tan(delta1)+tan(delta2))
    term3=((amp[2]*np.e**complex(0,(phase_20[2]+inp)*twoPi)))/cos(delta2)
    chiTOT=(term1+term2+term3)

    chiAMP=abs(chiTOT)

    chiAMPi=chiTOT.imag
    chiAMPr=chiTOT.real
    chiPHASE=(((np.arctan2(chiTOT.imag,chiTOT.real)))/twoPi)%1
    #chiPHASE=(0.5-chiPHASE)%1



    JX=J[0]**(2.*(terms[0]+terms[1]-2.)/2.)
    JY=J[1]**(2.*(terms[2]+terms[3])/2.)


    Invariance=JX*JY
    Facot4AMP=Invariance*4/2 # to for conversion complex, invariance = ((2*JX)^(j+k-2)/2)*((2*JY)^(l+m)/2)


    chiAMP=chiAMP/Facot4AMP
    chiAMPi=chiAMPi/Facot4AMP
    chiAMPr=chiAMPr/Facot4AMP

    return [chiAMP,chiAMPi,chiAMPr,chiPHASE]


def get_chi_terms(MADTwiss,filesF,plane,name,ListOfZeroDPPX,ListOfZeroDPPY):

    # bmps
    files = filesF[0]

    dbpms = utils.bpm.intersect(files)
    dbpms = utils.bpm.model_intersect(dbpms, MADTwiss)


    # initiliasing variables
    XIT=[]
    XITi=[]
    XITr=[]
    XIrmsT=[]
    XI_phase_T=[]
    XI_phaseRMS_T=[]

    POS1=[]
    POS2=[]
    POS3=[]

    XITMODEL=[]
    XITMODELi=[]
    XITMODELr=[]
    XITMODEL_phase=[]

    BPMS=[]
    invarianceJx=[]
    invarianceJy=[]

    for i in range(0,len(dbpms)): # ask rogelio
        bn1=str.upper(dbpms[i][1])
        BPMS.append(bn1)

    #### invariance
    for j in range(0,len(ListOfZeroDPPX)):
        # Since betax,rmsbbx,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJX = ( beta.beta_from_amplitude(MADTwiss,ListOfZeroDPPX,'H') )[3]
        # Since betay,rmsbby,bpms(return_value[0:3]) are not used, slice the return value([3]) (vimaier)
        invariantJY= ( beta.beta_from_amplitude(MADTwiss,ListOfZeroDPPY,'V') )[3]
        invarianceJx.append(invariantJX[0])
        invarianceJy.append(invariantJY[0])

    if DEBUG:
        print "invarianceJX:",invarianceJx
    #### model chi
    MADTwiss.chiterms(BPMS)
    if name=='chi3000':
        MODEL=MADTwiss.chi
    elif name=='chi4000':
        MODEL=MADTwiss.chi4000


    for i in range(0,len(MODEL)):

        MODEL[i]=MODEL[i]
        amp=abs(MODEL[i])
        ampi=MODEL[i].imag
        ampr=MODEL[i].real

        if(MODEL[i].real==0. ):

            phase=0

        else:

            phase=np.arctan2(MODEL[i].imag,MODEL[i].real)%1


        XITMODEL.append(amp)
        XITMODELi.append(ampi)
        XITMODELr.append(ampr)
        XITMODEL_phase.append(phase)

    XIMODEl=[XITMODEL,XITMODELi,XITMODELr,XITMODEL_phase]

    for i in range(0,len(dbpms)-2):

        XI=[]
        XIi=[]
        XIr=[]
        XIrms=[]
        XI_phase=[]
        XI_phaseRMS=[]

        bn1=str.upper(dbpms[i][1])
        bn2=str.upper(dbpms[i+1][1])
        bn3=str.upper(dbpms[i+2][1])

        filej=ListOfZeroDPPX[0]

        pos1=filej.S[filej.indx[bn1]]
        pos2=filej.S[filej.indx[bn2]]
        pos3=filej.S[filej.indx[bn3]]


        POS1.append(pos1)
        POS2.append(pos2)
        POS3.append(pos3)

        imaM=XITMODELi[i]
        realM=XITMODELr[i]

        for j in range(0,len(files)):
            jx=files[j]




            # for chi3000
            if name=='chi3000':
                phase1=jx.PHASE_20[jx.indx[bn1]]
                phase2=jx.PHASE_20[jx.indx[bn2]]
                phase3=jx.PHASE_20[jx.indx[bn3]]
                phase_SL=[phase1,phase2,phase3]

                terms=[3,0,0,0]
                amp1=jx.AMP_20[jx.indx[bn1]]
                amp2=jx.AMP_20[jx.indx[bn2]]
                amp3=jx.AMP_20[jx.indx[bn3]]
                amp=[amp1,amp2,amp3]

            # for chi4000
            elif name=='chi4000':
                phase1=jx.PHASE_30[jx.indx[bn1]]
                phase2=jx.PHASE_30[jx.indx[bn2]]
                phase3=jx.PHASE_30[jx.indx[bn3]]
                phase_SL=[phase1,phase2,phase3]

                terms=[4,0,0,0]
                amp1=jx.AMP_30[jx.indx[bn1]]
                amp2=jx.AMP_30[jx.indx[bn2]]
                amp3=jx.AMP_30[jx.indx[bn3]]
                amp=[amp1,amp2,amp3]

            phase11=jx.MUX[jx.indx[bn1]]
            phase12=jx.MUX[jx.indx[bn2]]
            phase13=jx.MUX[jx.indx[bn3]]
            phase=[phase11,phase12,phase13]


            J=[ invarianceJx[j],invarianceJy[j]]



            chi = _compute_chi_terms(amp,phase_SL,phase,terms,J,'H',imaM,realM)



            XI.append(chi[0])
            XIi.append(chi[1])
            XIr.append(chi[2])
            XI_phase.append(chi[3])



        XI=np.array(XI)
        XIi=np.array(XIi)
        XIr=np.array(XIr)
        try:
            XIrms=math.sqrt(np.average(XI*XI)-np.average(XI)**2+2.2e-16)
        except:
            XIrms=0
        XI_phase=np.array(XI_phase)
        try:
            XI_phaseRMS=math.sqrt(np.average(XI_phase*XI_phase)-np.average(XI_phase)**2+2.2e-16)
        except:
            XI_phaseRMS=0


        XIT.append(np.average(XI))
        XITi.append(np.average(XIi))
        XITr.append(np.average(XIr))
        XIrmsT.append(XIrms)
        XI_phase_T.append(np.average(XI_phase))
        XI_phaseRMS_T.append(XI_phaseRMS)

        POS=[POS1,POS2,POS3]

        XItot=[XIT,XITi,XITr,XIrmsT,XI_phase_T,XI_phaseRMS_T]

    return [dbpms,POS,XItot,XIMODEl]

def getchi1010(MADTwiss, plane, name, files_zero_dpp_x, files_zero_dpp_y):
    if len(files_zero_dpp_x) != len(files_zero_dpp_y):
        print "Different length of x, y files. Leaving getchi1010 with empty values."
        return [[], []]
    
    dbpms=utils.bpm.intersect(files_zero_dpp_x + files_zero_dpp_y)
    dbpms=utils.bpm.model_intersect(dbpms, MADTwiss)

    dbpmsy=utils.bpm.intersect(files_zero_dpp_y + files_zero_dpp_x)
    dbpmsy=utils.bpm.model_intersect(dbpmsy, MADTwiss)


    # initiliasing variables
    XIT=[]
    XIrmsT=[]
    XI_phase_T=[]
    XI_phaseRMS_T=[]

    XItot = []
    
    for i in range(0,len(dbpms)):

        XI=[]
        XIrms=[]
        XI_phase=[]
        XI_phaseRMS=[]

        bn=str.upper(dbpms[i][1])
        bny=str.upper(dbpmsy[i][1])

        for j in range(0,len(files_zero_dpp_x)):

            jx=files_zero_dpp_x[j]

            jy=files_zero_dpp_y[j]

            amp10x=jx.AMP01[jx.indx[bn]]
            amp10y=jy.AMP10[jy.indx[bny]]
            phase10x=jx.PHASE01[jx.indx[bn]]
            phasex=jx.MUX[jx.indx[bn]]

            XI1010=0.25*math.sqrt(amp10x*amp10y)
            phase1010=phase10x+phasex

            XI.append(XI1010)
            XI_phase.append(phase1010)

        XI=np.array(XI)
        XIrms=math.sqrt(np.average(XI*XI)-np.average(XI)**2+2.2e-16)
        XI_phase=np.array(XI_phase)
        XI_phaseRMS=math.sqrt(np.average(XI_phase*XI_phase)-np.average(XI_phase)**2+2.2e-16)


        XIT.append(np.average(XI))
        XIrmsT.append(XIrms)
        XI_phase_T.append(np.average(XI_phase))
        XI_phaseRMS_T.append(XI_phaseRMS)

        XItot=[XIT,XIrmsT,XI_phase_T,XI_phaseRMS_T]

    return [dbpms,XItot]





