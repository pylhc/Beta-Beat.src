'''
Created on 27 May 2013

@author: vimaier

@version: 0.0.1

GetLLM.algorithms.phase.py stores helper functions for phase calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>
'''

import sys
import math
import traceback

import numpy as np
from numpy import sin, cos, tan

import Utilities.bpm
import compensate_ac_effect


DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# main part
#===================================================================================================

class BetaData(object):
    """ File for storing results from beta computations. """

    def __init__(self):
        self.x_phase = None # beta x from phase
        self.x_phase_f = None # beta x from phase free
        self.y_phase = None # beta y from phase
        self.y_phase_f = None # beta y from phase free
        
        self.x_amp = None # beta x from amplitude
        self.y_amp = None # beta y from amplitude
        
        self.x_ratio = None # beta x ratio
        self.x_ratio_f = None # beta x ratio free
        self.y_ratio = None # beta x ratio
        self.y_ratio_f = None # beta x ratio free

def calculate_beta_from_phase(getllm_d, twiss_d, tune_d, phase_d, mad_twiss, mad_ac, files_dict):
    '''
    Calculates beta and fills the following TfsFiles:
        getbetax.out        getbetax_free.out        getbetax_free2.out
        getbetay.out        getbetay_free.out        getbetay_free2.out
        
    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            lhc_phase, accel and beam_direction are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances
        'phase_d': _PhaseData (In-param, values will only be read)
            Holds results from get_phases
    '''
    beta_d = BetaData()
    
    print 'Calculating beta'
    #---- H plane
    if twiss_d.has_zero_dpp_x():
        [beta_d.x_phase, rmsbbx, alfax, bpms] = beta_from_phase(mad_ac, twiss_d.zero_dpp_x, phase_d.ph_x, 'H')
        beta_d.x_phase['DPP'] = 0
        tfs_file = files_dict['getbetax.out']
        tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1))
        tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbx))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), beta_d.x_phase[bn1][0], beta_d.x_phase[bn1][1], beta_d.x_phase[bn1][2], alfax[bn1][0], alfax[bn1][1], alfax[bn1][2], mad_ac.BETX[mad_ac.indx[bn1]], mad_ac.ALFX[mad_ac.indx[bn1]], mad_ac.MUX[mad_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free beta
        if getllm_d.with_ac_calc:
            #-- from eq
            try:
                [beta_d.x_phase_f, rmsbbxf, alfaxf, bpmsf] = beta_from_phase(mad_twiss, twiss_d.zero_dpp_x, phase_d.x_f, 'H')
                tfs_file = files_dict['getbetax_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
                tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = str.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), beta_d.x_phase_f[bn1][0], beta_d.x_phase_f[bn1][1], beta_d.x_phase_f[bn1][2], alfaxf[bn1][0], alfaxf[bn1][1], alfaxf[bn1][2], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.ALFX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except:
                traceback.print_exc()
                
            #-- from the model
            [betaxf2, rmsbbxf2, alfaxf2, bpmsf2] = _get_free_beta(mad_ac, mad_twiss, beta_d.x_phase, rmsbbx, alfax, bpms, 'H')
            tfs_file = files_dict['getbetax_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
            tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf2))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "ERRBETX", "STDBETX", "ALFX", "ERRALFX", "STDALFX", "BETXMDL", "ALFXMDL", "MUXMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = str.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), betaxf2[bn1][0], betaxf2[bn1][1], betaxf2[bn1][2], alfaxf2[bn1][0], alfaxf2[bn1][1], alfaxf2[bn1][2], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.ALFX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
    #---- V plane
    if twiss_d.has_zero_dpp_y():
        [beta_d.y_phase, rmsbby, alfay, bpms] = beta_from_phase(mad_ac, twiss_d.zero_dpp_y, phase_d.ph_y, 'V')
        beta_d.y_phase['DPP'] = 0
        tfs_file = files_dict['getbetay.out']
        tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1))
        tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbby))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "ERRBETY", "STDBETY", "ALFY", "ERRALFY", "STDALFY", "BETYMDL", "ALFYMDL", "MUYMDL"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), beta_d.y_phase[bn1][0], beta_d.y_phase[bn1][1], beta_d.y_phase[bn1][2], alfay[bn1][0], alfay[bn1][1], alfay[bn1][2], mad_ac.BETY[mad_ac.indx[bn1]], mad_ac.ALFY[mad_ac.indx[bn1]], mad_ac.MUY[mad_ac.indx[bn1]]]
            tfs_file.add_table_row(list_row_entries)
        
        #-- ac to free beta
        if getllm_d.with_ac_calc:
            #-- from eq
            try:
                [beta_d.y_phase_f, rmsbbyf, alfayf, bpmsf] = beta_from_phase(mad_twiss, twiss_d.zero_dpp_y, phase_d.y_f, 'V')
                tfs_file = files_dict['getbetay_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
                tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "ERRBETY", "STDBETY", "ALFY", "ERRALFY", "STDALFY", "BETYMDL", "ALFYMDL", "MUYMDL"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = str.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), beta_d.y_phase_f[bn1][0], beta_d.y_phase_f[bn1][1], beta_d.y_phase_f[bn1][2], alfayf[bn1][0], alfayf[bn1][1], alfayf[bn1][2], mad_twiss.BETY[mad_twiss.indx[bn1]], mad_twiss.ALFY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]]]
                    tfs_file.add_table_row(list_row_entries)
            except:
                traceback.print_exc()
                
            #-- from the model
            [betayf2, rmsbbyf2, alfayf2, bpmsf2] = _get_free_beta(mad_ac, mad_twiss, beta_d.y_phase, rmsbby, alfay, bpms, 'V')
            tfs_file = files_dict['getbetay_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
            tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf2))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "ERRBETY", "STDBETY", "ALFY", "ERRALFY", "STDALFY", "BETYMDL", "ALFYMDL", "MUYMDL"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = str.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), betayf2[bn1][0], betayf2[bn1][1], betayf2[bn1][2], alfayf2[bn1][0], alfayf2[bn1][1], alfayf2[bn1][2], mad_twiss.BETY[mad_twiss.indx[bn1]], mad_twiss.ALFY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]]]
                tfs_file.add_table_row(list_row_entries)
    
    return beta_d
# END calculate_beta_from_phase -------------------------------------------------------------------------------


def calculate_beta_from_amplitude(getllm_d, twiss_d, tune_d, phase_d, beta_d, mad_twiss, mad_ac, files_dict):
    '''
    Calculates beta and fills the following TfsFiles:
        getampbetax.out        getampbetax_free.out        getampbetax_free2.out
        getampbetay.out        getampbetay_free.out        getampbetay_free2.out
        
    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            accel and beam_direction are used.
        'twiss_d': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In-param, values will only be read)
            Holds tunes and phase advances
        'phase_d': _PhaseData (In-param, values will only be read)
            Holds results from get_phases
        'beta_d': _BetaData (In/Out-param, values will be read and set)
            Holds results from get_beta. Beta from amp and ratios will be set.
        
    :Return: _BetaData
        the same instance as param beta_d to indicate that x_amp,y_amp and ratios were set.
    '''
    print 'Calculating beta from amplitude'

    #---- H plane
    if twiss_d.has_zero_dpp_x():
        [beta_d.x_amp, rmsbbx, bpms, inv_jx] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_x, 'H')
        beta_d.x_amp['DPP'] = 0
        #-- Rescaling
        beta_d.x_ratio = 0
        skipped_bpmx = []
        arcbpms = Utilities.bpm.filterbpm(bpms)
        for bpm in arcbpms:
            name = str.upper(bpm[1]) # second entry is the name
        #Skip BPM with strange data
            if abs(beta_d.x_phase[name][0] / beta_d.x_amp[name][0]) > 100:
                skipped_bpmx.append(name)
            elif (beta_d.x_amp[name][0] < 0 or beta_d.x_phase[name][0] < 0):
                skipped_bpmx.append(name)
            else:
                beta_d.x_ratio = beta_d.x_ratio + (beta_d.x_phase[name][0] / beta_d.x_amp[name][0])
        
        try:
            beta_d.x_ratio = beta_d.x_ratio / (len(arcbpms) - len(skipped_bpmx))
        except:
            traceback.print_exc()
            beta_d.x_ratio = 1
        betax_rescale = {}
        
        for bpm in bpms:
            name = str.upper(bpm[1])
            betax_rescale[name] = [beta_d.x_ratio * beta_d.x_amp[name][0], beta_d.x_ratio * beta_d.x_amp[name][1], beta_d.x_amp[name][2]]
        
        tfs_file = files_dict['getampbetax.out']
        tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1))
        tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbx))
        tfs_file.add_descriptor("RescalingFactor", "%le", str(beta_d.x_ratio))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), beta_d.x_amp[bn1][0], beta_d.x_amp[bn1][1], mad_ac.BETX[mad_ac.indx[bn1]], mad_ac.MUX[mad_ac.indx[bn1]], betax_rescale[bn1][0], betax_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries) 
            
        #-- ac to free amp beta
        if getllm_d.with_ac_calc: 
            #-- from eq
            try:
                [betaxf, rmsbbxf, bpmsf] = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, twiss_d.zero_dpp_x, tune_d.q1, tune_d.q1f, phase_d.acphasex_ac2bpmac, 'H', getllm_d.beam_direction, getllm_d.lhc_phase)[:3]
                #-- Rescaling
                beta_d.x_ratio_f = 0
                skipped_bpmxf = []
                arcbpms = Utilities.bpm.filterbpm(bpmsf)
                for bpm in arcbpms:
                    name = str.upper(bpm[1]) # second entry is the name
                #Skip BPM with strange data
                    if abs(beta_d.x_phase_f[name][0] / betaxf[name][0]) > 10:
                        skipped_bpmxf.append(name)
                    elif abs(beta_d.x_phase_f[name][0] / betaxf[name][0]) < 0.1:
                        skipped_bpmxf.append(name)
                    elif (betaxf[name][0] < 0 or beta_d.x_phase_f[name][0] < 0):
                        skipped_bpmxf.append(name)
                    else:
                        beta_d.x_ratio_f = beta_d.x_ratio_f + (beta_d.x_phase_f[name][0] / betaxf[name][0])
                
                try:
                    beta_d.x_ratio_f = beta_d.x_ratio_f / (len(arcbpms) - len(skipped_bpmxf))
                except:
                    traceback.print_exc()
                    beta_d.x_ratio_f = 1
                tfs_file = files_dict['getampbetax_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
                tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf))
                tfs_file.add_descriptor("RescalingFactor", "%le", str(beta_d.x_ratio_f))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = str.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), betaxf[bn1][0], betaxf[bn1][1], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]], beta_d.x_ratio_f * betaxf[bn1][0], beta_d.x_ratio_f * betaxf[bn1][1]]
                    tfs_file.add_table_row(list_row_entries) # Since invJxf(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            
            except:
                traceback.print_exc()
            #-- from the model
            # Since invJxf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betaxf2, rmsbbxf2, bpmsf2] = _get_free_amp_beta(beta_d.x_amp, rmsbbx, bpms, inv_jx, mad_ac, mad_twiss, 'H')[:3]
            betaxf2_rescale = _get_free_amp_beta(betax_rescale, rmsbbx, bpms, inv_jx, mad_ac, mad_twiss, 'H')[0]
            tfs_file = files_dict['getampbetax_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
            tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbxf2))
            tfs_file.add_descriptor("RescalingFactor", "%le", str(beta_d.x_ratio))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = str.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), betaxf2[bn1][0], betaxf2[bn1][1], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]], betaxf2_rescale[bn1][0], betaxf2_rescale[bn1][1]]
                tfs_file.add_table_row(list_row_entries) #---- V plane
    
    if twiss_d.has_zero_dpp_y():
        [beta_d.y_amp, rmsbby, bpms, inv_jy] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_y, 'V')
        beta_d.y_amp['DPP'] = 0
        #-- Rescaling
        beta_d.y_ratio = 0
        skipped_bpmy = []
        arcbpms = Utilities.bpm.filterbpm(bpms)
        for bpm in arcbpms:
            name = str.upper(bpm[1]) # second entry is the name
        #Skip BPM with strange data
            if abs(beta_d.y_phase[name][0] / beta_d.y_amp[name][0]) > 100:
                skipped_bpmy.append(name)
            elif (beta_d.y_amp[name][0] < 0 or beta_d.y_phase[name][0] < 0):
                skipped_bpmy.append(name)
            else:
                beta_d.y_ratio = beta_d.y_ratio + (beta_d.y_phase[name][0] / beta_d.y_amp[name][0])
        
        try:
            beta_d.y_ratio = beta_d.y_ratio / (len(arcbpms) - len(skipped_bpmy))
        except ZeroDivisionError:
            beta_d.y_ratio = 1
        betay_rescale = {}
        
        for bpm in bpms:
            name = str.upper(bpm[1])
            betay_rescale[name] = [beta_d.y_ratio * beta_d.y_amp[name][0], beta_d.y_ratio * beta_d.y_amp[name][1], beta_d.y_amp[name][2]]
        
        tfs_file = files_dict['getampbetay.out']
        tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1))
        tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2))
        tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbby))
        tfs_file.add_descriptor("RescalingFactor", "%le", str(beta_d.y_ratio))
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), beta_d.y_amp[bn1][0], beta_d.y_amp[bn1][1], mad_ac.BETY[mad_ac.indx[bn1]], mad_ac.MUY[mad_ac.indx[bn1]], betay_rescale[bn1][0], betay_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries) #-- ac to free amp beta
        
        if getllm_d.with_ac_calc: #-- from eq
            try:
                # Since invJyf(return_value[3]) is not used, slice the return value([:3]) (vimaier)
                [betayf, rmsbbyf, bpmsf] = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, twiss_d.zero_dpp_y, tune_d.q2, tune_d.q2f, phase_d.acphasey_ac2bpmac, 'V', getllm_d.beam_direction, getllm_d.accel)[:3] #-- Rescaling
                beta_d.y_ratio_f = 0
                skipped_bpmyf = []
                arcbpms = Utilities.bpm.filterbpm(bpmsf)
                for bpm in arcbpms:
                    name = str.upper(bpm[1]) # second entry is the name
                    #Skip BPM with strange data
                    if abs(beta_d.y_phase_f[name][0] / betayf[name][0]) > 10:
                        skipped_bpmyf.append(name)
                    elif (betayf[name][0] < 0 or beta_d.y_phase_f[name][0] < 0):
                        skipped_bpmyf.append(name)
                    elif abs(beta_d.y_phase_f[name][0] / betayf[name][0]) < 0.1:
                        skipped_bpmyf.append(name)
                    else:
                        beta_d.y_ratio_f = beta_d.y_ratio_f + (beta_d.y_phase_f[name][0] / betayf[name][0])
                
                try:
                    beta_d.y_ratio_f = beta_d.y_ratio_f / (len(arcbpms) - len(skipped_bpmyf))
                except:
                    traceback.print_exc()
                    beta_d.y_ratio_f = 1
                tfs_file = files_dict['getampbetay_free.out']
                tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
                tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
                tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf))
                tfs_file.add_descriptor("RescalingFactor", "%le", str(beta_d.y_ratio_f))
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = str.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), betayf[bn1][0], betayf[bn1][1], mad_twiss.BETY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]], (beta_d.y_ratio_f * betayf[bn1][0]), (beta_d.y_ratio_f * betayf[bn1][1])]
                    tfs_file.add_table_row(list_row_entries) # 'except ALL' catched a SystemExit from filterbpm().(vimaier)
            
            except SystemExit:
                sys.exit(1)
            except:
                #-- from the model
                traceback.print_exc()
            # Since invJyf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betayf2, rmsbbyf2, bpmsf2] = _get_free_amp_beta(beta_d.y_amp, rmsbby, bpms, inv_jy, mad_ac, mad_twiss, 'V')[:3]
            betayf2_rescale = _get_free_amp_beta(betay_rescale, rmsbby, bpms, inv_jy, mad_ac, mad_twiss, 'V')[0]
            tfs_file = files_dict['getampbetay_free2.out']
            tfs_file.add_descriptor("Q1", "%le", str(tune_d.q1f))
            tfs_file.add_descriptor("Q2", "%le", str(tune_d.q2f))
            tfs_file.add_descriptor("RMSbetabeat", "%le", str(rmsbbyf2))
            tfs_file.add_descriptor("RescalingFactor", "%le", str(beta_d.y_ratio))
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = str.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), betayf2[bn1][0], betayf2[bn1][1], mad_twiss.BETY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]], betayf2_rescale[bn1][0], betayf2_rescale[bn1][1]]
                tfs_file.add_table_row(list_row_entries)
    
    return beta_d
# END calculate_beta_from_amplitude ----------------------------------------------------------------


#===================================================================================================
# helper-functions
#===================================================================================================

def beta_from_phase(MADTwiss,ListOfFiles,phase,plane):

    alfa={}
    beta={}
    commonbpms = Utilities.bpm.intersect(ListOfFiles)
    commonbpms = Utilities.bpm.model_intersect(commonbpms,MADTwiss)
    alfii=[]
    betii=[]
    delbeta=[]
    for i in range(0,len(commonbpms)):
        if i==len(commonbpms)-2: # The last monitor but one
            bn1=str.upper(commonbpms[i][1])
            bn2=str.upper(commonbpms[i+1][1])
            bn3=str.upper(commonbpms[0][1])
        elif i==len(commonbpms)-1: # The last monitor
            bn1=str.upper(commonbpms[i][1])
            bn2=str.upper(commonbpms[0][1])
            bn3=str.upper(commonbpms[1][1])
        else : # Others
            bn1=str.upper(commonbpms[i][1])
            bn2=str.upper(commonbpms[i+1][1])
            bn3=str.upper(commonbpms[i+2][1])
        ph2pi12=2.*np.pi*phase[bn1][0]
        ph2pi23=2.*np.pi*phase[bn2][0]
        ph2pi13=2.*np.pi*phase[bn1][2]
        # Find the model transfer matrices for beta1
        phmdl12=2.*np.pi*phase[bn1][4]
        phmdl13=2.*np.pi*phase[bn1][5]
        phmdl23=2.*np.pi*phase[bn2][4]
        if plane=='H':
            betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
            betmdl2=MADTwiss.BETX[MADTwiss.indx[bn2]]
            betmdl3=MADTwiss.BETX[MADTwiss.indx[bn3]]
            alpmdl1=MADTwiss.ALFX[MADTwiss.indx[bn1]]
            alpmdl2=MADTwiss.ALFX[MADTwiss.indx[bn2]]
            alpmdl3=MADTwiss.ALFX[MADTwiss.indx[bn3]]
        elif plane=='V':
            betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
            betmdl2=MADTwiss.BETY[MADTwiss.indx[bn2]]
            betmdl3=MADTwiss.BETY[MADTwiss.indx[bn3]]
            alpmdl1=MADTwiss.ALFY[MADTwiss.indx[bn1]]
            alpmdl2=MADTwiss.ALFY[MADTwiss.indx[bn2]]
            alpmdl3=MADTwiss.ALFY[MADTwiss.indx[bn3]]
        if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
            print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
            sys.exit(1)

        # Find beta1 and alpha1 from phases assuming model transfer matrix
        # Matrix M: BPM1-> BPM2
        # Matrix N: BPM1-> BPM3
        M11=math.sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
        M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        N11=math.sqrt(betmdl3/betmdl1)*(cos(phmdl13)+alpmdl1*sin(phmdl13))
        N12=math.sqrt(betmdl1*betmdl3)*sin(phmdl13)

        denom=M11/M12-N11/N12+1e-16
        numer=1/tan(ph2pi12)-1/tan(ph2pi13)
        beti1=numer/denom
        # beterr1=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # beterr1=beterr1+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # beterr1=beterr1/abs(denom)
        betstd1=        (2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        betstd1=betstd1+(2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        betstd1=math.sqrt(betstd1)/abs(denom)

        denom=M12/M11-N12/N11+1e-16
        numer=-M12/M11/tan(ph2pi12)+N12/N11/tan(ph2pi13)
        alfi1=numer/denom
        # alferr1=abs(M12/M11*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # alferr1=alferr1+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # alferr1=alferr1/abs(denom)
        alfstd1=        (M12/M11*2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        alfstd1=alfstd1+(N12/N11*2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        alfstd1=math.sqrt(alfstd1)/denom

        # Find beta2 and alpha2 from phases assuming model transfer matrix
        # Matrix M: BPM1-> BPM2
        # Matrix N: BPM2-> BPM3
        M22=math.sqrt(betmdl1/betmdl2)*(cos(phmdl12)-alpmdl2*sin(phmdl12))
        M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
        N11=math.sqrt(betmdl3/betmdl2)*(cos(phmdl23)+alpmdl2*sin(phmdl23))
        N12=math.sqrt(betmdl2*betmdl3)*sin(phmdl23)

        denom=M22/M12+N11/N12+1e-16
        numer=1/tan(ph2pi12)+1/tan(ph2pi23)
        beti2=numer/denom
        # beterr2=abs(phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # beterr2=beterr2+abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # beterr2=beterr2/abs(denom)
        betstd2=        (2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        betstd2=betstd2+(2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        betstd2=math.sqrt(betstd2)/abs(denom)

        denom=M12/M22+N12/N11+1e-16
        numer=M12/M22/tan(ph2pi12)-N12/N11/tan(ph2pi23)
        alfi2=numer/denom
        # alferr2=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi12)**2.)/tan(ph2pi12)**2.)
        # alferr2=alferr2+abs(N12/N11*phase[bn1][3]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # alferr2=alferr2/abs(denom)
        alfstd2=        (M12/M22*2*np.pi*phase[bn1][1]/sin(ph2pi12)**2)**2
        alfstd2=alfstd2+(N12/N11*2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        alfstd2=math.sqrt(alfstd2)/abs(denom)

        # Find beta3 and alpha3 from phases assuming model transfer matrix
        # Matrix M: BPM2-> BPM3
        # Matrix N: BPM1-> BPM3
        M22=math.sqrt(betmdl2/betmdl3)*(cos(phmdl23)-alpmdl3*sin(phmdl23))
        M12=math.sqrt(betmdl2*betmdl3)*sin(phmdl23)
        N22=math.sqrt(betmdl1/betmdl3)*(cos(phmdl13)-alpmdl3*sin(phmdl13))
        N12=math.sqrt(betmdl1*betmdl3)*sin(phmdl13)

        denom=M22/M12-N22/N12+1e-16
        numer=1/tan(ph2pi23)-1/tan(ph2pi13)
        beti3=numer/denom
        # beterr3=abs(phase[bn2][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # beterr3=beterr3+abs(phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # beterr3=beterr3/abs(denom)
        betstd3=        (2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        betstd3=betstd3+(2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        betstd3=math.sqrt(betstd3)/abs(denom)

        denom=M12/M22-N12/N22+1e-16
        numer=M12/M22/tan(ph2pi23)-N12/N22/tan(ph2pi13)
        alfi3=numer/denom
        # alferr3=abs(M12/M22*phase[bn1][1]*(1.-tan(ph2pi23)**2.)/tan(ph2pi23)**2.)
        # alferr3=alferr3+abs(N22/N12*phase[bn1][3]*(1.-tan(ph2pi13)**2.)/tan(ph2pi13)**2.)
        # alferr3=alferr3/abs(denom)
        alfstd3=        (M12/M22*2*np.pi*phase[bn2][1]/sin(ph2pi23)**2)**2
        alfstd3=alfstd3+(N12/N22*2*np.pi*phase[bn1][3]/sin(ph2pi13)**2)**2
        alfstd3=math.sqrt(alfstd3)/abs(denom)

        betii.append([beti1,betstd1,beti2,betstd2,beti3,betstd3])
        alfii.append([alfi1,alfstd1,alfi2,alfstd2,alfi3,alfstd3])

    # Output the beta at all monitors even for the first, second, last and last but one using beta1,2 and 3!
    for i in range(0,len(commonbpms)):
        bn1=str.upper(commonbpms[i][1])
        if i==0: # The first monitor
            ib1=0
            ib2=len(commonbpms)-1
            ib3=len(commonbpms)-2

        elif i==1: # The second monitor
            ib1=i
            ib2=0
            ib3=len(commonbpms)-1
        else: # Others
            ib1=i
            ib2=i-1
            ib3=i-2
        # Averaging over beta/alpha1,2 and 3. Find beta/alpha errors
        beti=(betii[ib1][0]+betii[ib2][2]+betii[ib3][4])/3.
        betstd=math.sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/math.sqrt(3.)
        # betstd=math.sqrt(betii[ib1][1]**2.+betii[ib2][3]**2.+betii[ib3][5]**2.)/math.sqrt(3.*len(ListOfFiles))  # If we want to plot std "over files" !!!!!
        try:
            beterr=math.sqrt((betii[ib1][0]**2.+betii[ib2][2]**2.+betii[ib3][4]**2.)/3.-beti**2.)
        except:
            beterr=0
        alfi=(alfii[ib1][0]+alfii[ib2][2]+alfii[ib3][4])/3.
        alfstd=math.sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/math.sqrt(3.)
        # alfstd=math.sqrt(alfii[ib1][1]**2.+alfii[ib2][3]**2.+alfii[ib3][5]**2.)/math.sqrt(3.*len(ListOfFiles))  # If we want to plot std "over files" !!!!!
        try:
            alferr=math.sqrt((alfii[ib1][0]**2.+alfii[ib2][2]**2.+alfii[ib3][4]**2.)/3.-alfi**2.)
        except:
            alferr=0

        beta[bn1]=(beti,beterr,betstd)
        alfa[bn1]=(alfi,alferr,alfstd)
        if plane=='H':
            betmdl1=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            betmdl1=MADTwiss.BETY[MADTwiss.indx[bn1]]
        delbeta.append((beti-betmdl1)/betmdl1)


    delbeta=np.array(delbeta)
    rmsbb=math.sqrt(np.average(delbeta*delbeta))

    return [beta,rmsbb,alfa,commonbpms]


def beta_from_amplitude(MADTwiss,ListOfFiles,plane):

    beta={}
    root2J=[]
    commonbpms=Utilities.bpm.intersect(ListOfFiles)
    commonbpms=Utilities.bpm.model_intersect(commonbpms,MADTwiss)
    SumA=0.0
    Amp=[]
    Amp2=[]
    Kick2=[]
    for i in range(0,len(commonbpms)): # this loop have become complicated after modifications... anybody simplify?
        bn1=str.upper(commonbpms[i][1])
        if plane=='H':
            tembeta=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            tembeta=MADTwiss.BETY[MADTwiss.indx[bn1]]
        Ampi=0.0
        Ampj2=[]
        root2Ji=0.0
        jj=0
        for j in ListOfFiles:
            if i==0:
                Kick2.append(0)
            if plane=='H':
                Ampi+=j.AMPX[j.indx[bn1]]
                Ampj2.append(j.AMPX[j.indx[bn1]]**2)
                root2Ji+=j.PK2PK[j.indx[bn1]]/2.

            elif plane=='V':
                Ampi+=j.AMPY[j.indx[bn1]]
                Ampj2.append(j.AMPY[j.indx[bn1]]**2)
                root2Ji+=j.PK2PK[j.indx[bn1]]/2.
            Kick2[jj]+=Ampj2[jj]/tembeta
            jj+=1
        Ampi=Ampi/len(ListOfFiles)
        root2Ji=root2Ji/len(ListOfFiles)
        Amp.append(Ampi)
        Amp2.append(Ampj2)



        SumA+=Ampi**2/tembeta
        root2J.append(root2Ji/math.sqrt(tembeta))


    Kick=SumA/len(commonbpms) # Assuming the average of beta is constant
    Kick2=np.array(Kick2)
    Kick2=Kick2/len(commonbpms)
    Amp2=np.array(Amp2)
    root2J=np.array(root2J)
    root2Jave=np.average(root2J)
    root2Jrms=math.sqrt(np.average(root2J*root2J)-root2Jave**2+2.2e-16)

    delbeta=[]
    for i in range(0,len(commonbpms)):
        bn1=str.upper(commonbpms[i][1])
        location=commonbpms[i][0]
        for j in range(0,len(ListOfFiles)):
            Amp2[i][j]=Amp2[i][j]/Kick2[j]
        #print np.average(Amp2[i]*Amp2[i]),np.average(Amp2[i])**2
        try:
            betstd=math.sqrt(np.average(Amp2[i]*Amp2[i])-np.average(Amp2[i])**2+2.2e-16)
        except:
            betstd=0

        beta[bn1]=[Amp[i]**2/Kick,betstd,location]
        if plane=='H':
            betmdl=MADTwiss.BETX[MADTwiss.indx[bn1]]
        elif plane=='V':
            betmdl=MADTwiss.BETY[MADTwiss.indx[bn1]]
        delbeta.append((beta[bn1][0]-betmdl)/betmdl)

    invariantJ=[root2Jave,root2Jrms]

    delbeta=np.array(delbeta)
    rmsbb=math.sqrt(np.average(delbeta*delbeta))
    return [beta,rmsbb,commonbpms,invariantJ]

#===================================================================================================
# ac-dipole stuff
#===================================================================================================

def _get_free_beta(modelfree,modelac,betal,rmsbb,alfal,bpms,plane): # to check "+"
    if DEBUG:
        print "Calculating free beta using model"
    bpms=Utilities.bpm.model_intersect(bpms, modelfree)
    bpms=Utilities.bpm.model_intersect(bpms, modelac)
    betan={}
    alfan={}
    for bpma in bpms:

        bpm=bpma[1].upper()
        beta,beterr,betstd=betal[bpm]
        alfa,alferr,alfstd=alfal[bpm]

        if plane=="H":
            betmf=modelfree.BETX[modelfree.indx[bpm]]
            betma=modelac.BETX[modelac.indx[bpm]]
            bb=(betma-betmf)/betmf
            alfmf=modelfree.ALFX[modelfree.indx[bpm]]
            alfma=modelac.ALFX[modelac.indx[bpm]]
            aa=(alfma-alfmf)/alfmf
        else:
            betmf=modelfree.BETY[modelfree.indx[bpm]]
            betma=modelac.BETY[modelac.indx[bpm]]
            alfmf=modelfree.ALFY[modelfree.indx[bpm]]
            alfma=modelac.ALFY[modelac.indx[bpm]]
            bb=(betma-betmf)/betmf
            aa=(alfma-alfmf)/alfmf

        betan[bpm]=beta*(1+bb),beterr,betstd # has to be plus!
        alfan[bpm]=alfa*(1+aa),alferr,alfstd

    return betan,rmsbb,alfan,bpms


def _get_free_amp_beta(betai,rmsbb,bpms,invJ,modelac,modelfree,plane): # "-"

    #
    # Why difference in betabeta calculation ??
    #
    #

    betas={}

    if DEBUG:
        print "Calculating free beta from amplitude using model"

    for bpm in bpms:

        bpmm=bpm[1].upper()
        beta=betai[bpmm][0]

        if plane=="H":
            betmf=modelfree.BETX[modelfree.indx[bpmm]]
            betma=modelac.BETX[modelac.indx[bpmm]]
            bb=(betmf-betma)/betma

        else:
            betmf=modelfree.BETY[modelfree.indx[bpmm]]
            betma=modelac.BETY[modelac.indx[bpmm]]
            bb=(betmf-betma)/betma

        betas[bpmm]=[beta*(1+bb),betai[bpmm][1],betai[bpmm][2]]

    return betas,rmsbb,bpms,invJ

