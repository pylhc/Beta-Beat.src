'''
Created on 27 May 2013

@author: awegsche, vimaier

@version: 2016.11.p2

GetLLM.algorithms.beta.py stores helper functions for phase calculations for GetLLM.
This module is not intended to be executed. It stores only functions.


'''
import sys
import math
import traceback

import numpy as np
from numpy import sin, cos, tan

import Python_Classes4MAD.metaclass
import tfs_files
import utils.bpm
import compensate_ac_effect
import os
import re
import multiprocessing
import time
from constants import PI, TWOPI

__version__ = "2017.3.2"

DEBUG = sys.flags.debug  # True with python option -d! ("python -d GetLLM.py...") (vimaier)
PRINTTIMES = False

if False:
    from utils.progressbar import startProgress, progress, endProgress
else:
    def startProgress(name):
        print_("START " + name)
        
    def progress(i):
        return
    
    def endProgress():
        return

#--- Constants



DEFAULT_WRONG_BETA      = 1000                      #@IgnorePep8
EPSILON                 = 0#1.0E-16                 #@IgnorePep8
SEXT_FACT               = 2.0                       #@IgnorePep8
A_FACT                  = -.5                       #@IgnorePep8
BADPHASE                = .5
BETA_THRESHOLD          = 1e3                       #@IgnorePep8
ZERO_THRESHOLD          = 1e-2                      #@IgnorePep8
PHASE_THRESHOLD         = 1e-2                      #@IgnorePep8
MOD_POINTFIVE_LOWER     = PHASE_THRESHOLD           #@IgnorePep8
MOD_POINTFIVE_UPPER     = (BADPHASE - PHASE_THRESHOLD)    #@IgnorePep8
RCOND                   = 1.0e-14                    #@IgnorePep8

BOXLENGTH               = 50                        #@IgnorePep8
BOXINDENT               =  4                        #@IgnorePep8

#=======================================================================================================================
#--- classes
#=======================================================================================================================


class MeasuredValues:
    def __init__(self, alfa, beta, string="", use_it=False):
        self.beta = beta
        self.alfa = alfa
        self.patternstring = string
        self.use_it = use_it
        
        
class BetaData(object):
    """ File for storing results from beta computations. """

    def __init__(self):
        self.x_phase = None  # beta x from phase
        self.x_phase_f = None  # beta x from phase free
        self.y_phase = None  # beta y from phase
        self.y_phase_f = None  # beta y from phase free

        self.x_amp = None  # beta x from amplitude
        self.y_amp = None  # beta y from amplitude

        self.x_ratio = 0  # beta x ratio
        self.x_ratio_f = 0  # beta x ratio free
        self.y_ratio = 0  # beta x ratio
        self.y_ratio_f = 0  # beta x ratio free

#===================================================================================================
# main part
#===================================================================================================


def _write_getbeta_out(twiss_d_zero_dpp, q1, q2, mad_ac, number_of_bpms, range_of_bpms, beta_d_col,
                       data, rmsbbx, error_method, bpms,
                       tfs_file, mod_BET, mod_ALF, mod_MU, _plane_char,
                       dpp=0, dppq1=0):

    tfs_file.add_float_descriptor("Q1", q1)
    tfs_file.add_float_descriptor("Q2", q2)
    tfs_file.add_float_descriptor("RMSbetabeat", rmsbbx)
    tfs_file.add_float_descriptor("DPP", dpp)
        
    tfs_file.add_string_descriptor("BetaAlgorithmVersion", __version__)
    if error_method == "Estimated by std (3 BPM method), no bet_deviations.npy file found":
        tfs_file.add_float_descriptor("NumberOfBPMs", 3)
        tfs_file.add_float_descriptor("RangeOfBPMs", 5)
    else:
        tfs_file.add_float_descriptor("NumberOfBPMs", number_of_bpms)
        tfs_file.add_float_descriptor("RangeOfBPMs", range_of_bpms)
    tfs_file.add_string_descriptor("ErrorsFrom", error_method)
    tfs_file.add_float_descriptor("PhaseTheshold", PHASE_THRESHOLD)
    tfs_file.add_float_descriptor("RCond", RCOND)
    tfs_file.add_column_names(["NAME", "S", "COUNT",
                               "BET" + _plane_char, "SYSBET" + _plane_char, "STATBET" + _plane_char, "ERRBET" + _plane_char,
                               "CORR_ALFABETA",
                               "ALF" + _plane_char, "SYSALF" + _plane_char, "STATALF" + _plane_char, "ERRALF" + _plane_char,
                               "BET" + _plane_char + "MDL", "ALF" + _plane_char + "MDL", "MU" + _plane_char + "MDL",
                               "NCOMBINATIONS"])
    tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
    for bpm in bpms:
        name = bpm[1]
        try:  # TODO: An error was happening here, this is a fast fix, remove!!
            row = data[name]
        except KeyError:
            continue
        beta_d_col[name] = [row[0], row[1], row[2], row[3]]
        model_ac_index = mad_ac.indx[name]
        list_row_entries = ['"' + name + '"', bpm[0], len(twiss_d_zero_dpp),
                            row[0], row[1], row[2], row[3],
                            row[8],
                            row[4], row[5], row[6], row[7],
                            mod_BET[model_ac_index], mod_ALF[model_ac_index], mod_MU[model_ac_index],
                            row[10]]
        # list_row_entries = ['"' + name + '"', mad_ac.S[model_ac_index], len(twiss_d_zero_dpp),
        #                     row[0], row[1], row[2], row[3],
        #                     row[8],
        #                     row[4], row[5], row[6], row[7],
        #                     mod_BET[model_ac_index], mod_ALF[model_ac_index], mod_MU[model_ac_index],
        #                     row[10]]
        tfs_file.add_table_row(list_row_entries)


def calculate_beta_from_phase(getllm_d, twiss_d, tune_d, phase_d,
                              mad_twiss, mad_ac, mad_elem, mad_elem_centre, mad_best_knowledge, mad_ac_best_knowledge,
                              files_dict):
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
    
    debugfile = None
    if getllm_d.nprocesses == -1:
        getllm_d.nprocesses = multiprocessing.cpu_count()
    getllm_d.parallel = (getllm_d.nprocesses > 0)
    #---- H plane
    _plane_char = "X"
    print "\n"
    print_box_edge()
    print_box("Calculating beta from phase")
    print_box("Version: {0:5s}".format(__version__))
    
    if DEBUG:
        debugfile = open(files_dict['getbetax.out'].s_output_path + "/getbetax.debug", "w+")
        print_box_edge()
        print_("ATTENTION: DEBUG is set to true, calculation of beta functions will be done serially")
    elif getllm_d.parallel:
        print_box("parallel : TRUE")
        print_box("number of processes : {0:2d}".format(getllm_d.nprocesses))
        print_box_edge()
    else:
        print_box("parallel : FALSE")
        print_box_edge()
        
    starttime = time.time()

    if twiss_d.has_zero_dpp_x():
        print ""
        print_("Calculate beta from phase for plane " + _plane_char, ">")
        data, rmsbbx, bpms, error_method = beta_from_phase(mad_ac_best_knowledge, mad_elem, mad_elem_centre,
                                                           twiss_d.zero_dpp_x, phase_d.ph_x, 'H',
                                                           getllm_d, debugfile)
        beta_d.x_phase = {}
        beta_d.x_phase['DPP'] = 0
        tfs_file = files_dict['getbetax.out']
        
        _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1, tune_d.q2, mad_ac, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.x_phase,
                           data, rmsbbx, error_method, bpms,
                           tfs_file, mad_ac.BETX, mad_ac.ALFX, mad_ac.MUX, _plane_char)

        #-- ac to free beta
        if getllm_d.with_ac_calc:
            #-- from eq
            print ""
            print_("Calculate beta from phase for plane " + _plane_char + " with AC dipole (_free.out)", ">")
            try:
                if DEBUG:
                    debugfile = open(files_dict['getbetax_free.out'].s_output_path + "/getbetax_free.debug", "w+")

                dataf, rmsbbxf, bpmsf, error_method = beta_from_phase(mad_best_knowledge, mad_elem, mad_elem_centre,
                                                                      twiss_d.zero_dpp_x, phase_d.x_f, 'H',
                                                                      getllm_d, debugfile)
                tfs_file = files_dict['getbetax_free.out']
                beta_d.x_phase_f = {}
                _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, mad_twiss, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.x_phase_f,
                                   dataf, rmsbbxf, error_method, bpmsf,
                                   tfs_file, mad_twiss.BETX, mad_twiss.ALFX, mad_twiss.MUX, _plane_char)
                
            except:
                traceback.print_exc()

            #-- from the model
            print ""
            print_("Calculate beta from phase for plane " + _plane_char + " with AC dipole (_free2.out)", ">")
            dataf2, bpmsf2 = _get_free_beta(mad_twiss, mad_ac, data, bpms, 'H')
            tfs_file = files_dict['getbetax_free2.out']
            betaxf2 = {}
            rmsbbxf2 = 0
            
            _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, mad_twiss, getllm_d.number_of_bpms, getllm_d.range_of_bpms, betaxf2,
                               dataf2, rmsbbxf2, error_method, bpmsf2,
                               tfs_file, mad_twiss.BETX, mad_twiss.ALFX, mad_twiss.MUX, _plane_char)
            
    #---- V plane
    _plane_char = "Y"
    
    if twiss_d.has_zero_dpp_y():
        print ""
        print_("Calculate beta from phase for plane " + _plane_char, ">")

        if DEBUG:
            debugfile = open(files_dict['getbetay.out'].s_output_path + "/getbetay.debug", "w+")

        data, rmsbby, bpms, error_method = beta_from_phase(mad_ac_best_knowledge, mad_elem, mad_elem_centre,
                                                           twiss_d.zero_dpp_y, phase_d.ph_y, 'V',
                                                           getllm_d, debugfile)
        beta_d.y_phase = {}
        beta_d.y_phase['DPP'] = 0
        tfs_file = files_dict['getbetay.out']
        
        _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1, tune_d.q2, mad_ac, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.y_phase,
                           data, rmsbby, error_method, bpms,
                           tfs_file, mad_ac.BETY, mad_ac.ALFY, mad_ac.MUY, _plane_char)
        
        #-- ac to free beta
        if getllm_d.with_ac_calc:
            #-- from eq
            print ""
            print_("Calculate beta from phase for plane " + _plane_char + " with AC dipole (_free.out)", ">")

            if DEBUG:
                debugfile = open(files_dict['getbetay_free.out'].s_output_path + "/getbetay_free.debug", "w+")

            try:
                dataf, rmsbbyf, bpmsf, error_method = beta_from_phase(mad_best_knowledge, mad_elem, mad_elem_centre,
                                                                      twiss_d.zero_dpp_y, phase_d.y_f, 'V',
                                                                      getllm_d, debugfile)
                
                tfs_file = files_dict['getbetay_free.out']
                beta_d.y_phase_f = {}
                _write_getbeta_out(twiss_d.zero_dpp_y, tune_d.q1f, tune_d.q2f, mad_twiss, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.y_phase_f,
                                   dataf, rmsbbyf, error_method, bpmsf,
                                   tfs_file, mad_twiss.BETY, mad_twiss.ALFY, mad_twiss.MUY, _plane_char)
            except:
                traceback.print_exc()

            #-- from the model
            print ""
            print_("Calculate beta from phase for plane " + _plane_char + " with AC dipole (_free2.out)", ">")

            [datayf2, bpmsf2] = _get_free_beta(mad_twiss, mad_ac, data, bpms, 'V')
            tfs_file = files_dict['getbetay_free2.out']
            
            betayf2 = {}
            rmsbbyf2 = 0
            
            _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, mad_twiss, getllm_d.number_of_bpms, getllm_d.range_of_bpms, betayf2,
                               datayf2, rmsbbyf2, error_method, bpmsf2,
                               tfs_file, mad_twiss.BETY, mad_twiss.ALFY, mad_twiss.MUY, _plane_char)
            
    elapsed = time.time() - starttime
    
    print_box_edge()
    print_box("beta from phase finished")
    print_box("")
    print_box("elapsed time: {0:3.3f}s".format(elapsed))
    print_box_edge()
    if PRINTTIMES:
        timesfile = open("times.dat", "a")
        if getllm_d.parallel:
            timesfile.write("PARALLEL {0:f} {1:d} {2:d}\n".format(elapsed, 0, getllm_d.nprocesses))
        else:
            timesfile.write("SERIAL {0:f} {1:d} {2:d}\n".format(elapsed, 0, getllm_d.nprocesses))

    print "\n"
    
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
        arcbpms = utils.bpm.filterbpm(bpms)
        if len(arcbpms) == 0:
           arcbpms = bpms
        for bpm in arcbpms:
        #for bpm in bpms:
            name = str.upper(bpm[1])  # second entry is the name
        #Skip BPM with strange data
            if abs(beta_d.x_phase[name][0] / beta_d.x_amp[name][0]) > 100:
                skipped_bpmx.append(name)
            elif (beta_d.x_amp[name][0] < 0 or beta_d.x_phase[name][0] < 0):
                skipped_bpmx.append(name)
            else:
                beta_d.x_ratio = beta_d.x_ratio + (beta_d.x_phase[name][0] / beta_d.x_amp[name][0])

        try: 
            beta_d.x_ratio = beta_d.x_ratio / (len(arcbpms) - len(skipped_bpmx))
        except ZeroDivisionError:
            beta_d.x_ratio = 1
        except:
            traceback.print_exc()
            beta_d.x_ratio = 1

        betax_rescale = {}

        for bpm in bpms:
            name = str.upper(bpm[1])
            betax_rescale[name] = [beta_d.x_ratio * beta_d.x_amp[name][0], beta_d.x_ratio * beta_d.x_amp[name][1], beta_d.x_amp[name][2]]

        tfs_file = files_dict['getampbetax.out']
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_float_descriptor("RMSbetabeat", rmsbbx)
        tfs_file.add_float_descriptor("RescalingFactor", beta_d.x_ratio)
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
                betaxf, rmsbbxf, bpmsf = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, 
                                                                                        twiss_d.zero_dpp_x, 
                                                                                        tune_d.q1, 
                                                                                        tune_d.q1f, 
                                                                                        phase_d.acphasex_ac2bpmac, 
                                                                                        'H', 
                                                                                        getllm_d.beam_direction, 
                                                                                        getllm_d.lhc_phase) # AnIssue 
                #-- Rescaling
                beta_d.x_ratio_f = 0
                skipped_bpmxf = []
                arcbpms = utils.bpm.filterbpm(bpmsf)
                if len(arcbpms) == 0:
                   arcbpms = bpmsf
                #arcbpms = utils.bpm.filterbpm(bpmsf)
                for bpm in arcbpms:
                    name = str.upper(bpm[1])  # second entry is the name
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
                tfs_file.add_float_descriptor("Q1", tune_d.q1f)
                tfs_file.add_float_descriptor("Q2", tune_d.q2f)
                tfs_file.add_float_descriptor("RMSbetabeat", rmsbbxf)
                tfs_file.add_float_descriptor("RescalingFactor", beta_d.x_ratio_f)
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = str.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), betaxf[bn1][0], betaxf[bn1][1], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]], beta_d.x_ratio_f * betaxf[bn1][0], beta_d.x_ratio_f * betaxf[bn1][1]]
                    tfs_file.add_table_row(list_row_entries)

            except:
                traceback.print_exc()
            #-- from the model
            # Since invJxf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betaxf2, rmsbbxf2, bpmsf2] = _get_free_amp_beta(beta_d.x_amp, rmsbbx, bpms, inv_jx, mad_ac, mad_twiss, 'H')[:3]
            betaxf2_rescale = _get_free_amp_beta(betax_rescale, rmsbbx, bpms, inv_jx, mad_ac, mad_twiss, 'H')[0]
            tfs_file = files_dict['getampbetax_free2.out']
            tfs_file.add_float_descriptor("Q1", tune_d.q1f)
            tfs_file.add_float_descriptor("Q2", tune_d.q2f)
            tfs_file.add_float_descriptor("RMSbetabeat", rmsbbxf2)
            tfs_file.add_float_descriptor("RescalingFactor", beta_d.x_ratio)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = str.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), betaxf2[bn1][0], betaxf2[bn1][1], mad_twiss.BETX[mad_twiss.indx[bn1]], mad_twiss.MUX[mad_twiss.indx[bn1]], betaxf2_rescale[bn1][0], betaxf2_rescale[bn1][1]]
                tfs_file.add_table_row(list_row_entries)  # V plane

    if twiss_d.has_zero_dpp_y():
        [beta_d.y_amp, rmsbby, bpms, inv_jy] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_y, 'V')
        beta_d.y_amp['DPP'] = 0
        #-- Rescaling
        beta_d.y_ratio = 0
        skipped_bpmy = []
        arcbpms = utils.bpm.filterbpm(bpms)
        if len(arcbpms) == 0:
          arcbpms = bpms
        #arcbpms = utils.bpm.filterbpm(bpms)
        for bpm in arcbpms:
            name = str.upper(bpm[1])  # second entry is the name
            #Skip BPM with strange data
            if name in beta_d.y_phase:
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
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_float_descriptor("RMSbetabeat", rmsbby)
        tfs_file.add_float_descriptor("RescalingFactor", beta_d.y_ratio)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for i in range(0, len(bpms)):
            bn1 = str.upper(bpms[i][1])
            bns1 = bpms[i][0]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), beta_d.y_amp[bn1][0], beta_d.y_amp[bn1][1], mad_ac.BETY[mad_ac.indx[bn1]], mad_ac.MUY[mad_ac.indx[bn1]], betay_rescale[bn1][0], betay_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries)  # ac to free amp beta

        if getllm_d.with_ac_calc:  # from eq
            try:
                betayf, rmsbbyf, bpmsf = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, 
                                                                                        twiss_d.zero_dpp_y, 
                                                                                        tune_d.q2, 
                                                                                        tune_d.q2f, 
                                                                                        phase_d.acphasey_ac2bpmac, 
                                                                                        'V', 
                                                                                        getllm_d.beam_direction, 
                                                                                        getllm_d.lhc_phase)  # AnIssue
                beta_d.y_ratio_f = 0
                skipped_bpmyf = []
                arcbpms = utils.bpm.filterbpm(bpmsf)
                if len(arcbpms) == 0:
                   arcbpms = bpms
                #arcbpms = utils.bpm.filterbpm(bpmsf)
                for bpm in arcbpms:
                    name = str.upper(bpm[1])  # second entry is the name
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
                except ZeroDivisionError:
                    beta_d.y_ratio_f = 1
                tfs_file = files_dict['getampbetay_free.out']
                tfs_file.add_float_descriptor("Q1", tune_d.q1f)
                tfs_file.add_float_descriptor("Q2", tune_d.q2f)
                tfs_file.add_float_descriptor("RMSbetabeat", rmsbbyf)
                tfs_file.add_float_descriptor("RescalingFactor", beta_d.y_ratio_f)
                tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
                tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
                for i in range(0, len(bpmsf)):
                    bn1 = str.upper(bpmsf[i][1])
                    bns1 = bpmsf[i][0]
                    list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), betayf[bn1][0], betayf[bn1][1], mad_twiss.BETY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]], (beta_d.y_ratio_f * betayf[bn1][0]), (beta_d.y_ratio_f * betayf[bn1][1])]
                    tfs_file.add_table_row(list_row_entries)  # 'except ALL' catched a SystemExit from filterbpm().(vimaier)

            except SystemExit:
                traceback.print_exc()
                sys.exit(1)
            except:
                #-- from the model
                traceback.print_exc()
            # Since invJyf2(return_value[3]) is not used, slice the return value([:3]) (vimaier)
            [betayf2, rmsbbyf2, bpmsf2] = _get_free_amp_beta(beta_d.y_amp, rmsbby, bpms, inv_jy, mad_ac, mad_twiss, 'V')[:3]
            betayf2_rescale = _get_free_amp_beta(betay_rescale, rmsbby, bpms, inv_jy, mad_ac, mad_twiss, 'V')[0]
            tfs_file = files_dict['getampbetay_free2.out']
            tfs_file.add_float_descriptor("Q1", tune_d.q1f)
            tfs_file.add_float_descriptor("Q2", tune_d.q2f)
            tfs_file.add_float_descriptor("RMSbetabeat", rmsbbyf2)
            tfs_file.add_float_descriptor("RescalingFactor", beta_d.y_ratio)
            tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
            tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
            for i in range(0, len(bpmsf2)):
                bn1 = str.upper(bpmsf2[i][1])
                bns1 = bpmsf2[i][0]
                list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), betayf2[bn1][0], betayf2[bn1][1], mad_twiss.BETY[mad_twiss.indx[bn1]], mad_twiss.MUY[mad_twiss.indx[bn1]], betayf2_rescale[bn1][0], betayf2_rescale[bn1][1]]
                tfs_file.add_table_row(list_row_entries)

    return beta_d
# END calculate_beta_from_amplitude ----------------------------------------------------------------


def beta_from_phase(madTwiss, madElements, madElementsCentre, ListOfFiles, phase, plane, getllm_d, debugfile):
    '''
    Calculate the beta function from phase advances
    If range of BPMs is sufficiently large use averaging with weighted mean. The weights are determinde using either the
    output of Monte Carlo simulations or using analytical formulas derived from exacter formulas to calculate the beta
    function.
    

    :Return:
        'data':table
            columns:
             0: beti
             1: betstat
             2: betsys
             3: beterr
             4: alfi
             5: alfstat
             6: alfsys
             7: alferr
             8: corr
             9: delbeta
            10: n_combinations
        'rmsbb':float
            rms beta-beating
        'commonbpms':list
            intersection of common BPMs in measurement files and model

    
    :Parameters:
        'madTwiss':twiss
            model twiss file
        'ListOfFiles':twiss
            measurement files
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'use_only_three_bpms_for_beta_from_phase':Boolean
            self-explanatory
        'number_of_bpms':int
            use the number_of_bpm betas with the lowest phase advance (only used when monte carlo simulations are used)
        'range_of_bpms':int
            length of range of bpms used to calculate the beta function. Has to be an odd number 1 < r < 21
        'errordefspath':String
            the paths to the error definition file. If not specified, Monte Carlo simulation is used
    '''
    
    errors_method = "NOT DEFINED"
    if phase == {}:
        return [{}, 0.0, {}, errors_method]

    commonbpms = utils.bpm.intersect(ListOfFiles)
    commonbpms = utils.bpm.model_intersect(commonbpms, madTwiss)
    commonbpms = JPARC_intersect(plane, getllm_d, commonbpms)

    
    errorfile = None
    if not getllm_d.use_only_three_bpms_for_beta_from_phase:
        errorfile = create_errorfile(getllm_d.errordefspath, madTwiss, madElements, madElementsCentre, commonbpms, plane, getllm_d.accel)
   
    if 3 > len(commonbpms):
        print_("beta_from_phase: Less than three BPMs for plane " + plane + ". Returning empty values.")
        return ({}, 0.0, {}, errors_method)

    if 7 > len(commonbpms) and errorfile is None:
        print_("beta_from_phase: Less than seven BPMs for plane " + plane + ". Can not use optimized algorithm.")

    errors_method = "Covariance matrix"
    rmsbb = 0.0
    
    #---- Error definitions given and we decided to not use the simulation => use analytical formulas to calculate the
    # systematic errors
    if errorfile is not None and not getllm_d.use_only_three_bpms_for_beta_from_phase:
        
        rmsbb, errors_method, data = scan_all_BPMs_withsystematicerrors(madTwiss, errorfile, phase, plane, getllm_d, commonbpms, debugfile)
        return data, rmsbb, commonbpms, errors_method
    #---- use the simulations
    else:
        
        rmsbb, errors_method, data = scan_all_BPMs_sim_3bpm(madTwiss, phase, plane, getllm_d, commonbpms, debugfile)

    return data, rmsbb, commonbpms, errors_method


def beta_from_amplitude(mad_twiss, list_of_files, plane):

    beta = {}
    root2j = []
    commonbpms = utils.bpm.intersect(list_of_files)
    commonbpms = utils.bpm.model_intersect(commonbpms, mad_twiss)
    sum_a = 0.0
    amp = []
    amp2 = []
    kick2 = []
    for i in range(0, len(commonbpms)):  # this loop have become complicated after modifications... anybody simplify?
        bn1 = str.upper(commonbpms[i][1])
        if plane == 'H':
            tembeta = mad_twiss.BETX[mad_twiss.indx[bn1]]
        elif plane == 'V':
            tembeta = mad_twiss.BETY[mad_twiss.indx[bn1]]
        amp_i = 0.0
        amp_j2 = []
        root2j_i = 0.0
        counter = 0
        for tw_file in list_of_files:
            if i == 0:
                kick2.append(0)
            if plane == 'H':
                amp_i += tw_file.AMPX[tw_file.indx[bn1]]
                amp_j2.append(tw_file.AMPX[tw_file.indx[bn1]] ** 2)
                root2j_i += tw_file.PK2PK[tw_file.indx[bn1]] / 2.
            elif plane == 'V':
                amp_i += tw_file.AMPY[tw_file.indx[bn1]]
                amp_j2.append(tw_file.AMPY[tw_file.indx[bn1]] ** 2)
                root2j_i += tw_file.PK2PK[tw_file.indx[bn1]] / 2.

            kick2[counter] += amp_j2[counter] / tembeta
            counter += 1

        amp_i = amp_i / len(list_of_files)
        root2j_i = root2j_i / len(list_of_files)
        amp.append(amp_i)
        amp2.append(amp_j2)

        sum_a += amp_i ** 2 / tembeta
        root2j.append(root2j_i / math.sqrt(tembeta))

    kick = sum_a / len(commonbpms)  # Assuming the average of beta is constant
    kick2 = np.array(kick2)
    kick2 = kick2 / len(commonbpms)
    amp2 = np.array(amp2)
    root2j = np.array(root2j)
    root2j_ave = np.average(root2j)
    root2j_rms = math.sqrt(np.average(root2j * root2j) - root2j_ave**2 + 2.2e-16)

    delbeta = []
    for i in range(0, len(commonbpms)):
        bn1 = str.upper(commonbpms[i][1])
        location = commonbpms[i][0]
        for j in range(0, len(list_of_files)):
            amp2[i][j] = amp2[i][j] / kick2[j]
        #print np.average(amp2[i]*amp2[i]),np.average(amp2[i])**2
        try:
            betstd = math.sqrt(np.average(amp2[i] * amp2[i]) - np.average(amp2[i])**2 + 2.2e-16)
        except:
            betstd = 0

        beta[bn1] = [amp[i] ** 2 / kick, betstd, location]
        if plane == 'H':
            betmdl = mad_twiss.BETX[mad_twiss.indx[bn1]]
        elif plane == 'V':
            betmdl = mad_twiss.BETY[mad_twiss.indx[bn1]]
        delbeta.append((beta[bn1][0] - betmdl) / betmdl)

    invariant_j = [root2j_ave, root2j_rms]

    delbeta = np.array(delbeta)
    rmsbb = math.sqrt(np.average(delbeta * delbeta))
    return [beta, rmsbb, commonbpms, invariant_j]


#=======================================================================================================================
#---============== using the simulations to calculate the beta function and error bars =================================
#=======================================================================================================================

def scan_all_BPMs_sim_3bpm(madTwiss, phase, plane, getllm_d, commonbpms, debugfile):
    systematics_error_path = os.path.join(os.path.dirname(os.path.abspath(madTwiss.filename)), "bet_deviations.npy")
    systematic_errors = None
    
    montecarlo = True
    errors_method = "Monte-Carlo Simulations"
    
    use_only_three_bpms_for_beta_from_phase = getllm_d.use_only_three_bpms_for_beta_from_phase
    
    if use_only_three_bpms_for_beta_from_phase:
        montecarlo = False
        errors_method = "Standard (3 BPM method)"
    elif not os.path.isfile(systematics_error_path):
        montecarlo = False
        errors_method = "Stdandard (3 BPM method) because no bet_deviations.npy could be found"
        use_only_three_bpms_for_beta_from_phase = True
    else:
        systematic_errors = np.load(systematics_error_path)
        
    data = {}
    used_bpms = 0

    print_("Errors from " + errors_method)
    for i in range(0, len(commonbpms)):
        alfa_beta, probed_bpm_name, M = get_best_three_bpms_with_beta_and_alfa(madTwiss, phase, plane, commonbpms, i, use_only_three_bpms_for_beta_from_phase, getllm_d.number_of_bpms, getllm_d.range_of_bpms)
        alfi = sum([alfa_beta[i][3] for i in range(len(alfa_beta))]) / len(alfa_beta)
        alfstd = math.sqrt(sum([alfa_beta[i][2] ** 2 for i in range(len(alfa_beta))])) / math.sqrt(len(alfa_beta))
        try:
            alferr = math.sqrt(sum([alfa_beta[i][3] ** 2 for i in range(len(alfa_beta))]) / len(alfa_beta) - alfi ** 2.)
        except ValueError:
            alferr = 0
        if plane == 'H':
            betmdl1 = madTwiss.BETX[madTwiss.indx[probed_bpm_name]]
        elif plane == 'V':
            betmdl1 = madTwiss.BETY[madTwiss.indx[probed_bpm_name]]
        beti = DEFAULT_WRONG_BETA
        if montecarlo:
            T = np.transpose(np.matrix([alfa_beta[i][7] for i in range(len(alfa_beta))]))
            V1 = M * T
            V_stat = np.transpose(T) * V1
            V = np.zeros([len(alfa_beta), len(alfa_beta)])
            V_syst = np.zeros([len(alfa_beta), len(alfa_beta)])
            np.fill_diagonal(V_syst, 100000)
            #all_comb = [[alfa_beta[i][5], alfa_beta[i][6]] for i in range(len(alfa_beta))]
            if plane == 'H':
                sindex = 0
            elif plane == 'V':
                sindex = 1
            for comb1 in range(len(alfa_beta)):
                for comb2 in range(len(alfa_beta)):
                    possible_dict_key = [''.join([probed_bpm_name, alfa_beta[comb1][i], alfa_beta[comb1][j], alfa_beta[comb2][k], alfa_beta[comb2][l]]) for i in [5, 6] for j in [5, 6] if i != j for k in [5, 6] for l in [5, 6] if k != l]
                    for keyname in possible_dict_key:
                        if keyname in systematic_errors[sindex]:
                            V_syst[comb1][comb2] = systematic_errors[sindex][keyname] * betmdl1 ** 2
                            V_syst[comb2][comb1] = systematic_errors[sindex][keyname] * betmdl1 ** 2
            
            for k in range(len(alfa_beta)):
                for l in range(len(alfa_beta)):
                    V[k][l] = V_stat.item(k, l) + V_syst[k][l]
            
            try:
                V_inv = np.linalg.pinv(V)
            except:
                V_inv = np.diag(1.0, V.shape[0])
            if DEBUG:
                debugfile.write("\n\nbegin BPM " + probed_bpm_name + " :\n")
                printMatrix(debugfile, V_stat, "V_stat")
                printMatrix(debugfile, np.matrix(V_syst), "V_syst")
                printMatrix(debugfile, V, "V")
                printMatrix(debugfile, V_inv, "V_inv")
                printMatrix(debugfile, np.transpose(T), "T_trans")
                printMatrix(debugfile, M, "M")
            w = np.zeros(len(alfa_beta))
            V_inv_row_sum = V_inv.sum(axis=1, dtype='float')
            V_inv_sum = V_inv.sum(dtype='float')
            betstd = 0
            beterr = 0
            if V_inv_sum != 0:
                for i in range(len(w)):
                    w[i] = V_inv_row_sum[i] / V_inv_sum
                
                beti = float(sum([w[i] * alfa_beta[i][1] for i in range(len(alfa_beta))]))
                for i in range(len(alfa_beta)):
                    for j in range(len(alfa_beta)):
                        betstd = betstd + w[i] * w[j] * V_stat.item(i, j)
                
                betstd = np.sqrt(float(betstd))
                for i in range(len(alfa_beta)):
                    for j in range(len(alfa_beta)):
                        beterr = beterr + w[i] * w[j] * V_syst.item(i, j)
                if beterr < 0:
                    beterr = DEFAULT_WRONG_BETA
                else:
                    beterr = np.sqrt(float(beterr))
                used_bpms = len(w)
                if DEBUG:
                    debugfile.write("\ncombinations:\t")
                    for i in range(len(w)):
                        debugfile.write(alfa_beta[i][8] + "\t")
                    
                    debugfile.write("\n")
                    debugfile.write("\n")
                    debugfile.write("\nweights:\t")
                    #print "beti =", beti
                    for i in range(len(w)):
                        debugfile.write("{:.7f}".format(w[i]) + "\t")
                    
                    debugfile.write("\nbeta_values:\t")
                    for i in range(len(w)):
                        debugfile.write(str(alfa_beta[i][1]) + "\t")
                    
                    debugfile.write("\n")
                    debugfile.write("averaged beta: " + str(beti) + "\n")
                    debugfile.write("\nalfa_values:\t")
                    for i in range(len(w)):
                        debugfile.write(str(alfa_beta[i][3]) + "\t")
                    
                    debugfile.write("\n")
            else:
                betstd = DEFAULT_WRONG_BETA
                beterr = DEFAULT_WRONG_BETA
        else:
            used_bpms = -1
            beti = sum([alfa_beta[i][1] for i in range(len(alfa_beta))]) / len(alfa_beta)
            betstd = math.sqrt(sum([alfa_beta[i][0] ** 2 for i in range(len(alfa_beta))])) / math.sqrt(len(alfa_beta))
            try:
                beterr = math.sqrt(sum([alfa_beta[i][1] ** 2 for i in range(len(alfa_beta))]) / len(alfa_beta) - beti ** 2.)
            except ValueError:
                beterr = 0


        data[probed_bpm_name] = [beti, betstd, beterr, math.sqrt(beterr ** 2 + betstd ** 2),
                                 alfi, alfstd, alferr, math.sqrt(alferr ** 2 + alfstd ** 2),
                                 0.0,
                                 (beti - betmdl1) / betmdl1,
                                 used_bpms]
       
        if DEBUG:
            debugfile.write("end\n")

    if DEBUG:
        debugfile.close()
    return 0, errors_method, data


def get_best_three_bpms_with_beta_and_alfa(madTwiss, phase, plane, commonbpms, i, use_only_three_bpms_for_beta_from_phase, number_of_bpms, range_of_bpms):
    '''
    Sorts BPM sets for the beta calculation based on their phase advance.
    If less than 7 BPMs are available it will fall back to using only next neighbours.
    :Parameters:
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'commonbpms':list
            intersection of common BPMs in measurement files and model
        'i': integer
            current iterator in the loop of all BPMs
    :Return: tupel(candidate1, candidate2, candidate3, bn4)
        'candidate1-3':list
            contains calculated beta and alfa
        'bn4':string
            name of the probed BPM
    '''

    RANGE = int(range_of_bpms)
    probed_index = int((RANGE - 1) / 2.)
 
    if 7 > len(commonbpms) or use_only_three_bpms_for_beta_from_phase:
        bn1 = str.upper(commonbpms[(probed_index + i - 2) % len(commonbpms)][1])
        bn2 = str.upper(commonbpms[(probed_index + i - 1) % len(commonbpms)][1])
        bn3 = str.upper(commonbpms[(probed_index + i) % len(commonbpms)][1])
        bn4 = str.upper(commonbpms[(probed_index + i + 1) % len(commonbpms)][1])
        bn5 = str.upper(commonbpms[(probed_index + i + 2) % len(commonbpms)][1])
        candidates = []
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, _ = _beta_from_phase_BPM_right(bn1, bn2, bn3, madTwiss, phase, plane, 0, 0, 0, True)
        candidates.append([tbetstd, tbet, talfstd, talf])
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, _ = _beta_from_phase_BPM_mid(bn2, bn3, bn4, madTwiss, phase, plane, 0, 0, 0, True)
        candidates.append([tbetstd, tbet, talfstd, talf])
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, _ = _beta_from_phase_BPM_left(bn3, bn4, bn5, madTwiss, phase, plane, 0, 0, 0, True)
        candidates.append([tbetstd, tbet, talfstd, talf])
        return candidates, bn3, []

    bpm_name = {}
    for n in range(RANGE):
        bpm_name[n] = str.upper(commonbpms[(i + n) % len(commonbpms)][1])
        phase_err = {}
    if plane == 'H':
        for i in range(RANGE):
            if i < probed_index:
                phase_err[i] = phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[i]]] / madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]])
            elif i > probed_index:
                phase_err[i] = phase["".join([plane, bpm_name[probed_index], bpm_name[i]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[i]]] / madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]])
        phase_err[probed_index] = min([phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETX[madTwiss.indx[bpm_name[i]]]) for i in range(probed_index)] + [phase["".join([plane, bpm_name[probed_index], bpm_name[probed_index + 1 + i]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETX[madTwiss.indx[bpm_name[probed_index + 1 + i]]]) for i in range(probed_index)])
    if plane == 'V':
        for i in range(RANGE):
            if i < probed_index:
                phase_err[i] = phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[i]]] / madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]])
            if i > probed_index:
                phase_err[i] = phase["".join([plane, bpm_name[probed_index], bpm_name[i]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[i]]] / madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]])
        phase_err[probed_index] = min([phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETY[madTwiss.indx[bpm_name[i]]]) for i in range(probed_index)] + [phase["".join([plane, bpm_name[probed_index], bpm_name[probed_index + 1 + i]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETY[madTwiss.indx[bpm_name[probed_index + 1 + i]]]) for i in range(probed_index)])
   
    M = np.zeros([RANGE - 1, RANGE - 1])
    for k in range(RANGE - 1):
        for l in range(RANGE - 1):
            if k == l and k < probed_index:
                M[k][l] = (2*np.pi)**2 * phase["".join([plane, bpm_name[probed_index], bpm_name[probed_index + k + 1]])][1]**2
            elif k == l and k >= probed_index:
                M[k][l] = (2*np.pi)**2 * phase["".join([plane, bpm_name[RANGE - 2 - k], bpm_name[probed_index]])][1]**2
            elif (k < probed_index and l >= probed_index) or (k >= probed_index and l < probed_index):
                M[k][l] = -(2*np.pi)**2 * phase_err[probed_index]**2
            else:
                M[k][l] = (2*np.pi)**2 * phase_err[probed_index]**2
    candidates = []
    left_bpm = range(probed_index)
    right_bpm = range(probed_index + 1, RANGE)
    left_combo = [[x, y] for x in left_bpm for y in left_bpm if x < y]
    right_combo = [[x, y] for x in right_bpm for y in right_bpm if x < y]
    mid_combo = [[x, y] for x in left_bpm for y in right_bpm]
    
    #right_combo = []    #ABB
    #mid_combo = []      # BAB
    #left_combo = []     # BBA
 
    for n in left_combo:
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, use_it = _beta_from_phase_BPM_right(bpm_name[n[0]], bpm_name[n[1]], bpm_name[probed_index], madTwiss, phase, plane, phase_err[n[0]], phase_err[n[1]], phase_err[probed_index])
        if use_it:
            t_matrix_row = [0] * (RANGE-1)
            t_matrix_row[RANGE-2 - n[0]] = t1
            t_matrix_row[RANGE-2 - n[1]] = t2
            patternstr = ["x"] * RANGE
            patternstr[n[0]] = "B"
            patternstr[n[1]] = "B"
            patternstr[probed_index] = "A"
            candidates.append([tbetstd, tbet, talfstd, talf, mdlerr, bpm_name[n[0]], bpm_name[n[1]], t_matrix_row, "".join(patternstr)])

    for n in mid_combo:
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, use_it = _beta_from_phase_BPM_mid(bpm_name[n[0]], bpm_name[probed_index], bpm_name[n[1]], madTwiss, phase, plane, phase_err[n[0]], phase_err[probed_index], phase_err[n[1]])
        if use_it:
            t_matrix_row = [0] * (RANGE - 1)
            t_matrix_row[RANGE - 2 - n[0]] = t1
            t_matrix_row[n[1] - 1 - probed_index] = t2
            patternstr = ["x"] * RANGE
            patternstr[n[0]] = "B"
            patternstr[n[1]] = "B"
            patternstr[probed_index] = "A"
            candidates.append([tbetstd, tbet, talfstd, talf, mdlerr, bpm_name[n[0]], bpm_name[n[1]], t_matrix_row, "".join(patternstr)])

    for n in right_combo:
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, use_it = _beta_from_phase_BPM_left(bpm_name[probed_index], bpm_name[n[0]], bpm_name[n[1]], madTwiss, phase, plane, phase_err[probed_index], phase_err[n[0]], phase_err[n[1]])
        if use_it:
            t_matrix_row = [0] * (RANGE-1)
            t_matrix_row[n[0] - 1 - probed_index] = t1
            t_matrix_row[n[1] - 1 - probed_index] = t2
            patternstr = ["x"] * RANGE
            patternstr[n[0]] = "B"
            patternstr[n[1]] = "B"
            patternstr[probed_index] = "A"
            candidates.append([tbetstd, tbet, talfstd, talf, mdlerr, bpm_name[n[0]], bpm_name[n[1]], t_matrix_row, "".join(patternstr)])

    #sort_cand = sorted(candidates, key=lambda x: x[4])
    #return [sort_cand[i] for i in range(NUM_BPM_COMBOS)], bpm_name[probed_index], M
    return candidates, bpm_name[probed_index], M


def _beta_from_phase_BPM_left(bn1, bn2, bn3, madTwiss, phase, plane, p1, p2, p3, is3BPM=False):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is left of the other two BPMs (ABB)
    :Parameters:
        'bn1':string
            Name of probed BPM
        'bn2':string
            Name of first BPM right of the probed BPM
        'bn3':string
            Name of second BPM right of the probed BPM
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'betstd':float
            calculated error on beta function at probed BPM
        'alf':float
            calculated alfa function at probed BPM
        'alfstd':float
            calculated error on alfa function at probed BPM
    '''
    ph2pi12=2.*np.pi*phase["".join([plane,bn1,bn2])][0]
    ph2pi13=2.*np.pi*phase["".join([plane,bn1,bn3])][0]

    # Find the model transfer matrices for beta1
    phmdl12 = 2.*np.pi*phase["".join([plane,bn1,bn2])][2]
    phmdl13=2.*np.pi*phase["".join([plane,bn1,bn3])][2]
    if plane=='H':
        betmdl1=madTwiss.BETX[madTwiss.indx[bn1]]
        betmdl2=madTwiss.BETX[madTwiss.indx[bn2]]
        betmdl3=madTwiss.BETX[madTwiss.indx[bn3]]
        alpmdl1=madTwiss.ALFX[madTwiss.indx[bn1]]
    elif plane=='V':
        betmdl1=madTwiss.BETY[madTwiss.indx[bn1]]
        betmdl2=madTwiss.BETY[madTwiss.indx[bn2]]
        betmdl3=madTwiss.BETY[madTwiss.indx[bn3]]
        alpmdl1=madTwiss.ALFY[madTwiss.indx[bn1]]
    if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)
        
    if not is3BPM and (
        bad_phase(phmdl12) or bad_phase(phmdl13) or bad_phase(phmdl12 - phmdl13) or bad_phase(ph2pi12) or bad_phase(ph2pi13) or bad_phase(ph2pi12 - ph2pi13)):
        return 0, 0, 0, 0, 0, 0, 0, False

    # Find beta1 and alpha1 from phases assuming model transfer matrix
    # Matrix M: BPM1-> BPM2
    # Matrix N: BPM1-> BPM3
    M11=math.sqrt(betmdl2/betmdl1)*(cos(phmdl12)+alpmdl1*sin(phmdl12))
    M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
    N11=math.sqrt(betmdl3/betmdl1)*(cos(phmdl13)+alpmdl1*sin(phmdl13))
    N12=math.sqrt(betmdl1*betmdl3)*sin(phmdl13)

    denom=M11/M12-N11/N12+1e-16
    numer=1/tan(ph2pi12)-1/tan(ph2pi13)
    bet=numer/denom

    betstd=        (2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
    betstd=betstd+(2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
    betstd=math.sqrt(betstd)/abs(denom)

    mdlerr=        (2*np.pi*0.001/sin(phmdl12)**2)**2
    mdlerr=mdlerr+(2*np.pi*0.001/sin(phmdl13)**2)**2
    mdlerr=math.sqrt(mdlerr)/abs(denom)    

    term1 = 1/sin(phmdl12)**2/denom
    term2 = -1/sin(phmdl13)**2/denom

    denom=M12/M11-N12/N11+1e-16
    numer=-M12/M11/tan(ph2pi12)+N12/N11/tan(ph2pi13)
    alf=numer/denom

    alfstd=        (M12/M11*2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
    alfstd=alfstd+(N12/N11*2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
    alfstd=math.sqrt(alfstd)/denom

    return bet, betstd, alf, alfstd, mdlerr, term1, term2, True

def _beta_from_phase_BPM_mid(bn1,bn2,bn3,madTwiss,phase,plane,p1,p2,p3, is3BPM=False):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is between the other two BPMs
    :Parameters:
        'bn1':string
            Name of BPM left of the probed BPM
        'bn2':string
            Name of probed BPM
        'bn3':string
            Name of BPM right of the probed BPM
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'betstd':float
            calculated error on beta function at probed BPM
        'alf':float
            calculated alfa function at probed BPM
        'alfstd':float
            calculated error on alfa function at probed BPM
    '''
    ph2pi12=2.*np.pi*phase["".join([plane,bn1,bn2])][0]
    ph2pi23=2.*np.pi*phase["".join([plane,bn2,bn3])][0]

    # Find the model transfer matrices for beta1
    phmdl12=2.*np.pi*phase["".join([plane,bn1,bn2])][2]
    phmdl23=2.*np.pi*phase["".join([plane,bn2,bn3])][2]
    if plane=='H':
        betmdl1=madTwiss.BETX[madTwiss.indx[bn1]]
        betmdl2=madTwiss.BETX[madTwiss.indx[bn2]]
        betmdl3=madTwiss.BETX[madTwiss.indx[bn3]]
        alpmdl2=madTwiss.ALFX[madTwiss.indx[bn2]]
    elif plane=='V':
        betmdl1=madTwiss.BETY[madTwiss.indx[bn1]]
        betmdl2=madTwiss.BETY[madTwiss.indx[bn2]]
        betmdl3=madTwiss.BETY[madTwiss.indx[bn3]]
        alpmdl2=madTwiss.ALFY[madTwiss.indx[bn2]]
    if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)
        
    if not is3BPM and (
        bad_phase(phmdl12) or bad_phase(phmdl23) or bad_phase(phmdl23 - phmdl12) or 
        bad_phase(ph2pi12) or bad_phase(ph2pi23) or bad_phase(ph2pi12 - ph2pi23)):
        return 0, 0, 0, 0, 0, 0, 0, False


    # Find beta2 and alpha2 from phases assuming model transfer matrix
    # Matrix M: BPM1-> BPM2
    # Matrix N: BPM2-> BPM3
    M22=math.sqrt(betmdl1/betmdl2)*(cos(phmdl12)-alpmdl2*sin(phmdl12))
    M12=math.sqrt(betmdl1*betmdl2)*sin(phmdl12)
    N11=math.sqrt(betmdl3/betmdl2)*(cos(phmdl23)+alpmdl2*sin(phmdl23))
    N12=math.sqrt(betmdl2*betmdl3)*sin(phmdl23)

    denom=M22/M12+N11/N12+EPSILON
    #denom = (1.0/tan(phmdl12) + 1.0/tan(phmdl23) +EPSILON)/betmdl2
    numer=1/tan(ph2pi12)+1/tan(ph2pi23)
    bet=numer/denom

    betstd=        (2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
    betstd=betstd+(2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
    betstd=math.sqrt(betstd)/abs(denom)

    mdlerr=        (2*np.pi*0.001/sin(phmdl12)**2)**2
    mdlerr=mdlerr+(2*np.pi*0.001/sin(phmdl23)**2)**2
    mdlerr=math.sqrt(mdlerr)/abs(denom)

    term2 = 1/sin(phmdl23)**2/denom  #sign
    term1 = -1/sin(phmdl12)**2/denom  #sign

    denom=M12/M22+N12/N11+1e-16
    numer=M12/M22/tan(ph2pi12)-N12/N11/tan(ph2pi23)
    alf=numer/denom

    alfstd=        (M12/M22*2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
    alfstd=alfstd+(N12/N11*2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
    alfstd=math.sqrt(alfstd)/abs(denom)

    return bet, betstd, alf, alfstd, mdlerr, term1, term2, True

def _beta_from_phase_BPM_right(bn1,bn2,bn3,madTwiss,phase,plane,p1,p2,p3, is3BPM=False):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is right the other two BPMs
    :Parameters:
        'bn1':string
            Name of second BPM left of the probed BPM
        'bn2':string
            Name of first BPM left of the probed BPM
        'bn3':string
            Name of probed BPM
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'betstd':float
            calculated error on beta function at probed BPM
        'alf':float
            calculated alfa function at probed BPM
        'alfstd':float
            calculated error on alfa function at probed BPM
    '''
    ph2pi23=2.*np.pi*phase["".join([plane,bn2,bn3])][0]
    ph2pi13=2.*np.pi*phase["".join([plane,bn1,bn3])][0]

    # Find the model transfer matrices for beta1
    phmdl13=2.*np.pi*phase["".join([plane,bn1,bn3])][2]
    phmdl23=2.*np.pi*phase["".join([plane,bn2,bn3])][2]
    if plane=='H':
        betmdl1=madTwiss.BETX[madTwiss.indx[bn1]]
        betmdl2=madTwiss.BETX[madTwiss.indx[bn2]]
        betmdl3=madTwiss.BETX[madTwiss.indx[bn3]]
        alpmdl3=madTwiss.ALFX[madTwiss.indx[bn3]]
    elif plane=='V':
        betmdl1=madTwiss.BETY[madTwiss.indx[bn1]]
        betmdl2=madTwiss.BETY[madTwiss.indx[bn2]]
        betmdl3=madTwiss.BETY[madTwiss.indx[bn3]]
        alpmdl3=madTwiss.ALFY[madTwiss.indx[bn3]]
    if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)

    if not is3BPM and (
        bad_phase(phmdl23) or bad_phase(phmdl13) or bad_phase(phmdl13 - phmdl23) or
        bad_phase(ph2pi13) or bad_phase(ph2pi23) or bad_phase(ph2pi13 - ph2pi23)):
        return 0, 0, 0, 0, 0, 0, 0, False
    # Find beta3 and alpha3 from phases assuming model transfer matrix
    # Matrix M: BPM2-> BPM3
    # Matrix N: BPM1-> BPM3
    M22=math.sqrt(betmdl2/betmdl3)*(cos(phmdl23)-alpmdl3*sin(phmdl23))
    M12=math.sqrt(betmdl2*betmdl3)*sin(phmdl23)
    N22=math.sqrt(betmdl1/betmdl3)*(cos(phmdl13)-alpmdl3*sin(phmdl13))
    N12=math.sqrt(betmdl1*betmdl3)*sin(phmdl13)

    denom=M22/M12-N22/N12+1e-16
    numer=1/tan(ph2pi23)-1/tan(ph2pi13)
    bet=numer/denom

    betstd=        (2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
    betstd=betstd+(2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
    betstd=math.sqrt(betstd)/abs(denom)

    mdlerr=        (2*np.pi*0.001/sin(phmdl23)**2)**2
    mdlerr=mdlerr+(2*np.pi*0.001/sin(phmdl13)**2)**2
    mdlerr=math.sqrt(mdlerr)/abs(denom)

    term2 = -1/sin(phmdl23)**2/denom  #sign
    term1 = 1/sin(phmdl13)**2/denom  #sign

    denom=M12/M22-N12/N22+1e-16
    numer=M12/M22/tan(ph2pi23)-N12/N22/tan(ph2pi13)
    alf=numer/denom

    alfstd=        (M12/M22*2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
    alfstd=alfstd+(N12/N22*2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
    alfstd=math.sqrt(alfstd)/abs(denom)


    return bet, betstd, alf, alfstd, mdlerr, term1, term2, True

#=======================================================================================================================
#---============== using analytical formula ============================================================================
#=======================================================================================================================



def scan_all_BPMs_withsystematicerrors(madTwiss, errorfile, phase, plane, getllm_d, commonbpms, debugfile):
    '''
    If errorfile is given (!= "0") this function calculates the beta function for each BPM using the analytic expression
    for the estimation of the error matrix.
    
    :Parameters:
        'madTwiss':twiss
            model twiss file
        'errorfile':twiss
            final errorfile (not error definitions)
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'alfa, beta':vector
            out vectors which will be filled with the beta nad alfa functions
        'commonbpms':list
            intersection of common BPMs in measurement files and model
    :Returns:
        'rmsbb'
        'errors_method':string
        'result':dictionary
            columns:
           key: probed_bpm_name
             0: beti
             1: betstat
             2: betsys
             3: beterr
             4: alfi
             5: alfstat
             6: alfsys
             7: alferr
             8: corr
             9: delbeta
    '''
    
    print_("INFO: errorfile given. Create list_of_Ks")
    list_of_Ks = []
    
    errors_method = "Analytical Formula"
    print_("Errors from " + errors_method)
   
    #---- create list of Ks, i.e. assign to each BPM a vector with all the errore lements that come after the bpm
    # and their respective errors
    # and model phases so that list_of_Ks[n] yields the error elements between BPM[n] and BPM[n+1]
    # update 2016-07-28: list_of_Ks[n][k], n: BPM number, k=0: quadrupole field errors,
    # k=1: transversal sextupole missalignments
    # k=2: longitudinal quadrupole missalignments
    for n in range(len(commonbpms) + getllm_d.range_of_bpms + 1):
        index_n = errorfile.indx[commonbpms[n % len(commonbpms)][1]]
        index_nplus1 = errorfile.indx[commonbpms[(n + 1) % len(commonbpms)][1]]
              
        quad_fields = []
        sext_trans = []
        quad_missal = []
               
        if index_n < index_nplus1:
            for i in range(index_n + 1, index_nplus1):
                if errorfile.dK1[i] != 0:
                    quad_fields.append(i)
                if errorfile.dX[i] != 0:
                    sext_trans.append(i)
                if errorfile.dS[i] != 0:
                    quad_missal.append(i)
                    
        else:
            for i in range(index_n + 1, len(errorfile.NAME)):
                if errorfile.dK1[i] != 0:
                    quad_fields.append(i)
                if errorfile.dX[i] != 0:
                    sext_trans.append(i)
                if errorfile.dS[i] != 0:
                    quad_missal.append(i)
                    
            for i in range(index_nplus1):  # ums Eck
                if errorfile.dK1[i] != 0:
                    quad_fields.append(i)
                if errorfile.dX[i] != 0:
                    sext_trans.append(i)
                if errorfile.dS[i] != 0:
                    quad_missal.append(i)
       
        list_of_Ks.append([quad_fields, sext_trans, quad_missal])
              
    width = getllm_d.range_of_bpms / 2
    left_bpm = range(-width, 0)
    right_bpm = range(0 + 1, width + 1)
    BBA_combo = [[x, y] for x in left_bpm for y in left_bpm if x < y]
    ABB_combo = [[x, y] for x in right_bpm for y in right_bpm if x < y]
    BAB_combo = [[x, y] for x in left_bpm for y in right_bpm]
    result = {}
    
    def collect(row):
        if row[11]:
            result[row[0]] = row[1:]
        
    def collectblock(block):
        for row in block:
            if row[11]:
                result[row[0]] = row[1:]
                
    st = time.time()
    if getllm_d.parallel and not DEBUG:
        
        chunksize = int(len(commonbpms) / getllm_d.nprocesses) + 1
        pool = multiprocessing.Pool()
        n = int(len(commonbpms) / chunksize)
        
        for i in range(n):
            pool.apply_async(scan_several_BPMs_withsystematicerrors,
                             (madTwiss, errorfile,
                              phase, plane, getllm_d.range_of_bpms, commonbpms, debugfile, list_of_Ks,
                              i * chunksize, (i + 1) * chunksize, BBA_combo, ABB_combo, BAB_combo),
                             callback=collectblock)
        pool.apply_async(scan_several_BPMs_withsystematicerrors,
                         (madTwiss, errorfile,
                          phase, plane, getllm_d.range_of_bpms, commonbpms, debugfile, list_of_Ks,
                          n * chunksize, len(commonbpms), BBA_combo, ABB_combo, BAB_combo),
                         callback=collectblock)
        pool.close()
        pool.join()
    else:
        startProgress("Scan all BPMs")
        for i in range(0, len(commonbpms)):
            if (i % 20 == 0):
                progress(float(i) * 100.0 / len(commonbpms))
            row = scan_one_BPM_withsystematicerrors(madTwiss, errorfile, phase, plane, getllm_d.range_of_bpms, commonbpms,
                                                    debugfile, list_of_Ks, i,
                                                    BBA_combo, ABB_combo, BAB_combo)
            collect(row)
        endProgress()
    et = time.time()
    
    print_("time elapsed = {0:3.3f}".format(et - st))
    if DEBUG:
        debugfile.close()
    rmsbb = -1
    return rmsbb, errors_method, result


def scan_several_BPMs_withsystematicerrors(madTwiss, errorfile,
                                           phase, plane, range_of_bpms, commonbpms, debugfile, list_of_Ks,
                                           begin, end, BBA_combo, ABB_combo, BAB_combo):
    block = []
    for i in range(begin, end):
        block.append(scan_one_BPM_withsystematicerrors(madTwiss, errorfile, phase, plane, range_of_bpms, commonbpms,
                                                       debugfile, list_of_Ks, i,
                                                       BBA_combo, ABB_combo, BAB_combo))
    return block
    

def get_beta_from_phase_3bpm(madTwiss, phase, plane, range_of_bpms, commonbpms, debugfile, i, probed_bpm_name_):
    alfa_beta, probed_bpm_name, _ = get_best_three_bpms_with_beta_and_alfa(madTwiss, phase, plane, commonbpms, i, True, 3, range_of_bpms)
    beti = sum([alfa_beta[i][1] for i in range(len(alfa_beta))]) / len(alfa_beta)
    betstat = math.sqrt(sum([alfa_beta[i][0] ** 2 for i in range(len(alfa_beta))])) / math.sqrt(len(alfa_beta))
    try:
        betsys = math.sqrt(sum([alfa_beta[i][1] ** 2 for i in range(len(alfa_beta))]) / len(alfa_beta) - beti ** 2.)
    except ValueError:
        betsys = DEFAULT_WRONG_BETA
    alfi = sum([alfa_beta[i][3] for i in range(len(alfa_beta))]) / len(alfa_beta)
    alfstat = math.sqrt(sum([alfa_beta[i][2] ** 2 for i in range(len(alfa_beta))])) / math.sqrt(len(alfa_beta))
    try:
        alfsys = math.sqrt(sum([alfa_beta[i][3] ** 2 for i in range(len(alfa_beta))]) / len(alfa_beta) - alfi ** 2.)
    except ValueError:
        alfsys = DEFAULT_WRONG_BETA
    if DEBUG:
        debugfile.write("\n\nbegin BPM " + probed_bpm_name + " 3BPM\n")
        debugfile.write("end\n")
    
    if probed_bpm_name != probed_bpm_name_:
        print "\33[31m ------------------------------------ PROBLEM ------------------\33[0m"
    return probed_bpm_name, beti, betstat, betsys, math.sqrt(betstat ** 2 + betsys ** 2), alfi, alfstat, alfsys, math.sqrt(alfstat ** 2 + alfsys ** 2)


def scan_one_BPM_withsystematicerrors(madTwiss, errorfile,
                                      phase, plane, range_of_bpms, commonbpms, debugfile, list_of_Ks,
                                      Index, BBA_combo, ABB_combo, BAB_combo):
    
    T_Alfa, T_Beta, alfa_beta, probed_bpm_name, M = get_beta_from_phase_systematic_errors(madTwiss, errorfile,
                                                                                          phase, plane, commonbpms, list_of_Ks, Index,
                                                                                          range_of_bpms,
                                                                                          ABB_combo, BAB_combo, BBA_combo)
    if plane == 'H':
        betmdl1 = madTwiss.BETX[madTwiss.indx[probed_bpm_name]]
    elif plane == 'V':
        betmdl1 = madTwiss.BETY[madTwiss.indx[probed_bpm_name]]
    
    beti        = DEFAULT_WRONG_BETA    #@IgnorePep8
    betstat     = .0                    #@IgnorePep8
    betsys      = .0                    #@IgnorePep8
    beterr      = DEFAULT_WRONG_BETA    #@IgnorePep8
    alfi        = DEFAULT_WRONG_BETA    #@IgnorePep8
    alfstat     = .0                    #@IgnorePep8
    alfsys      = .0                    #@IgnorePep8
    alferr      = DEFAULT_WRONG_BETA    #@IgnorePep8
    corr        = .0                    #@IgnorePep8
    used_bpms   = 0                     #@IgnorePep8
    
    if len(alfa_beta) == 0:
        print "\33[35mWARNING For bpm " + probed_bpm_name + " are 0 combinations with good angle left, FALLBACK to 3BPM.\33[0m"
   
        if DEBUG:
            debugfile.write("\n\nbegin BPM " + probed_bpm_name + " :\n")
            debugfile.write("ERROR: <= 3 BPM combinations left. USING 3BPM")
            
            printMatrix(debugfile, T_Beta, "T_b")
            printMatrix(debugfile, T_Alfa, "T_Alfa")
            printMatrix(debugfile, M, "M")
            debugfile.write("MUX: {0:.10e}\n".format(phase[probed_bpm_name][0]))
            debugfile.write("end\n")
            
        probed_bpm_name, beti, betstat, betsys, beterr, alfi, alfstat, alfsys, alferr = get_beta_from_phase_3bpm(madTwiss, phase, plane, range_of_bpms, commonbpms, debugfile, Index, probed_bpm_name)
        return [probed_bpm_name,
                beti, betstat, betsys, beterr,
                alfi, alfstat, alfsys, beterr,
                .0, (beti - betmdl1) / betmdl1,
                -1]
    try:
        betas = np.array([x.beta for x in alfa_beta])
        alfas = np.array([x.alfa for x in alfa_beta])
        len_w = len(betas)
    
        #--- calculate V and its inverse
        V_Beta = T_Beta * M * np.transpose(T_Beta)
        V_Alfa = T_Alfa * M * np.transpose(T_Alfa)
    except:
        print "\33[31;1m SOMETHING WENT WRONG\33[0m Use 3 BPM method"
        probed_bpm_name, beti, betstat, betsys, beterr, alfi, alfstat, alfsys, alferr = get_beta_from_phase_3bpm(madTwiss, phase, plane, range_of_bpms, commonbpms, debugfile, Index, probed_bpm_name)
        
        if DEBUG:
            
            # something went wrong, we have to log precisely what happened:
            
            debugfile.write("\n\nbegin BPM " + probed_bpm_name + " :\n")
            debugfile.write("ERROR: SOMETHING WENT WRONG in building V_Beta and V_alfa")
            
            printMatrix(debugfile, T_Beta, "T_b")
            printMatrix(debugfile, T_Alfa, "T_Alfa")
            printMatrix(debugfile, M, "M")
            printMatrix(debugfile, V_Beta, "Vb")
            printMatrix(debugfile, V_Alfa, "Va")
            debugfile.write("\ncombinations:\t")
            for i in range(len_w):
                debugfile.write(alfa_beta[i].patternstring + "\t")
            
            debugfile.write("\n")
            debugfile.write("\nweights:\t")
            for i in range(len_w):
                debugfile.write("{0:.7e}".format(0.0) + "\t")
            
            debugfile.write("\nbeta_values:\t")
            for i in range(len_w):
                debugfile.write(str(betas[i]) + "\t")
            
            debugfile.write("\n")
            debugfile.write("averaged beta: " + str(beti) + "\n")
            debugfile.write("\nalfa_values:\t")
            for i in range(len_w):
                debugfile.write(str(alfa_beta[i].alfa) + "\t")
            
            debugfile.write("\n")
            debugfile.write("MUX: {0:.10e}\n".format(phase[probed_bpm_name][0]))
            debugfile.write("end\n")
        
        return [probed_bpm_name,
                beti, betstat, betsys, beterr,
                alfi, alfstat, alfsys, beterr,
                .0, (beti - betmdl1) / betmdl1,
                -1]

    # TODO The LinalgError doesnt seem to be defined in 2.6 version of python, this should be checked.
    # TODO Should we really output a zero matrix if the pinv fails? Play with rcond in pinv.
    try:
        V_Beta_inv = np.linalg.pinv(V_Beta, rcond=RCOND)
        w = np.sum(V_Beta_inv, axis=1)
        VBeta_inv_sum = np.sum(w)
        beterr = math.sqrt(float(np.dot(np.transpose(w), np.dot(V_Beta, w)) / VBeta_inv_sum ** 2))
        beti = float(np.dot(np.transpose(w), betas) / VBeta_inv_sum)
        used_bpms = len(w)
    except:
        _, beti, betstat, betsys, beterr, _, _, _, _ = get_beta_from_phase_3bpm(madTwiss, phase, plane, range_of_bpms, commonbpms, debugfile, Index, probed_bpm_name, beti, alfi)

        used_bpms = -1
        print "WARN: LinAlgEror in V_Beta_inv for " + probed_bpm_name
        
    try:
        V_Alfa_inv = np.linalg.pinv(V_Alfa, rcond=RCOND)
        walfa = np.sum(V_Alfa_inv, axis=1)
        VAlfa_inv_sum = np.sum(walfa)
        
        alfi = float(np.dot(np.transpose(walfa), alfas) / VAlfa_inv_sum)
        alferr = math.sqrt(float(np.dot(np.transpose(walfa), np.dot(V_Alfa, walfa)) / VAlfa_inv_sum ** 2))
    except:
        _, _, _, _, _, alfi, alfstat, alfsys, alferr = get_beta_from_phase_3bpm(madTwiss, phase, plane, range_of_bpms, commonbpms, debugfile, Index, probed_bpm_name, beti, alfi)

        print "WARN: LinAlgEror in V_Alfa_inv for " + probed_bpm_name
    #--- calculate weights, weighted beta and beta_error

    #--- calculate correlation coefficient
                                                                                            
    #         for i in range(M.shape[0]):
    #             betaterm = 0.0
    #             alfaterm = 0.0
    #             for k in range(len_w):
    #                 betaterm += w[k] * t_rows_beta_reduced[k][i]
    #                 alfaterm -= walfa[k] * t_rows_alfa_reduced[k][i]
    #             rho_alfa_beta += betaterm * alfaterm * M[i][i]
    #         rho_alfa_beta /= (beterr * alferr)
    # so far, beta_syst and beta_stat are not separated
    #--- DEBUG
    if DEBUG:
        debugfile.write("\n\nbegin BPM " + probed_bpm_name + " :\n")
        printMatrix(debugfile, T_Beta, "T_b")
        printMatrix(debugfile, T_Alfa, "T_Alfa")
        printMatrix(debugfile, M, "M")
        printMatrix(debugfile, V_Beta, "Vb")
        printMatrix(debugfile, V_Beta_inv, "Vb_inv")
        printMatrix(debugfile, V_Alfa, "Va")
        printMatrix(debugfile, V_Alfa_inv, "Va_inv")
        Zero = V_Beta * V_Beta_inv * V_Beta - V_Beta
        sumnorm = Zero.sum(dtype='float') / (Zero.shape[0] * Zero.shape[1])
        debugfile.write("Zeros sum norm: " + str(sumnorm) + "\n")
        printMatrix(debugfile, Zero, "Zero")
        Einheit = np.dot(V_Beta, V_Beta_inv)
        printMatrix(debugfile, Einheit, "V_inv*V")
        debugfile.write("\ncombinations:\t")
        for i in range(len_w):
            debugfile.write(alfa_beta[i].patternstring + "\t")
        
        debugfile.write("\n")
        debugfile.write("\nweights:\t")
        for i in range(len_w):
            debugfile.write("{0:.7e}".format(float(w[i])) + "\t")
        
        debugfile.write("\nbeta_values:\t")
        for i in range(len_w):
            debugfile.write(str(betas[i]) + "\t")
        
        debugfile.write("\n")
        debugfile.write("averaged beta: " + str(beti) + "\n")
        debugfile.write("\nalfa_values:\t")
        for i in range(len_w):
            debugfile.write(str(alfa_beta[i].alfa) + "\t")
        
        debugfile.write("\n")
        debugfile.write("MUX: {0:.10e}\n".format(phase[probed_bpm_name][0]))
        debugfile.write("end\n")
        
    return [probed_bpm_name,
            beti, betstat, betsys, beterr,
            alfi, alfstat, alfsys, alferr,
            corr, (beti - betmdl1) / betmdl1,
            used_bpms]


def get_beta_from_phase_systematic_errors(madTwiss, errorfile, phase, plane, commonbpms,
                                          list_of_Ks, CurrentIndex, range_of_bpms,
                                          ABB_combo, BAB_combo, BBA_combo):
    '''
    Calculates the beta function at one BPM for all combinations.
    If less than 7 BPMs are available it will fall back to using only next neighbours.
    :Parameters:
        'madTwiss':twiss
            model twiss file
        'madTwiss':twiss
            final error twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'commonbpms':list
            intersection of common BPMs in measurement files and model
        'list_of_Ks':vector (of vectors)
            contains information about the errors the i-th entry in list_of_Ks yields the errors of all the lattice
            elements that come after the i-th BPM
        'CurrentIndex': integer
            current iterator in the loop of all BPMs
        'xyz_combo':vector
            contain all the BPM combinations for the three different cases ABB, BAB, BBA
    :Return: tupel(T_trans, betadata, bpm_name[probed_index], M)
        'T_Trans':np.matrix
            the T matrix to transform the variance matrix M into the covariance matrix V.
            In comparision with Andy's paper, the transposed is calculated, thus the name T_Trans
        'betadata':list
            contains calculated betas with information about the combination that was used
            consists of: betafunction, name of first used bpm, name of second used bpm, patternstring
            the pattern string is primarily for debugging purposes, to see which combination got which weight
            should be discarded afterwards
        'M':np-matrix
            variance matrix for all error kinds.
    '''
       
    RANGE = int(range_of_bpms)
    probed_index = int((RANGE - 1) / 2.)
    
    sizeOfMatrix = RANGE
    
        #---- add the block M_K
    for k in range(RANGE):
        sizeOfMatrix += len(list_of_Ks[(CurrentIndex + k) % len(list_of_Ks)][0])
        
        #---- add the block M_X, the sextupole transversal missalignments
    for k in range(RANGE):
        sizeOfMatrix += len(list_of_Ks[(CurrentIndex + k) % len(list_of_Ks)][1])
    
        #---- add the block M_S, the quadrupole longitudinal missalignments
    for k in range(RANGE):
        sizeOfMatrix += len(list_of_Ks[(CurrentIndex + k) % len(list_of_Ks)][2])
        
        #---- add the block M_BPM, the BPM missalignment errors
    sizeOfMatrix += RANGE
      
    #--- Make a dictionary<int, string> to get the names of the BPMs
    bpm_name = [""] * RANGE
    err_diagonal = np.zeros(sizeOfMatrix)
    for n in range(RANGE):
        bpm_name[n] = str.upper(commonbpms[(CurrentIndex + n) % len(commonbpms)][1])
    if plane == 'H':
        for i in range(RANGE):
            if i < probed_index:
                err_diagonal[i] = phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[i]]] / madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]])
            elif i > probed_index:
                err_diagonal[i] = phase["".join([plane, bpm_name[probed_index], bpm_name[i]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[i]]] / madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]])
        err_diagonal[probed_index] = min([phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETX[madTwiss.indx[bpm_name[i]]]) for i in range(probed_index)] + [phase["".join([plane, bpm_name[probed_index], bpm_name[probed_index + 1 + i]])][1] / np.sqrt(1 + madTwiss.BETX[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETX[madTwiss.indx[bpm_name[probed_index + 1 + i]]]) for i in range(probed_index)])
    if plane == 'V':
        for i in range(RANGE):
            if i < probed_index:
                err_diagonal[i] = phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[i]]] / madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]])
            if i > probed_index:
                err_diagonal[i] = phase["".join([plane, bpm_name[probed_index], bpm_name[i]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[i]]] / madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]])
        err_diagonal[probed_index] = min([phase["".join([plane, bpm_name[i], bpm_name[probed_index]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETY[madTwiss.indx[bpm_name[i]]]) for i in range(probed_index)] + [phase["".join([plane, bpm_name[probed_index], bpm_name[probed_index + 1 + i]])][1] / np.sqrt(1 + madTwiss.BETY[madTwiss.indx[bpm_name[probed_index]]] / madTwiss.BETY[madTwiss.indx[bpm_name[probed_index + 1 + i]]]) for i in range(probed_index)])
    err_diagonal = (TWOPI * err_diagonal) ** 2

    position = RANGE

        #---- assign field errors
    for j in range(RANGE):
        index = (CurrentIndex + j) % len(list_of_Ks)
        for k in range(len(list_of_Ks[index][0])):
            err_diagonal[k + position] = errorfile.dK1[list_of_Ks[index][0][k]] ** 2
        position += len(list_of_Ks[index][0])
        
        #---- assign sextupole transversal missalignments
    for j in range(RANGE):
        index = (CurrentIndex + j) % len(list_of_Ks)
        for k in range(len(list_of_Ks[index][1])):
            err_diagonal[k + position] = errorfile.dX[list_of_Ks[index][1][k]] ** 2
        position += len(list_of_Ks[index][1])
        
        #---- assign longitudinal missalignments
    for j in range(RANGE):
        index = (CurrentIndex + j) % len(list_of_Ks)
        for k in range(len(list_of_Ks[index][2])):
            err_diagonal[k + position] = errorfile.dS[list_of_Ks[index][2][k]] ** 2
        position += len(list_of_Ks[index][2])
        
        #---- assign BPM missalignments
    for j in range(RANGE):
        err_diagonal[j + position] = errorfile.dS[errorfile.indx[bpm_name[j]]] ** 2
        
    M = np.diag(err_diagonal)
               
    #---- calculate betas_from_phase for the three cases.
    #     and add matrix_row for the given combination
    matrix_rows_Beta = []
    matrix_rows_Alfa = []
    
    beta_alfa = []
    
    for n in BBA_combo:
        n0 = probed_index + n[0]
        n1 = probed_index + n[1]

        measured, TrowBeta, TrowAlfa = _beta_from_phase_BPM_BBA_with_systematicerrors(CurrentIndex,
                                                                                      bpm_name[n0], bpm_name[n1], bpm_name[probed_index], n0, n1, probed_index,
                                                                                      madTwiss, errorfile, phase, plane,
                                                                                      list_of_Ks, sizeOfMatrix, RANGE)
                            
        if measured.use_it:
            beta_alfa.append(measured)
            matrix_rows_Beta.append(TrowBeta)
            matrix_rows_Alfa.append(TrowAlfa)

    for n in BAB_combo:
        n0 = probed_index + n[0]
        n1 = probed_index + n[1]

        measured, TrowBeta, TrowAlfa = _beta_from_phase_BPM_BAB_with_systematicerrors(CurrentIndex,
                                                                                      bpm_name[n0], bpm_name[probed_index], bpm_name[n1], n0, probed_index, n1,
                                                                                      madTwiss, errorfile, phase, plane,
                                                                                      list_of_Ks, sizeOfMatrix, RANGE)
        if measured.use_it:
            beta_alfa.append(measured)
            matrix_rows_Alfa.append(TrowAlfa)
            matrix_rows_Beta.append(TrowBeta)
             
    for n in ABB_combo:
        n0 = probed_index + n[0]
        n1 = probed_index + n[1]
        
        measured, TrowBeta, TrowAlfa = _beta_from_phase_BPM_ABB_with_systematicerrors(CurrentIndex,
                                                                                      bpm_name[probed_index], bpm_name[n0], bpm_name[n1], probed_index, n0, n1,
                                                                                      madTwiss, errorfile, phase, plane,
                                                                                      list_of_Ks, sizeOfMatrix, RANGE)
        if measured.use_it:
            beta_alfa.append(measured)
            matrix_rows_Alfa.append(TrowAlfa)
            matrix_rows_Beta.append(TrowBeta)
    
    return np.matrix(matrix_rows_Alfa), np.matrix(matrix_rows_Beta), beta_alfa, bpm_name[probed_index], M


def _beta_from_phase_BPM_ABB_with_systematicerrors(I, bn1, bn2, bn3, bi1, bi2, bi3, madTwiss, errorfile, phase, plane, list_of_Ks, matrixSize, RANGE):
    '''
       Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is left of the other two BPMs (case ABB)
    Calculates also the corresponding column of the T-matrix. (awegsche June 2016)
    
    :Parameters:
        'bn1,bn2,bn3':string
            the names of the three BPMs
        'bi1,bi2,bi3':string
            the indices of the three BPMs (important to find the errors)
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'list_of_Ks':vector (of vectors)
            contains information about the errors the i-th entry in list_of_Ks yields the errors of all the lattice
            elements that come after the i-th BPM
        'matrixSize':int
            the size of the T-matrix
        'RANGE':int
            the range of BPMs
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'alf':float
            calculated error on beta function at probed BPM
        '0':float
            0
        'T':numpy matrix
            column of the T-matrix
        'T_Alf':numpy matrix
            column of the T-matrix for alfa
        'patternstring':string
            [For Debugging] string which represents the combination used
            For example AxxBxB
    '''
    I1 = madTwiss.indx[bn1]
    I2 = madTwiss.indx[bn2]
    I3 = madTwiss.indx[bn3]
    
    ph2pi12 = phase["".join([plane, bn1, bn2])][0]
    ph2pi13 = phase["".join([plane, bn1, bn3])][0]

    phmdl12 = phase["".join([plane, bn1, bn2])][2]
    phmdl13 = phase["".join([plane, bn1, bn3])][2]

    if plane == 'H':
        betmdl1 = madTwiss.BETX[I1]
        betmdl2 = madTwiss.BETX[I2]
        betmdl3 = madTwiss.BETX[I3]
        alfmdl1 = madTwiss.ALFX[I1]

    elif plane == 'V':
        betmdl1 = madTwiss.BETY[madTwiss.indx[bn1]]
        betmdl2 = madTwiss.BETY[madTwiss.indx[bn2]]
        betmdl3 = madTwiss.BETY[madTwiss.indx[bn3]]
        alfmdl1 = madTwiss.ALFY[madTwiss.indx[bn1]]

    if betmdl3 < 0 or betmdl2 < 0 or betmdl1 < 0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)
    if (bad_phase(phmdl12) or bad_phase(phmdl13) or bad_phase(phmdl12 - phmdl13)) or bad_phase(ph2pi12) or bad_phase(ph2pi13) or bad_phase(ph2pi12 - ph2pi13):
        return MeasuredValues(0, 0), [], []
   
    phmdl12 *= TWOPI
    phmdl13 *= TWOPI
    ph2pi12 *= TWOPI
    ph2pi13 *= TWOPI
      
    #--- Calculate beta
    cotphmdl12 = 1.0 / tan(phmdl12)
    cotphmdl13 = 1.0 / tan(phmdl13)
#     if is_small(cotphmdl12) or is_small(cotphmdl13):
#         return 0, 0, 0, [], [], "", False
    
    denom = (cotphmdl12 - cotphmdl13 + EPSILON) / betmdl1
    numer = 1.0 / tan(ph2pi12) - 1.0 / tan(ph2pi13)
    
#     if denom > BETA_THRESHOLD or denom < ZERO_THRESHOLD or numer > BETA_THRESHOLD or numer < ZERO_THRESHOLD:
#         return 0, 0, 0, [], [], "", False
    
    bet = numer / denom
    
    #--- Calculate alfa
    mdlterm = 1.0 / tan(phmdl12) + 1.0 / tan(phmdl13)
    denomalf = (mdlterm + 2.0 * alfmdl1) / betmdl1
    numeralf = 1.0 / tan(ph2pi12) + 1.0 / tan(ph2pi13)
    
    alf = 0.5 * (denomalf * bet - numeralf)
    
    T = [0] * matrixSize
    T_Alf = [0] * matrixSize
    
    phi_2 = madTwiss.MUX[I2] * TWOPI
    phi_3 = madTwiss.MUX[I3] * TWOPI
    
    frac = 1.0 / denom
    if plane == 'V':
        frac *= -1.0
        phi_2 = madTwiss.MUY[I2] * TWOPI
        phi_3 = madTwiss.MUY[I3] * TWOPI
     
    s_i2 = sin(phmdl12) ** 2
    s_i3 = sin(phmdl13) ** 2
          
    #--- Phase Advance
    phi_err1 = (1.0 / s_i2 - 1.0 / s_i3) / denom
    phi_err2 = (-1.0 / s_i2) / denom
    phi_err3 = 1.0 / s_i3 / denom
    
    T[bi1] = phi_err1
    T[bi2] = phi_err2
    T[bi3] = phi_err3
    
    T_Alf[bi1] = A_FACT * (-1.0 / s_i2 - 1.0 / s_i3 + phi_err1 * denomalf)
    T_Alf[bi2] = A_FACT * (1.0 / s_i2 + phi_err2 * denomalf)
    T_Alf[bi3] = A_FACT * (1.0 / s_i3 + phi_err3 * denomalf)
    
    # the first columns belong to the phase errors:
    K_offset = RANGE
        
    # then we jump to the first non-zero column in the K1 errors
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
       
    #--- Quad Fielderrors
    _assign_quaderrors(I, bi1, bi2,
                       errorfile, list_of_Ks,
                       denomalf, s_i2, T, T_Alf, phi_2, -frac, K_offset,
                       .5, .5)
       
    K_offset = _assign_quaderrors(I, bi1, bi3,
                                  errorfile, list_of_Ks,
                                  denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                  .5, .5)
            
    #--- Sext Transverse Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
        
    _assign_sext_errors(I, bi1, bi2,
                        errorfile, list_of_Ks,
                        denomalf, s_i2, T, T_Alf, phi_2, frac, K_offset,
                        -.5, .5)
    
    K_offset = _assign_sext_errors(I, bi1, bi3,
                                   errorfile, list_of_Ks,
                                   denomalf, s_i3, T, T_Alf, phi_3, -frac, K_offset,
                                   -.5, .5)

    #--- Quad Longitudinal Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])

    _assign_quadlongmissal(I, bi1, bi2,
                           errorfile, list_of_Ks,
                           denomalf, s_i2, T, T_Alf, phi_2, -frac, K_offset,
                           .5, .5)
       
    K_offset = _assign_quadlongmissal(I, bi1, bi3,
                                      errorfile, list_of_Ks,
                                      denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                      .5, .5)
    
    #--- BPM Missalignments
    # jump to end of RANGE
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
          
    errindx1 = errorfile.indx[bn1]
    errindx2 = errorfile.indx[bn2]
    errindx3 = errorfile.indx[bn3]
      
    if errorfile.dS[errindx1] != 0:
        numerphi = -1.0 / (betmdl2 * sin(phmdl12) ** 2) + 1.0 / (betmdl3 * sin(phmdl13) ** 2)
        T[K_offset + bi1] = numerphi / denom
          
    if errorfile.dS[errindx2] != 0:
        numerphi = 1.0 / (betmdl2 * sin(phmdl12) ** 2)
        T[K_offset + bi2] = numerphi / denom
       
    if errorfile.dS[errindx3] != 0:
        numerphi = -1.0 / (betmdl3 * sin(phmdl13) ** 2)
        T[K_offset + bi3] = numerphi / denom
         
    patternstr = ["x"] * RANGE
    patternstr[bi1] = "A"
    patternstr[bi2] = "B"
    patternstr[bi3] = "B"
 
    return MeasuredValues(alf, bet, "".join(patternstr), True), T, T_Alf


def _beta_from_phase_BPM_BAB_with_systematicerrors(I, bn1, bn2, bn3, bi1, bi2, bi3, madTwiss, errorfile, phase, plane, list_of_Ks, matrixSize, RANGE):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is betweeb the other two BPMs (case BAB)
    Calculates also the corresponding column of the T-matrix. (awegsche June 2016)
    
    :Parameters:
        'bn1,bn2,bn3':string
            the names of the three BPMs
        'bi1,bi2,bi3':string
            the indices of the three BPMs (important to find the errors)
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'matrixSize':int
            the size of the T-matrix
        'RANGE':int
            the range of BPMs
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'alf':float
            calculated error on beta function at probed BPM
        '0':float
            0
        'T':numpy matrix
            column of the T-matrix
        'T_Alf':numpy matrix
            column of the T-matrix for alfa
        'patternstring':string
            [For Debugging] string which represents the combination used
            For example AxxBxB
    '''
    
    I1 = madTwiss.indx[bn1]
    I2 = madTwiss.indx[bn2]
    I3 = madTwiss.indx[bn3]
    
    ph2pi21 = -phase["".join([plane, bn1, bn2])][0]
    ph2pi23 =  phase["".join([plane, bn2, bn3])][0]

    # Find the model transfer matrices for beta1
    phmdl21 = -phase["".join([plane, bn1, bn2])][2]
    phmdl23 =  phase["".join([plane, bn2, bn3])][2]
    if plane == 'H':
        betmdl1 = madTwiss.BETX[I1]
        betmdl2 = madTwiss.BETX[I2]
        betmdl3 = madTwiss.BETX[I3]
        alpmdl2 = madTwiss.ALFX[I2]
    elif plane == 'V':
        betmdl1 = madTwiss.BETY[madTwiss.indx[bn1]]
        betmdl2 = madTwiss.BETY[madTwiss.indx[bn2]]
        betmdl3 = madTwiss.BETY[madTwiss.indx[bn3]]
        alpmdl2 = madTwiss.ALFY[I2]

    if betmdl3 < 0 or betmdl2 < 0 or betmdl1 < 0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)
    if bad_phase(phmdl21) or bad_phase(phmdl23) or bad_phase(phmdl23 - phmdl21) or bad_phase(ph2pi21) or bad_phase(ph2pi23) or bad_phase(ph2pi21 - ph2pi23):
        return MeasuredValues(0, 0), [], []
    
    phmdl21 *= TWOPI
    phmdl23 *= TWOPI
    ph2pi21 *= TWOPI
    ph2pi23 *= TWOPI

    #--- Calculate beta
    cotphmdl21 = 1.0 / tan(phmdl21)
    cotphmdl23 = 1.0 / tan(phmdl23)
     
    denom = (cotphmdl21 - cotphmdl23 + EPSILON) / betmdl2
    numer = 1.0 / tan(ph2pi21) - 1.0 / tan(ph2pi23)
    bet = numer / denom
    
    #--- Calculate alfa
    mdlterm = 1.0 / tan(phmdl21) + 1.0 / tan(phmdl23)
    denomalf = (mdlterm + 2.0 * alpmdl2) / betmdl2
    numeralf = 1.0 / tan(ph2pi21) + 1.0 / tan(ph2pi23)
    alf = 0.5 * (denomalf * bet - numeralf)
      
    s_i1 = sin(phmdl21) ** 2
    s_i3 = sin(phmdl23) ** 2
     
    T = [0] * (matrixSize)
    T_Alf = [0] * matrixSize
    
    #--- Phase Advance
    phi_err1 = (-1.0 / s_i1) / denom
    phi_err2 = (1.0 / s_i1 - 1.0 / s_i3) / denom
    phi_err3 = 1.0 / s_i3 / denom
    
    T[bi1] = phi_err1
    T[bi2] = phi_err2
    T[bi3] = phi_err3
    
    T_Alf[bi1] = A_FACT * (1.0 / s_i1 + phi_err1 * denomalf)
    T_Alf[bi2] = A_FACT * (-1.0 / s_i1 - 1.0 / s_i3 + phi_err2 * denomalf)
    T_Alf[bi3] = A_FACT * (1.0 / s_i3 + phi_err3 * denomalf)
        
    phi_1 = madTwiss.MUX[I1] * TWOPI
    phi_3 = madTwiss.MUX[I3] * TWOPI
    
    frac = 1.0 / denom
    if plane == 'V':
        frac *= -1.0
        phi_1 = madTwiss.MUY[I1] * TWOPI
        phi_3 = madTwiss.MUY[I3] * TWOPI
    # the first columns belong to the phase errors:
    
    K_offset = RANGE
    
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
        
    #--- Quad Fielderrors
    K_offset = _assign_quaderrors(I, bi1, bi2,
                                  errorfile, list_of_Ks,
                                  denomalf, s_i1, T, T_Alf, phi_1, frac,
                                  K_offset,
                                  -.5, .5)
     
    K_offset = _assign_quaderrors(I, bi2, bi3,
                                  errorfile, list_of_Ks,
                                  denomalf, s_i3, T, T_Alf, phi_3, frac,
                                  K_offset,
                                  .5, .5)
     
    #--- Sext Transverse Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
        
    K_offset = _assign_sext_errors(I, bi1, bi2,
                                   errorfile, list_of_Ks,
                                   denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset,
                                   .5, .5)
                
    K_offset = _assign_sext_errors(I, bi2, bi3,
                                   errorfile, list_of_Ks,
                                   denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                   -.5, .5)
    #--- Quad Longitudinal Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])

    K_offset = _assign_quadlongmissal(I, bi1, bi2,
                                      errorfile, list_of_Ks,
                                      denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset,
                                      -.5, .5)
    
    K_offset = _assign_quadlongmissal(I, bi2, bi3,
                                      errorfile, list_of_Ks,
                                      denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                      .5, .5)
    #--- BPM Missalignments
    # jump to end of RANGE
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
          
    errindx1 = errorfile.indx[bn1]
    errindx2 = errorfile.indx[bn2]
    errindx3 = errorfile.indx[bn3]
      
    if errorfile.dS[errindx2] != 0:
        numerphi = -1.0 / (betmdl1 * sin(phmdl21) ** 2) + 1.0 / (betmdl3 * sin(phmdl23) ** 2)
        T[K_offset + bi2] = numerphi / denom
          
    if errorfile.dS[errindx1] != 0:
        numerphi = 1.0 / (betmdl1 * sin(phmdl21) ** 2)
        T[K_offset + bi1] = numerphi / denom
       
    if errorfile.dS[errindx3] != 0:
        numerphi = -1.0 / (betmdl3 * sin(phmdl23) ** 2)
        T[K_offset + bi3] = numerphi / denom
     
    patternstr = ["x"] * RANGE
    patternstr[bi1] = "B"
    patternstr[bi2] = "A"
    patternstr[bi3] = "B"
    
    return MeasuredValues(alf, bet, "".join(patternstr), True), T, T_Alf


def _beta_from_phase_BPM_BBA_with_systematicerrors(I, bn1, bn2, bn3, bi1, bi2, bi3, madTwiss, errorfile, phase, plane, list_of_Ks, matrixSize, RANGE):
    '''
        Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is right of the other two BPMs (case BBA)
    Calculates also the corresponding column of the T-matrix. (awegsche June 2016)
    
    :Parameters:
        'bn1,bn2,bn3':string
            the names of the three BPMs
        'bi1,bi2,bi3':string
            the indices of the three BPMs (important to find the errors)
        'madTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
        'matrixSize':int
            the size of the T-matrix
        'RANGE':int
            the range of BPMs
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'alf':float
            calculated error on beta function at probed BPM
        '0':float
            0
        'T':numpy matrix
            column of the T-matrix
        'T_Alf':numpy matrix
            column of the T-matrix for alfa
        'patternstring':string
            [For Debugging] string which represents the combination used
            For example AxxBxB
    '''
    
    I1 = madTwiss.indx[bn1]
    I2 = madTwiss.indx[bn2]
    I3 = madTwiss.indx[bn3]
   
    ph2pi32 = -phase["".join([plane, bn2, bn3])][0]
    ph2pi31 = -phase["".join([plane, bn1, bn3])][0]

    # Find the model transfer matrices for beta1
    phmdl32 = -phase["".join([plane, bn2, bn3])][2]
    phmdl31 = -phase["".join([plane, bn1, bn3])][2]

    if plane == 'H':
        betmdl1 = madTwiss.BETX[I1]
        betmdl2 = madTwiss.BETX[I2]
        betmdl3 = madTwiss.BETX[I3]
        alpmdl3 = madTwiss.ALFX[I3]
    elif plane == 'V':
        betmdl1 = madTwiss.BETY[I1]
        betmdl2 = madTwiss.BETY[I2]
        betmdl3 = madTwiss.BETY[I3]
        alpmdl3 = madTwiss.ALFY[I3]

    if betmdl3 < 0 or betmdl2 < 0 or betmdl1 < 0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)
    if bad_phase(phmdl32) or bad_phase(phmdl31) or bad_phase(phmdl31 - phmdl32) or bad_phase(ph2pi31) or bad_phase(ph2pi32) or bad_phase(ph2pi31 - ph2pi32):
        return MeasuredValues(0, 0), [], []
    
    phmdl31 *= TWOPI
    phmdl32 *= TWOPI
    ph2pi31 *= TWOPI
    ph2pi32 *= TWOPI
    
    cotphmdl32 = 1.0 / tan(phmdl32)
    cotphmdl31 = 1.0 / tan(phmdl31)
    
    #--- Calculate beta
    denom = (cotphmdl32 - cotphmdl31 + EPSILON) / betmdl3
    numer = 1.0 / tan(ph2pi32) - 1.0 / tan(ph2pi31)
    bet = numer / denom
  
    #--- Calculate alfa
    mdlterm = 1.0 / tan(phmdl32) + 1.0 / tan(phmdl31)
    denomalf = (mdlterm + 2.0 * alpmdl3) / betmdl3
    numeralf = 1.0 / tan(ph2pi32) + 1.0 / tan(ph2pi31)
    alf = 0.5 * (denomalf * bet - numeralf)
    
    s_i2 = sin(phmdl32) ** 2
    s_i1 = sin(phmdl31) ** 2
    
    T = [0] * matrixSize
    T_Alf = [0] * matrixSize
        
    #--- Phase Advance
    phi_err1 = 1.0 / s_i1 / denom
    phi_err2 = -1.0 / s_i2 / denom
    phi_err3 = -(1.0 / s_i1 - 1.0 / s_i2) / denom
   
    T[bi1] = phi_err1
    T[bi2] = phi_err2
    T[bi3] = phi_err3
    
    T_Alf[bi1] = A_FACT * (1.0 / s_i1 + phi_err1 * denomalf)
    T_Alf[bi2] = A_FACT * (1.0 / s_i2 + phi_err2 * denomalf)
    T_Alf[bi3] = A_FACT * (-1.0 / s_i1 - 1.0 / s_i2 + phi_err3 * denomalf)
                               
    phi_2 = madTwiss.MUX[I2] * TWOPI
    phi_1 = madTwiss.MUX[I1] * TWOPI
    
    frac = 1.0 / denom
    if plane == 'V':
        frac *= -1.0
        phi_2 = madTwiss.MUY[I2] * TWOPI
        phi_1 = madTwiss.MUY[I1] * TWOPI
           
    # the first columns belong to the phase errors:
    
    K_offset = RANGE
        
    #--- Quad Fielderrors
    # then we jump to the first non-zero column in the K1 errors
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    K_begin = K_offset
    #K_offset = RANGE
    
    # assign T matrix elements for phase errors between BPM 1 and 3
    for k in range(bi1, bi3):
        which_k = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[which_k][0])):
            idx_k = list_of_Ks[which_k][0][w]
            err_beta = -frac * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += -.5 * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += .5 * err_beta * denomalf
        K_offset += len(list_of_Ks[which_k][0])
    
    # go back because the second h_ij begins at 2
    # so we have to find the position of BPM2 in the matrix
    K_offset = K_begin
    for k in range(bi1, bi2):
        which_k = (k + I) % len(list_of_Ks)
        K_offset += len(list_of_Ks[which_k][0])
    # and assign the T matrix alements between BPM 1 and 3
    for k in range(bi2, bi3):
        which_k = (k + I) % len(list_of_Ks)
        
        for w in range(len(list_of_Ks[which_k][0])):
            idx_k = list_of_Ks[which_k][0][w]
            err_beta = frac * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_2) ** 2 / s_i2)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += -.5 * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_2) ** 2 / s_i2)
            T_Alf[K_offset + w] += .5 * err_beta * denomalf
        K_offset += len(list_of_Ks[which_k][0])
    
    #--- Sext Trasverse Missalignments
    # jump to the end of RANGE then to b1
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    K_begin = K_offset
    
    for k in range(bi1, bi3):
        which_k = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[which_k][1])):
            idx_k = list_of_Ks[which_k][1][w]
            err_beta = SEXT_FACT * frac * errorfile.K2L[idx_k] * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += .5 * SEXT_FACT * errorfile.K2L[idx_k] * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += .5 * err_beta * denomalf
        K_offset += len(list_of_Ks[which_k][1])
        
    K_offset = K_begin
    for k in range(bi1, bi2):
        which_k = (k + I) % len(list_of_Ks)
        K_offset += len(list_of_Ks[which_k][1])
    # and assign the T matrix alements between BPM 1 and 3
    for k in range(bi2, bi3):
        which_k = (k + I) % len(list_of_Ks)
        
        for w in range(len(list_of_Ks[which_k][1])):
            idx_k = list_of_Ks[which_k][1][w]
            err_beta = -SEXT_FACT * frac * errorfile.K2L[idx_k] * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_2) ** 2 / s_i2)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += .5 * SEXT_FACT * errorfile.K2L[idx_k] * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_2) ** 2 / s_i2)
            T_Alf[K_offset + w] += .5 * err_beta * denomalf
        K_offset += len(list_of_Ks[which_k][1])
    
    #--- Quad Longitudinal Missalignments
    # jump to the end of RANGE then to b1
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
    K_begin = K_offset
      
    for k in range(bi1, bi3):
        which_k = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[which_k][2])):
            idx_k = list_of_Ks[which_k][2][w]
            err_beta = -frac * errorfile.K1LEND[idx_k] * errorfile.BETEND[idx_k] * (sin(errorfile.MUEND[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += -.5 * errorfile.K1LEND[idx_k] * errorfile.BETEND[idx_k] * (sin(errorfile.MUEND[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += .5 * err_beta * denomalf
        K_offset += len(list_of_Ks[which_k][2])
        
    K_offset = K_begin
    for k in range(bi1, bi2):
        which_k = (k + I) % len(list_of_Ks)
        K_offset += len(list_of_Ks[which_k][2])
    # and assign the T matrix alements between BPM 1 and 3
    for k in range(bi2, bi3):
        which_k = (k + I) % len(list_of_Ks)
        
        for w in range(len(list_of_Ks[which_k][2])):
            idx_k = list_of_Ks[which_k][2][w]
            err_beta = frac * errorfile.K1LEND[idx_k] * errorfile.BETEND[idx_k] * (sin(errorfile.MUEND[idx_k] * TWOPI - phi_2) ** 2 / s_i2)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += -.5 * errorfile.K1LEND[idx_k] * errorfile.BETEND[idx_k] * (sin(errorfile.MUEND[idx_k] * TWOPI - phi_2) ** 2 / s_i2)
            T_Alf[K_offset + w] += .5 * err_beta * denomalf
        K_offset += len(list_of_Ks[which_k][2])
        
    #--- BPM Missalignments
    # jump to end of RANGE
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
          
    errindx1 = errorfile.indx[bn1]
    errindx2 = errorfile.indx[bn2]
    errindx3 = errorfile.indx[bn3]
      
    if errorfile.dS[errindx3] != 0:
        numerphi = -1.0 / (betmdl2 * sin(phmdl32) ** 2) + 1.0 / (betmdl1 * sin(phmdl31) ** 2)
        T[K_offset + bi3] = numerphi / denom
          
    if errorfile.dS[errindx2] != 0:
        numerphi = 1.0 / (betmdl2 * sin(phmdl32) ** 2)
        T[K_offset + bi2] = numerphi / denom
       
    if errorfile.dS[errindx1] != 0:
        numerphi = -1.0 / (betmdl1 * sin(phmdl31) ** 2)
        T[K_offset + bi1] = numerphi / denom
            
    patternstr = ["x"] * RANGE
    patternstr[bi1] = "B"
    patternstr[bi2] = "B"
    patternstr[bi3] = "A"

    return MeasuredValues(alf, bet, "".join(patternstr), True), T, T_Alf


def _assign_quaderrors(I, bi1, bi2, errorfile, list_of_Ks, denomalf, sinus_ij_squared, T_Bet, T_Alf, reference_phi, frac, K_offset, Alf_fact_1, Alf_fact_2):
    for k in range(bi1, bi2):
        whichK = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[whichK][0])):
            idx_k = list_of_Ks[whichK][0][w]
            err_beta = frac * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - reference_phi) ** 2 / sinus_ij_squared)
            T_Bet[K_offset + w] += err_beta
            T_Alf[K_offset + w] += Alf_fact_1 * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - reference_phi) ** 2 / sinus_ij_squared)
            T_Alf[K_offset + w] += Alf_fact_2 * err_beta * denomalf
        
        K_offset += len(list_of_Ks[whichK][0])
    
    return K_offset


def _assign_sext_errors(I, bi1, bi2, errorfile, list_of_Ks, denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset, Alf_fact_1, Alf_fact_2):
    for k in range(bi1, bi2):
        whichK = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[whichK][1])):
            idx_k = list_of_Ks[whichK][1][w]
            err_beta = -SEXT_FACT * frac * errorfile.K2L[idx_k] * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += Alf_fact_1 * SEXT_FACT * errorfile.K2L[idx_k] * errorfile.BET[idx_k] * (sin(errorfile.MU[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += Alf_fact_2 * err_beta * denomalf
        
        K_offset += len(list_of_Ks[whichK][1])
    
    return K_offset


def _assign_quadlongmissal(I, bi1, bi2, errorfile, list_of_Ks, denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset, Alf_fact_1, Alf_fact_2):
    for k in range(bi1, bi2):
        whichK = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[whichK][2])):
            idx_k = list_of_Ks[whichK][2][w]
            err_beta = frac * errorfile.K1LEND[idx_k] * errorfile.BETEND[idx_k] * (sin(errorfile.MUEND[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += Alf_fact_1 * errorfile.K1LEND[idx_k] * errorfile.BETEND[idx_k] * (sin(errorfile.MUEND[idx_k] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += Alf_fact_2 * err_beta * denomalf
        
        K_offset += len(list_of_Ks[whichK][2])
    
    return K_offset

#===================================================================================================
#--- ac-dipole stuff
#===================================================================================================


def _get_free_beta(modelfree, modelac, data, bpms, plane):  # to check "+"
    data2 = {}
    
    bpms = utils.bpm.model_intersect(bpms, modelfree)
    bpms = utils.bpm.model_intersect(bpms, modelac)
    for bpma in bpms:
        bpm = bpma[1].upper()
        beta, betsys, betstat, beterr = data[bpm][0:4]
        alfa, alfsys, alfstat, alferr = data[bpm][4:8]

        if plane == "H":
            betmf = modelfree.BETX[modelfree.indx[bpm]]
            betma = modelac.BETX[modelac.indx[bpm]]
            bb = betma / betmf
            alfmf = modelfree.ALFX[modelfree.indx[bpm]]
            alfma = modelac.ALFX[modelac.indx[bpm]]
            aa = alfma / alfmf
        else:
            betmf = modelfree.BETY[modelfree.indx[bpm]]
            betma = modelac.BETY[modelac.indx[bpm]]
            alfmf = modelfree.ALFY[modelfree.indx[bpm]]
            alfma = modelac.ALFY[modelac.indx[bpm]]
            bb = betma / betmf
            aa = alfma / alfmf
        data2[bpm] = beta / bb, betsys, betstat, beterr, alfa * aa, alfsys, alfstat, alferr, data[bpm][8], data[bpm][9], data[bpm][10]
    return data2, bpms


def _get_free_amp_beta(betai, rmsbb, bpms, inv_j, mad_ac, mad_twiss, plane):
    #
    # Why difference in betabeta calculation ??
    #
    #
    betas = {}

    if DEBUG:
        print "Calculating free beta from amplitude using model"

    for bpm in bpms:
        bpmm = bpm[1].upper()
        beta = betai[bpmm][0]

        if plane == "H":
            betmf = mad_twiss.BETX[mad_twiss.indx[bpmm]]
            betma = mad_ac.BETX[mad_ac.indx[bpmm]]
            bb = (betmf - betma) / betma
        else:
            betmf = mad_twiss.BETY[mad_twiss.indx[bpmm]]
            betma = mad_ac.BETY[mad_ac.indx[bpmm]]
            bb = (betmf - betma) / betma

        betas[bpmm] = [beta * (1.0 + bb), betai[bpmm][1], betai[bpmm][2]]

    return betas, rmsbb, bpms, inv_j

#=======================================================================================================================
#--- Helper / Debug Functions
#=======================================================================================================================


def create_errorfile(errordefspath, model, twiss_full, twiss_full_centre, commonbpms, plane, accel):
    '''
    Creates a file in Twiss format that contains the information about the expected errors for each element .
    
    There has to be an error definition file called "errordefs".
    
    There has to be a twiss model file called "twiss_full.dat" with all the elments in the lattice (also drift spaces) which contains the
    following columns:
    NAME, S, BETX, BETY, MUX, MUY, K1L, K2L
    '''
    
    if errordefspath is None:
        return None
    
    #bpms = []
    #for bpm in commonbpms:
    #    bpms.append(bpm[1])
    if accel == "JPARC":
        bpmre = re.compile("^MO[HV]\\.")
    elif accel == "PETRA":
        bpmre = re.compile("^BPM")
    else:
        bpmre = re.compile("^BPM.*B[12]$")
    

    print_("Create errorfile")
    print_("")
    
    # if something in loading / writing the files goes wrong, return None
    # which forces the script to fall back to 3bpm
    try:
        definitions = Python_Classes4MAD.metaclass.twiss(errordefspath)
        filename = "error_elements_" + plane + ".dat"
        errorfile = tfs_files.tfs_file_writer.TfsFileWriter(filename)
    except:
        print >> sys.stderr, "loading errorfile didnt work"
        print >> sys.stderr, "errordefspath = {0:s}".format(errordefspath)
        return None
     
    errorfile.add_column_names(     ["NAME",    "BET",  "BETEND",   "MU",   "MUEND",    "dK1",  "K1L",  "K1LEND",   "K2L",  "dX",   "dS", "DEBUG"])  #@IgnorePep8
    errorfile.add_column_datatypes( ["%s",      "%le",  "%le",      "%le",  "%le",      "%le",  "%le",  "%le",      "%le",  "%le",  "%le", "%s"])  #@IgnorePep8
    
    mainfield = definitions.RELATIVE == "MAINFIELD"
    
    regex_list = []
    for pattern in definitions.PATTERN:
        regex_list.append(re.compile(pattern))

    # OLD:
    for index_twissfull in range(len(twiss_full_centre.NAME)):
        BET = twiss_full_centre.BETX[index_twissfull]
        MU = twiss_full_centre.MUX[index_twissfull]

        if plane == 'V':
            BET = twiss_full_centre.BETY[index_twissfull]
            MU = twiss_full_centre.MUY[index_twissfull]
            
        BET_end = twiss_full.BETX[index_twissfull]
        MU_end = twiss_full.MUX[index_twissfull]
        BETminus1_end = twiss_full.BETX[index_twissfull - 1]
        MUminus1_end = twiss_full.MUX[index_twissfull - 1]

        if plane == 'V':
            BET_end = twiss_full.BETY[index_twissfull]
            MU_end = twiss_full.MUY[index_twissfull]
            BETminus1_end = twiss_full.BETY[index_twissfull - 1]
            MUminus1_end = twiss_full.MUY[index_twissfull - 1]

        found = False
        for index_defs in range(len(definitions.PATTERN)):
            regex = regex_list[index_defs]
            if regex.match(twiss_full.NAME[index_twissfull]):

                found = True
                isQuad = False
                MF = 1000
                if mainfield:
                    if definitions.MAINFIELD[index_defs] == "QUAD":
                        MF = twiss_full_centre.K1L[index_twissfull]
                        isQuad = True
                    elif definitions.MAINFIELD[index_defs] == "SEXT":
                        MF = twiss_full_centre.K2L[index_twissfull]
                    elif definitions.MAINFIELD[index_defs] == "DIPL":
                        MF = twiss_full_centre.K0L[index_twissfull]
                else:
                    MF = twiss_full.K1L[index_twissfull]
               
                errorfile.add_table_row([
                                        twiss_full.NAME[index_twissfull],
                                        BET,
                                        BET_end,
                                        MU,
                                        MU_end,
                                        definitions.dK1[index_defs] * MF,
                                        twiss_full_centre.K1L[index_twissfull],
                                        twiss_full.K1L[index_twissfull],
                                        twiss_full_centre.K2L[index_twissfull],
                                        definitions.dX[index_defs],
                                        definitions.dS[index_defs],
                                        "OK"
                                        ])
                if definitions.dS[index_defs] != 0 and isQuad:
                    errorfile.add_table_row([
                                            twiss_full.NAME[index_twissfull - 1],
                                            0,     # BET
                                            BETminus1_end,
                                            0,     # MU
                                            MUminus1_end,
                                            0,     # dK1,
                                            0,     # K1L,
                                            - twiss_full.K1L[index_twissfull],  # same index
                                            0,     # K2L
                                            0,     # dX
                                            definitions.dS[index_defs],  # here no -1 because the same dS applies
                                            "OK"])

        if not found:  # if element doesn't have any error add it nevertheless if it is a BPM
            if bpmre.match(twiss_full_centre.NAME[index_twissfull]):
                index_model = model.indx[twiss_full.NAME[index_twissfull]]
                errorfile.add_table_row([
                                        model.NAME[index_model],
                                        BET,
                                        0,  # BETEND
                                        MU,
                                        0,  # MUEND
                                        0,  # dK1
                                        0,  # K1L
                                        0,  # K1LEND
                                        0,  # K2L
                                        0,  # dX
                                        0,  # dS
                                        "NOT_FOUND"])
    errorfile.write_to_file(True)

    print_("DONE creating errofile.")

    return Python_Classes4MAD.metaclass.twiss(filename)


def printMatrix(debugfile, M, name):
    debugfile.write("begin Matrix " + name + "\n" + str(M.shape[0]) + " " + str(M.shape[1]) + "\n")

    np.savetxt(debugfile, M, fmt="%18.10e")
    debugfile.write("\nend\n")


def bad_phase(phi):
    modphi = phi % BADPHASE
    return (modphi < MOD_POINTFIVE_LOWER or modphi > MOD_POINTFIVE_UPPER)


def JPARC_intersect(plane, getllm_d, commonbpms):
    if getllm_d.accel == 'JPARC': # if JPARC, H-BPMs and V-BPMs are at different positions, so take the right ones
        bpm_regex = re.compile("^MOH")
        if plane == 'V':
            bpm_regex = re.compile("^MOV")

        commonbpms = [bpm for bpm in commonbpms if bpm_regex.match(bpm[1])]


    return commonbpms


def is_small(x):
    return abs(x) < ZERO_THRESHOLD


def print_box(string):
    print "=" + " " * BOXINDENT + string + " " * (BOXLENGTH - 3 - BOXINDENT - len(string)) + "="
    
    
def print_(string, prefix=" "):
    print " " * (BOXINDENT + 1) + prefix + " " + string
    
    
def print_box_edge():
    print "= " * (BOXLENGTH / 2)

