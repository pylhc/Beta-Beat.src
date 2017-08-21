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
import pandas as pd

from scipy.linalg import circulant
import Python_Classes4MAD.metaclass
import Utilities.bpm
from GetLLMError import GetLLMError
import compensate_ac_effect
import os
import re
import multiprocessing
import time
from constants import PI, TWOPI
from numpy.linalg.linalg import LinAlgError
from model.accelerators.accelerator import AccExcitationMode
from Utilities import tfs_pandas
__version__ = "2017.8.1"

DEBUG = sys.flags.debug  # True with python option -d! ("python -d GetLLM.py...") (vimaier)
PRINTTIMES = False

if False:
    from Utilities.progressbar import startProgress, progress, endProgress
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
PHASE_THRESHOLD         = 1.0e-2                      #@IgnorePep8
COT_THRESHOLD           = 16.0
MOD_POINTFIVE_LOWER     = PHASE_THRESHOLD           #@IgnorePep8
MOD_POINTFIVE_UPPER     = (BADPHASE - PHASE_THRESHOLD)    #@IgnorePep8
RCOND                   = 1.0e-14                    #@IgnorePep8

BOXLENGTH               = 50                        #@IgnorePep8
BOXINDENT               =  4                        #@IgnorePep8
CALCULATE_BETA_HOR = True
CALCULATE_BETA_VER = True

# ================ Errorfile indices:
BET_INDEX       = 1
BETEND_INDEX    = 2
MU_INDEX        = 3
MUEND_INDEX     = 4
DK1_INDEX       = 5
K1L_INDEX       = 6
K1LEND_INDEX    = 7
K2L_INDEX       = 8
DX_INDEX        = 9
DS_INDEX        = 10

# ================ Multiprocessing indices:
BETI_MP     = 0
BETSTAT_MP  = 1
BETSYST_MP  = 2
BETERR_MP   = 3
ALFI_MP     = 4
ALFSTAT_MP  = 5
ALFSYS_MP   = 6
ALFERR_MP   = 7
CORR_MP     = 8
DELBETA_MP  = 9
NCOMB_MP    = 10
# ================ Errors method
METH_IND        = -1
METH_3BPM       = 0            
METH_A_NBPM     = 1
METH_MC_NBPM    = 2

ID_TO_METHOD = {
        METH_3BPM:"3BPM method",
        METH_A_NBPM:"Analytical N-BPM method"}

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


ID_INVALID = 0
IDQUAD = 1
IDSEXT = 2
IDBPM = 3
IDDIPL = 4

MAINFIELD = {
        IDQUAD:"K1L",
        IDSEXT:"K2L",
        IDDIPL:"K0L"}

def gettype(_type):
    if _type == "QUAD":
        return IDQUAD
    elif _type == "SEXT":
        return IDSEXT
    elif _type == "BPM":
        return IDBPM
    elif _type == "DIPL":
        return IDDIPL
    print "-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<- INVALID"
    return ID_INVALID

class UncertaintyDefinition:
    def __init__(self, _pattern, _dk1=0, _ds=0, _dx=0, _type="QUAD"):
        self.pattern = _pattern
        self.dX = _dx
        self.dS = _ds
        self.dK1 = _dk1
        self.type = gettype(_type)
        
        def settype(self, _type):
         self.type = gettype(_type)
           
   
class UncertaintyDefinitionRE:
    def __init__(self, _pattern, _dk1=0, _ds=0, _dx=0, _type="QUAD"):
        self.pattern = re.compile(_pattern)
        self.dX = _dx
        self.dS = _ds
        self.dK1 = _dk1
        self.type = gettype(_type)
        
    def settype(self, _type):
        self.type = gettype(_type)
        
    def match(self, string):
        return self.pattern.match(string)
             
        
class UncertaintyInformation:
    def __init__(self, _name, _bet, _betend, _mu, _muend, _dk1, _k1l, _k1lend, _k2l, _dx, _ds, _debug):
        self.name = _name
        self.bet = _bet
        self.betend = _betend
        self.mu = _mu
        self.muend = _muend
        self.dk1 = _dk1
        self.k1l = _k1l
        self.k1lend = _k1lend
        self.k2l = _k2l
        self.dx = _dx
        self.ds = _ds
        self.debug = _debug


class ErrorFile:
    def __init__(self):
       
        self.indx = {}
        self.elements = None
        self.size = 0
        self.names = []
        
    def add(self, uni):
        if self.elements is None:
            self.elements =  [0.0, uni.bet, uni.betend, uni.mu, uni.muend, uni.dk1, uni.k1l, uni.k1lend, uni.k2l, uni.dx, uni.ds],
        else:
            self.elements = np.concatenate(
                (self.elements,
                 ([0.0, uni.bet, uni.betend, uni.mu, uni.muend, uni.dk1, uni.k1l, uni.k1lend, uni.k2l, uni.dx, uni.ds],))
                                           )
        self.size = len(self.elements)
        self.indx[uni.name] = self.size - 1
        self.names.append(uni.name)
#         print "{:12s}, at index {:d} {:6.1f} {:6.1f} {:5.1f} {:5.1f} {:12.3e} {:12.3e} {:9.4f} {:9.4f} {:9.4f} {:9.4f}".format(
#             uni.name, self.size, uni.bet, uni.betend, uni.mu, uni.muend, uni.dk1, uni.k1l, uni.k1lend, uni.k2l, uni.dx, uni.ds
#                                     )
        
        
class Uncertainties:  # error definition file
    
    def __init__(self):
        self.version = "1"
        self.keys = []
        self.regex = []
        self.properties = {}
        
    
    def open(self, filename):
        with open(filename, "r") as F:
        
            line = F.readline()
            
            if line.strip() == "version 2":
                line = F.readline()
                for line in F:
                    
                    line = re.split("//", line)[0].strip()  # separating code from comments
                    if len(line) == 0:
                        continue
                    match = re.search(r'prop (\w+)\s+=\s(\w[\w\s]*)', line)
                    if match is not None:
                        print "adding {} = {} to properties".format(match.group(1), match.group(2))
                        self.properties[match.group(1)] = match.group(2).strip()
                    else:
                        words = re.split('\s', line)
                        
                        if words[0].startswith("re:"):
                            ud = UncertaintyDefinitionRE(words[0].split(":")[1])
                            
                            for word in words:
                                kv = word.split("=")
                                if kv[0] == "dK1":
                                    ud.dK1 = float(kv[1])
                                elif kv[0] == "dS":
                                    ud.dS = float(kv[1])
                                elif kv[0] == "dX":
                                    ud.dX = float(kv[1])
                                elif kv[0] == "Type":
                                    ud.settype(kv[1])
                            self.regex.append(ud)
                        else:
                            ud = UncertaintyDefinition(words[0])
                            
                            for word in words:
                                kv = word.split("=")
                                if kv[0] == "dK1":
                                    ud.dK1 = float(kv[1])
                                elif kv[0] == "dS":
                                    ud.dS = float(kv[1])
                                elif kv[0] == "dX":
                                    ud.dX = float(kv[1])
                                elif kv[0] == "Type":
                                    ud.settype(kv[1])
                            self.keys.append(ud)
            
                return True
            else:
                
                try:
                    definitions = Python_Classes4MAD.metaclass.twiss(filename)
                    
                except:
                    print >> sys.stderr, "loading errorfile didn't work"
                    print >> sys.stderr, "errordefspath = {0:s}".format(filename)
                    return False
                
                print_("error definitions file version 1")
                self.properties["RELATIVE"] = definitions.RELATIVE
                self.properties["RADIUS"] = definitions.RADIUS
                 
                for index in range(len(definitions.PATTERN)):
                    pattern = definitions.PATTERN[index]
                    self.regex.append(UncertaintyDefinitionRE(
                        pattern,
                        definitions.dK1[index],
                        definitions.dS[index],
                        definitions.dX[index],
                        definitions.MAINFIELD[index]))
                return True
    def create_errorfile(self, twiss_full, twiss_full_centre):
        starttime = time.time()
        twiss_full.loc[:]["MUX_END"] = np.roll(twiss_full.loc[:]["MUX"], 1)
        twiss_full.loc[:]["MUY_END"] = np.roll(twiss_full.loc[:]["MUY"], 1)
        twiss_full.loc[:]["BETX_END"] = np.roll(twiss_full.loc[:]["BETX"], 1)
        twiss_full.loc[:]["BETY_END"] = np.roll(twiss_full.loc[:]["BETY"], 1)
        twiss_full.loc[:]["UNC"] = False
        twiss_full.loc[:]["dK1"] = 0
        twiss_full.loc[:]["dS"] = 0
        twiss_full.loc[:]["dX"] = 0
        twiss_full.loc[:]["BPMdS"] = 0
        for reg in self.regex:
            reg_mask = twiss_full.index.str.match(reg.pattern)
            twiss_full.loc[reg_mask, "dK1"] = (reg.dK1 * twiss_full.loc[reg_mask, "K1L"]) **2 # TODO change K1L --> mainfield if necessary
            twiss_full.loc[reg_mask, "dX"] = reg.dX**2
            if reg.type == IDBPM:
                twiss_full.loc[reg_mask, "BPMdS"] = reg.dS**2
            else:
                twiss_full.loc[reg_mask, "dS"] = reg.dS**2
            twiss_full.loc[reg_mask, "UNC"] = True
        twiss_full.loc[:]["dS"] -= np.roll(twiss_full.loc[:]["dS"], 1)
        twiss_full.loc[:]["dK1_END"] = -np.roll(twiss_full.loc[:]["dK1"], 1)
        twiss_full.loc[:]["UNC"] |= np.roll(twiss_full.loc[:]["UNC"], 1)
        print "==============================", time.time() - starttime
        print_("DONE creating errofile.")
        
        tfs_pandas.write_tfs(twiss_full, {}, "dump_twiss_full")
        
        return twiss_full[twiss_full["UNC"] == True]

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
    if error_method == METH_3BPM:
        tfs_file.add_float_descriptor("NumberOfBPMs", 3)
        tfs_file.add_float_descriptor("RangeOfBPMs", 5)
    else:
        tfs_file.add_float_descriptor("NumberOfBPMs", number_of_bpms)
        tfs_file.add_float_descriptor("RangeOfBPMs", range_of_bpms)
    tfs_file.add_string_descriptor("ErrorsFrom", ID_TO_METHOD[error_method])
    tfs_file.add_float_descriptor("PhaseTheshold", PHASE_THRESHOLD)
    tfs_file.add_float_descriptor("RCond", RCOND)
    
    tfs_file.add_column_names(["NAME", "S",
                               "BET" + _plane_char,
                               "STATBET" + _plane_char,
                               "SYSBET" + _plane_char,
                               "ERRBET" + _plane_char,
                               "ALF" + _plane_char,
                               "STATALF" + _plane_char,
                               "SYSALF" + _plane_char,
                               "ERRALF" + _plane_char,
                               "CORR",
                               "BET" + _plane_char + "MDL",
                               "BBEAT",
                               "NCOMB"])
    tfs_file.add_column_datatypes(["%s", "%le",
                               "%le" ,
                               "%le" ,
                               "%le" ,
                               "%le" ,
                               "%le" ,
                               "%le",
                               "%le",
                               "%le",
                               "%le",
                               "%le",
                               "%le",
                               "%le"])
    for row in data:
        tfs_file.add_table_row(row)


def calculate_beta_from_phase(getllm_d, twiss_d, tune_d, phase_d,
                              model, model_driven, elements, elements_centre,
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
    
    commonbpms_x = twiss_d.zero_dpp_commonbpms_x
    commonbpms_y = twiss_d.zero_dpp_commonbpms_y
    
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
        debugfilename = files_dict['getbetax.out'].s_output_path + "/getbetax.debug"
        print debugfilename
        debugfile = open(debugfilename, "w+")
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
    
    # check whether analytical N-BPM method should be used    
    # if yes, uncertainty estimates will be distributed to the elements
    # do this only once for both planes
    error_method = METH_IND
    unc_elements = None

    if getllm_d.use_only_three_bpms_for_beta_from_phase:
        error_method = METH_3BPM
    elif getllm_d.errordefspath is not None:
        unc = Uncertainties()
        unc.open(getllm_d.errordefspath)
        unc_elements = unc.create_errorfile(elements, elements_centre)
        error_method = METH_A_NBPM
    else:
        error_method = METH_3BPM  # fall back to three BPM method. MC is not supported in this version
    
    #--- =========== HORIZONTAL
    if CALCULATE_BETA_HOR:
        if twiss_d.has_zero_dpp_x():
            
            print ""
            print_("Calculate free beta from phase for plane " + _plane_char + " (_free.out)", ">")
            if DEBUG:
                debugfile = open(files_dict['getbetax_free.out'].s_output_path + "/getbetax_free.debug", "w+")

            dataf, rmsbbxf, bpmsf, error_method_x = beta_from_phase(model, unc_elements, elements_centre,
                                                                  twiss_d.zero_dpp_x, commonbpms_x, phase_d.phase_advances_free_x, 'H',
                                                                  getllm_d, debugfile, error_method, tune_d.q1f, tune_d.q1mdl)
            beta_d.x_phase_f = {}
            _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, model, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.x_phase_f,
                               dataf, rmsbbxf, error_method_x, bpmsf,
                               files_dict['getbetax_free.out'], model.BETX, model.ALFX, model.MUX, _plane_char)
            
            
            if getllm_d.accelerator.excitation is not AccExcitationMode.FREE: 
                print ""
                print_("Calculate beta from phase for plane " + _plane_char, ">")
                data, rmsbbx, bpms, error_method_x = beta_from_phase(model_driven, unc_elements, elements_centre,
                                                                   twiss_d.zero_dpp_x, commonbpms_x, phase_d.phase_advances_x, 'H',
                                                                   getllm_d, debugfile, error_method, tune_d.q1, tune_d.q1mdl)
                beta_d.x_phase = {}
                beta_d.x_phase['DPP'] = 0
                tfs_file = files_dict['getbetax.out']
                
                _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1, tune_d.q2, model_driven, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.x_phase,
                                   data, rmsbbx, error_method_x, bpms,
                                   tfs_file, model_driven.BETX, model_driven.ALFX, model_driven.MUX, _plane_char)

                print_("Skip free2 calculation")
#                print_("Calculate beta from phase for plane " + _plane_char + " with AC dipole (_free2.out)", ">")
#                dataf2, bpmsf2 = _get_free_beta(model, model_driven, data, bpms, 'H')
##                tfs_file = files_dict['getbetax_free2.out']
#                betaxf2 = {}
#                rmsbbxf2 = 0
#                tfs_pandas.write_tfs(dataf2, headers, os.path.join( getllm_d.outputpath, "NEWgetbetax_free2.out"))
##                _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, model, getllm_d.number_of_bpms, getllm_d.range_of_bpms, betaxf2,
##                                   dataf2, rmsbbxf2, error_method, bpmsf2,
##                                   tfs_file, model.BETX, model.ALFX, model.MUX, _plane_char)
                
    #--- =========== VERTICAL
    if CALCULATE_BETA_VER:
        _plane_char = "Y"
        if twiss_d.has_zero_dpp_y():
            print ""
            print_("Calculate free beta from phase for plane " + _plane_char + " (_free.out)", ">")
    
            if DEBUG:
                debugfile = open(files_dict['getbetay_free.out'].s_output_path + "/getbetay_free.debug", "w+")
    

            dataf, rmsbby, bpms, error_method_y = beta_from_phase(model, unc_elements, elements_centre,
                                                               twiss_d.zero_dpp_y, commonbpms_y, phase_d.phase_advances_free_y, 'V',
                                                               getllm_d, debugfile, error_method, tune_d.q2f, tune_d.q2mdl)
            beta_d.y_phase = {}
            beta_d.y_phase['DPP'] = 0
            tfs_file = files_dict['getbetay_free.out']
            
            _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, model, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.y_phase,
                               dataf, rmsbby, error_method_y, bpms,
                               tfs_file, model.BETY, model.ALFY, model.MUY, _plane_char)
            
            #-- ac to free beta
            if getllm_d.accelerator.excitation is not AccExcitationMode.FREE:
                #-- from eq
                print ""
                print_("Calculate beta from phase for plane " + _plane_char, ">")
    
                data, rmsbbyf, bpmsf, error_method_y = beta_from_phase(model_driven, unc_elements, elements_centre,
                                                                      twiss_d.zero_dpp_y, commonbpms_y, phase_d.phase_advances_y, 'V',
                                                                      getllm_d, debugfile, error_method, tune_d.q2, tune_d.q2mdl)
   
             
                tfs_file = files_dict['getbetay.out']
                beta_d.y_phase_f = {}
                _write_getbeta_out(twiss_d.zero_dpp_y, tune_d.q1, tune_d.q2, model_driven, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.y_phase_f,
                                   dataf, rmsbbyf, error_method_y, bpmsf,
                                   tfs_file, model_driven.BETY, model_driven.ALFY, model_driven.MUY, _plane_char)

#                #-- from the model
                print_("Skip free2 calculation")
#                print_("Calculate beta from phase for plane " + _plane_char + " with AC dipole (_free2.out)", ">")
#                [datayf2, bpmsf2] = _get_free_beta(model, model_driven, data, bpms, 'V')
##                tfs_file = files_dict['getbetay_free2.out']
#                betayf2 = {}
#                rmsbbyf2 = 0
#                tfs_pandas.write_tfs(datayf2, headers, os.path.join( getllm_d.outputpath, "NEWgetbetay_free2.out"))
##                _write_getbeta_out(twiss_d.zero_dpp_x, tune_d.q1f, tune_d.q2f, model, getllm_d.number_of_bpms, getllm_d.range_of_bpms, betayf2,
##                                   datayf2, rmsbbyf2, error_method, bpmsf2,
##                                   tfs_file, model.BETY, model.ALFY, model.MUY, _plane_char)
#                
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
    
    commonbpms_x = twiss_d.zero_dpp_commonbpms_x
    commonbpms_y = twiss_d.zero_dpp_commonbpms_y


    #---- H plane
    if twiss_d.has_zero_dpp_x():
        [beta_d.x_amp, rmsbbx, bpms, inv_jx] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_x, 'H', commonbpms_x)
        beta_d.x_amp['DPP'] = 0
        #-- Rescaling
        beta_d.x_ratio = 0
        skipped_bpmx = []
        arcbpms= mad_twiss.index[getllm_d.accelerator.get_arc_bpms_mask(mad_twiss.index)]
        for bpm in arcbpms:
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
                betaxf, rmsbbxf, bpmsf = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, twiss_d.zero_dpp_x, tune_d.q1, tune_d.q1f, phase_d.acphasex_ac2bpmac, 'H', getllm_d)
                #-- Rescaling
                beta_d.x_ratio_f = 0
                skipped_bpmxf = []
                arcbpms = Utilities.bpm.filterbpm(bpmsf)
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
        [beta_d.y_amp, rmsbby, bpms, inv_jy] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_y, 'V', commonbpms_y)
        beta_d.y_amp['DPP'] = 0
        #-- Rescaling
        beta_d.y_ratio = 0
        skipped_bpmy = []
        arcbpms = Utilities.bpm.filterbpm(bpms)
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
                betayf, rmsbbyf, bpmsf = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, twiss_d.zero_dpp_y, tune_d.q2, tune_d.q2f, phase_d.acphasey_ac2bpmac, 'V', getllm_d)  # Rescaling
                beta_d.y_ratio_f = 0
                skipped_bpmyf = []
                arcbpms = Utilities.bpm.filterbpm(bpmsf)
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


def beta_from_phase(madTwiss, madElements, madElementsCentre, ListOfFiles, commonbpms, phase, plane, getllm_d, debugfile, errors_method, tune, mdltune):
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
    
    rmsbb = 0.0
    
    #---- Error definitions given and we decided to use  => use analytical formulas to calculate the
    # systematic errors
    if errors_method == METH_A_NBPM:
        
        
        rmsbb, errors_method, data = scan_all_BPMs_withsystematicerrors(madTwiss, madElements,
                                                                        phase, plane, getllm_d, commonbpms, debugfile, errors_method,
                                                                        tune, mdltune)
        return data, rmsbb, commonbpms, errors_method
    #---- use the simulations
    else:
        
        rmsbb, errors_method, data = scan_all_BPMs_sim_3bpm(madTwiss,
                                                            phase, plane, getllm_d, commonbpms, debugfile, errors_method,
                                                            tune, mdltune)

    return data, rmsbb, commonbpms, errors_method


def beta_from_amplitude(mad_twiss, list_of_files, plane, commonbpms):

    beta = {}
    root2j = []
    sum_a = 0.0
    amp = []
    amp2 = []
    kick2 = []
    for i in range(0, len(commonbpms)):  # this loop have become complicated after modifications... anybody simplify?
        bn1 = str.upper(commonbpms.index[i])
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
        bn1 = str.upper(commonbpms.index[i])
        location = commonbpms.iloc[i]
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
def scan_all_BPMs_sim_3bpm(madTwiss, phase, plane, getllm_d, commonbpms, debugfile, errors_method, tune, mdltune):
    number_commonbpms = commonbpms.shape[0]
    plane_bet = "BETX" if plane == "H" else "BETY"
    plane_alf = "ALFX" if plane == "H" else "ALFY"
#    np.set_printoptions(edgeitems=7, precision=3, linewidth=150)
    print ID_TO_METHOD[errors_method]
    if errors_method == METH_3BPM:
        madTwiss_intersected = madTwiss.loc[commonbpms.index]
        starttime = time.time()
        
        # setup the used variables
        # tilt phase advances in order to have the phase advances in a neighbourhood
        tilted_meas = tilt_slice_matrix(phase["MEAS"].as_matrix(), 2, 5, tune)
        tilted_model = tilt_slice_matrix(phase["MODEL"].as_matrix(), 2, 5, mdltune)
        print tilted_meas
        betmdl = madTwiss_intersected.loc[:][plane_bet]
        alfmdl = madTwiss_intersected.loc[:][plane_alf]

        
        # calculate cotangens of all the phase advances in the neighbourhood
        cot_phase_meas = 1.0 / tan(tilted_meas * TWOPI)
        cot_phase_model = 1.0 / tan(tilted_model * TWOPI)
       
        # calculate enumerators and denominators for far more cases than needed
        # shift1 are the cases BBA, ABB, AxBB, AxxBB etc. (the used BPMs are adjacent)
        # shift2 are the cases where the used BPMs are separated by one. only BAB is used for  3-BPM
        cot_phase_meas_shift1 = cot_phase_meas - np.roll(cot_phase_meas, -1, axis=0)
        cot_phase_model_shift1 = cot_phase_model - np.roll(cot_phase_model, -1, axis=0) + 1.0e-16
        cot_phase_meas_shift2 = cot_phase_meas - np.roll(cot_phase_meas, -2, axis=0)
        cot_phase_model_shift2 = cot_phase_model - np.roll(cot_phase_model, -2, axis=0)+ 1.0e-16
        
        # calculate the sum of the fractions
        bet_frac = (cot_phase_meas_shift1[0]/cot_phase_model_shift1[0] +
                    cot_phase_meas_shift1[3]/cot_phase_model_shift1[3] +
                    cot_phase_meas_shift2[1]/cot_phase_model_shift2[1]) / 3.0
        
        # multiply the fractions by betmdl and calculate the arithmetic mean
        beti = bet_frac * betmdl
        
        alfi = (bet_frac * ( cot_phase_model[1] + cot_phase_model[3] + 2.0 * alfmdl)
                        - ( cot_phase_meas[1] + cot_phase_meas[3])) / 2.0
        betstd=np.zeros(number_commonbpms)
        beterr=np.zeros(number_commonbpms)
        alfstd=np.zeros(number_commonbpms)
        alferr=np.zeros(number_commonbpms)

        
        print "===================================", time.time() - starttime
        starttime = time.time()
        
        
#        data_ = pd.DataFrame(columns=commonbpms.index,
#                             index=["NAME", "S", plane_bet, "BETSTD", "BETSYS", "BETERR", plane_alf, "ALFSTD", "ALFSYS", "ALFERR", "CORR", plane_bet + "MDL", "BBEAT", "NUM"],
#                             data=[
#                                     commonbpms.index, madTwiss_intersected.loc[:]["S"], beti, betstd, beterr, np.sqrt(beterr ** 2 + betstd ** 2),
#                                     alfi, alfstd, alferr, np.sqrt(alferr ** 2 + alfstd ** 2),
#                                     alfi,
#                                     betmdl,
#                                     bet_frac - 1.0,
#                                     alfi]
#                             )
        
        print "===================================", time.time() - starttime
        print_("Errors from " + ID_TO_METHOD[errors_method])


        return 0, errors_method, np.transpose([
                                     commonbpms.index, madTwiss_intersected.loc[:]["S"],
                                     beti, betstd, beterr, np.sqrt(beterr ** 2 + betstd ** 2),
                                     alfi, alfstd, alferr, np.sqrt(alferr ** 2 + alfstd ** 2),
                                     alfi,
                                     betmdl,
                                     bet_frac - 1.0,
                                     alfi])

    raise GetLLMError("Monte Carlo N-BPM is not implemented. Please contact awegsche")
    

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
        bn1 = str.upper(commonbpms.index[(probed_index + i - 2) % len(commonbpms)])
        bn2 = str.upper(commonbpms.index[(probed_index + i - 1) % len(commonbpms)])
        bn3 = str.upper(commonbpms.index[(probed_index + i) % len(commonbpms)])
        bn4 = str.upper(commonbpms.index[(probed_index + i + 1) % len(commonbpms)])
        bn5 = str.upper(commonbpms.index[(probed_index + i + 2) % len(commonbpms)])
        candidates = []
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, _ = _beta_from_phase_BPM_right(bn1, bn2, bn3, madTwiss, phase, plane, 0, 0, 0, True)
        candidates.append([tbetstd, tbet, talfstd, talf])
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, _ = _beta_from_phase_BPM_mid(bn2, bn3, bn4, madTwiss, phase, plane, 0, 0, 0, True)
        candidates.append([tbetstd, tbet, talfstd, talf])
        tbet, tbetstd, talf, talfstd, mdlerr, t1, t2, _ = _beta_from_phase_BPM_left(bn3, bn4, bn5, madTwiss, phase, plane, 0, 0, 0, True)
        candidates.append([tbetstd, tbet, talfstd, talf])
        return candidates, bn3, []
    
    # TODO remove the code hereafter or implement Monte Cralo N-BPM for accelerator classes

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
      
        t_matrix_row = [0] * (RANGE-1)
        t_matrix_row[n[0] - 1 - probed_index] = t1
        t_matrix_row[n[1] - 1 - probed_index] = t2
        patternstr = ["x"] * RANGE
        patternstr[n[0]] = "B"
        patternstr[n[1]] = "B"
        patternstr[probed_index] = "A"
        candidates.append([tbetstd, tbet, talfstd, talf, mdlerr, bpm_name[n[0]], bpm_name[n[1]], t_matrix_row, "".join(patternstr)])

    sort_cand = sorted(candidates, key=lambda x: x[4])
    return [sort_cand[i] for i in range(number_of_bpms)], bpm_name[probed_index], M
#     return candidates, bpm_name[probed_index], M


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
    ph2pi12=phase.loc["MEAS", bn1, bn2]
    ph2pi13=phase.loc["MEAS", bn1, bn3]
    betstd = 0
    alfstd = 0

    # Find the model transfer matrices for beta1
    phmdl12 = phase.loc["MODEL", bn1, bn2]
    phmdl13= phase.loc["MODEL", bn1, bn3]
    if plane=='H':
        betmdl1=madTwiss.loc[bn1, "BETX"]
        betmdl2=madTwiss.loc[bn2, "BETX"]
        betmdl3=madTwiss.loc[bn3, "BETX"]
        alpmdl1=madTwiss.loc[bn1, "ALFX"]
    elif plane=='V':
        betmdl1=madTwiss.loc[bn1, "BETY"]
        betmdl2=madTwiss.loc[bn2, "BETY"]
        betmdl3=madTwiss.loc[bn3, "BETY"]
        alpmdl1=madTwiss.loc[bn1, "ALFY"]
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
    bet=numer/denom

#    betstd=        (2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
#    betstd=betstd+(2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
#    betstd=math.sqrt(betstd)/abs(denom)

    mdlerr=        (2*np.pi*0.001/sin(phmdl12)**2)**2
    mdlerr=mdlerr+(2*np.pi*0.001/sin(phmdl13)**2)**2
    mdlerr=math.sqrt(mdlerr)/abs(denom)    

    term1 = 1/sin(phmdl12)**2/denom
    term2 = -1/sin(phmdl13)**2/denom

    denom=M12/M11-N12/N11+1e-16
    numer=-M12/M11/tan(ph2pi12)+N12/N11/tan(ph2pi13)
    alf=numer/denom

#    alfstd=        (M12/M11*2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
#    alfstd=alfstd+(N12/N11*2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
#    alfstd=math.sqrt(alfstd)/denom

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
    ph2pi12=phase.loc["MEAS", bn1, bn2]
    ph2pi23=phase.loc["MEAS", bn2, bn3]
    
    betstd = 0
    alfstd = 0

    # Find the model transfer matrices for beta1
    phmdl12=phase.loc["MODEL", bn1, bn2]
    phmdl23=phase.loc["MODEL", bn2, bn3]
    if plane=='H':
        betmdl1=madTwiss.loc[bn1, "BETX"]
        betmdl2=madTwiss.loc[bn2, "BETX"]
        betmdl3=madTwiss.loc[bn3, "BETX"]
        alpmdl2=madTwiss.loc[bn2, "ALFX"]
    elif plane=='V':
        betmdl1=madTwiss.loc[bn1, "BETY"]
        betmdl2=madTwiss.loc[bn2, "BETY"]
        betmdl3=madTwiss.loc[bn3, "BETY"]
        alpmdl2=madTwiss.loc[bn2, "ALFY"]
    if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)
        

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

#    betstd=        (2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
#    betstd=betstd+(2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
#    betstd=math.sqrt(betstd)/abs(denom)

    mdlerr=        (2*np.pi*0.001/sin(phmdl12)**2)**2
    mdlerr=mdlerr+(2*np.pi*0.001/sin(phmdl23)**2)**2
    mdlerr=math.sqrt(mdlerr)/abs(denom)

    term2 = 1/sin(phmdl23)**2/denom  #sign
    term1 = -1/sin(phmdl12)**2/denom  #sign

    denom=M12/M22+N12/N11+1e-16
    numer=M12/M22/tan(ph2pi12)-N12/N11/tan(ph2pi23)
    alf=numer/denom

#    alfstd=        (M12/M22*2*np.pi*phase["".join([plane,bn1,bn2])][1]/sin(ph2pi12)**2)**2
#    alfstd=alfstd+(N12/N11*2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
#    alfstd=math.sqrt(alfstd)/abs(denom)

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
    ph2pi23=phase.loc["MEAS", bn2, bn3]
    ph2pi13=phase.loc["MEAS", bn1, bn3]

    err23=phase.loc["ERRMEAS", bn2, bn3]
    err13=phase.loc["ERRMEAS", bn1, bn3]
    
    mdlerr = 0
    alfstd = 0

    # Find the model transfer matrices for beta1
    phmdl23=phase.loc["MODEL", bn2, bn3]
    phmdl13=phase.loc["MODEL", bn2, bn3]
    if plane=='H':
        betmdl1=madTwiss.loc[bn1, "BETX"]
        betmdl2=madTwiss.loc[bn2, "BETX"]
        betmdl3=madTwiss.loc[bn3, "BETX"]
        alpmdl3=madTwiss.loc[bn3, "ALFX"]
    elif plane=='V':
        betmdl1=madTwiss.loc[bn1, "BETY"]
        betmdl2=madTwiss.loc[bn2, "BETY"]
        betmdl3=madTwiss.loc[bn3, "BETY"]
        alpmdl3=madTwiss.loc[bn3, "ALFY"]
    if betmdl3 < 0 or betmdl2<0 or betmdl1<0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)

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

    betstd =          (err23 / sin(ph2pi23)**2)**2
    betstd = betstd + (err13 / sin(ph2pi13)**2)**2
    betstd = math.sqrt(betstd) / abs(denom)

#    mdlerr=        (2*np.pi*0.001/sin(phmdl23)**2)**2
#    mdlerr=mdlerr+(2*np.pi*0.001/sin(phmdl13)**2)**2
#    mdlerr=math.sqrt(mdlerr)/abs(denom)

    term2 = -1/sin(phmdl23)**2/denom  #sign
    term1 = 1/sin(phmdl13)**2/denom  #sign

    denom=M12/M22-N12/N22+1e-16
    numer=M12/M22/tan(ph2pi23)-N12/N22/tan(ph2pi13)
    alf=numer/denom

#    alfstd=        (M12/M22*2*np.pi*phase["".join([plane,bn2,bn3])][1]/sin(ph2pi23)**2)**2
#    alfstd=alfstd+(N12/N22*2*np.pi*phase["".join([plane,bn1,bn3])][1]/sin(ph2pi13)**2)**2
#    alfstd=math.sqrt(alfstd)/abs(denom)


    return bet, betstd, alf, alfstd, mdlerr, term1, term2, True

#=======================================================================================================================
#---============== using analytical formula ============================================================================
#=======================================================================================================================



def scan_all_BPMs_withsystematicerrors(madTwiss, madElements,
                                       phase, plane, getllm_d, commonbpms, debugfile, errors_method,
                                       tune, mdltune):
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
    
    print_("Errors from " + ID_TO_METHOD[errors_method])
    
    
    width = getllm_d.range_of_bpms / 2
    left_bpm = range(-width, 0)
    right_bpm = range(0 + 1, width)
    BBA_combo = [[x, y] for x in left_bpm for y in left_bpm if x < y]
    ABB_combo = [[x, y] for x in right_bpm for y in right_bpm if x < y]
    BAB_combo = [[x, y] for x in left_bpm for y in right_bpm]
    
    # get the model values only for used elements, so that commonbps[i] = masTwiss[i]
    madTwiss_intersected = madTwiss.loc[commonbpms.index]
    mu = "MUX" if plane == "H" else "MUY"
    mu_elements = madElements.loc[:][mu].values
#    print mu_elements
#    elements_phases = (madTwiss_intersected.loc[:][mu].values[:, np.newaxis] - 
#                       mu_elements[np.newaxis, :]) * TWOPI
#    print elements_phases.shape
    
#    cot_meas = 1.0 / tan(phase["MEAS"] * TWOPI)
#    cot_model = 1.0 / tan(phase["MODEL"] * TWOPI)
    phases_meas = phase["MEAS"] * TWOPI
    phases_model = phase["MODEL"] * TWOPI
    phases_err = phase["ERRMEAS"] * TWOPI
#    cot_meas = 1.0 / tan(tilt_slice_matrix(phase["MEAS"], 20, 40, tune) * TWOPI)
#    cot_model = 1.0 / tan(tilt_slice_matrix(phase["MODEL"], 20, 40, mdltune) * TWOPI)
#
#    mdltune = mdltune % 1.0

    result = np.array(np.empty(madTwiss_intersected.shape[0]), 
                      dtype = "S24, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, i4")
#                      dtype=[("NAME", str), ("S", np.float64),
#                             ("BET", np.float64),("BETSTAT", np.float64),("BETSYS", np.float64),("BETERR", np.float64),
#                             ("ALF", np.float64),("ALFSTAT", np.float64),("ALFSYS", np.float64),("ALFERR", np.float64),
#                             ("CORR", np.float64),("BETMDL", np.float64),("BB", np.float64),("NCOMB", int)])
#    
    
    print "skuhfrghrjghj------------------------------"


    def collect(row):
#        if row[11]:
            result[row[0]]= row[1:]
        
    def collectblock(block):
        for row in block:
#            if row[11]:
            result[row[0]] = row[1:]
                
    st = time.time()
    if getllm_d.parallel and not DEBUG:
        
        chunksize = int(len(commonbpms) / getllm_d.nprocesses) + 1
        pool = multiprocessing.Pool()
        n = int(len(commonbpms) / chunksize)
        
        for i in range(n):
            pool.apply_async(scan_several_BPMs_withsystematicerrors,
                             (madTwiss_intersected, madElements,
                              phases_meas, phases_err,
                              plane, getllm_d.range_of_bpms, commonbpms, debugfile,
                              i * chunksize, (i + 1) * chunksize, BBA_combo, ABB_combo, BAB_combo,
                              tune, mdltune),
                             callback=collectblock)
        pool.apply_async(scan_several_BPMs_withsystematicerrors,
                         (madTwiss_intersected, madElements,
                          phases_meas, phases_model,
                          plane, getllm_d.range_of_bpms, commonbpms, debugfile,
                          n * chunksize, len(commonbpms), BBA_combo, ABB_combo, BAB_combo,
                          tune, mdltune),
                         callback=collectblock)
        pool.close()
        pool.join()
    else:
        startProgress("Scan all BPMs")
        for i in range(0, len(commonbpms)):
            if (i % 20 == 0):
                progress(float(i) * 100.0 / len(commonbpms))
            row = scan_one_BPM_withsystematicerrors(madTwiss_intersected, madElements,
                                                    phases_meas, phases_err,
                                                    plane, getllm_d.range_of_bpms, commonbpms,
                                                    debugfile, i,
                                                    BBA_combo, ABB_combo, BAB_combo,
                                                    tune, mdltune)
            collect(row)
        endProgress()
    et = time.time()
    
    print_("time elapsed = {0:3.3f}".format(et - st))
    
    if DEBUG:
        debugfile.close()
    rmsbb = -1
    return rmsbb, errors_method, result


def scan_several_BPMs_withsystematicerrors(madTwiss, madElements,
                                           cot_meas, phases_err,
                                           plane, range_of_bpms, commonbpms, debugfile,
                                           begin, end, BBA_combo, ABB_combo, BAB_combo,
                                           tune, mdltune):
    block = []
    for i in range(begin, end):
        block.append(scan_one_BPM_withsystematicerrors(madTwiss, madElements,
                                                       cot_meas, phases_err,
                                                       plane, range_of_bpms, commonbpms,
                                                       debugfile, i,
                                                       BBA_combo, ABB_combo, BAB_combo,
                                                       tune, mdltune))
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


def scan_one_BPM_withsystematicerrors(madTwiss, madElements,
                                      phases_meas, phases_err,
                                      plane, range_of_bpms, commonbpms, debugfile,
                                      Index, BBA_combo, ABB_combo, BAB_combo,
                                      tune, mdltune):
    '''
    Scans the range of BPMs in order to get the final value for one BPM in the lattice
    
    :Parameters:
        'madTwiss':tfs_pandas
            The model twiss table, contains all the BPMs. Has to be already intersected with the common BPMs.
            
                  | S  | BETX  | MUX | ... | BETY  | ...
            ------|----|-------|-----|-----|-------|-----
             BPM1 | s1 | beta1 | mu1 | ... | bety1 |
             BPM2 | s2 | beta2 | mu2 | ... | bety2 |
             ...  |    |       |     |     |       |
             BPMn | sn | betan | mun | ... | betyn |
             
        'madElements':tfs_pandas
            Twiss table of all elements with known uncertainties. The keys are already intersected with the list of
            common BPMs.
            
                  | S      | BETX      | MUX     | ... | dK1    | dS      | BPMdS
            ------|--------|-----------|---------|-----|--------|---------|--------
             BPM1 | s1     | beta1     | mu1     | ... | 0      | 0       | 1e-3
             MQ1  | s(El1) | beta(El1) | mu(El1) | ... | 1e-4   | 0       | 0
             MQ2  | s(El2) | beta(El2) | mu(El2) | ... | 2e-4   | 0       | 0
             MS1  | s(El3) | beta(El3) | mu(El3) | ... | 0      | 1.0e-3  | 0
             BPM2 | s2     | beta2     | mu2     | ... | 0      | 0       | 1e-3
             ...  | ...    | ...       | ...     | ... | ...    | ...     | ...
             BPMn | sn     | betan     | mun     | ... | 0      | 0       | 1e-3
             
             for later distinction we denote the row index of madElements by <l>:
                 s<1> := s1
                 s<2> := s(El1)
                 s<3> := s(El2)
                 etc.
             
        'phases_meas/phases_model':numpy.ndarray
            matrix of the cotanges of the meas/model phase advances.
            
                  | BPM1   | BPM2   | ... | BPMn 
            ------|--------|--------|-----|-----
             BPM1 | 0      | phi_21 | ... | phi_n1 
             BPM2 | phi_12 | 0      | ... | phi_n2
             ...  | ...    | ...    | ... | ...   
             BPMn | phi_1n | phi_2n | ... | 0 
             
             index and columns: madTwiss.index
             
             cotangens of shifted and sliced:
                 
                      | BPM1 | BPM2 | ...  | BPMn 
            ----------|------|------|------|-----
             BPM(i-2) | *    | *    | ...  | * 
             BPM(i-1) | *    | *    | ...  | * 
             BPMi     | 0    | 0    | 0..0 | 0   
             BPM(i+1) | *    | *    | ...  | * 
             BPM(i+2) | *    | *    | ...  | * 
             
             where * = cot(phi_j - phi_i) =: cot(phi_ij)

           
        'phases_elements':numpy.ndarray
            
                 | BPM1 | MQ1 | MQ2 | MS1 | BPM2 | ... | BPMn
            -----|------|-----|-----|-----|------|-----|------
            BPM1 | 0    | *   | *   | *   | *    | ... | *
            BPM2 | *    | *   | *   | *   | 0    | ... | *
            BPM3 | *    | *   | *   | *   | *    | ... | *
            ...  | ...  | ... | ... | ... | ...  | ... | ...
            BPMn | *    | *   | *   | *   | *    | ... | 0

            where * = phi_j - phi<i>
            
            columns: madElements.index
            rows: madTwiss.index
            
    IMPORTANT NOTE: from the above only madTwiss and madElements are pandas.DataFrames, cot_meas and sin_squared_elements
                    are numpy arrays and thus are not equipped with an index or column headers.
            
    The heavy fancy indexing part is explained in detail:
        In the following the probed BPM has the name probed_bpm and the index i. The range of BPMs has the length J and
        m = floor(J/2).
        
        The range of BPMs is then [BPM(i-m), ..., BPM(i-1), BPMi, BPM(i+1), ... BPM(i+m)]
        The combinations are split into the three cases BBA, BAB, ABB
        BBA_combo = [(-m, -m+1), ... (-2, -1)]
        BAB_combo = [(-m, 1), ... (-m, m), (-m+1, 1), ... (-m+1, m), ... (-1, m)]
        ABB_combo = [(1,2), ... (m-1, m)]
        
        for each (j,k) in combos: calculate beta and row of the Jacobian T:
        
            beta_i(combo) = (cot(phi_ij) - cot(phi_ik)) / (cot(phimdl_ij) - cot(phimdl_ik)) * betmdl_i for (j,k) in combo
            
        So, before we start the loop, let's slice out the interval of all BPMs and all elements that are reached by the
        combinations (outer interval) at this time we can apply the tune jump, wrap the interval and calculate the 
        trigonometric functions of the phase advances
            
         [] for this we need cot(phi_(i)(i-m) ... cot(phi_(i)(i+m)) and idem for the model phases 
            
            K-part of the Jacobian
            T_k(combo) = betmdl_i * betmdl<l> / (cot(phimdl_ij) - cot(phimdl_ik)) *
                            (
                            sin^2(phimdl_j - phimdl<l>) / sin^2(phimdl_j - phimdl_i) * A(i,j) - 
                            sin^2(phimdl_k - phimdl<l>) / sin^2(phimdl_k - phimdl_i) * A(i,k)
                            )
            where A(i,j) = 1 if i < j, else -1 and it is 0 if <l> is not between BPMi and BPMj
            => BBA_combo: A(i,j) = -1, A(i,k) = -1
               BAB_combo: A(i,j) = -1, A(i,k) = 1
               ABB_combo: A(i,j) = A(i,k) = 1
         [] for this we need sin^2(phimdl_(i-m) - phimdl_(i-m+<1>)), sin^2(phimdl_(i-m) - phimdl_(i-m+<2>)), ...
                              sin^2(phimdl_(i+m) - phimdl_(i+m-<2>)), sin^2(phimdl_(i+m) - phimdl_(i+m-<1>))
            which is again the range of BPMs but with all other Elements lying between them.
    
        Take the left-most combination, that is (i-m, i-m+1). By construction of the combo sets this is the first
        combination in BBA_combo. Analogously the right-most combination is the last element of ABB_combo.
        
        IDEA: Use the index of madTwiss and madElements to get the location of the start and end elements and then slice 
        intelligently. 
        
        indx_first = indx_(i-m) = i-m
        indx_last = indx(i+m) = i+m
        name_first = madTwiss.index[indx_first], same for last
        indx<first> = madElements.index.get_loc(name_first), same for last
        
        if indx_(i-m) < 0: 
            have to wrap around and add the tune
            outer_interval_meas = [phases_meas[index%length] ...] + tune concat [phases_meas[0] ...]
                same for outer_interval_model
            outer_interval_elements = [phases_elements[indx<first>] ... phases_elements[length]] + tune
                                        concat  [phases_elements[0] ... phases_elements[indx<last>]]
        else if indx_(i+m) > length: analogously
        else: simple
        
        
        
        
        
    '''
    probed_bpm_name = madTwiss.index[Index]
    s = madTwiss.get_value(probed_bpm_name, "S")


    if plane == 'H':
        betmdl1 = madTwiss.get_value(probed_bpm_name, "BETX")
        mu_column = "MUX"
        bet_column = "BETX"
    elif plane == 'V':
        betmdl1 = madTwiss.get_value(probed_bpm_name, "BETX")
        mu_column = "MUY"
        bet_column = "BETY"
    
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
    
    m = range_of_bpms / 2
    indx_first = Index - m
    indx_last = Index + m
    name_first = madTwiss.index[indx_first]
    name_last = madTwiss.index[indx_last% len(madTwiss.index)]
    probed_bpm_name = madTwiss.index[Index]
    len_bpms_total = phases_meas.shape[0]
    
    indx_el_first = madElements.index.get_loc(name_first)
    indx_el_last= madElements.index.get_loc(name_last )
    
    
    if indx_first < 0:
        outerMeasPhaseAdv = pd.concat((
                phases_meas.iloc[Index, indx_first % len_bpms_total:] - tune * TWOPI,
                phases_meas.iloc[Index, :indx_last + 1]))
        outerMeasErr = pd.concat((
                phases_err.iloc[Index, indx_first % len_bpms_total:],
                phases_err.iloc[Index, :indx_last + 1]))
        outerMdlPh = np.concatenate((
                madTwiss.iloc[indx_first % len_bpms_total:][mu_column] - mdltune,
                madTwiss.iloc[:indx_last + 1][mu_column])) * TWOPI
        outerElmts = pd.concat((
                madElements.iloc[indx_el_first:],
                madElements.iloc[:indx_el_last + 1]))
        outerElmtsPh = np.concatenate((
                madElements.iloc[indx_el_first:][mu_column] - mdltune,
                madElements.iloc[:indx_el_last + 1][mu_column])) * TWOPI
    elif indx_last >= len_bpms_total:
        outerMeasPhaseAdv = pd.concat((
                phases_meas.iloc[Index, indx_first:],
                phases_meas.iloc[Index, :indx_last % len_bpms_total] + tune * TWOPI))
        outerMeasErr = pd.concat((
                phases_err.iloc[Index, indx_first:],
                phases_err.iloc[Index, :indx_last % len_bpms_total]))
        outerMdlPh = np.concatenate((
                madTwiss.iloc[indx_first:][mu_column],
                madTwiss.iloc[:indx_last % len_bpms_total][mu_column]  + mdltune)) * TWOPI
        outerElmts = pd.concat((
                madElements.iloc[indx_el_first:],
                madElements.iloc[:indx_el_last]))
        outerElmtsPh = np.concatenate((
                madElements.iloc[indx_el_first:][mu_column],
                madElements.iloc[:indx_el_last + 1][mu_column] + mdltune)) * TWOPI
        
    else:
        outerMeasPhaseAdv = phases_meas.iloc[Index, indx_first : indx_last]
        outerMeasErr = phases_err.iloc[Index, indx_first : indx_last]
        outerMdlPh = madTwiss.iloc[indx_first:indx_last][mu_column].as_matrix() * TWOPI
        outerElmts = madElements.iloc[indx_el_first:indx_el_last]
        outerElmtsPh = madElements.iloc[indx_el_first:indx_el_last][mu_column] * TWOPI
    outerElPhAdv = (outerElmtsPh[:, np.newaxis] - outerMdlPh[np.newaxis, :])
    indx_el_probed = outerElmts.index.get_loc(probed_bpm_name)
    outerElmtsBet = outerElmts.loc[:][bet_column].as_matrix()
        
    cot_meas = 1.0 / tan(outerMeasPhaseAdv.as_matrix())
    cot_model = 1.0 / tan((outerMdlPh - outerMdlPh[m])) 
    outerElPhAdv = sin(outerElPhAdv)
    sin_squared_elements = np.multiply(outerElPhAdv, outerElPhAdv)
    
#    print "outerMdlPh"
#    print outerMdlPh
#    print "outerELmtsPh"
#    print outerElmtsPh
#    print "outerElPhAdv"
#    print outerElPhAdv
#        
#    print interval_all
    betas = np.empty(len(BBA_combo) + len(BAB_combo) + len(ABB_combo))
    beta_mask = np.empty(len(BBA_combo) + len(BAB_combo) + len(ABB_combo), dtype=bool)
#    print betmdl1
#    M = np.zeros((len(sin_squared_elements) * 3,
#                  len(sin_squared_elements) * 3))
    
#    diag = np.concatenate((outerElmts.loc[:]["dK1"],outerElmts.loc[:]["dS"],outerElmts.loc[:]["dX"]))
    diag = np.concatenate((outerMeasErr.as_matrix(), outerElmts.loc[:]["dK1"]))
    mask = diag != 0
    
    T_Beta = np.zeros((len(betas),
                   len(diag) ))

    M = np.diag(diag[mask])
    
#    print cot_model[m]
#    print cot_model
    
    for i, combo in enumerate(BBA_combo):
        ix = combo[0] + m
        iy = combo[1] + m
        
        if (abs(cot_model[ix]) > COT_THRESHOLD or 
            abs(cot_model[iy]) > COT_THRESHOLD or
            abs(cot_meas[ix]) > COT_THRESHOLD or
            abs(cot_meas[iy]) > COT_THRESHOLD or
            (cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
            beta_mask[i] = False
            continue
        beta_mask[i] = True
        
        denom = (cot_model[ix] - cot_model[iy]) / betmdl1
        betas[i] = (cot_meas[ix] - cot_meas[iy]) / denom
        
        
        
        xloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[ix])
        yloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[iy])
        
        elementPh_XA = sin_squared_elements[xloc:indx_el_probed, ix]
        elementPh_YA = sin_squared_elements[yloc:indx_el_probed, iy]
        elementBet_XA = outerElmtsBet[xloc:indx_el_probed]
        elementBet_YA = outerElmtsBet[yloc:indx_el_probed]
        denom_sinx = sin_squared_elements[xloc, m]
        denom_siny = sin_squared_elements[yloc, m]
        
        T_Beta[i][ix] = -1.0 / (denom_sinx * denom)
        T_Beta[i][iy] = 1.0 / (denom_siny * denom)

        T_Beta[i][xloc+range_of_bpms:indx_el_probed+range_of_bpms] += elementPh_XA * elementBet_XA / (denom_sinx * denom)
        T_Beta[i][yloc+range_of_bpms:indx_el_probed+range_of_bpms] -= elementPh_YA * elementBet_YA / (denom_siny * denom)
    offset_i = len(BBA_combo)
        
    for i, combo in enumerate(BAB_combo):
        ix = combo[0] + m
        iy = combo[1] + m
        
        if (abs(cot_model[ix]) > COT_THRESHOLD or 
            abs(cot_model[iy]) > COT_THRESHOLD or
            abs(cot_meas[ix]) > COT_THRESHOLD or
            abs(cot_meas[iy]) > COT_THRESHOLD or
            (cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
            beta_mask[i + offset_i] = False
            continue
        beta_mask[i + offset_i] = True

        
        denom = (cot_model[ix] - cot_model[iy]) / betmdl1
        betas[i + offset_i] = (cot_meas[ix] - cot_meas[iy]) / denom
        
        xloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[ix])
        yloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[iy])
        
        elementPh_XA = sin_squared_elements[xloc:indx_el_probed, ix]
        elementPh_YA = sin_squared_elements[indx_el_probed:yloc, iy]
        elementBet_XA = outerElmtsBet[xloc:indx_el_probed]
        elementBet_YA = outerElmtsBet[indx_el_probed:yloc]
        denom_sinx = sin_squared_elements[xloc, m]
        denom_siny = sin_squared_elements[yloc, m]
        
        T_Beta[i + offset_i][ix] = -1.0 / (denom_sinx * denom)
        T_Beta[i + offset_i][iy] = 1.0 / (denom_siny * denom)

        
        T_Beta[i + offset_i, xloc+range_of_bpms:indx_el_probed+range_of_bpms] += elementPh_XA * elementBet_XA / (denom_sinx * denom)
        T_Beta[i + offset_i, indx_el_probed+range_of_bpms:yloc+range_of_bpms] += elementPh_YA * elementBet_YA / (denom_siny * denom)
    offset_i += len(BAB_combo)
        
    for i, combo in enumerate(ABB_combo):
        ix = combo[0] + m 
        iy = combo[1] + range_of_bpms / 2
        
        if (abs(cot_model[ix]) > COT_THRESHOLD or 
            abs(cot_model[iy]) > COT_THRESHOLD or
            abs(cot_meas[ix]) > COT_THRESHOLD or
            abs(cot_meas[iy]) > COT_THRESHOLD or
            (cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
            beta_mask[i + offset_i] = False
            continue
        beta_mask[i + offset_i] = True

        
        denom = (cot_model[ix] - cot_model[iy]) / betmdl1
        betas[i + offset_i] = (cot_meas[ix] - cot_meas[iy]) / denom
        
        xloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[ix])
        yloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[iy])
        
        elementPh_XA = sin_squared_elements[indx_el_probed:xloc, ix]
        elementPh_YA = sin_squared_elements[indx_el_probed:yloc, iy]
        elementBet_XA = outerElmtsBet[indx_el_probed:xloc]
        elementBet_YA = outerElmtsBet[indx_el_probed:yloc]
        denom_sinx = sin_squared_elements[xloc, m]
        denom_siny = sin_squared_elements[yloc, m]
        
        T_Beta[i + offset_i][ix] = -1.0 / (denom_sinx * denom)
        T_Beta[i + offset_i][iy] = 1.0 / (denom_siny * denom)
        
#        print elementPh_XA * elementBet_XA / (denom_sinx * denom)

        
        T_Beta[i + offset_i][indx_el_probed+range_of_bpms:xloc+range_of_bpms] -= elementPh_XA * elementBet_XA / (denom_sinx * denom)
        T_Beta[i + offset_i][indx_el_probed+range_of_bpms:yloc+range_of_bpms] += elementPh_YA * elementBet_YA / (denom_siny * denom)

    T_Beta = T_Beta[:, mask]
    T_Beta = T_Beta[beta_mask]
    betas = betas[beta_mask]
    
  
    
    V_Beta = np.dot(T_Beta, np.dot(M,np.transpose(T_Beta)))
    try:
        V_Beta_inv = np.linalg.pinv(V_Beta, rcond=RCOND)
        w = np.sum(V_Beta_inv, axis=1)
#        print w
#        raw_input()
        VBeta_inv_sum = np.sum(w)
        beterr = math.sqrt(float(np.dot(np.transpose(w), np.dot(V_Beta, w)) / VBeta_inv_sum ** 2))
        beti = float(np.dot(np.transpose(w), betas) / VBeta_inv_sum)
        used_bpms = len(w)
    except ValueError:
        print "ValueError"
    
    return (Index, probed_bpm_name, s,
                beti, betstat, betsys, beterr,
                alfi, alfstat, alfsys, alferr,
                .0, betmdl1, (beti - betmdl1) / betmdl1,
                len(betas))
    
    
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
        
        return (probed_bpm_name,
                beti, betstat, betsys, beterr,
                alfi, alfstat, alfsys, beterr,
                .0, (beti - betmdl1) / betmdl1,
                -1)

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


def get_beta_from_phase_systematic_errors(madTwiss, madElements, phase, plane, commonbpms,
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

            err_diagonal[k + position] = errorfile.elements[list_of_Ks[index][0][k], DK1_INDEX] ** 2

        position += len(list_of_Ks[index][0])
        
        #---- assign sextupole transversal missalignments
    for j in range(RANGE):
        index = (CurrentIndex + j) % len(list_of_Ks)
        for k in range(len(list_of_Ks[index][1])):

            err_diagonal[k + position] = errorfile.elements[list_of_Ks[index][1][k], DX_INDEX] ** 2

        position += len(list_of_Ks[index][1])
        
        #---- assign longitudinal missalignments
    for j in range(RANGE):
        index = (CurrentIndex + j) % len(list_of_Ks)
        for k in range(len(list_of_Ks[index][2])):

            err_diagonal[k + position] = errorfile.elements[list_of_Ks[index][2][k], DS_INDEX] ** 2

        position += len(list_of_Ks[index][2])
        
        #---- assign BPM missalignments
    for j in range(RANGE):
        err_diagonal[j + position] = errorfile.elements[errorfile.indx[bpm_name[j]], DS_INDEX] ** 2
        
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
    
    elements = errorfile.elements
   
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
                       elements, list_of_Ks,
                       denomalf, s_i2, T, T_Alf, phi_2, -frac, K_offset,
                       .5, .5)
       
    K_offset = _assign_quaderrors(I, bi1, bi3,
                                  elements, list_of_Ks,
                                  denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                  .5, .5)
            
    #--- Sext Transverse Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
        
    _assign_sext_errors(I, bi1, bi2,
                        elements, list_of_Ks,
                        denomalf, s_i2, T, T_Alf, phi_2, frac, K_offset,
                        -.5, .5)
    
    K_offset = _assign_sext_errors(I, bi1, bi3,
                                   elements, list_of_Ks,
                                   denomalf, s_i3, T, T_Alf, phi_3, -frac, K_offset,
                                   -.5, .5)

    #--- Quad Longitudinal Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])

    _assign_quadlongmissal(I, bi1, bi2,
                           elements, list_of_Ks,
                           denomalf, s_i2, T, T_Alf, phi_2, -frac, K_offset,
                           .5, .5)
       
    K_offset = _assign_quadlongmissal(I, bi1, bi3,
                                      elements, list_of_Ks,
                                      denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                      .5, .5)
    
    #--- BPM Missalignments
    # jump to end of RANGE
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
          
    errindx1 = errorfile.indx[bn1]
    errindx2 = errorfile.indx[bn2]
    errindx3 = errorfile.indx[bn3]
      
    if elements[errindx1, DS_INDEX] != 0:
        numerphi = -1.0 / (betmdl1 * sin(phmdl12) ** 2) + 1.0 / (betmdl1 * sin(phmdl13) ** 2)
        T[K_offset + bi1] = numerphi / denom - 2 * alfmdl1
          
    if elements[errindx2, DS_INDEX] != 0:
        numerphi = 1.0 / (betmdl2 * sin(phmdl12) ** 2)
        T[K_offset + bi2] = numerphi / denom
       
    if elements[errindx3, DS_INDEX] != 0:
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
    
    elements = errorfile.elements
    
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
                                  elements, list_of_Ks,
                                  denomalf, s_i1, T, T_Alf, phi_1, frac,
                                  K_offset,
                                  -.5, .5)
     
    K_offset = _assign_quaderrors(I, bi2, bi3,
                                  elements, list_of_Ks,
                                  denomalf, s_i3, T, T_Alf, phi_3, frac,
                                  K_offset,
                                  .5, .5)
     
    #--- Sext Transverse Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
        
    K_offset = _assign_sext_errors(I, bi1, bi2,
                                   elements, list_of_Ks,
                                   denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset,
                                   .5, .5)
                
    K_offset = _assign_sext_errors(I, bi2, bi3,
                                   elements, list_of_Ks,
                                   denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                   -.5, .5)
    #--- Quad Longitudinal Missalignments
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])

    K_offset = _assign_quadlongmissal(I, bi1, bi2,
                                      elements, list_of_Ks,
                                      denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset,
                                      -.5, .5)
    
    K_offset = _assign_quadlongmissal(I, bi2, bi3,
                                      elements, list_of_Ks,
                                      denomalf, s_i3, T, T_Alf, phi_3, frac, K_offset,
                                      .5, .5)
    #--- BPM Missalignments
    # jump to end of RANGE
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
          
    errindx1 = errorfile.indx[bn1]
    errindx2 = errorfile.indx[bn2]
    errindx3 = errorfile.indx[bn3]
      
    if elements[errindx2, DS_INDEX] != 0:
        numerphi = -1.0 / (betmdl2 * sin(phmdl21) ** 2) + 1.0 / (betmdl2 * sin(phmdl23) ** 2)
        T[K_offset + bi2] = numerphi / denom - 2.0 * alpmdl2
          
    if elements[errindx1, DS_INDEX] != 0:
        numerphi = 1.0 / (betmdl1 * sin(phmdl21) ** 2)
        T[K_offset + bi1] = numerphi / denom
       
    if elements[errindx3, DS_INDEX] != 0:
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
    
    elements = errorfile.elements
    
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
    
    K_offset = _assign_quaderrors(I, bi1, bi3,
                                  elements, list_of_Ks,
                                  denomalf, s_i1, T, T_Alf, phi_1, -frac,
                                  K_offset,
                                  -.5, .5)
    
    # go back because the second h_ij begins at 2
    # so we have to find the position of BPM2 in the matrix
    K_offset = K_begin
    for k in range(bi1, bi2):
        which_k = (k + I) % len(list_of_Ks)
        K_offset += len(list_of_Ks[which_k][0])
     
    K_offset = _assign_quaderrors(I, bi2, bi3,
                                  elements, list_of_Ks,
                                  denomalf, s_i2, T, T_Alf, phi_2, frac,
                                  K_offset,
                                  -.5, .5)
    
    #--- Sext Trasverse Missalignments
    # jump to the end of RANGE then to b1
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][0])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    K_begin = K_offset
    
    K_offset = _assign_sext_errors(I, bi1, bi3,
                                   elements, list_of_Ks,
                                   denomalf, s_i1, T, T_Alf, phi_1, -frac, K_offset,
                                   .5, .5)
    
    K_offset = K_begin
    for k in range(bi1, bi2):
        which_k = (k + I) % len(list_of_Ks)
        K_offset += len(list_of_Ks[which_k][1])
                
    K_offset = _assign_sext_errors(I, bi2, bi3,
                                   elements, list_of_Ks,
                                   denomalf, s_i2, T, T_Alf, phi_2, frac, K_offset,
                                   .5, .5)
    
    #--- Quad Longitudinal Missalignments
    # jump to the end of RANGE then to b1
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][1])
    for k in range(bi1):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
    K_begin = K_offset
    
    K_offset = _assign_quadlongmissal(I, bi1, bi3,
                                      elements, list_of_Ks,
                                      denomalf, s_i1, T, T_Alf, phi_1, -frac, K_offset,
                                      -.5, .5)
    
    K_offset = K_begin
    for k in range(bi1, bi2):
        which_k = (k + I) % len(list_of_Ks)
        K_offset += len(list_of_Ks[which_k][2])
    
    K_offset = _assign_quadlongmissal(I, bi2, bi3,
                                      elements, list_of_Ks,
                                      denomalf, s_i2, T, T_Alf, phi_2, frac, K_offset,
                                      -.5, .5)
      
   
    #--- BPM Missalignments
    # jump to end of RANGE
    for k in range(bi3, RANGE):
        K_offset += len(list_of_Ks[(k + I) % len(list_of_Ks)][2])
          
    errindx1 = errorfile.indx[bn1]
    errindx2 = errorfile.indx[bn2]
    errindx3 = errorfile.indx[bn3]
      
    if elements[errindx3, DS_INDEX] != 0:
        numerphi = -1.0 / (betmdl3 * sin(phmdl32) ** 2) + 1.0 / (betmdl3 * sin(phmdl31) ** 2)
        T[K_offset + bi3] = numerphi / denom - 2.0 * alpmdl3
          
    if elements[errindx2, DS_INDEX] != 0:
        numerphi = 1.0 / (betmdl2 * sin(phmdl32) ** 2)
        T[K_offset + bi2] = numerphi / denom
       
    if elements[errindx1, DS_INDEX] != 0:
        numerphi = -1.0 / (betmdl1 * sin(phmdl31) ** 2)
        T[K_offset + bi1] = numerphi / denom
            
    patternstr = ["x"] * RANGE
    patternstr[bi1] = "B"
    patternstr[bi2] = "B"
    patternstr[bi3] = "A"

    return MeasuredValues(alf, bet, "".join(patternstr), True), T, T_Alf


def _assign_quaderrors(I, bi1, bi2, elements, list_of_Ks, denomalf, sinus_ij_squared, T_Bet, T_Alf, reference_phi, frac, K_offset, Alf_fact_1, Alf_fact_2):
    for k in range(bi1, bi2):
        whichK = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[whichK][0])):
            idx_k = list_of_Ks[whichK][0][w]
            err_beta = frac * elements[idx_k, BET_INDEX] * (sin(elements[idx_k, MU_INDEX] * TWOPI - reference_phi) ** 2 / sinus_ij_squared)
            T_Bet[K_offset + w] += err_beta
            T_Alf[K_offset + w] += Alf_fact_1 * elements[idx_k, BET_INDEX] * (sin(elements[idx_k, MU_INDEX] * TWOPI - reference_phi) ** 2 / sinus_ij_squared)
            T_Alf[K_offset + w] += Alf_fact_2 * err_beta * denomalf
        
        K_offset += len(list_of_Ks[whichK][0])
    
    return K_offset


def _assign_sext_errors(I, bi1, bi2, elements, list_of_Ks, denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset, Alf_fact_1, Alf_fact_2):
    for k in range(bi1, bi2):
        whichK = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[whichK][1])):
            idx_k = list_of_Ks[whichK][1][w]
            err_beta = -SEXT_FACT * frac * elements[idx_k, K2L_INDEX] * elements[idx_k, BET_INDEX] * (sin(elements[idx_k, MU_INDEX] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += Alf_fact_1 * SEXT_FACT * elements[idx_k, K2L_INDEX] * elements[idx_k, BET_INDEX] * (sin(elements[idx_k, MU_INDEX] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += Alf_fact_2 * err_beta * denomalf
        
        K_offset += len(list_of_Ks[whichK][1])
    
    return K_offset


def _assign_quadlongmissal(I, bi1, bi2, elements, list_of_Ks, denomalf, s_i1, T, T_Alf, phi_1, frac, K_offset, Alf_fact_1, Alf_fact_2):
    for k in range(bi1, bi2):
        whichK = (k + I) % len(list_of_Ks)
        for w in range(len(list_of_Ks[whichK][2])):
            idx_k = list_of_Ks[whichK][2][w]
            err_beta = frac * elements[idx_k, K1LEND_INDEX] * elements[idx_k, BETEND_INDEX] * (sin(elements[idx_k, MUEND_INDEX] * TWOPI - phi_1) ** 2 / s_i1)
            T[K_offset + w] += err_beta
            T_Alf[K_offset + w] += Alf_fact_1 * elements[idx_k, K1LEND_INDEX] * elements[idx_k, BETEND_INDEX] * (sin(elements[idx_k, MUEND_INDEX] * TWOPI - phi_1) ** 2 / s_i1)
            T_Alf[K_offset + w] += Alf_fact_2 * err_beta * denomalf
        
        K_offset += len(list_of_Ks[whichK][2])
    
    return K_offset

#===================================================================================================
#--- ac-dipole stuff
#===================================================================================================


def _get_free_beta(modelfree, modelac, data, bpms, plane):  # to check "+"

    print data    
    beta = data["BET"]
    betsys = data["BETSYS"]
    betstd = data["BETSTD"]
    beterr = data["BETERR"]
    alfa = data["ALF"]
    alfsys = data["ALFSYS"]
    alfstd = data["ALFSTD"]
    alferr = data["ALFERR"]

    if plane == "H":
        betmf = modelfree.loc[bpms.index, "BETX"]
        betma = modelac.loc[bpms.index, "BETX"]
        bb = betma / betmf
        alfmf = modelfree.loc[bpms.index, "ALFX"]
        alfma = modelac.loc[bpms.index, "ALFX"]
        aa = alfma / alfmf
    else:
        betmf = modelfree.loc[bpms.index, "BETY"]
        betma = modelac.loc[bpms.index, "BETY"]
        bb = betma / betmf
        alfmf = modelfree.loc[bpms.index, "ALFY"]
        alfma = modelac.loc[bpms.index, "ALFY"]
        aa = alfma / alfmf
    data2 = beta / bb, betsys, betstat, beterr, alfa * aa, alfsys, alfstat, alferr, data[bpms.index][8], data[bpms.index][9], data[bpms.index][10]
#    for i,bpm in enumerate(bpms.index):
#        beta, betsys, betstat, beterr = data[bpm][["BET", "BETSYS", "BETSTD", "BETERR"]]
#        alfa, alfsys, alfstat, alferr = data[bpm][["ALF", "ALFSYS", "ALFSTD", "ALFERR"]]
#
#        if plane == "H":
#            betmf = modelfree.BETX[bpm]
#            betma = modelac.BETX[bpm]
#            bb = betma / betmf
#            alfmf = modelfree.ALFX[bpm]
#            alfma = modelac.ALFX[bpm]
#            aa = alfma / alfmf
#        else:
#            betmf = modelfree.BETY[bpm]
#            betma = modelac.BETY[bpm]
#            alfmf = modelfree.ALFY[bpm]
#            alfma = modelac.ALFY[bpm]
#            bb = betma / betmf
#            aa = alfma / alfmf
#        data2[bpm] = beta / bb, betsys, betstat, beterr, alfa * aa, alfsys, alfstat, alferr, data[bpm][8], data[bpm][9], data[bpm][10]
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

def create_listofKs(commonbpms, errorfile, el, getllm_d):
    list_of_Ks = []
    
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
                
                if el[i, DK1_INDEX] != 0:
                    quad_fields.append(i)
                if el[i, DX_INDEX] != 0:
                    sext_trans.append(i)
                if el[i, DS_INDEX] != 0:
                    quad_missal.append(i)
                    
        else:
            for i in range(index_n + 1, errorfile.size):
                if el[i, DK1_INDEX] != 0:
                    quad_fields.append(i)
                if el[i, DX_INDEX] != 0:
                    sext_trans.append(i)
                if el[i, DS_INDEX] != 0:
                    quad_missal.append(i)
                    
            for i in range(index_nplus1):  # ums Eck
                if el[i, DK1_INDEX] != 0:
                    quad_fields.append(i)
                if el[i, DX_INDEX] != 0:
                    sext_trans.append(i)
                if el[i, DS_INDEX] != 0:
                    quad_missal.append(i)
                    
        list_of_Ks.append([quad_fields, sext_trans, quad_missal])
        
    print_("done creating list of Ks")
    return list_of_Ks


def tilt_slice_matrix(matrix, slice_shift, slice_width, tune=0):
    print "tune=", tune
    invrange = matrix.shape[0] - 1 - np.arange(matrix.shape[0])
    matrix[matrix.shape[0] - slice_shift:,:slice_shift] += tune
    matrix[:slice_shift, matrix.shape[1] - slice_shift:] -= tune
    return np.roll(matrix[np.arange(matrix.shape[0]), circulant(invrange)[invrange]],
                          slice_shift, axis=0)[:slice_width]



def printMatrix(debugfile, M, name):
    debugfile.write("begin Matrix " + name + "\n" + str(M.shape[0]) + " " + str(M.shape[1]) + "\n")

    np.savetxt(debugfile, M, fmt="%18.10e")
    debugfile.write("\nend\n")


def bad_phase(phi):
    modphi = phi % BADPHASE
    return (modphi < MOD_POINTFIVE_LOWER or modphi > MOD_POINTFIVE_UPPER)

    
def is_small(x):
    return abs(x) < ZERO_THRESHOLD


def print_box(string):
    print "=" + " " * BOXINDENT + string + " " * (BOXLENGTH - 3 - BOXINDENT - len(string)) + "="
    
    
def print_(string, prefix=" "):
    print " " * (BOXINDENT + 1) + prefix + " " + string
    
    
def print_box_edge():
    print "= " * (BOXLENGTH / 2)

