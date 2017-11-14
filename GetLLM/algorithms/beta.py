'''
.. module: beta
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
from math import sqrt
import pandas as pd

from scipy.linalg import circulant
import Python_Classes4MAD.metaclass
import Utilities.bpm
#from GetLLM.GetLLMError import GetLLMError
import compensate_ac_effect
import os
import re
import multiprocessing
import time
from constants import PI, TWOPI
from model.accelerators.accelerator import AccExcitationMode
from Utilities import tfs_pandas
from Utilities import logging_tools

__version__ = "2017.10.b"

DEBUG = sys.flags.debug  # True with python option -d! ("python -d GetLLM.py...") (vimaier)
PRINTTIMES = False
LOGGER = logging_tools.get_logger(__name__)

#if False:
#    from Utilities.progressbar import startProgress, progress, endProgress
#else:
def startProgress(name):
    _info_("START " + name)
 
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
COT_THRESHOLD           = 1.0e6
MOD_POINTFIVE_LOWER     = PHASE_THRESHOLD           #@IgnorePep8
MOD_POINTFIVE_UPPER     = (BADPHASE - PHASE_THRESHOLD)    #@IgnorePep8
RCOND                   = 1.0e-10                    #@IgnorePep8

BOXLENGTH               = 50                        #@IgnorePep8
BOXINDENT               =  4                        #@IgnorePep8
CALCULATE_BETA_HOR = True
CALCULATE_BETA_VER = True

# ================ Column Indices :

INDX_NAME   = 0
INDX_S      = 1
INDX_BET    = 2
INDX_BETSTAT= 3
INDX_BETSYS = 4
INDX_BETERR = 5
INDX_ALF    = 6
INDX_ALFSTAT= 7
INDX_ALFSYS = 8
INDX_ALFERR = 9
INDX_CORR   = 10
INDX_BETMDL = 11
INDX_BETBEAT= 12
INDX_NCOMB  = 13

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
    '''
    class that stores information about the calculated alpha and beta values
    '''
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
    _warning_("-<-<-<-<-< INVALID UncertaintyDefinition type  '{:s}' -<-<-<-<-<-<-<-<".format(_type))
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
        self.tas = _type
        
    def settype(self, _type):
        self.type = gettype(_type)
        
    def match(self, string):
        return self.pattern.match(string)

    def to_string(self):
        return "/{:s}/ dK1={:g}, dS={:g}, dX={:g}, {:6s}".format(self.pattern.pattern, self.dK1, self.dS, self.dX,
                                                                 self.tas)
             
        
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
        _debug_("opening error definition file '{:s}'".format(filename))
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
                        _debug_("adding {} = {} to properties".format(match.group(1), match.group(2)))
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
                    _error_("loading errorfile didn't work")
                    _error_("errordefspath = {0:s}".format(filename))
                    return False
                
                _debug_("error definitions file version 1")
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
        '''
        Adds uncertainty information to twiss_full.
        
        :Sources of Errors:
            dK1:    quadrupolar field errors
            dS:     quadrupole longitudinal misalignments
            dX:     sextupole transverse misalignments
            BPMdS:  BPM longitudinal misalignments 
        '''
        
        _debug_("Start creating uncertainty information")

        # create new columns, fill MUX/Y_END and BETX/Y_END
        twiss_full.loc[:]["MUX_END"] = np.roll(twiss_full.loc[:]["MUX"], 1)
        twiss_full.loc[:]["MUY_END"] = np.roll(twiss_full.loc[:]["MUY"], 1)
        twiss_full.loc[:]["BETX_END"] = np.roll(twiss_full.loc[:]["BETX"], 1)
        twiss_full.loc[:]["BETY_END"] = np.roll(twiss_full.loc[:]["BETY"], 1)
        twiss_full.loc[:]["UNC"] = False
        twiss_full.loc[:]["dK1"] = 0
        twiss_full.loc[:]["dS"] = 0
        twiss_full.loc[:]["dX"] = 0
        twiss_full.loc[:]["BPMdS"] = 0
        
        # loop over uncertainty definitions, fill the respective columns, set UNC to true
        for reg in self.regex:
            _debug_("creating uncertainty information for RegEx {:s}".format(reg.to_string()))
            reg_mask = twiss_full.index.str.match(reg.pattern)
            twiss_full.loc[reg_mask, "dK1"] = (reg.dK1 * twiss_full.loc[reg_mask, "K1L"]) **2 # TODO change K1L --> mainfield if necessary
            twiss_full.loc[reg_mask, "dX"] = reg.dX**2
            if reg.type == IDBPM:
                twiss_full.loc[reg_mask, "BPMdS"] = reg.dS**2
            else:
                twiss_full.loc[reg_mask, "dS"] = reg.dS**2
            twiss_full.loc[reg_mask, "UNC"] = True
            
        # in case of quadrupole longitudinal misalignments, the element (DRIFT) in front of the misaligned quadrupole
        # will be used for the thin lens approximation of the misalignment
        twiss_full.loc[:]["dS"] -= np.roll(twiss_full.loc[:]["dS"], 1)
        twiss_full.loc[:]["dK1_END"] = -np.roll(twiss_full.loc[:]["dK1"], 1)
        twiss_full.loc[:]["UNC"] |= np.roll(twiss_full.loc[:]["UNC"], 1)
        
        # dump the modified twiss_full and return it to the beta calculation
        _info_("DONE creating uncertainty information")
        _debug_("dumping new twiss_full.dat")
        tfs_pandas.write_tfs(twiss_full, {}, "dump_twiss_full")
        return twiss_full[twiss_full["UNC"] == True]

#===================================================================================================
# main part
#===================================================================================================


def _write_getbeta_out(q1, q2, number_of_bpms, range_of_bpms, beta_d_phase,
                       data, rmsbbx, error_method, bpms, tfs_file, _plane_char,
                       dpp=0, dppq1=0):
    '''
    Writes the file ``getbeta<x/y>.out``. 
    
    :Parameters:
        q1, q2
            tunes
        number_of_bpms
            number of bpm combinations to keep for Monte Carlo N-BPM method
        range_of_bpms
            range of BPMs from which the combinations are taken
        beta_d_col
            no idea
        data 
            the result of the method
        rmsbbx
            RMS beta beating
        error_method
            the ID of the used error method
        tfs_file
            the tfs file to which the results will be written
            
    '''

    LOGGER.debug("Writing beta from phase results")

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
                               "NCOMBINATIONS",
                              "NFILES"])
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
                               "%le",
                                  "%le"])
    for i, row in enumerate(data):
        if not np.isnan(row[2]):
            beta_d_phase[row[0]] = [row[1], row[2], row[3], row[4]]

            tfs_file.add_table_row(list(row) + [int(bpms.iloc[i].loc["NFILES"])])


def calculate_beta_from_phase(getllm_d, twiss_d, tune_d, phase_d,
                              model, model_driven, elements, elements_centre,
                              files_dict):
    '''
    Calculates beta from phase using either the 3-BPM or N-BPM method.
    Fills the following TfsFiles:
        ``getbetax.out        getbetax_free.out        getbetax_free2.out``
        ``getbetay.out        getbetay_free.out        getbetay_free2.out``
        
    

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
    
    commonbpms_x = twiss_d.zero_dpp_unionbpms_x
    commonbpms_y = twiss_d.zero_dpp_unionbpms_y
    
    debugfile = None
    if getllm_d.nprocesses == -1:
        getllm_d.nprocesses = multiprocessing.cpu_count()
    getllm_d.parallel = (getllm_d.nprocesses > 0)
    #---- H plane
    _plane_char = "X"
    _box_edge_()
    _info_box_("Calculating beta from phase")
    _info_box_("Version: {0:5s}".format(__version__))

    _debug_value_box_("range of BPMs", str(getllm_d.range_of_bpms))
    _debug_value_box_("cot of phase threshold", "{:g}".format(COT_THRESHOLD))
    
    if DEBUG:
        debugfilename = files_dict['getbetax.out'].s_output_path + "/getbetax.debug"
        debugfile = open(debugfilename, "w+")
        _info_("ATTENTION: DEBUG is set to true, calculation of beta functions will be done serially")
    elif getllm_d.parallel:
        _info_value_box_("parallel", "TRUE")
        _info_value_box_("number of processes", "{0:2d}".format(getllm_d.nprocesses))
    else:
        _info_value_box_("parallel", "FALSE")
    
    _debug_value_box_("quad field errors", "[YES]")
    _debug_value_box_("quad long misalignments", "[YES]")
    _debug_value_box_("sext transverse misalignments", "[YES]")
    _debug_value_box_("BPM long misalignments", "[YES]")
    _debug_value_box_("dipole K1 errors", "[ NO]")
    _debug_value_box_("analytical alpha", "[ NO]")
    
    _box_edge_()
        
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
        _info_("")
        _info_("Accelerator Error Definition", ">")
        unc.open(getllm_d.errordefspath)
        unc_elements = unc.create_errorfile(elements, elements_centre)
        error_method = METH_A_NBPM
        _info_("")
    else:
        error_method = METH_3BPM  # fall back to three BPM method. MC is not supported in this version
    
    #--- =========== HORIZONTAL
    if CALCULATE_BETA_HOR:
        if twiss_d.has_zero_dpp_x():
            
            _info_("Calculate free beta from phase for plane " + _plane_char + " (_free.out)", ">")
            if DEBUG:
                debugfile = open(files_dict['getbetax_free.out'].s_output_path + "/getbetax_free.debug", "w+")

            dataf, rms_bb, bpmsf, error_method_x = beta_from_phase(model, unc_elements, elements_centre,
                                                                  twiss_d.zero_dpp_x, commonbpms_x, phase_d.phase_advances_free_x, 'H',
                                                                  getllm_d, debugfile, error_method, tune_d.q1f, tune_d.q1mdl)
        
            beta_d.x_phase = {}  # this is no typo, beta from amplitude still has the old convention that the free file is called not free
            beta_d.x_phase['DPP'] = 0
            
            _info_("RMS Betabeat: {:6.2f} ".format(rms_bb), ">")
            
            _write_getbeta_out(tune_d.q1f, tune_d.q2f, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.x_phase,
                               dataf, rms_bb, error_method_x, commonbpms_x,
                               files_dict['getbetax_free.out'], _plane_char)
            
            
            if getllm_d.accelerator.excitation is not AccExcitationMode.FREE: 
                _info_("Calculate beta from phase for plane " + _plane_char, ">")
                data, rms_bb, bpms, error_method_x = beta_from_phase(model_driven, unc_elements, elements_centre,
                                                                   twiss_d.zero_dpp_x, commonbpms_x, phase_d.phase_advances_x, 'H',
                                                                   getllm_d, debugfile, error_method, tune_d.q1, tune_d.q1mdl)
                beta_d.x_phase_f = {}
                beta_d.x_phase_f['DPP'] = 0
                _info_("RMS Betabeat: {:6.2f} ".format(rms_bb), ">")
                tfs_file = files_dict['getbetax.out']
                
                _write_getbeta_out(tune_d.q1, tune_d.q2, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.x_phase_f,
                                   data, rms_bb, error_method_x, commonbpms_x,
                                   tfs_file, model_driven.BETX, _plane_char)

                _debug_("Skip free2 calculation")
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
            _info_("Calculate free beta from phase for plane " + _plane_char + " (_free.out)", ">")
    
            if DEBUG:
                debugfile = open(files_dict['getbetay_free.out'].s_output_path + "/getbetay_free.debug", "w+")
    

            dataf, rms_bb, bpms, error_method_y = beta_from_phase(model, unc_elements, elements_centre,
                                                               twiss_d.zero_dpp_y, commonbpms_y, phase_d.phase_advances_free_y, 'V',
                                                               getllm_d, debugfile, error_method, tune_d.q2f, tune_d.q2mdl)
            beta_d.y_phase = {}
            beta_d.y_phase['DPP'] = 0
            _info_("RMS Betabeat: {:6.2f} ".format(rms_bb), ">")
            tfs_file = files_dict['getbetay_free.out']
            
            _write_getbeta_out(tune_d.q1f, tune_d.q2f, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.y_phase,
                               dataf, rms_bb, error_method_y, commonbpms_y,
                               tfs_file, _plane_char)
            
            #-- ac to free beta
            if getllm_d.accelerator.excitation is not AccExcitationMode.FREE:
                #-- from eq
                _info_("Calculate beta from phase for plane " + _plane_char, ">")
    
                data, rms_bb, bpmsf, error_method_y = beta_from_phase(model_driven, unc_elements, elements_centre,
                                                                      twiss_d.zero_dpp_y, commonbpms_y, phase_d.phase_advances_y, 'V',
                                                                      getllm_d, debugfile, error_method, tune_d.q2, tune_d.q2mdl)
   
             
                tfs_file = files_dict['getbetay.out']
                
                _info_("RMS Betabeat: {:6.2f} ".format(rms_bb), ">")
                beta_d.y_phase_f = {}
                beta_d.y_phase_f['DPP'] = 0

                _write_getbeta_out(tune_d.q1, tune_d.q2, getllm_d.number_of_bpms, getllm_d.range_of_bpms, beta_d.y_phase_f,
                                   dataf, rms_bb, error_method_y, commonbpms_y,
                                   tfs_file, _plane_char)

#                #-- from the model
                _info_("Skip free2 calculation")
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
        
        _box_edge_()
        _info_box_("beta from phase finished")
        _info_box_("")
        _info_box_("elapsed time: {0:3.3f}s".format(elapsed))
        _box_edge_()
        if PRINTTIMES:
            timesfile = open("times.dat", "a")
            if getllm_d.parallel:
                timesfile.write("PARALLEL {0:f} {1:d} {2:d}\n".format(elapsed, 0, getllm_d.nprocesses))
            else:
                timesfile.write("SERIAL {0:f} {1:d} {2:d}\n".format(elapsed, 0, getllm_d.nprocesses))
    
    
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
    _info_('Calculating beta from amplitude')
    
    commonbpms_x = twiss_d.zero_dpp_commonbpms_x
    commonbpms_y = twiss_d.zero_dpp_commonbpms_y


    #---- H plane
    if twiss_d.has_zero_dpp_x():
        # here was bpms instead of _ which makes the following code confusing because it looks like it is a modified
        # commonbpms
        [beta_d.x_amp, rmsbbx, _, inv_jx] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_x, 'H', commonbpms_x)
        beta_d.x_amp['DPP'] = 0
        #-- Rescaling
        beta_d.x_ratio = 0
        skipped_bpmx = []
        
        # I tried to make sense of the original line and understood it as follows: 
        # out of the used BPMs we take the ones that lie in the arcs.
        arcbpms = commonbpms_x.index[getllm_d.accelerator.get_arc_bpms_mask(commonbpms_x.index)]
        
        
        for bpm in arcbpms:
            name = str.upper(bpm)  # second entry is the name
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

        for bpm in commonbpms_x.index:
            name = str.upper(bpm)
            betax_rescale[name] = [beta_d.x_ratio * beta_d.x_amp[name][0], beta_d.x_ratio * beta_d.x_amp[name][1], beta_d.x_amp[name][2]]

        tfs_file = files_dict['getampbetax.out']
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_float_descriptor("RMSbetabeat", rmsbbx)
        tfs_file.add_float_descriptor("RescalingFactor", beta_d.x_ratio)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETX", "BETXSTD", "BETXMDL", "MUXMDL", "BETXRES", "BETXSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        for name in commonbpms_x.index:
            bn1 = str.upper(name)
            bns1 = commonbpms_x.loc[name, "S"]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_x), beta_d.x_amp[bn1][0], beta_d.x_amp[bn1][1], mad_ac.BETX[mad_ac.indx[bn1]], mad_ac.MUX[mad_ac.indx[bn1]], betax_rescale[bn1][0], betax_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries)

        #-- ac to free amp beta
        if getllm_d.accelerator.excitation is not AccExcitationMode.FREE:
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
        [beta_d.y_amp, rmsbby, _, inv_jy] = beta_from_amplitude(mad_ac, twiss_d.zero_dpp_y, 'V', commonbpms_y)
        beta_d.y_amp['DPP'] = 0
        #-- Rescaling
        beta_d.y_ratio = 0
        skipped_bpmy = []
        arcbpms = commonbpms_y.index[getllm_d.accelerator.get_arc_bpms_mask(commonbpms_y.index)]
        for bpm in arcbpms:
            name = str.upper(bpm)  # second entry is the name
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

        for bpm in commonbpms_y.index:
            name = str.upper(bpm)
            betay_rescale[name] = [beta_d.y_ratio * beta_d.y_amp[name][0], beta_d.y_ratio * beta_d.y_amp[name][1], beta_d.y_amp[name][2]]

        tfs_file = files_dict['getampbetay.out']
        tfs_file.add_float_descriptor("Q1", tune_d.q1)
        tfs_file.add_float_descriptor("Q2", tune_d.q2)
        tfs_file.add_float_descriptor("RMSbetabeat", rmsbby)
        tfs_file.add_float_descriptor("RescalingFactor", beta_d.y_ratio)
        tfs_file.add_column_names(["NAME", "S", "COUNT", "BETY", "BETYSTD", "BETYMDL", "MUYMDL", "BETYRES", "BETYSTDRES"])
        tfs_file.add_column_datatypes(["%s", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le"])
        
        for name in commonbpms_y.index:
            bn1 = str.upper(name)
            bns1 = commonbpms_y.loc[name, "S"]
            list_row_entries = ['"' + bn1 + '"', bns1, len(twiss_d.zero_dpp_y), beta_d.y_amp[bn1][0], beta_d.y_amp[bn1][1], mad_ac.BETY[mad_ac.indx[bn1]], mad_ac.MUY[mad_ac.indx[bn1]], betay_rescale[bn1][0], betay_rescale[bn1][1]]
            tfs_file.add_table_row(list_row_entries)  # ac to free amp beta

        if getllm_d.accelerator.excitation is not AccExcitationMode.FREE: # from eq
            try:
                betayf, rmsbbyf, bpmsf = compensate_ac_effect.get_free_beta_from_amp_eq(mad_ac, twiss_d.zero_dpp_y, tune_d.q2, tune_d.q2f, phase_d.acphasey_ac2bpmac, 'V', getllm_d)  # Rescaling
                # OK, bpmsf output of get_free_beta_from_amp_eq is probably similar to the old commonbpms and should be replaced
                
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
            length of range of bpms used to calculate the beta function. Has to be an odd number > 7
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
                amp_i += tw_file.loc[bn1, "AMPX"]
                amp_j2.append(tw_file.loc[bn1, "AMPX"] ** 2)
                root2j_i += tw_file.loc[bn1, "PK2PK"] / 2.
            elif plane == 'V':
                amp_i += tw_file.loc[bn1, "AMPY"]
                amp_j2.append(tw_file.loc[bn1, "AMPY"] ** 2)
                root2j_i += tw_file.loc[bn1, "PK2PK"] / 2.

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


#======================================================================================================================
#---============== calculate beta and alpha using the old 3 BPM method ================================================
#======================================================================================================================
def scan_all_BPMs_sim_3bpm(madTwiss, phase, plane, getllm_d, commonbpms, debugfile, errors_method, tune, mdltune):
    '''
    Calculates beta from phase using the old 3-BPM method
    Fast 'vectorized' pandas / numpy code
    
    ``phase["MEAS"]``, ``phase["MODEL"]``, ``phase["ERRMEAS"]`` (from ``get_phases``) are of the form:
    
    +----------+----------+----------+----------+----------+
    |          |   BPM1   |   BPM2   |   BPM3   |   BPM4   | 
    +----------+----------+----------+----------+----------+
    |   BPM1   |    0     |  phi_21  |  phi_31  |  phi_41  | 
    +----------+----------+----------+----------+----------+
    |   BPM2   |  phi_12  |     0    |  phi_32  |  phi_42  | 
    +----------+----------+----------+----------+----------+
    |   BPM3   |  phi_13  |  phi_23  |    0     |  phi_43  | 
    +----------+----------+----------+----------+----------+
    
    aa ``tilt_slice_matrix(matrix, shift, slice, tune)`` brings it into the form:
    
    +-----------+--------+--------+--------+--------+
    |           |  BPM1  |  BPM2  |  BPM3  |  BPM4  |                    
    +-----------+--------+--------+--------+--------+    
    | BPM_(i-1) | phi_1n | phi_21 | phi_32 | phi_43 |
    +-----------+--------+--------+--------+--------+
    | BPM_i     |    0   |    0   |    0   |    0   |
    +-----------+--------+--------+--------+--------+     
    | BPM_(i+1) | phi_12 | phi_23 | phi_34 | phi_45 |      
    +-----------+--------+--------+--------+--------+
        
    ``cot_phase_*_shift1``:
    
    +-----------------------------+-----------------------------+-----------------------------+
    | cot(phi_1n) - cot(phi_1n-1) |  cot(phi_21) - cot(phi_2n)  |   cot(phi_32) - cot(phi_31) | 
    +-----------------------------+-----------------------------+-----------------------------+
    |         NaN                 |         NaN                 |         NaN                 |
    +-----------------------------+-----------------------------+-----------------------------+
    |         NaN                 |         NaN                 |         NaN                 |
    +-----------------------------+-----------------------------+-----------------------------+
    |  cot(phi_13) - cot(phi_12)  |  cot(phi_24) - cot(phi_23)  |   cot(phi_35) - cot(phi_34) | 
    +-----------------------------+-----------------------------+-----------------------------+
    
    for the combination xxxABBx: first row
    for the combinstion xBBAxxx: fourth row and
    for the combination xxBABxx: second row of ``cot_phase_*_shift2``
    '''
    number_commonbpms = commonbpms.shape[0]
    plane_bet = "BETX" if plane == "H" else "BETY"
    plane_alf = "ALFX" if plane == "H" else "ALFY"
    if errors_method == METH_3BPM:
        madTwiss_intersected = madTwiss.loc[commonbpms.index]
        starttime = time.time()
        
        # ====== setup the used variables =============================================================================
        # tilt phase advances in order to have the phase advances in a neighbourhood
        tilted_meas = tilt_slice_matrix(phase["MEAS"].as_matrix(), 2, 5, tune) * TWOPI
        tilted_model = tilt_slice_matrix(phase["MODEL"].as_matrix(), 2, 5, mdltune) * TWOPI
        tilted_errmeas = tilt_slice_matrix(phase["ERRMEAS"].as_matrix(), 2, 5, mdltune) * TWOPI
        
        betmdl = madTwiss_intersected.loc[:][plane_bet].as_matrix()
        alfmdl = madTwiss_intersected.loc[:][plane_alf].as_matrix()

        # ======= main part, calculate the beta and alpha function ====================================================
        
        # calculate cotangens of all the phase advances in the neighbourhood
        with np.errstate(divide='ignore'):
            cot_phase_meas = 1.0 / tan(tilted_meas)
            cot_phase_model = 1.0 / tan(tilted_model)
       
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
        
        # alpha
        alfi = (bet_frac * ( cot_phase_model[1] + cot_phase_model[3] + 2.0 * alfmdl)
                        - ( cot_phase_meas[1] + cot_phase_meas[3])) / 2.0
        
        # ======= error propagation ===================================================================================
        
        # error = sqrt( errphi_ij^2 * (d beta / dphi_ij)^2 )
        # calculate sin(phimdl_ij)
        sin_model = sin(tilted_model) 
        # calculate errphi_ij^2 / sin^2 phimdl_ij * beta
        with np.errstate(divide='ignore', invalid='ignore'):
            sin_squared_model = tilted_errmeas / np.multiply(sin_model, sin_model) * betmdl
        # square it again beacause it's used in a vector length
        sin_squared_model = np.multiply(sin_squared_model, sin_squared_model)
        
        sin_squ_model_shift1 = sin_squared_model + np.roll(sin_squared_model, -1, axis=0) / np.multiply(cot_phase_model_shift1, cot_phase_model_shift1)
        sin_squ_model_shift2 = sin_squared_model + np.roll(sin_squared_model, -2, axis=0) / np.multiply(cot_phase_model_shift2, cot_phase_model_shift2)
        beterr = np.sqrt(sin_squ_model_shift1[0] + sin_squ_model_shift1[3] + sin_squ_model_shift2[1]) / 3.0
        
        betstd=np.zeros(number_commonbpms)
        alfstd=np.zeros(number_commonbpms)
        alferr=np.zeros(number_commonbpms)
        
        # ======= print error method and return the data rows for getbetax/y.out ======================================

        _info_("Errors from " + ID_TO_METHOD[errors_method])
        
        bb = (bet_frac - 1.0) * 100.0
        
        return sqrt(np.mean(np.multiply(bb,bb))), errors_method, np.transpose([
                                     commonbpms.index, madTwiss_intersected.loc[:]["S"],
                                     beti, betstd, betstd, beterr,
                                     alfi, alfstd, alfstd, alfstd,
                                     alfi,
                                     betmdl,
                                     bb,
                                     alfi])
    raise GetLLMError("Monte Carlo N-BPM is not implemented.")
    

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
    
    _info_("Errors from " + ID_TO_METHOD[errors_method])
    
    # =============== setup ===========================================================================================
    # setup combinations
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

    # for fast access
    phases_meas = phase["MEAS"] * TWOPI
    phases_model = phase["MODEL"] * TWOPI
    phases_err = phase["ERRMEAS"] * TWOPI

    # setup the results matrix
    result = np.array(np.empty(madTwiss_intersected.shape[0]), 
                      dtype = [("NAME","S24"), ("S", "f8"),
                               ("BET", "f8"), ("BETSTAT", "f8"), ("BETSYS", "f8"), ("BETERR", "f8"),
                               ("ALF", "f8"), ("ALFSTAT", "f8"), ("ALFSYS", "f8"), ("ALFERR", "f8"),
                               ("CORR", "f8"), ("BETMDL", "f8"), ("BETBEAT", "f8"), ("NCOMB", "i4")])
#                      dtype = "S24, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, i4")

    # ==========
    # define functions in a function -- python witchcraft, burn it!!!!! 
    def collect(row):
#        if row[11]:
            result[row[0]]= row[1:]
        
    def collectblock(block):
        for row in block:
#            if row[11]:
            result[row[0]] = row[1:]
            
     # =============== calculate the betas ============================================================================
           
    st = time.time()
    if getllm_d.parallel:
        
        # setup thread pool and data chunks
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
                             
        # calculate the last, incomplete chunk
        pool.apply_async(scan_several_BPMs_withsystematicerrors,
                         (madTwiss_intersected, madElements,
                          phases_meas, phases_model,
                          plane, getllm_d.range_of_bpms, commonbpms, debugfile,
                          n * chunksize, len(commonbpms), BBA_combo, ABB_combo, BAB_combo,
                          tune, mdltune),
                         callback=collectblock)
                         
        # wait for all the threads to finish and join the results
        pool.close()
        pool.join()
    else:  # not parallel
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
    
    _debug_("time elapsed = {0:3.3f}".format(et - st))
    
    rmsbb = sqrt(np.mean(np.multiply(result["BETBEAT"], result["BETBEAT"])))
    
    #result["BETERR"] *= np.sqrt(commonbpms.loc[:, "NFILES"])
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
    

def scan_one_BPM_withsystematicerrors(madTwiss, madElements,
                                      phases_meas, phases_err,
                                      plane, range_of_bpms, commonbpms, debugfile,
                                      Index, BBA_combo, ABB_combo, BAB_combo,
                                      tune, mdltune):
    '''
    Scans the range of BPMs in order to get the final value for one BPM in the lattice
    
    :Parameters:
        'madTwiss':tfs_pandas
            The model twiss table, contains all the BPMs. Has to be already intersected with the common BPMs:
            
           +-------+----+-------+-----+-----+-------+
           |       | S  | BETX  | MUX | ... | BETY  |
           +-------+----+-------+-----+-----+-------+
           |  BPM1 | s1 | beta1 | mu1 | ... | bety1 |
           +-------+----+-------+-----+-----+-------+
           |  BPM2 | s2 | beta2 | mu2 | ... | bety2 |
           +-------+----+-------+-----+-----+-------+
           |  ...  |    |       |     |     |       |
           +-------+----+-------+-----+-----+-------+
           |  BPMn | sn | betan | mun | ... | betyn |
           +-------+----+-------+-----+-----+-------+

        'madElements':tfs_pandas
            Twiss table of all elements with known uncertainties. The keys are already intersected with the list of
            common BPMs:
            
           +------+--------+-----------+---------+-----+--------+---------+--------+
           |      | S      | BETX      | MUX     | ... | dK1    | dS      | BPMdS  |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | BPM1 | s1     | beta1     | mu1     | ... | 0      | 0       | 1e-3   |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | MQ1  | s(El1) | beta(El1) | mu(El1) | ... | 1e-4   | 0       | 0      |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | MQ2  | s(El2) | beta(El2) | mu(El2) | ... | 2e-4   | 0       | 0      |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | MS1  | s(El3) | beta(El3) | mu(El3) | ... | 0      | 1.0e-3  | 0      |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | BPM2 | s2     | beta2     | mu2     | ... | 0      | 0       | 1e-3   |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | ...  | ...    | ...       | ...     | ... | ...    | ...     | ...    |
           +------+--------+-----------+---------+-----+--------+---------+--------+
           | BPMn | sn     | betan     | mun     | ... | 0      | 0       | 1e-3   |
           +------+--------+-----------+---------+-----+--------+---------+--------+
             
             for later distinction we denote the row index of madElements by <l>:
                 s<1> := s1
                 s<2> := s(El1)
                 s<3> := s(El2)
                 etc.
             
        'phases_meas/phases_model':numpy.ndarray
            matrix of the cotanges of the meas/model phase advances.
            
                  | BPM1   | BPM2   | ... | BPMn 
            ------+--------+--------+-----+-----
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
            
         - for this we need cot(phi_(i)(i-m) ... cot(phi_(i)(i+m)) and idem for the model phases 
            
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
         - for this we need sin^2(phimdl_(i-m) - phimdl_(i-m+<1>)), sin^2(phimdl_(i-m) - phimdl_(i-m+<2>)), ...
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
    used_bpms   = -2                     #@IgnorePep8
    
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
        if Index == len_bpms_total - 1:
            __index = 0
        else:
            __index = Index
        outerMeasPhaseAdv = pd.concat((
                phases_meas.iloc[__index, indx_first:],
                phases_meas.iloc[__index, :indx_last % len_bpms_total] + tune * TWOPI))
        outerMeasErr = pd.concat((
                phases_err.iloc[__index, indx_first:],
                phases_err.iloc[__index, :indx_last % len_bpms_total]))
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
        
    outerMeasErr = np.multiply(outerMeasErr, outerMeasErr)

    outerElPhAdv = (outerElmtsPh[:, np.newaxis] - outerMdlPh[np.newaxis, :])
    outerElK2 = outerElmts.loc[:, "K2L"].as_matrix()
    indx_el_probed = outerElmts.index.get_loc(probed_bpm_name)
    outerElmtsBet = outerElmts.loc[:][bet_column].as_matrix()
      
    with np.errstate(divide='ignore'):
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
    diag = np.concatenate((outerMeasErr.as_matrix(), outerElmts.loc[:]["dK1"], outerElmts.loc[:]["dX"]))
    mask = diag != 0
    
    T_Beta = np.zeros((len(betas),
                   len(diag) ))
    
    M = np.diag(diag[mask])
    
#    print cot_model[m]
#    print cot_model
    
    for i, combo in enumerate(BBA_combo):
        ix = combo[0] + m
        iy = combo[1] + m
        
        # remove bad combination
        if (abs(cot_model[ix]) > COT_THRESHOLD or 
            abs(cot_model[iy]) > COT_THRESHOLD or
            abs(cot_meas[ix]) > COT_THRESHOLD or
            abs(cot_meas[iy]) > COT_THRESHOLD ):
            #or abs(cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
            beta_mask[i] = False
            continue
        beta_mask[i] = True
        
        # calculate beta
        denom = (cot_model[ix] - cot_model[iy]) / betmdl1
        betas[i] = (cot_meas[ix] - cot_meas[iy]) / denom
        
        # slice       
        xloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[ix])
        yloc = outerElmts.index.get_loc(outerMeasPhaseAdv.index[iy])
        
        # get betas and sin for the elements in the slice
        elementPh_XA = sin_squared_elements[xloc:indx_el_probed, ix]
        elementPh_YA = sin_squared_elements[yloc:indx_el_probed, iy]
        elementBet_XA = outerElmtsBet[xloc:indx_el_probed]
        elementBet_YA = outerElmtsBet[yloc:indx_el_probed]
        elementK2_XA = outerElK2[xloc:indx_el_probed]
        elementK2_YA = outerElK2[yloc:indx_el_probed]
        denom_sinx = sin_squared_elements[xloc, m]
        denom_siny = sin_squared_elements[yloc, m]
        
        # apply phase uncertainty
        T_Beta[i][ix] = -1.0 / (denom_sinx * denom)
        T_Beta[i][iy] = 1.0 / (denom_siny * denom)

        # apply quadrupolar field uncertainty (quadrupole longitudinal misalignment already included)
        
        bet_sin_ix = elementPh_XA * elementBet_XA / (denom_sinx * denom)
        bet_sin_iy = elementPh_YA * elementBet_YA / (denom_siny * denom)
        
        T_Beta[i][xloc+range_of_bpms:indx_el_probed+range_of_bpms] += bet_sin_ix
        T_Beta[i][yloc+range_of_bpms:indx_el_probed+range_of_bpms] -= bet_sin_iy
        
        y_offset = range_of_bpms + len(outerElmts)
        
        # apply sextupole transverse misalignment
        T_Beta[i][xloc + y_offset : indx_el_probed + y_offset] += elementK2_XA * bet_sin_ix
        T_Beta[i][yloc + y_offset : indx_el_probed + y_offset] -= elementK2_YA * bet_sin_iy
        
    offset_i = len(BBA_combo)
        
    for i, combo in enumerate(BAB_combo):
        ix = combo[0] + m
        iy = combo[1] + m
        
        if (abs(cot_model[ix]) > COT_THRESHOLD or 
            abs(cot_model[iy]) > COT_THRESHOLD or
            abs(cot_meas[ix]) > COT_THRESHOLD or
            abs(cot_meas[iy]) > COT_THRESHOLD ):
            #or abs(cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
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
        elementK2_XA = outerElK2[xloc:indx_el_probed]
        elementK2_YA = outerElK2[indx_el_probed:yloc]
        denom_sinx = sin_squared_elements[xloc, m]
        denom_siny = sin_squared_elements[yloc, m]
        
        T_Beta[i + offset_i][ix] = -1.0 / (denom_sinx * denom)
        T_Beta[i + offset_i][iy] = 1.0 / (denom_siny * denom)

        
        T_Beta[i + offset_i, xloc+range_of_bpms:indx_el_probed+range_of_bpms] += elementPh_XA * elementBet_XA / (denom_sinx * denom)
        T_Beta[i + offset_i, indx_el_probed+range_of_bpms:yloc+range_of_bpms] += elementPh_YA * elementBet_YA / (denom_siny * denom)
        
        y_offset = range_of_bpms + len(outerElmts)
        
        # apply sextupole transverse misalignment
        T_Beta[i + offset_i][xloc + y_offset : indx_el_probed + y_offset] += elementPh_XA * elementBet_XA * elementK2_XA / (denom_sinx * denom)
        T_Beta[i + offset_i][indx_el_probed + y_offset : yloc + y_offset] -= elementPh_YA * elementBet_YA * elementK2_YA / (denom_siny * denom)


    offset_i += len(BAB_combo)
        
    for i, combo in enumerate(ABB_combo):
        ix = combo[0] + m 
        iy = combo[1] + range_of_bpms / 2
        
        if (abs(cot_model[ix]) > COT_THRESHOLD or 
            abs(cot_model[iy]) > COT_THRESHOLD or
            abs(cot_meas[ix]) > COT_THRESHOLD or
            abs(cot_meas[iy]) > COT_THRESHOLD ):
            #or abs(cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
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
        elementK2_XA = outerElK2[indx_el_probed:xloc]
        elementK2_YA = outerElK2[indx_el_probed:yloc]
        denom_sinx = sin_squared_elements[xloc, m]
        denom_siny = sin_squared_elements[yloc, m]
        
        T_Beta[i + offset_i][ix] = -1.0 / (denom_sinx * denom)
        T_Beta[i + offset_i][iy] = 1.0 / (denom_siny * denom)
      
        T_Beta[i + offset_i][indx_el_probed+range_of_bpms:xloc+range_of_bpms] -= elementPh_XA * elementBet_XA / (denom_sinx * denom)
        T_Beta[i + offset_i][indx_el_probed+range_of_bpms:yloc+range_of_bpms] += elementPh_YA * elementBet_YA / (denom_siny * denom)
        
        y_offset = range_of_bpms + len(outerElmts)
        
        # apply sextupole transverse misalignment
        T_Beta[i + offset_i][indx_el_probed + y_offset : xloc + y_offset] += 2.0 *elementPh_XA * elementBet_XA * elementK2_XA / (denom_sinx * denom)
        T_Beta[i + offset_i][indx_el_probed + y_offset : yloc + y_offset] -= 2.0 *elementPh_YA * elementBet_YA * elementK2_YA / (denom_siny * denom)
    
    
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
        if VBeta_inv_sum == 0:
            raise ValueError
        beterr = math.sqrt(float(np.dot(np.transpose(w), np.dot(V_Beta, w)) / VBeta_inv_sum ** 2))
        beti = float(np.dot(np.transpose(w), betas) / VBeta_inv_sum)
        used_bpms = len(w)
    except ValueError:
        _error_("ValueError")
    
    return (Index, probed_bpm_name, s,
                beti, betstat, betsys, beterr,
                alfi, alfstat, alfsys, alferr,
                .0, betmdl1, (beti - betmdl1) / betmdl1 * 100.0,
                len(betas))
    
 
#===================================================================================================
#--- ac-dipole stuff
#===================================================================================================


def _get_free_beta(modelfree, modelac, data, bpms, plane):  # to check "+"

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

    _debug_("Calculating free beta from amplitude using model")

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
#
#def create_listofKs(commonbpms, errorfile, el, getllm_d):
#    list_of_Ks = []
#    
#    #---- create list of Ks, i.e. assign to each BPM a vector with all the errore lements that come after the bpm
#    # and their respective errors
#    # and model phases so that list_of_Ks[n] yields the error elements between BPM[n] and BPM[n+1]
#    # update 2016-07-28: list_of_Ks[n][k], n: BPM number, k=0: quadrupole field errors,
#    # k=1: transversal sextupole missalignments
#    # k=2: longitudinal quadrupole missalignments
#   
#    for n in range(len(commonbpms) + getllm_d.range_of_bpms + 1):
#        index_n = errorfile.indx[commonbpms[n % len(commonbpms)][1]]
#        index_nplus1 = errorfile.indx[commonbpms[(n + 1) % len(commonbpms)][1]]
#              
#        quad_fields = []
#        sext_trans = []
#        quad_missal = []
#
#        if index_n < index_nplus1:
#            for i in range(index_n + 1, index_nplus1):
#                
#                if el[i, DK1_INDEX] != 0:
#                    quad_fields.append(i)
#                if el[i, DX_INDEX] != 0:
#                    sext_trans.append(i)
#                if el[i, DS_INDEX] != 0:
#                    quad_missal.append(i)
#                    
#        else:
#            for i in range(index_n + 1, errorfile.size):
#                if el[i, DK1_INDEX] != 0:
#                    quad_fields.append(i)
#                if el[i, DX_INDEX] != 0:
#                    sext_trans.append(i)
#                if el[i, DS_INDEX] != 0:
#                    quad_missal.append(i)
#                    
#            for i in range(index_nplus1):  # ums Eck
#                if el[i, DK1_INDEX] != 0:
#                    quad_fields.append(i)
#                if el[i, DX_INDEX] != 0:
#                    sext_trans.append(i)
#                if el[i, DS_INDEX] != 0:
#                    quad_missal.append(i)
#                    
#        list_of_Ks.append([quad_fields, sext_trans, quad_missal])
#        
#    print_("done creating list of Ks")
#    return list_of_Ks


def tilt_slice_matrix(matrix, slice_shift, slice_width, tune=0):
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

# ====================================================================================================================
# =================== LOGGING stuff ==================================================================================
# ====================================================================================================================

logger_box_value_format = "=    {:.<36s}{:<24s} ="
logger_box_format = "=    {:<60s} ="
logger_boxedge_format = "= " * 28
logger_info = "    {:s} {:<s}"


def _info_box_(string):
    LOGGER.info(" " + logger_box_format.format(string))

def _info_value_box_(key, value):
    LOGGER.info(" " + logger_box_value_format.format(key, value))

def _debug_value_box_(key, value):
    LOGGER.debug(logger_box_value_format.format(key, value))

def _info_(string, prefix=" "):
    LOGGER.info(logger_info.format(prefix, string))

def _debug_(string):
    LOGGER.debug(logger_info.format(" ", string))

def _box_edge_():
    LOGGER.info(" " + logger_boxedge_format)

def _error_(message):
    LOGGER.error(">>>" + message)

def _warning_(message):
    LOGGER.warning(" " + message)
