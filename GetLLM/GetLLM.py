import GetLLM
r'''
.. module: GetLLM.GetLLM

Created on 11/09/09

:author: Glenn Vanbavinckhove  (gvanbavi@cern.ch)

:version: 3.00dev


GetLLM calculates a large collection of optics functions and parameters at the BPMs using the output from DRIVE.

The GetLLM output is distributed over different files according to the observable and the transverse plane:
 - getCO*; files containing the closed orbit at all BPMs, one file per plane.
 - getIP*; files containing beta and phase advance functions at the IPs.
 - getampbeta*; files containing the beta function as compputed from amplitude, one file per plane (in case of ACdipole free and driven data is also computed).
 - getbeta*; files containing the beta function as computed from phase advances between 3 BPMs. There is one file per plane (in case of AC-dipole free and driven data is computed).
 - getchi*; files containing the $\chi$ terms.
 - getcouple*; files containing the coupling resonances $f_{1001}$  and $f_{0101}$ (in case of AC-dipole free and driven data is computed).
 - getsex*; files containing the sextupolar resonance driving terms.
 - getphase*; files containing the phase advances between BPMs, one file per plane (in case of AC-dipole free and driven data is computed).
 - getphasetot*; files containing the total phase-advance using the first BPM as reference, one file per plane (in case of ACdipole free and driven data is computed).
 - getkick*; files containing the kick amplitudes and the linear invariants as derived from peak-to-peak values and spectral lines amplitudes.
 - getD*; files containing the dispersion, one per plane ( if off-momentum data is acquired).
 - getNDx*; files containing the normalized horizontal dispersion ( if off-momentum data is acquired).


Usage1::

    >pythonafs ../GetLLM.py -m ../../MODEL/SPS/twiss.dat -f ../../MODEL/SPS/SimulatedData/ALLBPMs.3 -o ./

Usage2::

    >pythonafs ../GetLLM.py -m ../../MODEL/SPS/twiss.dat -d mydictionary.py -f 37gev270amp2_12.sdds.new -o ./


Change history::

        git log GetLLM.py


        --- STRUCTURE ---

_parse_args()-function
    _parse_args
main()-function
    main
helper-functions
    _intial_setup       note the missing i in the name!
    _create_tfs_files   
    _analyse_src_files
    _check_bpm_compatibility
    _calculate_orbit
    _phase_and_beta_for_non_zero_dpp
    _calculate_getsextupoles
    _calculate_kick
    _get_calibrated_amplitudes
    _copy_calibration_files
helper-classes
    _GetllmData
        __init__
        set_outputpath
        set_bpmu_and_cut_for_closed_orbit
    _TwissData
        __init__
        has_zero_dpp_x
        has_non_zero_dpp_x
        has_zero_dpp_y
        has_non_zero_dpp_y
        has_no_input_files
    _TuneData
        __init__
        initialize_tunes
main invocation
    _start
    call _start()
'''
import os
import sys
import traceback
import datetime
import math
import re
import argparse
from collections import OrderedDict
import __init__  # @UnusedImport init will include paths
import Python_Classes4MAD.metaclass
from tfs_utils_getllm import GetllmTfsFile, GetllmPandasTfs
from GetLLMError import GetLLMError, CriticalGetLLMError
import algorithms.helper
import algorithms.phase
import algorithms.beta
import algorithms.compensate_excitation
import algorithms.dispersion
import algorithms.coupling
import algorithms.resonant_driving_terms
import algorithms.interaction_point
import algorithms.chi_terms
import algorithms.lobster
import algorithms.kick
import utils.iotools
from model import manager, creator
from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas
from utils import logging_tools
from utils.entrypoint import ArgumentError
import pandas as pd
from time import time

import copy
import numpy as np


####
#######
#########
VERSION = 'V3.0.0 Dev'
#########
#######
####
DEBUG = sys.flags.debug  # True with python option -d! ("python -d GetLLM.py...") (vimaier)
LOGGER = logging_tools.get_logger(__name__)
BREAK = False

# default arguments:

ACCEL       = "LHCB1"   #@IgnorePep8
DICT        = "0"       #@IgnorePep8
MODELTWISS  = "0"       #@IgnorePep8
FILES       = "0"       #@IgnorePep8
COCUT       = 4000      #@IgnorePep8
OUTPATH     = "./"      #@IgnorePep8
NBCPL       = 2         #@IgnorePep8
NONLINEAR   = False     #@IgnorePep8
TBTANA      = "SUSSIX"  #@IgnorePep8
BPMUNIT     = "um"      #@IgnorePep8
LHCPHASE    = "0"       #@IgnorePep8
BBTHRESH    = "0.15"    #@IgnorePep8
ERRTHRESH   = "0.15"    #@IgnorePep8
NUMBER_OF_BPMS  = 10    #@IgnorePep8
RANGE_OF_BPMS   = 11    #@IgnorePep8
AVERAGE_TUNE    = 1     #@IgnorePep8
CALIBRATION     = None  #@IgnorePep8
ERRORDEFS       = None  #@IgnorePep8
NPROCESSES      = 16    #@IgnorePep8
ONLYCOUPLING    = 0     #@IgnorePep8
USE_ONLY_THREE_BPMS_FOR_BETA_FROM_PHASE   = 0    #@IgnorePep8
DPP_TOLERANCE = 0.0001
UNION       = 1

BAD_BPMS_hor = ["BPM.23L6.B1", "BPM.22R8.B1", "BPMYB.4L2.B1", "BPMSX.4L2.B1", "BPMYB.4L2.B1", "BPM.14L4.B1", "BPM23L6.B1",
                "BPM.16R3.B1", "BPM20L2.B2", "BPM6L1.B2", "BPM24R2.B2", "BPM.16L5.B2", "BPM.12R4.B1", "BPMSW.1L2.B1",
                "BPM.25L5.B1"]
BAD_BPMS_ver = ["BPMSX.4L2.B1", "BPM.22R8.B1", "BPM.23L6.B1", "BPMSX.4L2.B1", "BPM.31R2.B2", "BPM20L2.B2",
                "BPM.12R4.B1"]
BAD_BPMS_hor=[]
BAD_BPMS_ver=[]

# DEBUGGING
def print_time(index, t):
    LOGGER.debug(":::  GetLLM time  >>>>>>>>>> {:8.3f} s".format(t))


#---------------------------------------------------------------------------------------------------
# _parse_args()-function
#---------------------------------------------------------------------------------------------------
def _parse_args(start_args=sys.argv[1:]):
    ''' Parses command line arguments. '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files",
                    help="Files from analysis, separated by comma",
                    metavar="FILES", default=FILES, dest="files")
    parser.add_argument("-o", "--output",
                    help="Output Path",
                    metavar="OUT", default=OUTPATH, dest="output")
    parser.add_argument("-c", "--cocut",
                    help="Cut for closed orbit measurement [um]",
                    metavar="COCUT", default=COCUT, dest="COcut")
    parser.add_argument("-n", "--nbcpl",
                    help="Analysis option for coupling, 1 bpm or 2 bpms",
                    metavar="NBCPL", default=NBCPL, dest="NBcpl")
    parser.add_argument("--nonlinear",
                    help="Run the RDT analysis",
                    metavar="NONLINEAR", default=NONLINEAR, dest="nonlinear")
    parser.add_argument("-b", "--bpmu",
                    help="BPMunit: um, mm, cm, m (default um)",
                    metavar="BPMUNIT", default=BPMUNIT, dest="BPMUNIT")
    parser.add_argument("-p", "--lhcphase",
                    help="Compensate phase shifts by tunes for the LHC experiment data, off=0(default)/on=1",
                    metavar="LHCPHASE", default=LHCPHASE, dest="lhcphase")
    parser.add_argument("-k", "--bbthreshold",
                    help="Set beta-beating threshold for action calculations, default = 0.15",
                    metavar="BBTHRESH", default=BBTHRESH, dest="bbthreshold")
    parser.add_argument("-e", "--errthreshold",
                    help="Set beta relative uncertainty threshold for action calculations, default = 0.1",
                    metavar="ERRTHRESH", default=ERRTHRESH, dest="errthreshold")
    parser.add_argument("-g", "--threebpm",
                    help="Forces to use the 3 BPM method, yes=1/no=0, default = 0",
                    metavar="USE_ONLY_THREE_BPMS_FOR_BETA_FROM_PHASE", type=int,
                    default=USE_ONLY_THREE_BPMS_FOR_BETA_FROM_PHASE, dest="use_only_three_bpms_for_beta_from_phase")
    parser.add_argument("-j", "--numbpm",
                    help="Number of different BPM combinations for beta-calculation, default = 10",
                    metavar="NUMBER_OF_BPMS", type=int, default=NUMBER_OF_BPMS, dest="number_of_bpms")
    parser.add_argument("-i", "--range",
                    help="Range of BPM for beta-calculation (>=3 and odd), default = 11",
                    metavar="RANGE_OF_BPMS", type=int,  default=RANGE_OF_BPMS, dest="range_of_bpms")
    parser.add_argument("-r", "--average_tune",
                    help="Set to 1 to use average tune for all BPMs instead of specific for each one.",
                    metavar="AVERAGE_TUNE", type=int, default=AVERAGE_TUNE, dest="use_average")
    parser.add_argument("--calibration",
                    help="Path to the directory where the calibration files (calibration_x.out, calibration_y.out) are stored.",
                    metavar="CALIBRATION", default=CALIBRATION, dest="calibration_dir_path")
    parser.add_argument("-q", "--nprocesses", default=NPROCESSES, dest="nprocesses",
                      metavar="NPROCESSES", type=int,
                      help="""Sets the number of processes used. -1: take the number of CPUs 0: run serially >1: take the
                        specified number. default = {0:d}""".format(NPROCESSES))
    parser.add_argument("-u", "--union", type=int, dest="union", default=UNION,
                      help="""If given, the phase per BPM is calculated for each BPM with at least 3 valid measurements.
                        Otherwise (default) calculates the phase only for the intersection of all measurements.""")
    # The following is just for backwards compatitibility and should vanish once the GUI knows the new arguments
    parser.add_argument("-t", "--tbtana", # remove !
                    help="Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA",
                    metavar="TBTANA", default=TBTANA, dest="TBTana")
    parser.add_argument("--errordefs", # remove !
                    help="Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA",
                    metavar="TBTANA", default=None, dest="errordefspath")
    parser.add_argument("--coupling", default=ONLYCOUPLING, dest="onlycoupling",
                      metavar="ONLYCOUPLING", type=int,
                      help="When enabled only coupling is calculated. ")

    options, acc_args = parser.parse_known_args(args=start_args)


    try:  # for transition
        accelerator = manager.get_accel_instance(
            acc_args
        )
    except SystemExit:
        LOGGER.warning("Loading accelerator class failed.")
        LOGGER.warning("Suppose that the old command line arguments were used")
        LOGGER.warning("and guess what the right accelerator might be")
        accParser = argparse.ArgumentParser()
        accParser.add_argument("-a", "--accel", dest="accel")
        accParser.add_argument("-m", "--model", dest="modelpath")

        acc_opts, rest_args = accParser.parse_known_args(acc_args)
        if os.path.isdir(acc_opts.modelpath):
            model_dir = acc_opts.modelpath
        else:
            model_dir = os.path.dirname(acc_opts.modelpath)


        if acc_opts.accel == "LHCB1":
            accelerator = manager.get_accel_instance(
                accel="lhc",
                beam=1,
                lhc_mode="lhc_runII_2017",
                model_dir=model_dir
            )
        elif acc_opts.accel == "LHCB2":
            accelerator = manager.get_accel_instance(
                accel="lhc",
                beam=2,
                lhc_mode="lhc_runII_2017",
                model_dir=model_dir
            )
        else:
            LOGGER.error("Cannot understand accelerator arguments")
            LOGGER.error("Given are:")
            LOGGER.error(sys.argv)
            raise SyntaxError("Could not parse arguments")
    except IOError:
        accParser = argparse.ArgumentParser()
        accParser.add_argument("-a", "--accel", dest="accel")
        accParser.add_argument("-m", "--model", dest="modelpath")

        acc_opts, rest_args = accParser.parse_known_args(acc_args)
        if os.path.isdir(acc_opts.modelpath):
            model_dir = acc_opts.modelpath
        else:
            model_dir = os.path.dirname(acc_opts.modelpath)

        accelerator = manager.get_accel_instance(
            accel=acc_opts.accel,
            model_dir=model_dir
        )

    return options, accelerator


#---------------------------------------------------------------------------------------------------
# main()-function
#---------------------------------------------------------------------------------------------------
def main(accelerator,
         model_dir,
         outputpath,
         files_to_analyse,
         lhcphase=LHCPHASE,
         bpmu=BPMUNIT,
         cocut=COCUT,
         nbcpl=NBCPL,
         nonlinear=NONLINEAR,
         tbtana=TBTANA,
         bbthreshold=BBTHRESH,
         errthreshold=ERRTHRESH,
         use_only_three_bpms_for_beta_from_phase=USE_ONLY_THREE_BPMS_FOR_BETA_FROM_PHASE,
         number_of_bpms=NUMBER_OF_BPMS,
         range_of_bpms=RANGE_OF_BPMS,
         use_average=AVERAGE_TUNE,
         calibration_dir_path=CALIBRATION,
         nprocesses=NPROCESSES,
         onlycoupling=ONLYCOUPLING,
        union=False):
    '''
    GetLLM main function.

    :param string outputpath: The output path to store results
    :param string files_to_analyse: List of files, comma separated string.
    :param string model_filename: Path and filename of model to be used
    :param string dict_file: Name of the script which will be executed. Should store dictionary with
                    mappings of BPM names.
    :param string accel: Type of accelerator. LHCB1, LHCB2, LHCB4, RHIC, SPS
    :param string lhcphase: "0" or "1" -- Compensate phase shifts by tunes for the LHC experiment data,
                             off=0(default)/on=1
    :param string BPMU: BPMunit: "um", "mm", "cm", "m" (default "um")
    :param int COcut: Cut for closed orbit measurement [um]
    :param int NBcpl: For selecting the coupling measurement method 1 bpm or 2 bpms
    :param string TBTana: Turn-by-turn data analysis algorithm: SUSSIX, SVD or HA
    :param string bbthreshold: Threshold for _calculate_kick for (beta_d-beta_m)/beta_m
    :param string errthreshold: Threshold for calculate_kick for sigma(beta_d)/beta_d
    :param int use_only_three_bpms_for_beta_from_phase:
    :param int number_of_bpms: Number of BPM-combos for beta from phase
    :param int range_of_bpms: Range of BPMs for beta from phase
    :param int use_average: Uses AVG_MUX and AVG_MUY in _analyse_src_files if 1
    :returns: int  -- 0 if the function run successfully otherwise !=0.

    '''
    return_code = 0

    LOGGER.info("Starting GetLLM " + VERSION)
    global __getllm_starttime
    __getllm_starttime = time()


    use_average = (use_average == 1)
    use_only_three_bpms_for_beta_from_phase = (use_only_three_bpms_for_beta_from_phase == 1)
    union = (union == 1)

    # The following objects stores multiple variables for GetLLM to avoid having many local
    # variables. Identifiers supposed to be as short as possible.
    # --vimaier
    twiss_d = _TwissData()
    tune_d = _TuneData()
    getllm_d = _GetllmData()
    getllm_d.set_outputpath(outputpath)
    getllm_d.set_bpmu_and_cut_for_closed_orbit(cocut, bpmu)
    getllm_d.lhc_phase = lhcphase
    getllm_d.num_bpms_for_coupling = nbcpl
    getllm_d.number_of_bpms = number_of_bpms
    getllm_d.range_of_bpms = range_of_bpms
    getllm_d.use_only_three_bpms_for_beta_from_phase = use_only_three_bpms_for_beta_from_phase
    getllm_d.accelerator = accelerator
    getllm_d.nprocesses = nprocesses
    getllm_d.onlycoupling = onlycoupling
    getllm_d.union = union
    
    tune_d.q1mdl = accelerator.get_model_tfs().headers["Q1"]
    tune_d.q2mdl = accelerator.get_model_tfs().headers["Q2"]
    
    logging_tools.add_module_handler(logging_tools.file_handler(os.path.join(outputpath, "getllm.log")))
    
    # Setup

    if sys.flags.debug:
        LOGGER.info("     DEBUG ON")

    # Creates the output files dictionary
    files_dict, header_dict = _create_tfs_files(getllm_d, os.path.join(model_dir, "twiss.dat"), nonlinear)
    # Copy calibration files calibration_x/y.out from calibration_dir_path to outputpath
    calibration_twiss = _copy_calibration_files(outputpath, calibration_dir_path)
    print_time("BEFORE_ANALYSE_SRC", time() - __getllm_starttime)

    twisses_x, twisses_y = _analyse_src_files1(files_to_analyse)

    twiss_d = _analyse_src_files11(getllm_d, twiss_d, use_average,
                                    calibration_twiss, accelerator.get_model_tfs(), twisses_x, twisses_y)

    files_dict = _analyse_src_files2(twiss_d, nonlinear, files_dict, twisses_x, twisses_y)

    # Construct pseudo-double plane BPMs
    # TODO This should be in accelerator class
#    if (accelerator.__name__ == "SPS" or "RHIC" in accelerator.__name__) and twiss_d.has_zero_dpp_x() and twiss_d.has_zero_dpp_y():
#        [pseudo_list_x, pseudo_list_y] = algorithms.helper.pseudo_double_plane_monitors(accelerator.get_model_tfs(), twiss_d.zero_dpp_x, twiss_d.zero_dpp_y, bpm_dictionary)
#    else:
#        # Initialize variables otherwise calculate_coupling would raise an exception(vimaier)
#        pseudo_list_x = None
#        pseudo_list_y = None
    print_time("BEFORE_PHASE", time() - __getllm_starttime)

    #-------- START Phase for beta calculation with best knowledge model in ac phase compensation
    #temp_dict = copy.deepcopy(files_dict)

    try:
        phase_d_bk, tune_d = algorithms.phase.calculate_phase(
            getllm_d, twiss_d, tune_d, files_dict
        )
    except:
        _tb_()
        # if phase crashed, none of the subsequent algorithms can run. Thus
        raise CriticalGetLLMError(
            "get phase crashed. None of the following algorithms can work hence GetLLM will crash now. Good bye!"
        )
    print_time("AFTER_PHASE", time() - __getllm_starttime)

    #-------- START coupling.
    try:
        tune_d = algorithms.coupling.calculate_coupling(
            getllm_d, twiss_d, phase_d_bk, tune_d, accelerator, files_dict
        )
    except:
        _tb_()

    if getllm_d.onlycoupling == 0:
        #-------- START Beta

        try:
            beta_d, beta_driven_x, beta_free_x = algorithms.beta.calculate_beta_from_phase(
                getllm_d, twiss_d, tune_d, phase_d_bk, files_dict
            )
        except:
            _tb_()
        if use_only_three_bpms_for_beta_from_phase:
            print_time("AFTER_BETA_FROM_PHASE", time() - __getllm_starttime)
        else:
            print_time("AFTER_A_NBPM", time() - __getllm_starttime)

        # try:
        #     algorithms.lobster.get_local_observable(
        #         phase_d_bk, getllm_d.accelerator.get_model_tfs(), files_dict, tune_d.q1f)
        # except:
        #     _tb_()

        #------- START beta from amplitude
        try:
            beta_d = algorithms.beta.calculate_beta_from_amplitude(
                getllm_d, twiss_d, tune_d, phase_d_bk, beta_d,
                files_dict, accelerator
            )
        except:
            _tb_()

        # in the following functions, nothing should change, so we choose the models now
        mad_twiss = accelerator.get_model_tfs()
        mad_elements = accelerator.get_elements_tfs()
        if accelerator.excitation != AccExcitationMode.FREE:
            mad_ac = accelerator.get_driven_tfs()
        else:
            mad_ac = mad_twiss

        #-------- START IP
        try:
            algorithms.interaction_point.betastar_from_phase(getllm_d.accel, phase_d_bk, mad_twiss, files_dict)
        except:
            _tb_()

        #-------- START Dispersion
        try:
            algorithms.dispersion.calculate_orbit_and_dispersion(  # propagate here the unit
                twiss_d, tune_d, mad_twiss, header_dict, 'mm', getllm_d.cut_for_closed_orbit,
                beta_driven_x, getllm_d.outputpath)
        except:
            _tb_()

        #------ Start get Q,JX,delta
        try:
            files_dict, inv_x, inv_y = algorithms.kick.calculate_kick(
                mad_twiss, mad_ac, getllm_d, twisses_x, twisses_y, beta_d, phase_d_bk,
                getllm_d.outputpath, header_dict, files_dict)
        except:
            _tb_()
    else:
        LOGGER.info("GetLLM was only calculating coupling. Skipping the rest and returning ...")


    #-------- START RDTs
    if nonlinear:
        try:
            algorithms.resonant_driving_terms.calculate_RDTs(
                mad_twiss, getllm_d, twiss_d, phase_d_bk, tune_d, files_dict, inv_x, inv_y
            )
        except:
            _tb_()

# TODO: what does this?
#           if tbtana == "SUSSIX":
#               #------ Start getsextupoles @ Glenn Vanbavinckhove
#               files_dict = _calculate_getsextupoles(twiss_d, phase_d_bk, accelerator.get_model_tfs(), files_dict, tune_d.q1f)
#
#               #------ Start getchiterms @ Glenn Vanbavinckhove
#               files_dict = algorithms.chi_terms.calculate_chiterms(getllm_d, twiss_d, accelerator.get_model_tfs(), files_dict)

        # Write results to files in files_dict
    LOGGER.info("Writing files")
    for tfsfile in files_dict.itervalues():
        try:
            tfsfile.write_to_file(formatted=True)
        except AttributeError:
            LOGGER.error("tfsfile couldn't be written (AttributeError).")

    print_time("FINISH", time() - __getllm_starttime)
    return return_code
# END main() ---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# helper-functions
#---------------------------------------------------------------------------------------------------

def _create_tfs_files(getllm_d, model_filename, nonlinear):
    '''
    Creates the most tfs files and stores it in an dictionary whereby the key represents the file
    and the value is the corresponding GetllmTfsFile.

    :Return: dict: string --> GetllmTfsFile
            A dictionary of created GetllmTfsFile objects. Keys are the filenames and values are the
            GetllmTfsFile objects.
    '''
    # Static variable of GetllmTfsFile to save the outputfile, GetLLM version and model filename
    GetllmTfsFile.s_output_path = getllm_d.outputpath
    GetllmTfsFile.s_getllm_version = VERSION
    GetllmTfsFile.s_mad_filename = model_filename
    header_dict = OrderedDict([('GetLLMVersion', VERSION),
                               ('Command', sys.executable+" '"+"' '".join([]+sys.argv)+"'"),
                               ('CWD', os.getcwd()),
                               ('Date', datetime.datetime.today().strftime("%d. %B %Y, %H:%M:%S"))])
    files_dict = {}
    _fill_files_dict_plane(files_dict, "x", getllm_d)
    _fill_files_dict_plane(files_dict, "y", getllm_d)
    files_dict['getcouple.out'] = GetllmTfsFile('getcouple.out')
    if nonlinear:
        for rdt in algorithms.resonant_driving_terms.RDT_LIST:
            files_dict[rdt+'_line.out'] = GetllmTfsFile(rdt+'_line.out')
            files_dict[rdt+'.out'] = GetllmTfsFile(rdt+'.out')
    files_dict['getcouple_free.out'] = GetllmTfsFile('getcouple_free.out')
    files_dict['getcoupleterms.out'] = GetllmTfsFile('getcoupleterms.out')
    #if "LHC" in getllm_d.accelerator.__name__:
    files_dict['getIP.out'] = GetllmTfsFile('getIP.out')
    files_dict['getIPfromphase.out'] = GetllmTfsFile('getIPfromphase.out')
    files_dict['getIPfromphase_free.out'] = GetllmTfsFile('getIPfromphase_free.out')

    files_dict["getsex3000.out"] = GetllmTfsFile("getsex3000.out")
    files_dict['getchi3000.out'] = GetllmTfsFile('getchi3000.out')
    files_dict['getchi1010.out'] = GetllmTfsFile('getchi1010.out')
    files_dict['getkickac.out'] = GetllmTfsFile('getkickac.out')
    files_dict['getlobster.out'] = GetllmTfsFile('getlobster.out')
    files_dict['getrexter.out'] = GetllmPandasTfs('getrexter.out')
    return files_dict, header_dict
# END _create_tfs_files -----------------------------------------------------------------------------

def _fill_files_dict_plane(files_dict, plane, getllm_d):

    if getllm_d.accelerator.excitation == AccExcitationMode.FREE:
        files_dict['getphase{}_free.out'.format(plane)] = GetllmTfsFile('getphase{}.out'.format(plane))
        files_dict['getphasetot{}_free.out'.format(plane)] = GetllmTfsFile('getphasetot{}.out'.format(plane))
        files_dict['getbeta{}_free.out'.format(plane)] = GetllmPandasTfs('getbeta{}.out'.format(plane))
        files_dict['getampbeta{}.out'.format(plane)] = GetllmTfsFile('getampbeta{}.out'.format(plane))
        files_dict['getIP{}.out'.format(plane)] = GetllmTfsFile('getIP{}.out'.format(plane))
    else:
        files_dict['getphase{}.out'.format(plane)] = GetllmTfsFile('getphase{}.out'.format(plane))
        files_dict['getphasetot{}.out'.format(plane)] = GetllmTfsFile('getphasetot{}.out'.format(plane))
        files_dict['getphase{}_free.out'.format(plane)] = GetllmTfsFile('getphase{}_free.out'.format(plane))
        files_dict['getphasetot{}_free.out'.format(plane)] = GetllmTfsFile('getphasetot{}_free.out'.format(plane))
        files_dict['getbeta{}.out'.format(plane)] = GetllmPandasTfs('getbeta{}.out'.format(plane))
        files_dict['getbeta{}_free.out'.format(plane)] = GetllmPandasTfs('getbeta{}_free.out'.format(plane))
        files_dict['getampbeta{}.out'.format(plane)] = GetllmTfsFile('getampbeta{}.out'.format(plane))
        files_dict['getampbeta{}_free.out'.format(plane)] = GetllmTfsFile('getampbeta{}_free.out'.format(plane))
        files_dict['getIP{}.out'.format(plane)] = GetllmTfsFile('getIP{}.out'.format(plane))
        files_dict['getIP{}_free.out'.format(plane)] = GetllmTfsFile('getIP{}_free.out'.format(plane))


def _analyse_src_files1(files_to_analyse):

    LOGGER.debug("Start analysing source files")

    tfs_files_x, tfs_files_y = [], []
    for file_in in files_to_analyse.split(','):
        LOGGER.debug("> file: '{:s}'".format(file_in))
        # x file

        if os.path.isfile(file_in + ".linx"):
            file_x = file_in + ".linx"
        elif os.path.isfile(file_in + "_linx"):
            file_x = file_in + "_linx"

        if os.path.isfile(file_in + ".liny"):
            file_y = file_in + ".liny"
        elif os.path.isfile(file_in + "_liny"):
            file_y = file_in + "_liny"
        try:
            twiss_file_x = tfs_pandas.read_tfs(file_x).set_index("NAME")
            tfs_files_x.append(twiss_file_x)
        except IOError:
            LOGGER.warning("Cannot load file: " + file_x)
        except ValueError:
            pass  # Information printed by metaclass already

        try:
            twiss_file_y = tfs_pandas.read_tfs(file_y).set_index("NAME")
            tfs_files_y.append(twiss_file_y)
        except IOError:
            LOGGER.warning('Warning: There seems no ' + str(file_y) + ' file in the specified directory.')
        except ValueError:
            pass  # Information printed by metaclass already

    tfs_files_x = _arrange_dpp(tfs_files_x)
    tfs_files_y = _arrange_dpp(tfs_files_y)

    return tfs_files_x, tfs_files_y


def _analyse_src_files11(getllm_d, twiss_d, use_average, calibration_twiss, model, tfs_files_x,
                         tfs_files_y):

    union = getllm_d.union
    for twiss_file_x, twiss_file_y in zip(tfs_files_x, tfs_files_y):
        if use_average:
            twiss_file_x = twiss_file_x.rename(columns={"AVG_MUX": "MUX", "MUX": "MUX_OLD"})
        if calibration_twiss is not None:
            twiss_file_x["AMPX"], twiss_file_x["ERRAMPX"] = _get_calibrated_amplitudes(twiss_file_x, calibration_twiss, "X")
        try:
            dppi = float(twiss_file_x.headers["DPP"])
        except KeyError:
            dppi = 0.0
        if type(dppi) != float:
            print type(dppi)
            print >> sys.stderr, 'Warning: DPP may not be given as a number in ', twiss_file_x.filename, '...trying to forcibly cast it as a number'
            try:
                dppi = float(dppi)
                LOGGER.info('dppi = ' + dppi)
            except ValueError:
                print >> sys.stderr, 'but failing. DPP in ', twiss_file_x.filename, ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)
        if dppi == 0.0:  # abs(dppi) < DPP_THRESHOLD:
            twiss_d.zero_dpp_x.append(twiss_file_x)
        else:
            twiss_d.non_zero_dpp_x.append(twiss_file_x)

        # y file
        if use_average:
            twiss_file_y = twiss_file_y.rename(columns={"AVG_MUY": "MUY", "MUY": "MUY_OLD"})
        if calibration_twiss is not None:
            twiss_file_y.AMPY, twiss_file_y["ERRAMPY"] = _get_calibrated_amplitudes(twiss_file_y, calibration_twiss, "Y")
        try:
            dppi = float(twiss_file_y.headers["DPP"])
        except KeyError:
            dppi = 0.0
        if type(dppi) != float:
            print >> sys.stderr, 'Warning: DPP may not be given as a number in ', twiss_file_y.filename, '...trying to forcibly cast it as a number'
            try:
                dppi = float(dppi)
                print 'dppi= ', dppi
            except ValueError:
                print >> sys.stderr, 'but failing. DPP in ', twiss_file_y.filename, ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)
        if dppi == 0.0:  # abs(dppi) < DPP_THRESHOLD:
            twiss_d.zero_dpp_y.append(twiss_file_y)

        else:
            twiss_d.non_zero_dpp_y.append(twiss_file_y)

    if not twiss_d.has_zero_dpp_x():
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'

        if twiss_d.has_non_zero_dpp_x():
            twiss_d.zero_dpp_x = twiss_d.non_zero_dpp_x
            twiss_d.zero_dpp_y = twiss_d.non_zero_dpp_y
            twiss_d.non_zero_dpp_x = []
            twiss_d.non_zero_dpp_y = []

    if twiss_d.has_no_input_files():
        print >> sys.stderr, "No parsed input files"
        sys.exit(1)

    twiss_d.zero_dpp_commonbpms_x, twiss_d.zero_dpp_unionbpms_x = _get_commonbpms(twiss_d.zero_dpp_x, model, union, "H")
    twiss_d.zero_dpp_commonbpms_y, twiss_d.zero_dpp_unionbpms_y = _get_commonbpms(twiss_d.zero_dpp_y, model, union, "V")
    twiss_d.non_zero_dpp_commonbpms_x, twiss_d.non_zero_dpp_unionbpms_x = _get_commonbpms(twiss_d.non_zero_dpp_x, model,
                                                                                         union, "H")
    twiss_d.non_zero_dpp_commonbpms_y, twiss_d.non_zero_dpp_unionbpms_y = _get_commonbpms(twiss_d.non_zero_dpp_y, model,
                                                                                         union, "V")
    return twiss_d

def _analyse_src_files2(twiss_d, nonlinear, files_dict, tfs_files_x, tfs_files_y):

    for twiss_file_x, twiss_file_y in zip(tfs_files_x, tfs_files_y):
        try:
            dppi = float(twiss_file_x.headers["DPP"])
        except KeyError:
            dppi = 0.0
        if type(dppi) != float:
            print type(dppi)
            print >> sys.stderr, 'Warning: DPP may not be given as a number in ', twiss_file_x.filename, '...trying to forcibly cast it as a number'
            try:
                dppi = float(dppi)
                LOGGER.info('dppi = ' + dppi)
            except ValueError:
                print >> sys.stderr, 'but failing. DPP in ', twiss_file_x.filename, ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)
        if dppi == 0.0:  # abs(dppi) < DPP_THRESHOLD:
            files_dict['getcouple.out'].add_filename_to_getllm_header(twiss_file_x.filename[:-5])
            if nonlinear:
                for rdt in algorithms.resonant_driving_terms.RDT_LIST:
                    files_dict[rdt+'_line.out'].add_filename_to_getllm_header(twiss_file_x.filename[:-5])
                    files_dict[rdt+'.out'].add_filename_to_getllm_header(twiss_file_x.filename[:-5])
            _add_filename_to_header_for_files(
                files_dict, twiss_file_x.filename,
                [
                    'getphasex.out',
                    'getphasetotx.out',
                    'getbetax.out',
                    'getampbetax.out',
                    'getphasex_free.out',
                    'getphasex_free2.out',
                    'getphasetotx_free.out',
                    'getphasetotx_free2.out',
                    'getbetax_free.out',
                    'getbetax_free2.out',
                    'getampbetax_free.out',
                    'getampbetax_free2.out'
                ])
            _add_filename_to_header_for_files(
                files_dict, twiss_file_x.filename[:-5],
                [
                    'getIPx.out',
                    'getIPy.out',
                    'getIPfromphase.out',
                    'getIPx_free.out',
                    'getIPy_free.out',
                    'getIPx_free2.out',
                    'getIPy_free2.out',
                    'getIPfromphase_free.out',
                    'getIPfromphase_free2.out',
                    'getcouple_free.out',
                    'getcouple_free2.out'
                ])

        # y file
        try:
            dppi = float(twiss_file_y.headers["DPP"])
        except KeyError:
            dppi = 0.0
        if type(dppi) != float:
            print >> sys.stderr, 'Warning: DPP may not be given as a number in ', twiss_file_y.filename, '...trying to forcibly cast it as a number'
            try:
                dppi = float(dppi)
                print 'dppi= ', dppi
            except ValueError:
                print >> sys.stderr, 'but failing. DPP in ', twiss_file_y.filename, ' is something wrong. String? --- leaving GetLLM'
                print >> sys.stderr, traceback.format_exc()
                sys.exit(1)
        if dppi == 0.0:  # abs(dppi) < DPP_THRESHOLD:
            _add_filename_to_header_for_files(
                files_dict, twiss_file_y.filename,
                [
                    'getphasey.out',
                    'getphasetoty.out',
                    'getbetay.out',
                    'getampbetay.out',
                    'getphasey_free.out',
                    'getphasey_free2.out',
                    'getphasetoty_free.out',
                    'getphasetoty_free2.out',
                    'getbetay_free.out',
                    'getbetay_free2.out',
                    'getampbetay_free.out',
                    'getampbetay_free2.out'
                ])

    if not twiss_d.has_zero_dpp_x():
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'

        if twiss_d.has_non_zero_dpp_x():
            _add_filename_to_header_for_files(
                files_dict, "chrommode",
                [
                    'getphasex.out',
                    'getbetax.out',
                    'getampbetax.out',
                    'getcouple.out',
                    'getcouple_free.out',
                    'getcouple_free2.out',
                    'getphasey.out',
                    'getbetay_free.out',
                    'getampbetay.out'
                ])

    if twiss_d.has_no_input_files():
        print >> sys.stderr, "No parsed input files"
        sys.exit(1)

    return files_dict

def _add_filename_to_header_for_files(files_dict, filex, files_list):
    for key in files_list:
        if key in files_dict:
            files_dict[key].add_filename_to_getllm_header(filex)



def _arrange_dpp(list_of_tfs):
    '''
    Grouping of dpp-values in the given linx,liny-files and computing new values
    '''
    list_of_tfs_arranged = []
    for tfs_file in list_of_tfs:
        if "DPP" not in tfs_file.headers:
            tfs_file.headers["DPP"] = 0.0
    if len(list_of_tfs) == 1:
        only_dpp = list_of_tfs[0].headers["DPP"]
        if np.abs(only_dpp) > DPP_TOLERANCE:
            LOGGER.warn(
                'It looks like the file you are analyzing has too '
                'high momentum deviation {}. Optics parameters might '
                'be wrong.'.format(only_dpp)
            )
        list_of_tfs[0].headers["DPP"] = 0.0
        return list_of_tfs
    dpp_values = [tfs_file.DPP for tfs_file in list_of_tfs]
    if 0.0 in dpp_values:
        LOGGER.warn('Exact 0.0 found, the dp/p values are probably already grouped.')
        return list_of_tfs
    closest_to_zero = np.argmin(np.absolute(dpp_values))
    ordered_indices = np.argsort(dpp_values)
    ranges = _compute_ranges(list_of_tfs, ordered_indices)
    offset_range = _find_range_with_element(ranges, closest_to_zero)
    offset_dpps = _values_in_range(offset_range, dpp_values)
    LOGGER.debug("dp/p closest to zero is {}".format(dpp_values[closest_to_zero]))
    zero_offset = np.mean(offset_dpps)
    LOGGER.debug("Detected dpp differences, aranging as: {0}, zero offset: {1}."
                 .format(ranges, zero_offset))
    for idx in range(len(dpp_values)):
        range_to_use = _find_range_with_element(ranges, idx)
        dpps_from_range = _values_in_range(range_to_use, dpp_values)
        range_mean = np.mean(dpps_from_range)
        list_of_tfs[idx].headers["DPP"] = range_mean - zero_offset
        list_of_tfs_arranged.append(list_of_tfs[idx])
    return list_of_tfs_arranged


def _values_in_range(range_to_use, dpp_values):
    dpps_from_range = []
    for dpp_idx in range_to_use:
        dpps_from_range.append(dpp_values[dpp_idx])
    return dpps_from_range


def _find_range_with_element(ranges, element):
    range_with_element = None
    for dpp_range in ranges:
        if element in dpp_range:
            range_with_element = dpp_range
    return range_with_element


def _compute_ranges(list_of_tfs, ordered_indices):
    list_of_ranges = []
    last_range = None
    for idx in ordered_indices:
        if (list_of_ranges and
                _is_in_same_range(list_of_tfs[last_range[0]].DPP,
                                  list_of_tfs[idx].DPP)):
            last_range.append(idx)
        else:
            new_range = []
            new_range.append(idx)
            list_of_ranges.append(new_range)
            last_range = new_range
    return list_of_ranges


def _is_in_same_range(a, b):
    return b <= a + DPP_TOLERANCE and b >= a - DPP_TOLERANCE

# END _analyse_src_files ----------------------------------------------------------------------------


def _check_bpm_compatibility(twiss_d, mad_twiss):
    '''
    Checks the monitor compatibility between data and model. If a monitor will not be found in the
    model, a message will be print to sys.stderr.
    '''
    all_twiss_files = twiss_d.non_zero_dpp_x + twiss_d.zero_dpp_x + twiss_d.non_zero_dpp_y + twiss_d.zero_dpp_y
    for twiss_file in all_twiss_files:
        for bpm_name in twiss_file.NAME:
            try:
                mad_twiss.NAME[mad_twiss.indx[bpm_name]]
            except KeyError:
                try:
                    mad_twiss.NAME[mad_twiss.indx[str.upper(bpm_name)]]
                except KeyError:
                    print >> sys.stderr, 'Monitor ' + bpm_name + ' cannot be found in the model!'




def _calculate_getsextupoles(twiss_d, phase_d, mad_twiss, files_dict, q1f):
    '''
    Fills the following TfsFiles:
     - getsex3000.out

    :returns: dict string --> GetllmTfsFile -- The same instace of files_dict to indicate that the dict was extended.
    '''
    print "Calculating getsextupoles"
    # For getsex1200.out andgetsex2100.out take a look at older revisions. (vimaier)

    htot, afactor, pfactor = algorithms.helper.Getsextupole(mad_twiss, twiss_d.zero_dpp_x, phase_d.ph_x, q1f, 3, 0)

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
# END _calculate_getsextupoles ----------------------------------------------------------------------



def _get_calibrated_amplitudes(drive_file, calibration_twiss, plane):
    calibration_file = calibration_twiss[plane]
    cal_amplitudes = []
    err_cal_amplitudes = []
    for bpm_name in drive_file.index:
        drive_index = drive_file.indx[bpm_name]
        cal_amplitude = getattr(drive_file, "AMP" + plane)[drive_index]
        err_cal_amplitude = 0.
        if bpm_name in calibration_file.NAME:
            cal_index = calibration_file.indx[bpm_name]
            cal_amplitude = cal_amplitude * calibration_file.CALIBRATION[cal_index]
            err_cal_amplitude = calibration_file.ERROR_CALIBRATION[cal_index]
        cal_amplitudes.append(cal_amplitude)
        err_cal_amplitudes.append(err_cal_amplitude)
    return np.array(cal_amplitudes), np.array(err_cal_amplitudes)
# END _get_calibrated_amplitudes --------------------------------------------------------------------


def _copy_calibration_files(output_path, calibration_dir_path):
    calibration_twiss = {}
    if calibration_dir_path is not None:
        original_cal_file_path_x = os.path.join(calibration_dir_path, "calibration_x.out")
        original_cal_file_path_y = os.path.join(calibration_dir_path, "calibration_y.out")
        cal_file_path_x = os.path.join(output_path, "calibration_x.out")
        cal_file_path_y = os.path.join(output_path, "calibration_y.out")
        utils.iotools.copy_item(original_cal_file_path_x, cal_file_path_x)
        utils.iotools.copy_item(original_cal_file_path_y, cal_file_path_y)

        calibration_twiss["X"] = Python_Classes4MAD.metaclass.twiss(cal_file_path_x)
        calibration_twiss["Y"] = Python_Classes4MAD.metaclass.twiss(cal_file_path_y)
        return calibration_twiss
    else:
        return None
# END _copy_calibration_files --------------------------------------------------------------------


def _get_commonbpms(ListOfFiles, model, union, plane):
    """Returns a DataFrame with the list of common BPMs. 

    Hack for now: known bad BPMs are excluded at this point.
    """

    LOGGER.debug("Ignoring Bad BPMs")
    if plane == "H":
        good_bpms_index = model.index.drop(BAD_BPMS_hor, errors='ignore')
    else:
        good_bpms_index = model.index.drop(BAD_BPMS_ver, errors='ignore')

    if union:
        common = pd.DataFrame(model.loc[:, "S"])
        common["NFILES"] = 0
        for i in range(len(ListOfFiles)):
            common.loc[ListOfFiles[i].index, "NFILES"] += 1
        common = common.loc[good_bpms_index]
        union = common.drop(common[common.loc[:,"NFILES"] < min(len(ListOfFiles), 2)].index)
        inters = common.drop(common[common.loc[:,"NFILES"] < len(ListOfFiles)].index)

        return inters, union

    common_index = good_bpms_index
    for i in range(len(ListOfFiles)):
        common_index = common_index.intersection(ListOfFiles[i].index)

    commbpms = pd.DataFrame(model.loc[common_index, "S"])
    commbpms["NFILES"] = len(ListOfFiles)
    return commbpms, commbpms


def _tb_():
    if sys.stdout.isatty():
        err_exc = re.sub(r"line\s([0-9]+)", "\33[1mline \33[38;2;80;160;255m\\1\33[0m\33[21m", traceback.format_exc())
        err_exc = re.sub("File\\s\"([^\"]+)\",", "File \33[38;2;0;255;100m\\1\33[0m", err_exc)
        err_excs = err_exc.split("\n")
        for line in err_excs:
            LOGGER.error(line)
    else:
        LOGGER.error(traceback.format_exc())


#---------------------------------------------------------------------------------------------------
# helper classes for data structures
#---------------------------------------------------------------------------------------------------

class _GetllmData(object):
    ''' Holds some data from parameters of main function. '''

    def __init__(self):
        '''Constructor'''
        self.outputpath = ""
        self.list_of_input_files = []

        self.accelerator = None
        self.lhc_phase = ""
        self.bpm_unit = ""
        self.cut_for_closed_orbit = 0
        self.num_bpms_for_coupling = 0 
        self.use_only_three_bpms_for_beta_from_phase = True
        self.number_of_bpms = 0
        self.range_of_bpms = 0
        self.errordefspath = ""
        self.parallel = False
        self.nprocesses = 1
        self.important_pairs = {}
        self.onlycoupling = 0
        self.union = False

    def set_outputpath(self, outputpath):
        ''' Sets the outputpath and creates directories if they not exist.

        :param string outputpath: Path to output dir. If dir(s) to output do(es) not exist, it/they will be created.
        '''
        utils.iotools.create_dirs(outputpath)
        self.outputpath = outputpath

    def set_bpmu_and_cut_for_closed_orbit(self, cut_co, bpm_unit):
        ''' Calculates and sets the cut and bpm unit.
        :param int cut_co: Cut in um(micrometer).
        :param string bpm_unit: Indicates used unit. um, mm, cm or m
        '''
        self.bpm_unit = bpm_unit

        if bpm_unit == 'um':
            self.cut_for_closed_orbit = cut_co
        elif bpm_unit == 'mm':
            self.cut_for_closed_orbit = cut_co / 1.0e3
        elif bpm_unit == 'cm':
            self.cut_for_closed_orbit = cut_co / 1.0e4
        elif bpm_unit == 'm':
            self.cut_for_closed_orbit = cut_co / 1.0e6
        else:
            print >> sys.stderr, "Wrong BPM unit:", bpm_unit


class _TwissData(object):
    ''' Holds twiss instances of all src files. '''
    def __init__(self):
        '''Constructor'''
        self.zero_dpp_x = []  # List of src files which have dpp==0.0
        self.non_zero_dpp_x = []  # List of src files which have dpp!=0.0
        self.zero_dpp_y = []  # List of src files which have dpp==0.0
        self.non_zero_dpp_y = []  # List of src files which have dpp!=0.0    
        self.zero_dpp_commonbpms_x = []
        self.zero_dpp_commonbpms_y = []
        self.zero_dpp_unionbpms_x = None
        self.zero_dpp_unionbpms_y = None

    def has_zero_dpp_x(self):
        ''' Returns True if _linx file(s) exist(s) with dpp==0 '''
        return 0 != len(self.zero_dpp_x)

    def has_non_zero_dpp_x(self):
        ''' Returns True if _linx file(s) exist(s) with dpp!=0 '''
        return 0 != len(self.non_zero_dpp_x)

    def has_zero_dpp_y(self):
        ''' Returns True if _liny file(s) exist(s) with dpp==0 '''
        return 0 != len(self.zero_dpp_y)

    def has_non_zero_dpp_y(self):
        ''' Returns True if _liny file(s) exist(s) with dpp!=0 '''
        return 0 != len(self.non_zero_dpp_y)
    
    def has_no_input_files(self):
        return not self.has_zero_dpp_x() and not self.has_zero_dpp_y() and not self.has_non_zero_dpp_x() and not self.has_non_zero_dpp_y()

class _TuneData(object):
    ''' Used as data structure to hold tunes and phase advances. '''
    def __init__(self):
        '''Constructor'''
        self.q1 = 0.0  # Driven horizontal tune
        self.q2 = 0.0  # Driven vertical tune
        self.mux = 0.0  # Driven horizontal phase advance
        self.muy = 0.0  # Driven vertical phase advance

        # Free is from analytic equation
        self.q1f = 0.0  # Free horizontal tune
        self.q2f = 0.0  # Free vertical tune
        self.muxf = 0.0  # Free horizontal phase advance
        self.muyf = 0.0  # Free vertical phase advance
        self.q1mdl = 0.0
        self.q2mdl = 0.0

        self.q1mdlf = 0.0
        self.q2mdlf = 0.0

        # Free2 is using the effective model
        self.muxf2 = 0.0  # Free2 horizontal phase advance
        self.muyf2 = 0.0  # Free2 vertical phase advance

        self.delta1 = None  # Used later to calculate free Q1. Only if with ac calculation.
        self.delta2 = None  # Used later to calculate free Q2. Only if with ac calculation.

#---------------------------------------------------------------------------------------------------
# main invocation
#---------------------------------------------------------------------------------------------------

def _start():
    '''
    Starter function to avoid polluting global namespace with variables options,args.
    Before the following code was after 'if __name__=="__main__":'
    '''

    options, accelerator = _parse_args()

    if options.errordefspath is not None:
        accelerator.set_errordefspath(options.errordefspath)

    main(accelerator,
         accelerator.model_dir,
         outputpath=options.output,
         files_to_analyse=options.files,
         lhcphase=options.lhcphase,
         bpmu=options.BPMUNIT,
         cocut=float(options.COcut),
         nbcpl=int(options.NBcpl),
         nonlinear=options.nonlinear,
         bbthreshold=options.bbthreshold,
         errthreshold=options.errthreshold,
         use_only_three_bpms_for_beta_from_phase=options.use_only_three_bpms_for_beta_from_phase,
         number_of_bpms=options.number_of_bpms,
         range_of_bpms=options.range_of_bpms,
         use_average=options.use_average,
         calibration_dir_path=options.calibration_dir_path,
         nprocesses=options.nprocesses,
         union=options.union,
         onlycoupling=options.onlycoupling)

if __name__ == "__main__":
    _start()
