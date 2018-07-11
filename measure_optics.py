import os
from os.path import isfile
import sys
import traceback
import datetime
from time import time
import re
from collections import OrderedDict
import pandas as pd
import numpy as np
from GetLLM import algorithms
import Python_Classes4MAD.metaclass
from GetLLM.tfs_utils_getllm import GetllmTfsFile, GetllmPandasTfs
from GetLLM.GetLLMError import GetLLMError, CriticalGetLLMError
import utils.iotools
from GetLLM import optics_input
from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas
from utils import logging_tools
from utils.entrypoint import ArgumentError


VERSION = 'V3.0.0 Dev'
DEBUG = sys.flags.debug  # True with python option -d! ("python -d measure_optics.py...") (vimaier)
LOGGER = logging_tools.get_logger(__name__)
BREAK = False
PLANES = ('X', 'Y')

def print_time(index, t):
    LOGGER.debug(":::  GetLLM time  >>>>>>>>>> {:8.3f} s".format(t))


def main(optics_input):
    LOGGER.info("Starting GetLLM " + VERSION)
    global __getllm_starttime
    __getllm_starttime = time()
    twiss_d = _TwissData()
    tune_d = _TuneData()
    tune_d.q1mdl = optics_input.accelerator.get_model_tfs().headers["Q1"]
    tune_d.q2mdl = optics_input.accelerator.get_model_tfs().headers["Q2"]
    utils.iotools.create_dirs(optics_input.outputdir)
    logging_tools.add_module_handler(logging_tools.file_handler(os.path.join(optics_input.outputdir, "getllm.log")))
    header_dict = _get_header()
    if sys.flags.debug:
        LOGGER.info("     DEBUG ON")

    # Creates the output files dictionary
    files_dict = _create_tfs_files(optics_input, os.path.join(model_dir, "twiss.dat"))
        # Copy calibration files calibration_x/y.out from calibration_dir_path to outputpath
    calibration_twiss = _copy_calibration_files(optics_input.outputdir, calibration_dir_path)
    print_time("BEFORE_ANALYSE_SRC", time() - __getllm_starttime)
    input_files = InputFiles(optics_input.files)
    twisses_x, twisses_y = _analyse_src_files1(optics_input.files)

    twiss_d = _analyse_src_files11(optics_input, twiss_d, optics_input.no_averaged_tune,
                                    calibration_twiss, optics_input.accelerator.get_model_tfs(), twisses_x, twisses_y)

    files_dict = _analyse_src_files2(twiss_d, optics_input.nonlinear, files_dict, twisses_x, twisses_y)

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
            beta_d = algorithms.beta_from_amplitude.calculate_beta_from_amplitude(
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
    if InputFiles.nonlinear:
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

# END main() ---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# helper-functions
#---------------------------------------------------------------------------------------------------
def _get_header():
    return OrderedDict([('GetLLMVersion', VERSION),
                        ('Command', sys.executable + " '" + "' '".join([] + sys.argv) + "'"),
                        ('CWD', os.getcwd()),
                        ('Date', datetime.datetime.today().strftime("%d. %B %Y, %H:%M:%S"))])
    # TODO add model directory


def _create_tfs_files(getllm_d, model_filename):
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
    files_dict = {}
    _fill_files_dict_plane(files_dict, "x", getllm_d)
    _fill_files_dict_plane(files_dict, "y", getllm_d)
    files_dict['getcouple.out'] = GetllmTfsFile('getcouple.out')
    if getllm_d.nonlinear:
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
    return files_dict
# END _create_tfs_files -----------------------------------------------------------------------------

def _fill_files_dict_plane(files_dict, plane, getllm_d):

    if getllm_d.accelerator.excitation == AccExcitationMode.FREE:
        files_dict['getphase{}_free.out'.format(plane)] = GetllmTfsFile('getphase{}.out'.format(plane))
        files_dict['getphasetot{}_free.out'.format(plane)] = GetllmTfsFile('getphasetot{}.out'.format(plane))
        files_dict['getbeta{}_free.out'.format(plane)] = GetllmPandasTfs('getbeta{}.out'.format(plane))
        files_dict['getIP{}.out'.format(plane)] = GetllmTfsFile('getIP{}.out'.format(plane))
    else:
        files_dict['getphase{}.out'.format(plane)] = GetllmTfsFile('getphase{}.out'.format(plane))
        files_dict['getphasetot{}.out'.format(plane)] = GetllmTfsFile('getphasetot{}.out'.format(plane))
        files_dict['getphase{}_free.out'.format(plane)] = GetllmTfsFile('getphase{}_free.out'.format(plane))
        files_dict['getphasetot{}_free.out'.format(plane)] = GetllmTfsFile('getphasetot{}_free.out'.format(plane))
        files_dict['getbeta{}.out'.format(plane)] = GetllmPandasTfs('getbeta{}.out'.format(plane))
        files_dict['getbeta{}_free.out'.format(plane)] = GetllmPandasTfs('getbeta{}_free.out'.format(plane))
        files_dict['getIP{}.out'.format(plane)] = GetllmTfsFile('getIP{}.out'.format(plane))
        files_dict['getIP{}_free.out'.format(plane)] = GetllmTfsFile('getIP{}_free.out'.format(plane))


def _analyse_src_files1(files_to_analyse):
    LOGGER.debug("Start analysing source files")
    file_pandas = {'X': [], 'Y': []}
    for file_in in files_to_analyse.split(','):
        LOGGER.debug("> file: '{:s}'".format(file_in))
        for plane in ["X", "Y"]:
            if isfile(file_in + '.lin' + plane.lower()):
                file_to_load = file_in + '.lin' + plane.lower()
            else:
                file_to_load = file_in + '_lin' + plane.lower()
            try:
                file_pandas['X'].append(tfs_pandas.read_tfs(file_to_load).set_index("NAME"))
            except IOError:
                LOGGER.warning("Cannot load file: " + file_to_load)
            except ValueError:
                pass
    tfs_files_x = algorithms.dpp.arrange_dpp(file_pandas['X'])
    tfs_files_y = algorithms.dpp.arrange_dpp(file_pandas['Y'])

    return tfs_files_x, tfs_files_y


def _analyse_src_files11(getllm_d, twiss_d, no_average, calibration_twiss, model, tfs_files_x,
                         tfs_files_y):

    union = getllm_d.union
    for twiss_file_x, twiss_file_y in zip(tfs_files_x, tfs_files_y):
        if not no_average:
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
                    'getphasex_free.out',
                    'getphasex_free2.out',
                    'getphasetotx_free.out',
                    'getphasetotx_free2.out',
                    'getbetax_free.out',
                    'getbetax_free2.out'
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
                    'getphasey_free.out',
                    'getphasey_free2.out',
                    'getphasetoty_free.out',
                    'getphasetoty_free2.out',
                    'getbetay_free.out',
                    'getbetay_free2.out'
                ])

    if not twiss_d.has_zero_dpp_x():
        print 'Warning: you are running GetLLM without "linx of dp/p=0". Are you sure?'

        if twiss_d.has_non_zero_dpp_x():
            _add_filename_to_header_for_files(
                files_dict, "chrommode",
                [
                    'getphasex.out',
                    'getbetax.out',
                    'getcouple.out',
                    'getcouple_free.out',
                    'getcouple_free2.out',
                    'getphasey.out',
                    'getbetay_free.out'
                ])

    if twiss_d.has_no_input_files():
        print >> sys.stderr, "No parsed input files"
        sys.exit(1)

    return files_dict

def _add_filename_to_header_for_files(files_dict, filex, files_list):
    for key in files_list:
        if key in files_dict:
            files_dict[key].add_filename_to_getllm_header(filex)


# END _analyse_src_files ----------------------------------------------------------------------------

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


class InputFiles(dict):
    """
    Stores the input files, provides methods to gather quantity specific data
    Public methods:
        get_dpps(plane)
        get_joined_frame(plane, columns, zero_dpp=False, how='inner')
        get_columns(frame, column)
        get_data(frame, column)
    """
    def __init__(self, files_to_analyse):
        super(InputFiles, self).__init__(zip(PLANES, ([], [])))
        for file_in in files_to_analyse.split(','):
            for plane in PLANES:
                if isfile(file_in + '.lin' + plane.lower()):
                    file_to_load = file_in + '.lin' + plane.lower()
                else:
                    file_to_load = file_in + '_lin' + plane.lower()
                try:
                    self[plane].append(tfs_pandas.read_tfs(file_to_load).set_index("NAME"))
                except IOError:
                    LOGGER.warning("Cannot load file: " + file_to_load)
                except ValueError:
                    pass
            self[plane] = algorithms.dpp.arrange_dpp(self[plane])
        if len(self['X']) + len(self['Y']) == 0:
            raise IOError("No valid input files")

    def get_dpps(self, plane):
        """
        Gathers measured DPPs from input files corresponding to given plane
        Parameters:
            plane: "X" or "Y"

        Returns:
            numpy array of DPPs
        """
        return np.array([df.DPP for df in self[plane]])

    def _get_zero_dpp_frames(self, plane):
        return self[plane][np.logical_not(self.get_dpps(self, plane))]

    def _get_all_frames(self, plane):
        return self[plane]

    def get_joined_frame(self, plane, columns, zero_dpp=False, how='inner'):
        """
        Constructs merged DataFrame from InputFiles
        Parameters:
            plane:  "X" or "Y"
            columns: list of columns from input files
            zero_dpp: if True merges only zero-dpp files, default is False
            how: way of merging:  'inner' (intersection) or 'outer' (union), default is 'inner'
        Returns:
            merged DataFrame from InputFiles
        """
        if how not in ['inner', 'outer']:
            raise RuntimeWarning("'how' should be either 'inner' or 'outer', 'inner' will be used.")
        if zero_dpp:
            frames_to_join = self._get_zero_dpp_frames(plane)
        else:
            frames_to_join = self._get_all_frames(plane)
        if len(frames_to_join) == 0:
            raise ValueError("No data found")
        joined_frame = pd.DataFrame(self[plane][0]).loc[:, columns]
        for i, df in enumerate(self[plane][1:]):
            joined_frame = pd.merge(joined_frame, df.loc[:, columns], how=how, left_index=True,
                                    right_index=True, suffixes=('', '__' + str(i + 1)))
        for column in columns:
            joined_frame.rename(columns={column: column + '__0'}, inplace=True)
        return joined_frame

    @staticmethod
    def get_columns(frame, column):
        """
        Returns list of columns of frame corresponding to column in original files
        Parameters:
            frame:  joined frame
            column: name of column in original files
        Returns:
            list of columns
        """
        return frame.columns.str.startswith(column + '__').sort(key=lambda x: int(x.strip(column + '__')))

    @staticmethod
    def get_data(frame, column):
        """
        Returns data in columns of frame corresponding to column in original files
        Parameters:
            frame:  joined frame
            column: name of column in original files
        Returns:
            data in numpy array corresponding to column in original files
        """
        return frame.loc[:, frame.getcolumns(frame, column)].values


if __name__ == "__main__":
    main(*optics_input.parse_args())
