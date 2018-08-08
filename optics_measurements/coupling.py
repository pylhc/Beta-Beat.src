'''
Created on 27 May 2013

@author: ?, vimaier, dwierichs

@version: 0.0.1

GetLLM.algorithms.coupling.py stores helper functions for coupling calculations for GetLLM.
This module is not intended to be executed. It stores only functions.

Change history:
 - <version>, <author>, <date>:
    <description>

    - STRUCTURE - 
main part
    calculate_coupling
helper-functions
    GetCoupling1
    _get_coupling_2bpm
    getCandGammaQmin
    _find_sign_QxmQy
ac-dipol stuff
    getFreeCoupling
'''

import sys, os
import traceback
import math
import copy
import numpy as np
import pandas as pd

import utils.bpm
import phase
import helper
import compensate_excitation
from model.accelerators.accelerator import AccExcitationMode
from utils import logging_tools, stats, tfs_pandas
from twiss_optics.optics_class import TwissOptics

LOGGER = logging_tools.get_logger(__name__)
DEBUG = sys.flags.debug # True with python option -d! ("python -d GetLLM.py...") (vimaier)

#===================================================================================================
# main part
#===================================================================================================
def calculate_coupling(getllm_d, input_files, phase_d, tunes, header_dict):
    '''
    Calculates coupling and fills the following TfsFiles:
        getcouple.out        getcouple_free.out        getcouple_free2.out        getcoupleterms.out

    :Parameters:
        'getllm_d': _GetllmData (In-param, values will only be read)
            lhc_phase, accel, beam_direction and coupling_method are used.
        'input_files': _TwissData (In-param, values will only be read)
            Holds twiss instances of the src files.
        'tune_d': _TuneData (In/Out-param, values will be read and set)
            Holds tunes and phase advances. q1, mux, q2 and muy will be set if
            "coupling_method == 2" and accel is 'SPS' or 'RHIC'.

    :Return: _TuneData
        the same instance as param tune_d to indicate that tunes will be set.
    '''
    zero_dpp_frames_x = input_files.zero_dpp_frames('X')
    zero_dpp_frames_y = input_files.zero_dpp_frames('Y')
    LOGGER.info("Calculating coupling using the {0}-BPM-method and {1} file(s)"
                .format(getllm_d.coupling_method, len(zero_dpp_frames_x)))
    LOGGER.info("New coupling -- pandafied")

    # get the models -------------------------------------------------------------------------------
    accelerator = getllm_d.accelerator
    mad_twiss = accelerator.get_model_tfs()
    mad_elements = accelerator.get_elements_tfs()
    coupling_df = None
    free_coupling_df = None
    if accelerator.excitation is not AccExcitationMode.FREE:
        mad_ac = accelerator.get_driven_tfs()
    else:
        mad_ac = None

    if len(zero_dpp_frames_x) == len(zero_dpp_frames_y) and len(zero_dpp_frames_x) > 0:
        #-- Coupling in the model
        LOGGER.debug("Calculating coupling in the model.")
        optics_twiss = TwissOptics(mad_elements)
        optics_twiss.calc_cmatrix()
        optics_coupling = optics_twiss.get_coupling(method="cmatrix")  # this is not used!
        #-- Main part ------------------------------------------------------------------------------
        # 1-BPM method
        if getllm_d.coupling_method == 1:

            coupling_df = _get_coupling_1bpm(
                mad_twiss, zero_dpp_frames_x, zero_dpp_frames_y, tunes['X']["QF"], tunes['Y']["QF"],
                getllm_d.accelerator.get_beam_direction())
        # 2-BPM method
        elif getllm_d.coupling_method == 2:
            fwqw, second_ret = _get_coupling_2bpm(bpm
                mad_twiss, zero_dpp_frames_x, zero_dpp_frames_y, tunes['X']["QF"], tunes['Y']["QF"],
                phase_d['X']['F'], phase_d['Y']['F'],
                accelerator.get_beam_direction(), getllm_d.accelerator)
        else:
            raise ValueError('Number of monitors for coupling analysis should be 1 or 2 (option -n)')

        coupling_df["MDLF1001R"] = optics_coupling.loc[coupling_df.index, "F1001"].real
        coupling_df["MDLF1001I"] = optics_coupling.loc[coupling_df.index, "F1001"].imag
        coupling_df["MDLF1010R"] = optics_coupling.loc[coupling_df.index, "F1010"].real
        coupling_df["MDLF1010I"] = optics_coupling.loc[coupling_df.index, "F1010"].imag

        coupling_df["DELTAF1001"] = coupling_df["F1001W"] - \
            abs(optics_coupling.loc[coupling_df.index, "F1001"])

        tfs_pandas.write_tfs(os.path.join(getllm_d.outputdir, "getcoupling.out"),
                             coupling_df)
        #-- ac to free coupling
        if accelerator.excitation != AccExcitationMode.FREE:
            if getllm_d.coupling_method == 2:
                #-- analytic eqs
                [fwqwf, bpmsf] = compensate_excitation.GetFreeCoupling_Eq(
                    mad_twiss, input_files.zero_dpp_x,
                    input_files.zero_dpp_y, bpms, tune_d.q1, tune_d.q2,
                    tune_d.q1f, tune_d.q2f, accelerator)

                freecoupling_df.loc["COUNT"] = len(input_files.zero_dpp_x)
                freecoupling_df.loc["F1001W"] = abs(fwqwf[:][0][0])
                freecoupling_df.loc["F1001R"] = fwqwf[:][0][0].real
                freecoupling_df.loc["F1001I"] = fwqwf[:][0][0].imag
                freecoupling_df.loc["FWSTD1"] = fwqwf[:][0][1]
                freecoupling_df.loc["F1010W"] = abs(fwqwf[:][0][2])
                freecoupling_df.loc["F1010R"] = fwqwf[:][0][2].real
                freecoupling_df.loc["F1010I"] = fwqwf[:][0][2].imag
                freecoupling_df.loc["FWSTD2"] = fwqwf[:][0][3]
                freecoupling_df.loc["Q1001"] = fwqwf[:][1][0]
                freecoupling_df.loc["Q1001STD"] = fwqwf[:][1][1]
                freecoupling_df.loc["Q1010"] = fwqwf[:][1][2]
                freecoupling_df.loc["Q1010STD"] = fwqwf[:][1][3]
                freecoupling_df.loc["MDLF1001R"] = optics_coupling["F1001"].real
                freecoupling_df.loc["MDLF1001I"] = optics_coupling["F1001"].imag
                freecoupling_df.loc["MDLF1010R"] = optics_coupling["F1010"].real
                freecoupling_df.loc["MDLF1010I"] = optics_coupling["F1010"].imag

                freecoupling_df.headers = header_dict.copy()
                freecoupling_df.headers["CG"] = fwqwf['Global'][0]
                freecoupling_df.headers["QC"] = fwqwf['Global'][1]

            #-- global factor
            [fwqwf2, bpmsf2] = getFreeCoupling(tune_d.q1f, tune_d.q2f, tune_d.q1, tune_d.q2, fwqw, mad_twiss, bpms)

            coupleterms_df.loc["COUNT"] = len(input_files.zero_dpp_x)
            coupleterms_df.loc["F1001W"] = abs(fwqwf2[:][0][0])
            coupleterms_df.loc["F1001R"] = fwqwf2[:][0][0].real
            coupleterms_df.loc["F1001I"] = fwqwf2[:][0][0].imag
            coupleterms_df.loc["FWSTD1"] = fwqwf2[:][0][1]
            coupleterms_df.loc["F1010W"] = abs(fwqwf2[:][0][2])
            coupleterms_df.loc["F1010R"] = fwqwf2[:][0][2].real
            coupleterms_df.loc["F1010I"] = fwqwf2[:][0][2].imag
            coupleterms_df.loc["FWSTD2"] = fwqwf2[:][0][3]
            coupleterms_df.loc["Q1001"] = fwqwf2[:][1][0]
            coupleterms_df.loc["Q1001STD"] = fwqwf2[:][1][1]
            coupleterms_df.loc["Q1010"] = fwqwf2[:][1][2]
            coupleterms_df.loc["Q1010STD"] = fwqwf2[:][1][3]
            coupleterms_df.loc["MDLF1001R"] = optics_coupling["F1001"].real
            coupleterms_df.loc["MDLF1001I"] = optics_coupling["F1001"].imag
            coupleterms_df.loc["MDLF1010R"] = optics_coupling["F1010"].real
            coupleterms_df.loc["MDLF1010I"] = optics_coupling["F1010"].imag

            coupleterms_df.headers = header_dict.copy()
            coupleterms_df.headers["CG"] = fwqwf2['Global'][0]
            coupleterms_df.headers["QC"] = fwqwf2['Global'][1]


        #-- Convert to C-matrix:
        if accelerator.excitation != AccExcitationMode.FREE and (fwqwf is not None):
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

        return coupling_df, free_coupling_df, coupleterms_df
    # if numbers of files don't match or are 0:
    raise IOError("Input files unsuitable for coupling calculation. "
                  "Hor: {} files, Ver: {} files".format(len(zero_dpp_frames_x),
                                                        len(zero_dpp_frames_y)))


# --------------------------------------------------------------------------------------------------
# helper-functions
# --------------------------------------------------------------------------------------------------

def _get_coupling_1bpm(MADTwiss, list_zero_dpp_x, list_zero_dpp_y, tune_x, tune_y,
                       beam_direction):
    """Calculate coupling and phase with 1-BPM method for all BPMs and overall
    INPUT
     MADTwiss        - twiss instance of model from MAD
     list_zero_dpp_x - list with twiss objects for horizontal data
     list_zero_dpp_y - list with twiss objects for vertical data
     tune_x          - horizontal tune (use natural/free tunes!)
     tune_y          - vertical tune (use natural/free tunes!)
     outputpath      - directory which contains the Drive.inp file to determine operation point
    OUTPUT
     fwqw            - library with BPMs and corresponding results
     dbpms           - list of BPMs with correct phase
    Global: fwqw = [CG,QG,CG_std]
    """

    # Determine intersection of BPM-lists between measurement and model, refactor this maybe with
    # --- collect the needed columns of the Input files --------------------------------------------

    amp01 = _get_columns(list_zero_dpp_x, "AMP01")
    amp10 = _get_columns(list_zero_dpp_y, "AMP10")
    common_index = amp01.index.intersection(amp10.index)
    amp01 = amp01.loc[common_index]
    amp10 = amp10.loc[common_index]

    erramp01 = _get_columns(list_zero_dpp_x, "ERR_AMP01").loc[common_index]
    erramp10 = _get_columns(list_zero_dpp_y, "ERR_AMP10").loc[common_index]

    C01 = stats.weighted_mean(amp01.values,
                              errors=erramp01.values, axis=1)
    C10 = stats.weighted_mean(amp10.values,
                              errors=erramp10.values, axis=1)

    stdC01 = stats.weighted_error(amp01.values,
                                  errors=erramp01.values, axis=1)
    stdC10 = stats.weighted_error(amp10.values,
                                  errors=erramp10.values, axis=1)

    phase01 = stats.weighted_mean(
        _get_columns(list_zero_dpp_x, "PHASE01").loc[common_index],
        errors=_get_columns(list_zero_dpp_x, "PHASE01").loc[common_index], axis=1)
    phase10 = stats.weighted_mean(
        _get_columns(list_zero_dpp_y, "PHASE10").loc[common_index],
        errors=_get_columns(list_zero_dpp_y, "PHASE10").loc[common_index], axis=1)

    mux = np.mean(_get_columns(list_zero_dpp_x, "MUX").loc[common_index], axis=1)
    muy = np.mean(_get_columns(list_zero_dpp_y, "MUY").loc[common_index], axis=1)
    # phase of f1001 from spectral lines
    q1 = phase10 - mux
    q2 = muy - phase01

    # calculation of f1001 and error
    exp_q = (np.exp(2j * np.pi * q1) + np.exp(2j * np.pi * q2)) * 0.5
    f = 0.5 * np.arctan(np.sqrt(C01 * C10)) * exp_q
    f_err = np.sqrt((stdC01 * C10) ** 2 + (stdC10 * C01) ** 2)

    # select good BPMs and create output dataframe
    good_bpm = np.where(np.logical_and(abs(q1 - q2) >= 0.25, abs(q1 - q2) <= 0.75), True, False)
    common_index = common_index[good_bpm]
    f = f[good_bpm]
    f_err = f_err[good_bpm]
    coupling_df = tfs_pandas.TfsDataFrame(MADTwiss.loc[common_index, ["S"]])

    coupling_df = coupling_df.assign(F1001W=abs(f), F1001R=f.real, F1001I=f.imag, FWSTD1=f_err,
                                     F1010W=0.0, F1010R=0.0, F1010I=0.0, FWSTD2=0.0, Q1001=0.0,
                                     Q1001STD=0.0, Q1010=0.0, Q1010STD=0.0)

    # Calculate std of C- with error on weighted average
    coupling_df.headers["CG"] = 4.0 * \
        abs(tune_x - tune_y) * stats.weighted_mean(abs(f), errors=f_err)
    # ... error of C-
    coupling_df.headers["CG_std"] = 4.0 * \
        abs(tune_x - tune_y) * stats.weighted_error(abs(f), errors=f_err)
    # Determine global phase using average and operation point
    # TODO: lookup what is calculated here and implement
    LOGGER.warning("the calculation of global phase is not yet implemented.")
#    QG = (QG/len(dbpms)+0.5*(1.0-sign_QxmQy*0.5))%1.0
    # Cast determined results as global

    coupling_df.headers["Qx"] = tune_x
    coupling_df.headers["Qy"] = tune_y

    return coupling_df

def _get_coupling_2bpm(bpmMADTwiss, list_zero_dpp_x, list_zero_dpp_y, tune_x, tune_y, phasex, phasey, beam_direction, accel):
    """Calculate coupling and phase with 2-BPM method for all BPMs and overall
    INPUT
     MADTwiss        - twiss instance of model from MAD 
     list_zero_dpp_x - list with twiss objects for horizontal data
     list_zero_dpp_y - list with twiss objects for vertical data
     tune_x          - horizontal tune (use natural/free tunes!)
     tune_y          - vertical tune (use natural/free tunes!)
     phasex          - horizontal phase (usage with phase.py in algorithms)
     phasey          - vertical phase (usage with phase.py in algorithms)
     beam_direction  - direction of beam, LHCB1: 1, LHCB2: -1
     accel           - accelerator (LHCB1, LHCB2, LHCB4(?), SPS, RHIC,TEVATRON)
     outputpath      - directory which contains the Drive.inp file to determine operation point
    OUTPUT
     fwqw            - library with BPMs and corresponding results
     dbpms           - list of BPMs with correct phase
    """

    ### Prepare BPM lists ###

    # Check linx/liny files, if it's OK it is confirmed that ListofZeroDPPX[i] and ListofZeroDPPY[i]
    # come from the same (simultaneous) measurement. It might be redundant check.
    if len(list_zero_dpp_x) != len(list_zero_dpp_y):
        print >> sys.stderr, 'Leaving GetCoupling as linx and liny files seem not correctly paired...'
        dum0 = {"Global": [0.0, 0.0]}
        dum1 = []
        return [dum0, dum1]
    # Determine intersection of BPM-lists between measurement and model, create list dbpms

    modelphases_x = phasex["MODEL"]
    modelphases_y = phasey["MODEL"]
    LOGGER.debug(modelphases_x.index)

    LOGGER.debug("intersecting bpm lists for coupling calculation")
    XplusY = list_zero_dpp_x + list_zero_dpp_y
    dbpms = utils.bpm.intersect(XplusY)
    dbpms = utils.bpm.model_intersect(dbpms, MADTwiss)
    dbpms = dbpms.loc[dbpms.index.intersection(
        modelphases_x.index.intersection(
            modelphases_y.index))]
    LOGGER.debug("BPM list after intersections: {}, length = {}".format(type(dbpms), len(dbpms)))

    ### Calculate fw and qw, exclude BPMs having wrong phases ###

    # Initialize dictionary of BPMs with results
    fwqw = {}
    # Initialize dictionary of BPMs with results for old method
    f_old_out = {}
    # Initialize list of BPMs with correct phases
    dbpmt = []
    # Count number of BPM-pairs in intersection of model and measurement
    Numbpmpairs = len(dbpms) - 1
    bn0 = dbpms.index[0]
    # Loop through BPM-pairs
    for i in range(Numbpmpairs):
        # Get BPM names
        bn1 = dbpms.index[i]
        bn2 = dbpms.index[i + 1]
        delx = modelphases_x.loc[bn1, bn0] - 0.25  # Missprint in the coupling note
        dely = modelphases_y.loc[bn1, bn0] - 0.25

        # Initialize result-lists (coupling, error on coupling and phases)
        f1001ij = []
        f1010ij = []
        std_f1001ij = []
        std_f1010ij = []
        q1js = []
        q2js = []
        q1jd = []
        q2jd = []
        # Initialize as not bad BPM
        badbpm = 0
        # Loop through files to analyze
        for j in range(0, len(list_zero_dpp_x)):
            # Get twiss instance for current BPM and file
            tw_x = list_zero_dpp_x[j]
            tw_y = list_zero_dpp_y[j]
            # Get main amplitude
            [ampx_1, ampy_1] = [tw_x.loc[bn1, "AMPX"], tw_y.loc[bn1, "AMPY"]]
            [ampx_2, ampy_2] = [tw_x.loc[bn2, "AMPX"], tw_y.loc[bn2, "AMPY"]]
            # Exclude BPM if no main line was found
            if ampx_1 == 0 or ampy_1 == 0 or ampx_2 == 0 or ampy_2 == 0:
                badbpm = 1
                ampx_1 = 1  # Dummy value, badbpm variable makes sure the BPM is ignored
                ampy_1 = 1  # Dummy value, badbpm variable makes sure the BPM is ignored
                ampx_2 = 1  # Dummy value, badbpm variable makes sure the BPM is ignored
                ampy_2 = 1  # Dummy value, badbpm variable makes sure the BPM is ignored
            # Get coupled amplitude ratios
            [amp01_1, amp10_1] = [tw_x.loc[bn1, "AMP01"], tw_y.loc[bn1, "AMP10"]]
            [amp01_2, amp10_2] = [tw_x.loc[bn2, "AMP01"], tw_y.loc[bn2, "AMP10"]]
            # Replace secondary lines with amplitude infinity or 0 by noise average
            try:
                if amp01_1 == float("inf") or amp01_1 == 0:
                    amp01_1 = tw_x.loc[bn1, "AVG_NOISE"] / ampx_1
                if amp10_1 == float("inf") or amp10_1 == 0:
                    amp10_1 = tw_y.loc[bn1, "AVG_NOISE"] / ampy_1
                if amp01_2 == float("inf") or amp01_2 == 0:
                    amp01_2 = tw_x.loc[bn1, "AVG_NOISE"] / ampx_2
                if amp10_2 == float("inf") or amp10_2 == 0:
                    amp10_2 = tw_y.loc[bn1, "AVG_NOISE"] / ampy_2
            except AttributeError:
                print "AVG_NOISE column not found, cannot use noise floor."

            # Call routine in helper.py to get secondary lines for 2-BPM method
            [SA0p1ij,phi0p1ij] = helper.ComplexSecondaryLine(delx, amp01_1, amp01_2,
                    tw_x.loc[bn1, "PHASE01"], tw_x.loc[bn2, "PHASE01"])
            [SA0m1ij,phi0m1ij] = helper.ComplexSecondaryLine(delx, amp01_1, amp01_2,
                    -tw_x.loc[bn1, "PHASE01"], -tw_x.loc[bn2, "PHASE01"])
            [TBp10ij,phip10ij] = helper.ComplexSecondaryLine(dely, amp10_1, amp10_2,
                    tw_y.loc[bn1, "PHASE10"], tw_y.loc[bn2, "PHASE10"])
            [TBm10ij,phim10ij] = helper.ComplexSecondaryLine(dely, amp10_1, amp10_2,
                    -tw_y.loc[bn1, "PHASE10"], -tw_y.loc[bn2, "PHASE10"])

            # Get noise standard deviation and propagate to coupled amplitude ratio
            std_amp01_1 = tw_x.loc[bn1, "NOISE"]/ampx_1*math.sqrt(1+amp01_1**2)
            std_amp10_1 = tw_y.loc[bn1, "NOISE"]/ampy_1*math.sqrt(1+amp10_1**2)
            std_amp01_2 = tw_x.loc[bn2, "NOISE"]/ampx_2*math.sqrt(1+amp01_2**2)
            std_amp10_2 = tw_y.loc[bn2, "NOISE"]/ampy_2*math.sqrt(1+amp10_2**2)
            # Propagate to 2-BPM coupled amplitude ratio using a separate routine in helper.py
            std_SA0p1ij = helper.ComplexSecondaryLineSTD(delx, amp01_1, amp01_2,
                    tw_x.loc[bn1, "PHASE01"], tw_x.loc[bn2, "PHASE01"], std_amp01_1, std_amp01_2)
            std_SA0m1ij = helper.ComplexSecondaryLineSTD(delx, amp01_1, amp01_2,
                    -tw_x.loc[bn1, "PHASE01"], -tw_x.loc[bn2, "PHASE01"], std_amp01_1, std_amp01_2)
            std_TBp10ij = helper.ComplexSecondaryLineSTD(dely, amp10_1, amp10_2,
                    tw_y.loc[bn1, "PHASE10"], tw_y.loc[bn2, "PHASE10"], std_amp10_1, std_amp10_2)
            std_TBm10ij = helper.ComplexSecondaryLineSTD(dely, amp10_1, amp10_2,
                    -tw_y.loc[bn1, "PHASE10"], -tw_y.loc[bn2, "PHASE10"], std_amp10_1, std_amp10_2)

            # Append results for the coupling parameters
            f1001ij.append(0.5*math.sqrt(TBp10ij*SA0p1ij/2.0/2.0)) # division by 2 for each ratio as the scale of the #
            f1010ij.append(0.5*math.sqrt(TBm10ij*SA0m1ij/2.0/2.0)) # main lines is 2 (also see appendix of the note)  #
            # Propagate error to f1001 and f1010 if possible (no division by 0)
            if TBp10ij==0 or SA0p1ij==0:
                std_f1001ij.append(float("nan"))
            else:
                std_f1001ij.append(0.25*math.sqrt(4.0/TBp10ij/SA0p1ij)*math.sqrt((std_TBp10ij*SA0p1ij/4)**2+(TBp10ij*std_SA0p1ij/4)**2))
            if TBm10ij==0 or SA0m1ij==0:
                std_f1010ij.append(float("nan"))
            else:
                std_f1010ij.append(0.25*math.sqrt(4.0/TBm10ij/SA0m1ij)*math.sqrt((std_TBm10ij*SA0m1ij/4)**2+(TBm10ij*std_SA0m1ij/4)**2))

            if beam_direction == 1:
                q1jd.append((phi0p1ij-tw_y.loc[bn1, "MUY"]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append((-phip10ij+tw_x.loc[bn1, "MUX"]-0.25)%1.0)
            elif beam_direction == -1:
                q1jd.append((phi0p1ij-tw_y.loc[bn1, "MUY"]+0.25)%1.0) # note that phases are in units of 2pi
                q2jd.append(-(-phip10ij+tw_x.loc[bn1, "MUX"]-0.25)%1.0)

            # This sign change in the real part is to comply with MAD output
            q1jd[j] = (0.5-q1jd[j])%1.0 
            q2jd[j] = (0.5-q2jd[j])%1.0

            if beam_direction==1:
                q1js.append((phi0m1ij+tw_y.loc[bn1, "MUY"]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append((phim10ij+tw_x.loc[bn1, "MUX"]+0.25)%1.0)
            if beam_direction==-1:
                q1js.append((phi0m1ij+tw_y.loc[bn1, "MUY"]+0.25)%1.0) # note that phases are in units of 2pi
                q2js.append(-(phim10ij+tw_x.loc[bn1, "MUX"]+0.25)%1.0)
            # This sign change in the real part is to comply with MAD output
            q1js[j] = (0.5-q1js[j])%1.0
            q2js[j] = (0.5-q2js[j])%1.0

        q1jd = np.array(q1jd)
        q2jd = np.array(q2jd)
        q1d = stats.circular_mean(q1jd, period=1.0) % 1.0
        q2d = stats.circular_mean(q2jd, period=1.0) % 1.0

        q1js = np.array(q1js)
        q2js = np.array(q2js)
        q1s = stats.circular_mean(q1js, period=1.0) % 1.0
        q2s = stats.circular_mean(q2js, period=1.0) % 1.0

        if DEBUG:
            print("\n")
            print("i=",i," j=", j," ",bn1," <--> ",bn2)
            print("      q1s=",q1s," q2s=",q2s," q1d=",q1d," q2d=",q2d)

        # Take SPS and RHIC out of the badbpm procedure (badbpm stays 0 as initialized) and set phases 
        if (accel == "SPS" or accel == "RHIC"):
            print("accel is ", accel, " disabling wrong phase check")
            
            q1010i = q1d
            q1010i = q1s

        # Check phase and set badbpm for wrong phase (only for other accels than SPS and RHIC)
        else:
            if min(abs(q1d-q2d),1.0-abs(q1d-q2d))>0.25:
                if DEBUG:
                    print("Bad BPM because min(abs(q1d-q2d),1.0-abs(q1d-q2d))>0.25")
                badbpm = 1
            if min(abs(q1s-q2s),1.0-abs(q1s-q2s))>0.25:
                if DEBUG:
                    print("Bad BPM because min(abs(q1s-q2s),1.0-abs(q1s-q2s))>0.25")
                    print(" i j ",i, j," q1s, q2s",q1s,q2s) 
                badbpm = 1
           

        # If accel is SPS or RHIC or no no wrong phase was detected, process results
        if badbpm==0:
            # Cast coupling parameters and errors as np.arrays for later calculations
            f1001ij = np.array(f1001ij)
            f1010ij = np.array(f1010ij)
            std_f1001ij = np.array(std_f1001ij)
            std_f1010ij = np.array(std_f1010ij)

            # Old cminus method: averaging abs values
            if beam_direction==1:
                f_old_out[bn1] = np.average(abs(f1001ij), weights=1/std_f1001ij**2)
            elif beam_direction==-1:
                f_old_out[bn1] = np.average(abs(f1010ij), weights=1/std_f1010ij**2)

            # Use variance-weighted average to determine f and its std
            f1001i = np.average(f1001ij, weights=1/std_f1001ij**2)
            f1010i = np.average(f1010ij, weights=1/std_f1010ij**2)
            # Set std of f1001 and f1010 by calculating the std of the weighted average
            f1001istd = np.sqrt(1/sum(1/std_f1001ij**2))
            f1010istd = np.sqrt(1/sum(1/std_f1010ij**2))

            # Use routines in phase.py to get mean and std of the phase terms q1001 and q1010
            q1001i = stats.circular_mean(np.array([q1d,q2d]), period=1.0) % 1.0
            q1010i = stats.circular_mean(np.array([q1s,q2s]), period=1.0) % 1.0
            q1001istd = stats.circular_error(np.append(q1jd,q2jd), period=1.0, t_value_corr=False)
            q1010istd = stats.circular_error(np.append(q1js,q2js), period=1.0, t_value_corr=False)
            # Calculate complex coupling terms using phases from above
            f1001i = f1001i*complex(np.cos(2.0*np.pi*q1001i),np.sin(2.0*np.pi*q1001i))
            f1010i = f1010i*complex(np.cos(2.0*np.pi*q1010i),np.sin(2.0*np.pi*q1010i))
            # Add BPM to list of BPMs with correct phase
            dbpmt.append([dbpms.iloc[i], dbpms.index[i]])

            # Save results to BPM-results dictionary, sorted depending on beam_direction
            if beam_direction==1:
                fwqw[bn1] = [[f1001i,f1001istd,f1010i,f1010istd],[q1001i,q1001istd,q1010i,q1010istd]]
            elif beam_direction==-1:
                fwqw[bn1] = [[f1010i,f1010istd,f1001i,f1001istd],[q1010i,q1010istd,q1001i,q1001istd]]

    # Count number of skipped BPMs because of wrong phase
    Badbpms = len(dbpms)-len(dbpmt)
    # Rename list of BPMs with correct phase
    dbpms = dbpmt

    # Compute global values for coupling, error and phase
    # Initialize coupling term f
    f_new = complex(0,0)
    # Initialize denominator for weighted averaging
    denom = 0
    # Loop through BPMs with correct phase
    for i in range(0, len(dbpms) - 1):
        # Get BPM-names
        bn1 = str.upper(dbpms[i][1])
        bn2 = str.upper(dbpms[i+1][1])
        mux = MADTwiss.loc[bn1, "MUX"]
        muy = MADTwiss.loc[bn2, "MUY"]
        f_new += fwqw[bn2][0][0]*np.exp(complex(0,1)*2*np.pi*(mux-muy))/fwqw[bn2][0][1]**2 # Variance-weighted average for BPMs
        denom += 1/fwqw[bn2][0][1]**2 # denominator for weighted average
    
    N = len(dbpms)
    print("coupling.py: ",N, "denom = ",denom)
    if denom == 0:
        raise Exception("All BPMs were marked as bad")
     
    f_new_std = np.sqrt(1/denom)
    CG_new_abs = 4*abs(tune_x-tune_y)*abs(f_new)/denom
    CG_new_abs_std = 4*abs(tune_x-tune_y)*abs(f_new_std)
    CG_new_phase = np.angle(f_new)

    print('NewCMINUS: {0} +/- {1}'.format(CG_new_abs, CG_new_abs_std))
    print('Skipped BPMs: {0} (badbpm) of {1}'.format(Badbpms, Numbpmpairs))

    # Old formula
    # Initialize global values coupling CG and phse QG
    CG = 0.0
    QG = 0.0
    for i in range(0,len(dbpms)-1):
        bn1 = str.upper(dbpms[i][1])
        CG += abs(f_old_out[bn1])
        # For more than one file, this goes wrong, loops are mixed up for the phase calculation!
        tw_x = list_zero_dpp_x[0]
        tw_y = list_zero_dpp_y[0]
        QG += fwqw[bn1][1][0]-(tw_x.loc[bn1, "MUX"]-tw_y.loc[bn1, "MUY"])

    if len(dbpms)==0:
        print >> sys.stderr, 'Warning: There is no BPM to output linear coupling properly... leaving Getcoupling.'
        # Does this set a coupling without prefactor 4*(Qx-Qy) to global?
        fwqw['Global']=[CG,QG] #Quick fix Evian 2012
        return [fwqw,dbpms]
    else:
        CG_old = abs(4.0*(tune_x-tune_y)*CG/len(dbpms))
    print 'OldCMINUS' ,CG_old
    fwqw['Global'] = [CG_new_abs,CG_new_phase,CG_new_abs_std]

    return [fwqw,dbpms]

### END of _get_coupling_2bpm bpm###

def getCandGammaQmin(fqwq,bpms,tunex,tuney,twiss):
    # Cut the fractional part of Q1 and Q2
    QQ1 = float( int(twiss.Q1) )
    QQ2 = float( int(twiss.Q2) )

    tunex=float(tunex)+QQ1
    tuney=float(tuney)+QQ2

    tunefactor=(np.cos(2*np.pi*tunex)-np.cos(2*np.pi*tuney))/(np.pi*(np.sin(2*np.pi*tunex)+np.sin(2*np.pi*tuney)))

    coupleterms={}
    Qmin=[]

    if len(bpms)==0:
        print >> sys.stderr, "No bpms in getCandGammaQmin. Returning empty stuff"
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

### END of getCandGammaQmin ###

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

### END of _find_sign_QxmQy ###

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
    factor_top_diff=math.sqrt(np.sin(np.pi*(tunedrivenx-tunefreey))*np.sin(np.pi*(tunefreex-tunedriveny)))
    factor_bottom_diff=np.sin(np.pi*(tunefreex-tunefreey))

    factor_diff=abs((factor_top_diff/factor_bottom_diff))

    if DEBUG:
        print "Factor for coupling diff ",factor_diff

    # sum f1010
    factor_top_sum=math.sqrt(np.sin(np.pi*(tunedrivenx+tunefreey))*np.sin(np.pi*(tunefreex+tunedriveny)))
    factor_bottom_sum=np.sin(np.pi*(tunefreex+tunefreey))

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

### END of getFreeCoupling ###

def GetFreeCoupling_Eq(MADTwiss,FilesX,FilesY,Qh,Qv,Qx,Qy,psih_ac2bpmac,psiv_ac2bpmac,bd,acdipole,oa):

    #-- Details of this algorithms is in http://www.agsrhichome.bnl.gov/AP/ap_notes/ap_note_410.pdf

    #-- Check linx/liny files, may be redundant
    if len(FilesX)!=len(FilesY): return [{},[]]

    #-- Select common BPMs
    bpm=utils.bpm.model_intersect(utils.bpm.intersect(FilesX + FilesY), MADTwiss)
    bpm=[(b[0],str.upper(b[1])) for b in bpm]

    #-- Last BPM on the same turn to fix the phase shift by Q for exp data of LHC
    #if op=="1" and bd== 1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L2.B1']]
    #if op=="1" and bd==-1: s_lastbpm=MADTwiss.S[MADTwiss.indx['BPMSW.1L8.B2']]

    #-- Determine the BPM closest to the AC dipole and its position
    #BPMYB.6L4.B1 BPMYA.5L4.B1
    # BPMWA.B5L4.B1

    horBPMsCopensation =[]
    verBPMsCopensation = []
    #bpmac1_h=psih_ac2bpmac.keys()[0]
    #bpmac2_h=psih_ac2bpmac.keys()[1]

    #bpmac1_v = psiv_ac2bpmac.keys()[0]
    #bpmac2_v = psiv_ac2bpmac.keys()[1]
    for key in psih_ac2bpmac:
        if(key in list(zip(*bpm)[1])):
            horBPMsCopensation.append(key)
            verBPMsCopensation.append(key)
    fqwList = []
    for g in range(0, len(horBPMsCopensation)):
        k_bpmac_h =list(zip(*bpm)[1]).index(horBPMsCopensation[g])
        bpmac_h=horBPMsCopensation[g]

        k_bpmac_v=list(zip(*bpm)[1]).index(verBPMsCopensation[g])
        bpmac_v=verBPMsCopensation[g]
        '''
        try:
            k_bpmac_h=list(zip(*bpm)[1]).index(bpmac1_h)
            bpmac_h=bpmac1_h
        except:
            try:
                k_bpmac_h=list(zip(*bpm)[1]).index(bpmac2_h)
                bpmac_h=bpmac2_h
            except:
                print >> sys.stderr,'WARN: BPMs next to AC dipoles or ADT missing. AC or ADT dipole effects not calculated with analytic eqs for coupling'
                return [{},[]]
        #      if 'B5R4' in b: bpmac1=b
        #if 'A5R4' in b: bpmac2=b
        try:
            k_bpmac_v=list(zip(*bpm)[1]).index(bpmac1_v)
            bpmac_v=bpmac1_v
        except:
            try:
                k_bpmac_v=list(zip(*bpm)[1]).index(bpmac2_v)
                bpmac_v=bpmac2_v
            except:
                print >> sys.stderr,'WARN: BPMs next to AC dipoles or ADT missing. AC dipole or ADT effects not calculated with analytic eqs for coupling'
                return [{},[]]
        print k_bpmac_v, bpmac_v
        print k_bpmac_h, bpmac_h
        '''
       #-- Global parameters of the driven motion
        dh =Qh-Qx
        dv =Qv-Qy
        rh =sin(np.pi*(Qh-Qx))/sin(np.pi*(Qh+Qx))
        rv =sin(np.pi*(Qv-Qy))/sin(np.pi*(Qv+Qy))
        rch=sin(np.pi*(Qh-Qy))/sin(np.pi*(Qh+Qy))
        rcv=sin(np.pi*(Qx-Qv))/sin(np.pi*(Qx+Qv))

        #-- Loop for files
        f1001Abs =np.zeros((len(bpm),len(FilesX)))
        f1010Abs =np.zeros((len(bpm),len(FilesX)))
        f1001xArg=np.zeros((len(bpm),len(FilesX)))
        f1001yArg=np.zeros((len(bpm),len(FilesX)))
        f1010xArg=np.zeros((len(bpm),len(FilesX)))
        f1010yArg=np.zeros((len(bpm),len(FilesX)))
        for i in range(len(FilesX)):

            #-- Read amplitudes and phases
            amph  =     np.array([FilesX[i].AMPX[FilesX[i].indx[b[1]]]    for b in bpm])
            ampv  =     np.array([FilesY[i].AMPY[FilesY[i].indx[b[1]]]    for b in bpm])
            amph01=     np.array([FilesX[i].AMP01[FilesX[i].indx[b[1]]]   for b in bpm])
            ampv10=     np.array([FilesY[i].AMP10[FilesY[i].indx[b[1]]]   for b in bpm])
            psih  =2*np.pi*np.array([FilesX[i].MUX[FilesX[i].indx[b[1]]]     for b in bpm])
            psiv  =2*np.pi*np.array([FilesY[i].MUY[FilesY[i].indx[b[1]]]     for b in bpm])
            psih01=2*np.pi*np.array([FilesX[i].PHASE01[FilesX[i].indx[b[1]]] for b in bpm])
            psiv10=2*np.pi*np.array([FilesY[i].PHASE10[FilesY[i].indx[b[1]]] for b in bpm])
            #-- I'm not sure this is correct for the coupling so I comment out this part for now (by RM 9/30/11).
            #for k in range(len(bpm)):
            #       try:
            #               if bpm[k][0]>s_lastbpm:
            #                       psih[k]  +=bd*2*np.pi*Qh  #-- To fix the phase shift by Qh
            #                       psiv[k]  +=bd*2*np.pi*Qv  #-- To fix the phase shift by Qv
            #                       psih01[k]+=bd*2*np.pi*Qv  #-- To fix the phase shift by Qv
            #                       psiv10[k]+=bd*2*np.pi*Qh  #-- To fix the phase shift by Qh
            #       except: pass

            #-- Construct Fourier components
            #   * be careful for that the note is based on x+i(alf*x*bet*x')).
            #   * Calculating Eqs (87)-(92) by using Eqs (47) & (48) (but in the Fourier space) in the note.
            #   * Note that amph(v)01 is normalized by amph(v) and it is un-normalized in the following.
            dpsih  =np.append(psih[1:]  ,2*np.pi*Qh+psih[0]  )-psih
            dpsiv  =np.append(psiv[1:]  ,2*np.pi*Qv+psiv[0]  )-psiv
            dpsih01=np.append(psih01[1:],2*np.pi*Qv+psih01[0])-psih01
            dpsiv10=np.append(psiv10[1:],2*np.pi*Qh+psiv10[0])-psiv10

            X_m10=2*amph*np.exp(-1j*psih)
            Y_0m1=2*ampv*np.exp(-1j*psiv)
            X_0m1=amph*np.exp(-1j*psih01)/(1j*sin(dpsih))*(amph01*np.exp(1j*dpsih)-np.append(amph01[1:],amph01[0])*np.exp(-1j*dpsih01))
            X_0p1=amph*np.exp( 1j*psih01)/(1j*sin(dpsih))*(amph01*np.exp(1j*dpsih)-np.append(amph01[1:],amph01[0])*np.exp( 1j*dpsih01))
            Y_m10=ampv*np.exp(-1j*psiv10)/(1j*sin(dpsiv))*(ampv10*np.exp(1j*dpsiv)-np.append(ampv10[1:],ampv10[0])*np.exp(-1j*dpsiv10))
            Y_p10=ampv*np.exp( 1j*psiv10)/(1j*sin(dpsiv))*(ampv10*np.exp(1j*dpsiv)-np.append(ampv10[1:],ampv10[0])*np.exp( 1j*dpsiv10))

            #-- Construct f1001hv, f1001vh, f1010hv (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
            f1001hv=-np.conjugate(1/(2j)*Y_m10/X_m10)  #-- - sign from the different def
            f1001vh=-1/(2j)*X_0m1/Y_0m1             #-- - sign from the different def
            f1010hv=-1/(2j)*Y_p10/np.conjugate(X_m10)  #-- - sign from the different def
            f1010vh=-1/(2j)*X_0p1/np.conjugate(Y_0m1)  #-- - sign from the different def
    ##              f1001hv=conjugate(1/(2j)*Y_m10/X_m10)
    ##              f1001vh=1/(2j)*X_0m1/Y_0m1
    ##              f1010hv=1/(2j)*Y_p10/conjugate(X_m10)
    ##              f1010vh=1/(2j)*X_0p1/conjugate(Y_0m1)

            #-- Construct phases psih, psiv, Psih, Psiv w.r.t. the AC dipole
            psih=psih-(psih[k_bpmac_h]-psih_ac2bpmac[bpmac_h])
            psiv=psiv-(psiv[k_bpmac_v]-psiv_ac2bpmac[bpmac_v])
            print('the phase to the device', k_bpmac_h, psih[k_bpmac_h], bpmac_h, (psih[k_bpmac_h]-psih_ac2bpmac[bpmac_h]))
            Psih=psih-np.pi*Qh
            Psih[:k_bpmac_h]=Psih[:k_bpmac_h]+2*np.pi*Qh
            Psiv=psiv-np.pi*Qv
            Psiv[:k_bpmac_v]=Psiv[:k_bpmac_v]+2*np.pi*Qv

            Psix=np.arctan((1-rh)/(1+rh)*np.tan(Psih))%np.pi
            Psiy=np.arctan((1-rv)/(1+rv)*np.tan(Psiv))%np.pi
            for k in range(len(bpm)):
                if Psih[k]%(2*np.pi)>np.pi: Psix[k]=Psix[k]+np.pi
                if Psiv[k]%(2*np.pi)>np.pi: Psiy[k]=Psiy[k]+np.pi

            psix=Psix-np.pi*Qx
            psix[k_bpmac_h:]=psix[k_bpmac_h:]+2*np.pi*Qx
            psiy=Psiy-np.pi*Qy
            psiy[k_bpmac_v:]=psiy[k_bpmac_v:]+2*np.pi*Qy

            #-- Construct f1001h, f1001v, f1010h, f1010v (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
            f1001h=1/math.sqrt(1-rv**2)*(np.exp(-1j*(Psiv-Psiy))*f1001hv+rv*np.exp( 1j*(Psiv+Psiy))*f1010hv)
            f1010h=1/math.sqrt(1-rv**2)*(np.exp( 1j*(Psiv-Psiy))*f1010hv+rv*np.exp(-1j*(Psiv+Psiy))*f1001hv)
            f1001v=1/math.sqrt(1-rh**2)*(np.exp( 1j*(Psih-Psix))*f1001vh+rh*np.exp(-1j*(Psih+Psix))*np.conjugate(f1010vh))
            f1010v=1/math.sqrt(1-rh**2)*(np.exp( 1j*(Psih-Psix))*f1010vh+rh*np.exp(-1j*(Psih+Psix))*np.conjugate(f1001vh))

            #-- Construct f1001 and f1010 from h and v BPMs (these include math.sqrt(betv/beth) or math.sqrt(beth/betv))
            g1001h          =np.exp(-1j*((psih-psih[k_bpmac_h])-(psiy-psiy[k_bpmac_v])))*(ampv/amph*amph[k_bpmac_h]/ampv[k_bpmac_v])*f1001h[k_bpmac_h]
            g1001h[:k_bpmac_h]=1/(np.exp(2*np.pi*1j*(Qh-Qy))-1)*(f1001h-g1001h)[:k_bpmac_h]
            g1001h[k_bpmac_h:]=1/(1-np.exp(-2*np.pi*1j*(Qh-Qy)))*(f1001h-g1001h)[k_bpmac_h:]

            g1010h          =np.exp(-1j*((psih-psih[k_bpmac_h])+(psiy-psiy[k_bpmac_v])))*(ampv/amph*amph[k_bpmac_h]/ampv[k_bpmac_v])*f1010h[k_bpmac_h]
            g1010h[:k_bpmac_h]=1/(np.exp(2*np.pi*1j*(Qh+Qy))-1)*(f1010h-g1010h)[:k_bpmac_h]
            g1010h[k_bpmac_h:]=1/(1-np.exp(-2*np.pi*1j*(Qh+Qy)))*(f1010h-g1010h)[k_bpmac_h:]

            g1001v          =np.exp(-1j*((psix-psix[k_bpmac_h])-(psiv-psiv[k_bpmac_v])))*(amph/ampv*ampv[k_bpmac_v]/amph[k_bpmac_h])*f1001v[k_bpmac_v]
            g1001v[:k_bpmac_v]=1/(np.exp(2*np.pi*1j*(Qx-Qv))-1)*(f1001v-g1001v)[:k_bpmac_v]
            g1001v[k_bpmac_v:]=1/(1-np.exp(-2*np.pi*1j*(Qx-Qv)))*(f1001v-g1001v)[k_bpmac_v:]

            g1010v          =np.exp(-1j*((psix-psix[k_bpmac_h])+(psiv-psiv[k_bpmac_v])))*(amph/ampv*ampv[k_bpmac_v]/amph[k_bpmac_h])*f1010v[k_bpmac_v]
            g1010v[:k_bpmac_v]=1/(np.exp(2*np.pi*1j*(Qx+Qv))-1)*(f1010v-g1010v)[:k_bpmac_v]
            g1010v[k_bpmac_v:]=1/(1-np.exp(-2*np.pi*1j*(Qx+Qv)))*(f1010v-g1010v)[k_bpmac_v:]

            f1001x=np.exp(1j*(psih-psix))*f1001h
            f1001x=f1001x-rh*np.exp(-1j*(psih+psix))/rch*np.conjugate(f1010h)
            f1001x=f1001x-2j*sin(np.pi*dh)*np.exp(1j*(Psih-Psix))*g1001h
            f1001x=f1001x-2j*sin(np.pi*dh)*np.exp(-1j*(Psih+Psix))/rch*np.conjugate(g1010h)
            f1001x=1/math.sqrt(1-rh**2)*sin(np.pi*(Qh-Qy))/sin(np.pi*(Qx-Qy))*f1001x

            f1010x=np.exp(1j*(psih-psix))*f1010h
            f1010x=f1010x-rh*np.exp(-1j*(psih+psix))*rch*np.conjugate(f1001h)
            f1010x=f1010x-2j*sin(np.pi*dh)*np.exp(1j*(Psih-Psix))*g1010h
            f1010x=f1010x-2j*sin(np.pi*dh)*np.exp(-1j*(Psih+Psix))*rch*np.conjugate(g1001h)
            f1010x=1/math.sqrt(1-rh**2)*sin(np.pi*(Qh+Qy))/sin(np.pi*(Qx+Qy))*f1010x

            f1001y=np.exp(-1j*(psiv-psiy))*f1001v
            f1001y=f1001y+rv*np.exp(1j*(psiv+psiy))/rcv*f1010v
            f1001y=f1001y+2j*sin(np.pi*dv)*np.exp(-1j*(Psiv-Psiy))*g1001v
            f1001y=f1001y-2j*sin(np.pi*dv)*np.exp(1j*(Psiv+Psiy))/rcv*g1010v
            f1001y=1/math.sqrt(1-rv**2)*sin(np.pi*(Qx-Qv))/sin(np.pi*(Qx-Qy))*f1001y

            f1010y=np.exp(1j*(psiv-psiy))*f1010v
            f1010y=f1010y+rv*np.exp(-1j*(psiv+psiy))*rcv*f1001v
            f1010y=f1010y-2j*sin(np.pi*dv)*np.exp(1j*(Psiv-Psiy))*g1010v
            f1010y=f1010y+2j*sin(np.pi*dv)*np.exp(-1j*(Psiv+Psiy))*rcv*g1001v
            f1010y=1/math.sqrt(1-rv**2)*sin(np.pi*(Qx+Qv))/sin(np.pi*(Qx+Qy))*f1010y

            #-- For B2, must be double checked
            if bd == -1:
                f1001x=-np.conjugate(f1001x)
                f1001y=-np.conjugate(f1001y)
                f1010x=-np.conjugate(f1010x)
                f1010y=-np.conjugate(f1010y)

            #-- Separate to amplitudes and phases, amplitudes averaged to cancel math.sqrt(betv/beth) and math.sqrt(beth/betv)
            for k in range(len(bpm)):
                f1001Abs[k][i] =math.sqrt(abs(f1001x[k]*f1001y[k]))
                f1010Abs[k][i] =math.sqrt(abs(f1010x[k]*f1010y[k]))
                f1001xArg[k][i]=np.angle(f1001x[k])%(2*np.pi)
                f1001yArg[k][i]=np.angle(f1001y[k])%(2*np.pi)
                f1010xArg[k][i]=np.angle(f1010x[k])%(2*np.pi)
                f1010yArg[k][i]=np.angle(f1010y[k])%(2*np.pi)

        #-- Output
        fwqw={}
        goodbpm=[]
        for k in range(len(bpm)):

            #-- Bad BPM flag based on phase
            badbpm=0
            f1001xArgAve = phase.calc_phase_mean(f1001xArg[k],2*np.pi)
            f1001yArgAve = phase.calc_phase_mean(f1001yArg[k],2*np.pi)
            f1010xArgAve = phase.calc_phase_mean(f1010xArg[k],2*np.pi)
            f1010yArgAve = phase.calc_phase_mean(f1010yArg[k],2*np.pi)
            #This seems to be to conservative or somethings...
            if min(abs(f1001xArgAve-f1001yArgAve),2*np.pi-abs(f1001xArgAve-f1001yArgAve))>np.pi/2: badbpm=1
            if min(abs(f1010xArgAve-f1010yArgAve),2*np.pi-abs(f1010xArgAve-f1010yArgAve))>np.pi/2: badbpm=1


            #-- Output
	    badbpm=0
            if badbpm==0:
                f1001AbsAve = np.mean(f1001Abs[k])
                f1010AbsAve = np.mean(f1010Abs[k])
                f1001ArgAve = phase.calc_phase_mean(np.append(f1001xArg[k],f1001yArg[k]),2*np.pi)
                f1010ArgAve = phase.calc_phase_mean(np.append(f1010xArg[k],f1010yArg[k]),2*np.pi)
                f1001Ave = f1001AbsAve*np.exp(1j*f1001ArgAve)
                f1010Ave = f1010AbsAve*np.exp(1j*f1010ArgAve)
                f1001AbsStd = math.sqrt(np.mean((f1001Abs[k]-f1001AbsAve)**2))
                f1010AbsStd = math.sqrt(np.mean((f1010Abs[k]-f1010AbsAve)**2))
                f1001ArgStd = phase.calc_phase_std(np.append(f1001xArg[k],f1001yArg[k]),2*np.pi)
                f1010ArgStd = phase.calc_phase_std(np.append(f1010xArg[k],f1010yArg[k]),2*np.pi)
                fwqw[bpm[k][1]] = [[f1001Ave          ,f1001AbsStd       ,f1010Ave          ,f1010AbsStd       ],
                                 [f1001ArgAve/(2*np.pi),f1001ArgStd/(2*np.pi),f1010ArgAve/(2*np.pi),f1010ArgStd/(2*np.pi)]]  #-- Phases renormalized to [0,1)
                goodbpm.append(bpm[k])

    #-- Global parameters not implemented yet

        fqwList.append(fwqw)


    fwqw = copy.deepcopy(fqwList[0])
    for key in fwqw:
        for a in range(1, len(fqwList)):
            tmp = fqwList[a]
            fwqw[key][0][0] = fwqw[key][0][0] + tmp[key][0][0]
            fwqw[key][0][1] = fwqw[key][0][1] + tmp[key][0][1]
            fwqw[key][0][2] = fwqw[key][0][2] + tmp[key][0][2]
            fwqw[key][0][3] = fwqw[key][0][3] + tmp[key][0][3]
            if(key is 'BPMWB.4R5.B1'):
                print fwqw[key][1]
        fwqw[key][0][0]=fwqw[key][0][0]/len(fqwList)
        fwqw[key][0][1]=fwqw[key][0][1]/len(fqwList)
        fwqw[key][0][2]=fwqw[key][0][2]/len(fqwList)
        fwqw[key][0][3]=fwqw[key][0][3]/len(fqwList)
    fwqw['Global']=['"null"','"null"']
    return [fwqw,goodbpm]

def _get_columns(list_of_df, column):
    x_frame = list_of_df[0][[column]]
    for i, df in enumerate(list_of_df[1:]):
        x_frame = pd.merge(x_frame,
                        df[[column]],
                        how="inner", left_index=True,
                        right_index=True, suffixes=('', str(i + 1)))
    return x_frame
