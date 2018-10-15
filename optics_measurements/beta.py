'''
.. module: beta
Created on 27 May 2013

@author: awegsche, vimaier

@version: 2016.11.p2

GetLLM.algorithms.beta.py stores helper functions for phase calculations for GetLLM.
This module is not intended to be executed. It stores only functions.


'''
import sys
import os
import math
import time
import re

import numpy as np
from numpy import sin, tan
import pandas as pd

from scipy.linalg import circulant
from model.accelerators.accelerator import AccExcitationMode
from tfs_files import tfs_pandas
from utils import logging_tools

__version__ = "2018.7.a"

DEBUG = sys.flags.debug  # True with python option -d! ("python -d GetLLM.py...") (vimaier)
LOGGER = logging_tools.get_logger(__name__)

if DEBUG:
    from utils import debug_algorithms as DBG

# --- Constants

TWOPI = 2 * np.pi
DEFAULT_WRONG_BETA = 1000
EPSILON = 1.0E-16
ZERO_THRESHOLD = 1e-3
COT_THRESHOLD = 15.9
RCOND = 1.0e-10

BOXLENGTH = 50
BOXINDENT = 4
CALCULATE_BETA_HOR = True
CALCULATE_BETA_VER = True
APPLY = True

# --------------- Errors method
METH_IND = "Invalid"
METH_3BPM = "3BPM method"
METH_A_NBPM = "Analytical N-BPM method"
METH_NO_ERR = "No Errors"
METH_A_NBPM_PART = "Partly Analytical N-BPM method"


def rms(arr):  # TODO take this from the stats module
    return np.sqrt(np.mean(np.square(arr)))
# --------------------------------------------------------------------------------------------------
# main part
# --------------------------------------------------------------------------------------------------


def calculate_beta_from_phase(getllm_d, tune_d, phase_d, header_dict):
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
    # setting up
    accelerator = getllm_d.accelerator

    # selecting models -----------------------------------------------------------------------------

    free_model = accelerator.get_model_tfs()
    elements = accelerator.get_elements_tfs()

    # the following tries to get the best_knowledge model
    # if it doesn't find it, it takes the base model
    try:
        free_bk_model = accelerator.get_best_knowledge_model_tfs()
    except AttributeError:
        LOGGER.debug("No best knowledge model - using the normal one.")
        free_bk_model = free_model

    driven_model = None
    if accelerator.excitation != AccExcitationMode.FREE:
        # in the case of driven motion, we need the driven model as well
        driven_model = accelerator.get_driven_tfs()

    # print information to info and debug ----------------------------------------------------------
    LOGGER.info("Calculating beta from phase")
    LOGGER.info("Version: {0:5s}".format(__version__))

    LOGGER.info("range of BPMs: {}".format(getllm_d.range_of_bpms))
    LOGGER.info("cot of phase threshold: {:g}".format(COT_THRESHOLD))

    LOGGER.debug("quad field errors: [YES]")
    LOGGER.debug("quad long misalignments: [YES]")
    LOGGER.debug("sext transverse misalignments: [YES]")
    LOGGER.debug("BPM long misalignments: [YES]")
    LOGGER.debug("dipole K1 errors: [ NO]")
    LOGGER.debug("analytical alpha: [ YES]")

    # check whether analytical N-BPM method should be used
    # if yes, uncertainty estimates will be distributed to the elements
    error_method = METH_IND

    if getllm_d.three_bpm_method:
        error_method = METH_3BPM
    else:
        LOGGER.debug("Accelerator Error Definition")
        error_defs_path = getllm_d.accelerator.get_errordefspath()
        if error_defs_path is None:
            raise IOError("Error definition file '{}' could not be found"
                          .format(getllm_d.accelerator.get_errordefspath()))

        elements = _assign_uncertainties(elements, error_defs_path)
        error_method = METH_A_NBPM

    # start the calculation per plane --------------------------------------------------------------

    # ------------- HORIZONTAL
    if phase_d["X"]["F"]:
        beta_df_x, compensated_beta_df_x = _beta_from_phase_for_plane(
            free_model, driven_model, free_bk_model, elements,
            getllm_d.range_of_bpms, phase_d, error_method, tune_d, "X"
        )

    # ------------- VERTICAL
    if phase_d["Y"]["F"]:
        beta_df_y, compensated_beta_df_y = _beta_from_phase_for_plane(
            free_model, driven_model, free_bk_model, elements,
            getllm_d.range_of_bpms, phase_d, error_method, tune_d, "Y"
        )

    for df in [beta_df_x, compensated_beta_df_x, beta_df_y, compensated_beta_df_y]:
        if df is not None:
            _add_header(df, header_dict, error_method, getllm_d.range_of_bpms)
            LOGGER.debug("writing %s", df.headers["FILENAME"])
            tfs_pandas.write_tfs(os.path.join(getllm_d.outputdir, df.headers["FILENAME"]), df,
                                 save_index="NAME")

    return beta_df_x, compensated_beta_df_x, beta_df_y, compensated_beta_df_y


def _beta_from_phase_for_plane(free_model, driven_model, bk_model, elements, range_of_bpms,
                               phases, error_method, tunes, plane):
    """
    This function calculates and outputs the beta function measurement for the given plane.
    """
    plane_for_file = plane.lower()
    Q = tunes[plane]["Q"]
    Qf = tunes[plane]["QF"]
    Qmdl = tunes[plane]["QM"]
    Qmdlf = tunes[plane]["QFM"]
    phase_adv_free = phases[plane]["F"]
    phase_adv_driven = phases[plane]["D"]
    LOGGER.info("Beta {} free calculation".format(plane))
    # remove BPMs that are not in the input
    if phase_adv_driven is not None:
        comp_model = free_model.loc[phase_adv_driven["MEAS"].index]
        model = driven_model.loc[phase_adv_free["MEAS"].index]
    else:
        model = free_model.loc[phase_adv_free["MEAS"].index]

    bk_model = bk_model.loc[phase_adv_free["MEAS"].index]

    # if DEBUG create a binary debugfile where the algorithm is writing matrices, beta-values,
    # weights etc.
    if DEBUG:
        DBG.create_debugfile(
            "getbeta{}.bdebug".format(plane_for_file)  # TODO change working path
        )

    beta_df = _beta_from_phase(model, elements, phase_adv_free, plane, range_of_bpms,
                               error_method, Qf, Qmdlf % 1.0)

    beta_df.headers["FILENAME"] = "getbeta{}.out".format(plane_for_file)
    if DEBUG:
        DBG.close_file()

    compensated_beta_df = None

    if phase_adv_driven is not None:
        LOGGER.info("Beta {} compensated calculation".format(plane))
        if DEBUG:
            DBG.create_debugfile(
                "getbeta{}_free.bdebug".format(plane_for_file)  # TODO change working path
            )

        compensated_beta_df = _beta_from_phase(
            comp_model, elements,
            phase_adv_free, plane, range_of_bpms, error_method, Q, Qmdl % 1.0
        )
        compensated_beta_df.headers["FILENAME"] = "getbeta{}_free.out".format(plane_for_file)

        if DEBUG:
            DBG.close_file()

        LOGGER.warning("Skip free2 calculation")

    # add filename to header
    return beta_df, compensated_beta_df


def _beta_from_phase(madTwiss, madElements, phase, plane,
                     range_of_bpms, errors_method, tune, mdltune):
    '''
    Calculate the beta function from phase advances.

    Parameters:
        madTwiss: model tfs file
        madElements: model tfs file with all relevant elements (quadrupoles, sextupoles, drifts)
        phase: matrix with phase advances
        plane: 'X' or 'Y'
        getllm_d: GetLLM_Data
        errors_method: 3BPM or N-BPM method
        tune: measured tune
        mdltune: model tune
    '''
    plane_bet = "BET" + plane
    plane_alf = "ALF" + plane
    st = time.time()

    beta_df = tfs_pandas.TfsDataFrame(madTwiss).loc[phase["MEAS"].index, ["S", plane_bet, plane_alf]]

    beta_df = beta_df.rename(columns={plane_bet: plane_bet + "MDL", plane_alf: plane_alf + "MDL"})

    if errors_method == METH_A_NBPM:
        beta_df = _scan_all_BPMs_withsystematicerrors(madTwiss, madElements, phase, plane,
                                                      range_of_bpms,
                                                      tune, mdltune,
                                                      beta_df)
    # ---- use the simulations
    else:
        beta_df = _scan_all_BPMs_3bpm(phase, plane, tune, mdltune, beta_df)

    rmsbb = rms(beta_df["BBEAT" + plane]) * 100
    et = time.time() - st
    beta_df.headers["RMS_BETABEAT"] = "{:.3f} %".format(rmsbb)
    beta_df.headers["CALCULATION_TIME"] = "{:.2f} s".format(et)
    LOGGER.info(" - RMS beta beat: {:.2f}%".format(rmsbb))
    LOGGER.info(" - elapsed time: {:.2f}s".format(et))
    beta_df["DELTABET"+plane] = (beta_df.loc[:, 'BET' + plane] /
                                 beta_df.loc[:, 'BET' + plane + 'MDL'] - 1.0)

    # check if there were actually errors assigned
    if (len(madElements["dK1"].nonzero()[0]) + len(madElements["dX"].nonzero()[0]) +
        len(madElements["KdS"].nonzero()[0])) > 0:
        if -2 in beta_df["NCOMB"].values:
            LOGGER.warning("Analytical N-BPM method failed for at least one BPM. "
                           "Check debug log for more information.")
        if 0 in beta_df["NCOMB"].values:
            LOGGER.warning("No combinations left for at least one BPM. "
                           "Check debug log for more information.")
            errors_method = METH_A_NBPM
        beta_df = beta_df.loc[beta_df["NCOMB"] > 0]
    else:
        errors_method = METH_NO_ERR
        LOGGER.warning("No systematic errors were given or no element was found for the given "
                       "error definitions. The analytical N-BPM method was not used.")

    LOGGER.info("Errors from " + errors_method)
    beta_df.headers["ErrorsFrom:"] = errors_method
    return beta_df


# --------------------------------------------------------------------------------------------------
# ---------------- calculate beta and alpha using the old 3 BPM method -----------------------------
# --------------------------------------------------------------------------------------------------

def _scan_all_BPMs_3bpm(phase, plane, tune, mdltune,
                        beta_df):
    '''
    Calculates beta from phase using the old 3-BPM method

    ``phase["MEAS"]``, ``phase["MODEL"]``, ``phase["ERRMEAS"]`` (from ``get_phases``) are of the
    form:

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
    # tilt phase advances in order to have the phase advances in a neighbourhood
    tilted_meas = tilt_slice_matrix(phase["MEAS"].as_matrix(), 2, 5, tune) * TWOPI
    tilted_model = tilt_slice_matrix(phase["MODEL"].as_matrix(), 2, 5, mdltune) * TWOPI
    tilted_errmeas = tilt_slice_matrix(phase["ERRMEAS"].as_matrix(), 2, 5, mdltune) * TWOPI

    betmdl = beta_df.loc[:]["BET" + plane + "MDL"].as_matrix()
    alfmdl = beta_df.loc[:]["ALF" + plane + "MDL"].as_matrix()

    # ------ main part, calculate the beta and alpha function ----------------------------------

    # calculate cotangens of all the phase advances in the neighbourhood
    with np.errstate(divide='ignore'):
        cot_phase_meas = 1.0 / tan(tilted_meas)
        cot_phase_model = 1.0 / tan(tilted_model)

    # calculate enumerators and denominators for far more cases than needed
    # shift1 are the cases BBA, ABB, AxBB, AxxBB etc. (the used BPMs are adjacent)
    # shift2 are the cases where the used BPMs are separated by one. only BAB is used for  3-BPM
    cot_phase_meas_shift1 = cot_phase_meas - np.roll(cot_phase_meas, -1, axis=0)
    cot_phase_model_shift1 = cot_phase_model - np.roll(cot_phase_model, -1, axis=0) + EPSILON
    cot_phase_meas_shift2 = cot_phase_meas - np.roll(cot_phase_meas, -2, axis=0)
    cot_phase_model_shift2 = cot_phase_model - np.roll(cot_phase_model, -2, axis=0) + EPSILON

    # calculate the sum of the fractions
    bet_frac = (cot_phase_meas_shift1[0]/cot_phase_model_shift1[0] +
                cot_phase_meas_shift1[3]/cot_phase_model_shift1[3] +
                cot_phase_meas_shift2[1]/cot_phase_model_shift2[1]) / 3.0

    # multiply the fractions by betmdl and calculate the arithmetic mean
    beti = bet_frac * betmdl

    # alpha
    alfi = (bet_frac * (cot_phase_model[1] + cot_phase_model[3] + 2.0 * alfmdl) -
            (cot_phase_meas[1] + cot_phase_meas[3])) / 2.0

    # ------ error propagation -----------------------------------------------------------------

    # error = sqrt( errphi_ij^2 * (d beta / dphi_ij)^2 )
    # calculate sin(phimdl_ij)
    sin_model = sin(tilted_model)
    # calculate errphi_ij^2 / sin^2 phimdl_ij * beta
    with np.errstate(divide='ignore', invalid='ignore'):
        sin_squared_model = tilted_errmeas / np.multiply(sin_model, sin_model) * betmdl
    # square it again beacause it's used in a vector length
    sin_squared_model = np.multiply(sin_squared_model, sin_squared_model)

    sin_squ_model_shift1 = sin_squared_model + \
        np.roll(sin_squared_model, -1, axis=0) / \
        np.multiply(cot_phase_model_shift1, cot_phase_model_shift1)
    sin_squ_model_shift2 = sin_squared_model + \
        np.roll(sin_squared_model, -2, axis=0) / \
        np.multiply(cot_phase_model_shift2, cot_phase_model_shift2)
    beterr = np.sqrt(sin_squ_model_shift1[0] + sin_squ_model_shift1[3] +
                     sin_squ_model_shift2[1]) \
        / 3.0

    alferr = 0  # TODO calculate alferr

    bb = bet_frac - 1.0

    beta_df["BET" + plane] = beti
    beta_df["SYSBET" + plane] = 0
    beta_df["STATBET" + plane] = beterr
    beta_df["ERRBET" + plane] = beterr

    beta_df["ALF" + plane] = alfi
    beta_df["SYSALF" + plane] = 0
    beta_df["STATALF" + plane] = alferr
    beta_df["ERRALF" + plane] = alferr

    beta_df["BBEAT" + plane] = bb
    beta_df["COUNT"] = 0
    return beta_df


# --------------------------------------------------------------------------------------------------
# -------- using analytical formula ----------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

def _scan_all_BPMs_withsystematicerrors(madTwiss, madElements,
                                        phase, plane, range_of_bpms, tune, mdltune, beta_df):
    '''
    '''

    LOGGER.debug("starting scan_all_BPMs_withsystematicerrors")
    # ---------- setup -----------------------------------------------------------------------------
    # setup combinations
    width = range_of_bpms / 2
    left_bpm = range(-width, 0)
    right_bpm = range(0 + 1, width + 1)
    BBA_combo = [[x, y] for x in left_bpm for y in left_bpm if x < y]
    ABB_combo = [[x, y] for x in right_bpm for y in right_bpm if x < y]
    BAB_combo = [[x, y] for x in left_bpm for y in right_bpm]

    # for fast access
    phases_meas = phase["MEAS"] * TWOPI
    phases_err = phase["ERRMEAS"] * TWOPI
    phases_err.where(phases_err.notnull(), 1, inplace=True)

    result = np.ndarray(len(phases_meas.index), [('beti', float), ('betstat', float),
                                                 ('betsys', float), ('beterr', float),
                                                 ('alfi', float), ('alfstat', float),
                                                 ('alfsys', float), ('alferr', float),
                                                 ('corr', float), ('ncomb', int)])
    # ---------- calculate the betas --------------------------------------------------------------

    for i in range(0, len(phases_meas.index)):
        row = _scan_one_BPM_withsystematicerrors(madTwiss, madElements,
                                                 phases_meas, phases_err,
                                                 plane, range_of_bpms,
                                                 i,
                                                 BBA_combo, ABB_combo, BAB_combo,
                                                 tune, mdltune)
        result[row[0]] = row[1:]

    beta_df["BET" + plane] = result["beti"]
    beta_df["STATBET" + plane] = result["betstat"]
    beta_df["SYSBET" + plane] = result["betsys"]
    beta_df["ERRBET" + plane] = result["beterr"]
    beta_df["ALF" + plane] = result["alfi"]
    beta_df["STATALF" + plane] = result["alfstat"]
    beta_df["SYSALF" + plane] = result["alfsys"]
    beta_df["ERRALF" + plane] = result["alferr"]
    beta_df["CORR"] = result["corr"]
    beta_df["NCOMB"] = result["ncomb"]
    beta_df["BBEAT" + plane] = beta_df["BET" + plane] / beta_df["BET" + plane + "MDL"] - 1

    return beta_df


def scan_several_BPMs_withsystematicerrors(madTwiss, madElements, cot_meas, phases_err, plane,
                                           range_of_bpms, begin, end, BBA_combo, ABB_combo,
                                           BAB_combo, tune, mdltune):
    block = []
    for i in range(begin, end):
        block.append(_scan_one_BPM_withsystematicerrors(madTwiss, madElements,
                                                        cot_meas, phases_err,
                                                        plane, range_of_bpms,
                                                        i,
                                                        BBA_combo, ABB_combo, BAB_combo,
                                                        tune, mdltune))
    return block


def _scan_one_BPM_withsystematicerrors(madTwiss, madElements,
                                       phases_meas, phases_err,
                                       plane, range_of_bpms, Index, BBA_combo, ABB_combo, BAB_combo,
                                       tune, mdltune):
    '''
    Scans the range of BPMs in order to get the final value for one BPM in the lattice
    '''
    probed_bpm_name = madTwiss.index[Index]
    s = madTwiss.at[probed_bpm_name, "S"]

    betmdl1 = madTwiss.at[probed_bpm_name, "BET" + plane]
    alfmdl1 = madTwiss.at[probed_bpm_name, "ALF" + plane]
    mu_column = "MU" + plane
    bet_column = "BET" + plane

    beti = DEFAULT_WRONG_BETA
    betstat = .0
    betsys = .0
    beterr = DEFAULT_WRONG_BETA
    alfi = DEFAULT_WRONG_BETA
    alfstat = .0
    alfsys = .0
    alferr = DEFAULT_WRONG_BETA

    m = range_of_bpms / 2
    indx_first = Index - m
    indx_last = Index + m
    name_first = madTwiss.index[indx_first]
    name_last = madTwiss.index[indx_last % len(madTwiss.index)]
    probed_bpm_name = madTwiss.index[Index]
    len_bpms_total = phases_meas.shape[0]

    indx_el_first = madElements.index.get_loc(name_first)
    indx_el_last = madElements.index.get_loc(name_last)

    if indx_first < 0:
        outerMeasPhaseAdv = pd.concat((
            phases_meas.iloc[Index, indx_first % len_bpms_total:] - tune * TWOPI,
            phases_meas.iloc[Index, :indx_last+1]))
        outerMeasErr = pd.concat((
            phases_err.iloc[Index, indx_first % len_bpms_total:],
            phases_err.iloc[Index, :indx_last+1]))
        outerMdlPh = np.concatenate((
            madTwiss.iloc[indx_first % len_bpms_total:][mu_column] - mdltune,
            madTwiss.iloc[:indx_last+1][mu_column])) * TWOPI
        outerElmts = pd.concat((
            madElements.iloc[indx_el_first:],
            madElements.iloc[:indx_el_last + 1]))
        outerElmtsPh = np.concatenate((
            madElements.iloc[indx_el_first:][mu_column] - mdltune,
            madElements.iloc[:indx_el_last + 1][mu_column])) * TWOPI

    elif indx_last >= len_bpms_total:
        outerMeasPhaseAdv = pd.concat((
            phases_meas.iloc[Index, indx_first:],
            phases_meas.iloc[Index, :(indx_last + 1) % len_bpms_total] + tune * TWOPI))
        outerMeasErr = pd.concat((
            phases_err.iloc[Index, indx_first:],
            phases_err.iloc[Index, :(indx_last + 1) % len_bpms_total]))
        outerMdlPh = np.concatenate((
            madTwiss.iloc[indx_first:][mu_column],
            madTwiss.iloc[:(indx_last + 1) % len_bpms_total][mu_column] + mdltune)) * TWOPI
        outerElmts = pd.concat((
            madElements.iloc[indx_el_first:],
            madElements.iloc[:indx_el_last + 1]))
        outerElmtsPh = np.concatenate((
            madElements.iloc[indx_el_first:][mu_column],
            madElements.iloc[:indx_el_last + 1][mu_column] + mdltune)) * TWOPI

    else:
        outerMeasPhaseAdv = phases_meas.iloc[Index, indx_first: indx_last + 1]
        outerMeasErr = phases_err.iloc[Index, indx_first: indx_last + 1]
        outerMdlPh = madTwiss.iloc[indx_first:indx_last + 1][mu_column].as_matrix() * TWOPI
        outerElmts = madElements.iloc[indx_el_first:indx_el_last + 1]
        outerElmtsPh = madElements.iloc[indx_el_first:indx_el_last + 1][mu_column] * TWOPI

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

    betas = np.empty(len(BBA_combo) + len(BAB_combo) + len(ABB_combo))
    alfas = np.empty(len(BBA_combo) + len(BAB_combo) + len(ABB_combo))
    beta_mask = np.empty(len(BBA_combo) + len(BAB_combo) + len(ABB_combo), dtype=bool)

    diag = np.concatenate((outerMeasErr.as_matrix(), outerElmts.loc[:]["dK1"],
                           outerElmts.loc[:]["dX"], outerElmts.loc[:]["KdS"],
                           outerElmts.loc[:]["mKdS"]))
    mask = diag != 0

    T_Beta = np.zeros((len(betas),
                       len(diag)))
    T_Alfa = np.zeros((len(betas),
                       len(diag)))

    M = np.diag(diag[mask])
    line_length = len(diag)

    for i, combo in enumerate(BBA_combo):
        ix = combo[0] + m
        iy = combo[1] + m
        beta, alfa, betaline, alfaline = _calculate_beta_and_alfa_for_comb(
            ix, iy, sin_squared_elements, outerElmts, outerElmtsBet, outerElK2, cot_model, cot_meas,
            outerMeasPhaseAdv, combo, indx_el_probed, line_length, betmdl1, alfmdl1,
            range_of_bpms, m,
            1.0, -1.0, 1.0, -1.0)
        if beta > 0:
            T_Beta[i] = betaline
            T_Alfa[i] = alfaline
            betas[i] = beta
            alfas[i] = alfa
            beta_mask[i] = True
        else:
            beta_mask[i] = False

    for j, combo in enumerate(BAB_combo):
        ix = combo[0] + m
        iy = combo[1] + m
        i = j + len(BBA_combo)

        beta, alfa, betaline, alfaline = _calculate_beta_and_alfa_for_comb(
            ix, iy, sin_squared_elements, outerElmts, outerElmtsBet, outerElK2, cot_model, cot_meas,
            outerMeasPhaseAdv, combo, indx_el_probed, line_length, betmdl1, alfmdl1,
            range_of_bpms, m,
            1.0, 1.0, 1.0, 1.0)
        if beta > 0:
            T_Beta[i] = betaline
            T_Alfa[i] = alfaline
            betas[i] = beta
            alfas[i] = alfa
            beta_mask[i] = True
        else:
            beta_mask[i] = False

    for j, combo in enumerate(ABB_combo):
        ix = combo[0] + m
        iy = combo[1] + m

        i = j + len(BBA_combo) + len(BAB_combo)

        beta, alfa, betaline, alfaline = _calculate_beta_and_alfa_for_comb(
            ix, iy, sin_squared_elements, outerElmts, outerElmtsBet, outerElK2, cot_model, cot_meas,
            outerMeasPhaseAdv, combo, indx_el_probed, line_length, betmdl1, alfmdl1,
            range_of_bpms, m,
            -1.0, +1.0, -1.0, 1.0)

        if beta > 0:
            T_Beta[i] = betaline
            T_Alfa[i] = alfaline
            betas[i] = beta
            alfas[i] = alfa
            beta_mask[i] = True
        else:
            beta_mask[i] = False
    T_Beta = T_Beta[:, mask]
    T_Beta = T_Beta[beta_mask]
    betas = betas[beta_mask]
    V_Beta = np.dot(T_Beta, np.dot(M, np.transpose(T_Beta)))

    T_Alfa = T_Alfa[:, mask]
    T_Alfa = T_Alfa[beta_mask]
    alfas = alfas[beta_mask]
    V_Alfa = np.dot(T_Alfa, np.dot(M, np.transpose(T_Alfa)))
    
    if np.any(V_Beta):
        V_Beta_inv = np.linalg.pinv(V_Beta, rcond=RCOND)
        w = np.sum(V_Beta_inv, axis=1)
        VBeta_inv_sum = np.sum(w)
        if VBeta_inv_sum == 0:
            raise ValueError
        beterr = math.sqrt(float(np.dot(np.transpose(w), np.dot(V_Beta, w)) / VBeta_inv_sum ** 2))
        beti = float(np.dot(np.transpose(w), betas) / VBeta_inv_sum)
        ncomb = len(betas)
    elif len(betas) == 0:
        LOGGER.debug("No combinations left for {}.".format(probed_bpm_name))
        beti = DEFAULT_WRONG_BETA
        ncomb = 0
    else:
        LOGGER.debug("ValueError at {} in Beta calculation.".format(probed_bpm_name))
        LOGGER.debug("betas:\n" + str(betas))
        beti = np.average(betas)
        ncomb = -2
    if np.any(V_Alfa):
        V_Alfa_inv = np.linalg.pinv(V_Alfa, rcond=RCOND)
        w = np.sum(V_Alfa_inv, axis=1)
        VAlfa_inv_sum = np.sum(w)
        if VAlfa_inv_sum == 0:
            raise ValueError
        alferr = math.sqrt(float(np.dot(np.transpose(w), np.dot(V_Alfa, w)) / VAlfa_inv_sum ** 2))
        alfi = float(np.dot(np.transpose(w), alfas) / VAlfa_inv_sum)
        ncomb = len(betas)
    elif len(alfas) == 0:
        LOGGER.debug("No combinations left for {}.".format(probed_bpm_name))
        alfi = DEFAULT_WRONG_BETA
        ncomb = 0
    else:
        LOGGER.debug("ValueError at {} in Alfa calculation.".format(probed_bpm_name))
        LOGGER.debug("alfas:\n" + str(alfas))
        alfi = np.average(alfas)
        ncomb = -2

    # ----------------------------------------------------------------------------------------------
    # writing debug output
    # ----------------------------------------------------------------------------------------------
    if DEBUG:
        DBG.start_write_bpm(probed_bpm_name, s, beti, alfi, 0)
        DBG.write_matrix(T_Beta, "T_Beta")
        DBG.start_write_combinations(len(betas))
        combs = np.r_[BBA_combo, BAB_combo, ABB_combo]  # Stackexchange
        combs = combs[beta_mask]
        for n, beta_of_comb in enumerate(betas):
            DBG.write_bpm_combination(combs[n][0], combs[n][1], beta_of_comb, w[n] / VBeta_inv_sum)
        DBG.write_double("PHI{}MDL".format(plane), outerMdlPh[m])

        DBG.write_end()

    return (
        Index,
        beti, betstat, betsys, beterr,
        alfi, alfstat, alfsys, alferr,
        .0,
        ncomb
    )


def _calculate_beta_and_alfa_for_comb(ix, iy, sin_squared_elements, outerElmts, outerElmtsBet,
                                      outerElK2, cot_model, cot_meas, outerMeasPhaseAdv,
                                      combo, indx_el_probed, line_length, betmdl1, alfmdl1,
                                      range_of_bpms, m, fac1, fac2, sfac1, sfac2):
    """Calculates beta and alpha function as well as the respective covariance matrix lines for the
    given combination
    """
    betaline = np.zeros((line_length))
    alfaline = np.zeros((line_length))

    # remove bad combination
    if (abs(cot_model[ix]) > COT_THRESHOLD or
        abs(cot_model[iy]) > COT_THRESHOLD or
        abs(cot_meas[ix]) > COT_THRESHOLD or
        abs(cot_meas[iy]) > COT_THRESHOLD or
        abs(cot_model[ix] - cot_model[iy]) < ZERO_THRESHOLD):
        return -1.0, -1.0, None, None

    # calculate beta
    denom = (cot_model[ix] - cot_model[iy]) / betmdl1
    denomalf = denom * betmdl1 + 2 * alfmdl1
    beta_i = (cot_meas[ix] - cot_meas[iy]) / denom
    alfa_i = 0.5 * (denomalf * beta_i / betmdl1 -
                    (cot_meas[ix] + cot_meas[iy]))

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
    betaline[ix] = -1.0 / (denom_sinx * denom)
    betaline[iy] = 1.0 / (denom_siny * denom)

    alfaline[ix] = -1.0 / (denom_sinx * denom * betmdl1) * denomalf + 1.0 / denom_sinx
    alfaline[iy] = 1.0 / (denom_siny * denom * betmdl1) * denomalf + 1.0 / denom_siny

    # apply quadrupolar field uncertainty (quadrupole longitudinal misalignment already included)

    bet_sin_ix = elementPh_XA * elementBet_XA / (denom_sinx * denom)
    bet_sin_iy = elementPh_YA * elementBet_YA / (denom_siny * denom)

    betaline[xloc+range_of_bpms:indx_el_probed+range_of_bpms] += fac1 * bet_sin_ix
    betaline[yloc+range_of_bpms:indx_el_probed+range_of_bpms] += fac2 * bet_sin_iy

    alfaline[xloc+range_of_bpms:indx_el_probed+range_of_bpms] += fac1 * (
        .5 * (bet_sin_ix * denomalf + bet_sin_ix / betmdl1 * (cot_meas[ix] - cot_meas[iy])))

    alfaline[yloc+range_of_bpms:indx_el_probed+range_of_bpms] += fac2 * (
        .5 * (bet_sin_iy * denomalf + bet_sin_iy / betmdl1 * (cot_meas[ix] - cot_meas[iy])))

    y_offset = range_of_bpms + len(outerElmts)

    # apply sextupole transverse misalignment
    betaline[xloc + y_offset: indx_el_probed + y_offset] += fac1 * elementK2_XA * bet_sin_ix
    betaline[yloc + y_offset: indx_el_probed + y_offset] += fac2 * elementK2_YA * bet_sin_iy

    alfaline[xloc + y_offset: indx_el_probed + y_offset] += sfac1 * elementK2_XA * bet_sin_ix
    alfaline[yloc + y_offset: indx_el_probed + y_offset] += sfac2 * elementK2_YA * bet_sin_iy

    y_offset += len(outerElmts)

    # apply quadrupole longitudinal misalignments
    betaline[xloc + y_offset: indx_el_probed + y_offset] += fac1 * bet_sin_ix
    betaline[yloc + y_offset: indx_el_probed + y_offset] += fac2 * bet_sin_iy

    alfaline[xloc + y_offset: indx_el_probed + y_offset] += fac1 * (
        .5 * elementK2_XA * (bet_sin_ix * denomalf + bet_sin_ix / betmdl1 * (cot_meas[ix] -
                                                                             cot_meas[iy])))
    alfaline[yloc + y_offset: indx_el_probed + y_offset] += fac2 * (
        .5 * elementK2_YA * (bet_sin_iy * denomalf + bet_sin_iy / betmdl1 * (cot_meas[ix] -
                                                                             cot_meas[iy])))

    y_offset += len(outerElmts)

    betaline[xloc + y_offset: indx_el_probed + y_offset] -= fac1 * bet_sin_ix
    betaline[yloc + y_offset: indx_el_probed + y_offset] -= fac2 * bet_sin_iy

    alfaline[xloc + y_offset: indx_el_probed + y_offset] -= fac1 * (
        .5 * (bet_sin_ix * denomalf + bet_sin_ix / betmdl1 * (cot_meas[ix] - cot_meas[iy])))
    alfaline[yloc + y_offset: indx_el_probed + y_offset] -= fac2 * (
        .5 * (bet_sin_iy * denomalf + bet_sin_iy / betmdl1 * (cot_meas[ix] - cot_meas[iy])))

    return beta_i, alfa_i, betaline, alfaline


def _assign_uncertainties(twiss_full, errordefspath):
    '''
    Adds uncertainty information to twiss_full.

    :Sources of Errors:
        dK1:    quadrupolar field errors
        dS:     quadrupole longitudinal misalignments
        dX:     sextupole transverse misalignments
        BPMdS:  BPM longitudinal misalignments
    '''

    LOGGER.debug("Start creating uncertainty information")

    errdefs = tfs_pandas.read_tfs(errordefspath)

    # create new columns
    twiss_full = twiss_full.assign(UNC=False, dK1=0, KdS=0, mKdS=0, dX=0, BPMdS=0)

    # loop over uncertainty definitions, fill the respective columns, set UNC to true
    for indx in errdefs.index:
        patt = errdefs.loc[indx, "PATTERN"]
        if patt.startswith("key:"):
            LOGGER.debug("creating uncertainty information for {:s}".format(patt))
            mask = patt.split(":")[1]
        else:
            reg = re.compile(patt)
            LOGGER.debug("creating uncertainty information for RegEx {:s}".format(patt))
            mask = twiss_full.index.str.contains(reg)  # pandas chose to call it 'contains'.

        twiss_full.loc[mask, "dK1"] = (errdefs.loc[indx, "dK1"] * twiss_full.loc[mask, "K1L"]) ** 2
        twiss_full.loc[mask, "dX"] = errdefs.loc[indx, "dX"]*2
        if errdefs.loc[indx, "MAINFIELD"] == "BPM":
            twiss_full.loc[mask, "BPMdS"] = errdefs.loc[indx, "dS"]**2
        else:
            twiss_full.loc[mask, "KdS"] = (errdefs.loc[indx, "dS"] * twiss_full.loc[mask, "K1L"]) ** 2
        twiss_full.loc[mask, "UNC"] = True

    # in case of quadrupole longitudinal misalignments, the element (DRIFT) in front of the
    # misaligned quadrupole will be used for the thin lens approximation of the misalignment
    twiss_full["mKdS"] = np.roll(twiss_full.loc[:]["KdS"], 1)
    twiss_full.loc[:, "UNC"] = np.logical_or(abs(np.roll(twiss_full.loc[:, "dK1"], -1)) > 1.0e-12,
                                             twiss_full.loc[:, "UNC"])

    LOGGER.debug("DONE creating uncertainty information")
    return twiss_full[twiss_full["UNC"]]


def _add_header(df, header_dict, error_method, range_of_bpms):
    '''
    Adds common header elements to the headers of df
    df is an out parameter
    '''
    for key in header_dict.keys():
        df.headers[key] = header_dict[key]
    df.headers['BetaAlgorithmVersion'] = __version__
    df.headers['RCond'] = RCOND
    df.headers['RangeOfBPMs'] = "Adjacent" if error_method == METH_3BPM else range_of_bpms


# --------------------------------------------------------------------------------------------------
# -- Helper / Debug Functions
# --------------------------------------------------------------------------------------------------

def tilt_slice_matrix(matrix, slice_shift, slice_width, tune=0):
    """Tilts and slices the ``matrix``

    Tilting means shifting each column upwards one step more than the previous columnns, i.e.

    a a a a a       a b c d
    b b b b b       b c d e
    c c c c c  -->  c d e f
    ...             ...
    y y y y y       y z a b
    z z z z z       z a b c

    """
    invrange = matrix.shape[0] - 1 - np.arange(matrix.shape[0])
    matrix[matrix.shape[0] - slice_shift:, :slice_shift] += tune
    matrix[:slice_shift, matrix.shape[1] - slice_shift:] -= tune
    return np.roll(matrix[np.arange(matrix.shape[0]), circulant(invrange)[invrange]],
                   slice_shift, axis=0)[:slice_width]
