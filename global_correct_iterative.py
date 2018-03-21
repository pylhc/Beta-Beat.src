"""
Iterative Correction Scheme.

The response matrices :math:`R_{O}` for the observables :math:`O` (e.g. BBX, MUX, ...)
are loaded from a file and then the equation

.. math:: R_{O} \cdot \delta var = O_{meas} - O_{model}
    :label: eq1

is being solved for :math:`\delta var` via a chosen method (at the moment only numpys pinv,
which creates a pseudo-inverse via svd is used).

The response matrices are hereby merged into one matrix for all observables to solve vor all
:math:`\delta var` at the same time.

To normalize the observables to another ``weigths`` (W) can be applied.

Furthermore, an ``errorcut``, specifying the maximum errorbar for a BPM to be used, and
``modelcut``, specifying the maximum distance between measurement and model for a BPM to be used,
can be defined. Data from BPMs outside of those cut-values will be discarded.
These cuts are defined for each observable separately.

After each iteration the model variables are changed by :math:`-\delta var` and the
observables are recalculated by Mad-X.
:eq:`eq1` is then solved again.


:author: Lukas Malina, Joschua Dilly


Possible problems and notes:
 * do we need error cut, when we use error-based weights? probably not anymore
 * error-based weights default? likely - but be carefull with low tune errors vs
svd cut in pseudoinverse
 * manual creation of pd.DataFrame varslist, deltas? maybe
tunes in tfs_pandas single value or a column?
 * Minimal strength removed
 * Check the output files and when they are written
 * There should be some summation/renaming for iterations
 * For two beam correction
 * The two beams can be treated separately until the calcultation of correction
 * Values missing in the response (i.e. correctors of the other beam) shall be
treated as zeros
 * Missing a part that treats the output from LSA

"""
import cPickle
import datetime
import os
import pickle
import time

import numpy as np
import pandas as pd

import madx_wrapper
from correction.fullresponse.response_twiss import TwissResponse
from twiss_optics.sequence_parser import EXT as VARMAP_EXT
from model import manager
from segment_by_segment.segment_by_segment import GetLlmMeasurement
from utils import logging_tools
from utils import tfs_pandas as tfs, iotools
from utils.dict_tools import DotDict
from utils.entrypoint import entrypoint, EntryPointParameters
from utils.logging_tools import log_pandas_settings_with_copy
from twiss_optics.optics_class import TwissOptics

LOG = logging_tools.get_logger(__name__)

DEV_NULL = os.devnull


# Configuration ##################################################################

DEFAULT_ARGS = {
    "optics_file": None,
    "output_path": None,
    "svd_cut": 0.01,
    "optics_params": ['MUX', 'MUY', 'BBX', 'BBY', 'NDX', 'Q'],
    "variables": ["MQM", "MQT", "MQTL", "MQY"],
    "beta_file_name": "getbeta",
    "method": "pinv",
    "max_iter": 3,
    "eps": None,
}


# Define functions here, to new optics params
def _get_default_values():
    return {
        'modelcut': {
            'MUX': 0.05, 'MUY': 0.05,
            'BBX': 0.2, 'BBY': 0.2,
            'BETX': 0.2, 'BETY': 0.2,
            'DX': 0.2, 'DY': 0.2,
            'NDX': 0.2, 'Q': 0.1,
            'F1001R': 0.2, 'F1001I': 0.2,
            'F1010R': 0.2, 'F1010I': 0.2,
        },
        'errorcut': {
            'MUX': 0.035, 'MUY': 0.035,
            'BBX': 0.02, 'BBY': 0.02,
            'BETX': 0.02, 'BETY': 0.02,
            'DX': 0.02, 'DY': 0.02,
            'NDX': 0.02, 'Q': 0.027,
            'F1001R': 0.02, 'F1001I': 0.02,
            'F1010R': 0.02, 'F1010I': 0.02,
        },
        'weights': {
            'MUX': 1, 'MUY': 1,
            'BBX': 0, 'BBY': 0,
            'BETX': 0, 'BETY': 0,
            'DX': 0, 'DY': 0,
            'NDX': 0, 'Q': 10,
            'F1001R': 0, 'F1001I': 0,
            'F1010R': 0, 'F1010I': 0,
        },
    }


def _get_measurement_filters():
    return {
        'MUX': _get_filtered_phases, 'MUY': _get_filtered_phases,
        'BBX': _get_filtered_betabeat, 'BBY': _get_filtered_betabeat,
        'BETX': _get_filtered_generic, 'BETY': _get_filtered_generic,
        'DX': _get_filtered_generic, 'DY': _get_filtered_generic,
        'NDX': _get_filtered_generic, 'Q': _get_tunes,
        'F1001R': _get_filtered_generic, 'F1001I': _get_filtered_generic,
        'F1010R': _get_filtered_generic, 'F1010I': _get_filtered_generic,
        }


def _get_response_filters():
    return {
        'MUX': _get_phase_response, 'MUY': _get_phase_response,
        'BBX': _get_generic_response, 'BBY': _get_generic_response,
        'BETX': _get_generic_response, 'BETY': _get_generic_response,
        'DX': _get_generic_response, 'DY': _get_generic_response,
        'NDX': _get_generic_response, 'Q': _get_tune_response,
        'F1001R': _get_generic_response, 'F1001I': _get_generic_response,
        'F1010R': _get_generic_response, 'F1010I': _get_generic_response,
        }


def _get_model_appenders():
    return {
        'MUX': _get_model_phases, 'MUY': _get_model_phases,
        'BBX': _get_model_betabeat, 'BBY': _get_model_betabeat,
        'BETX': _get_model_generic, 'BETY': _get_model_generic,
        'DX': _get_model_generic, 'DY': _get_model_generic,
        'NDX': _get_model_norm_disp, 'Q': _get_model_tunes,
        'F1001R': _get_model_generic, 'F1001I': _get_model_generic,
        'F1010R': _get_model_generic, 'F1010I': _get_model_generic,
        }


def _get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--meas_dir",
        help="Path to the directory containing the measurement files.",
        name="meas_dir",
        required=True,
    )
    params.add_parameter(
        flags="--model_dir",
        help="Path to the model to use.",
        name="model_dir",
        required=True,
    )
    params.add_parameter(
        flags="--fullresponse",
        help=("Path to the fullresponse binary file."
             " If not given, calculates the response analytically."),
        name="fullresponse_path",
    )
    params.add_parameter(
        flags="--update_response",
        help="If True, it will update the (analytical) response per iteration.",
        name="update_response",
        action="store_true",
    )
    params.add_parameter(
        flags="--optics_params",
        help="List of parameters to correct upon (e.g. BBX BBY)",
        name="optics_params",
        type=str,
        nargs="+",
        default=DEFAULT_ARGS["optics_params"],
    )
    params.add_parameter(
        flags="--optics_file",
        help=("Path to the optics file to use, usually modifiers.madx. If "
              "not present will default to model_path/modifiers.madx"),
        name="optics_file",
    )
    params.add_parameter(
        flags="--output_dir",
        help=("Path to the directory where to write the output files, will "
              "default to the --meas input path."),
        name="output_path",
        default=DEFAULT_ARGS["output_path"],
    )
    params.add_parameter(
        flags="--svd_cut",
        help=("Cutoff for small singular values of the pseudo inverse. (Method: 'pinv')"
              "Singular values smaller than rcond*largest_singular_value are set to zero"),
        name="svd_cut",
        type=float,
        default=DEFAULT_ARGS["svd_cut"],
    )
    params.add_parameter(
        flags="--model_cut",
        help=("Reject BPMs whose deviation to the model is higher than the "
              "correspoding input. Input in order of optics_params."),
        name="modelcut",
        nargs="+",
        type=float,
    )
    params.add_parameter(
        flags="--error_cut",
        help=("Reject BPMs whose error bar is higher than the "
              "corresponding input. Input in order of optics_params."),
        name="errorcut",
        nargs="+",
        type=float,
    )
    params.add_parameter(
        flags="--weights",
        help=("Weight to apply to each measured quantity. "
              "Input in order of optics_params."),
        name="weights",
        nargs="+",
        type=float,
    )
    params.add_parameter(
        flags="--use_errorbars",
        help=("If True, it will take into account the measured errorbars "
              "in the correction."),
        name="use_errorbars",
        action="store_true",
    )
    params.add_parameter(
        flags="--variables",
        help="List of names of the variables classes to use.",
        name="variable_categories",
        nargs="+",
        default=DEFAULT_ARGS["variables"],
    )
    params.add_parameter(
        flags="--beta_file_name",
        help="Prefix of the beta file to use. E.g.: getkmodbeta",
        name="beta_file_name",
        default=DEFAULT_ARGS["beta_file_name"],
    )
    params.add_parameter(
        flags="--virt_flag",
        help="If true, it will use virtual correctors.",
        name="virt_flag",
        action="store_true",
    )
    params.add_parameter(
        flags="--method",
        help="Optimization method to use. (Not implemented yet)",
        name="method",
        type=str,
        default=DEFAULT_ARGS["method"],
        choices=["pinv"]
    )
    params.add_parameter(
        flags="--max_iter",
        help=("Maximum number of correction re-iterations to perform."
              "A value of `0` means the correction is calculated once (like in the old days)."),
        name="max_iter",
        type=int,
        default=DEFAULT_ARGS["max_iter"],
    )
    params.add_parameter(
        flags="--eps",
        help=("Convergence criterion." 
              "If <|delta(PARAM * WEIGHT)|> < eps, stop iteration.(Not implemented yet)"),
        name="eps",
        type=float,
        default=DEFAULT_ARGS["eps"],
    )
    params.add_parameter(
        flags="--debug",
        help="Print debug information.",
        name="debug",
        action="store_true",
    )
    return params


# Entry Point ##################################################################


@entrypoint(_get_params())
def global_correction(opt, accel_opt):
    """ Do the global correction. Iteratively.

    Keyword Args:
        Required
        meas_dir: Path to the directory containing the measurement files.
                       **Flags**: --meas_dir
        model_dir: Path to the dir containing the model (twiss.dat or twiss_elements.dat) to use.
                             **Flags**: --model_dir
        Optional
        beta_file_name: Prefix of the beta file to use. E.g.: getkmodbeta
                        **Flags**: --beta_file_name
                        **Default**: ``getbeta``
        debug: Print debug information.
               **Flags**: --debug
               **Action**: ``store_true``
        eps (float): (Not implemented yet) Convergence criterion.
                     If :math:`<|\Delta(PARAM \cdot WEIGHT)|> < \epsilon`, stop iteration.
                     **Flags**: --eps
                     **Default**: ``None``
        errorcut (float): Reject BPMs whose error bar is higher than the corresponding input.
                          Input in order of optics_params.
                          **Flags**: --error_cut
        fullresponse_path: Path to the fullresponse binary file.
                           If not given, calculates the response analytically.
                           **Flags**: --fullresponse
        max_iter (int): Maximum number of correction re-iterations to perform.
                        A value of `0` means the correction is calculated once
                        (like in the old days).
                        **Flags**: --max_iter
                        **Default**: ``3``
        method (str): Optimization method to use. (Not implemented yet)
                      **Flags**: --method
                      **Choices**: ['pinv']
                      **Default**: ``pinv``
        modelcut (float): Reject BPMs whose deviation to the model is higher than the
                          correspoding input. Input in order of optics_params.
                          **Flags**: --model_cut
        optics_file: Path to the optics file to use, usually modifiers.madx.
                     If not present will default to model_path/modifiers.madx
                     **Flags**: --optics_file
        optics_params (str): List of parameters to correct upon (e.g. BBX BBY)
                             **Flags**: --optics_params
                             **Default**: ``['MUX', 'MUY', 'BBX', 'BBY', 'NDX', 'Q']``
        output_path: Path to the directory where to write the output files,
                     will default to the --meas input path.
                     **Flags**: --output_dir
                     **Default**: ``None``
        svd_cut (float): Cutoff for small singular values of the pseudo inverse. (Method: 'pinv')
                         Singular values smaller than
                         :math:`rcond \cdot largest_singular_value` are set to zero
                         **Flags**: --svd_cut
                         **Default**: ``0.01``
        use_errorbars: If True, it will take into account the measured errorbars in the correction.
                       **Flags**: --use_errorbars
                       **Action**: ``store_true``
        variable_categories: List of names of the variables classes to use.
                             **Flags**: --variables
                             **Default**: ``['MQM', 'MQT', 'MQTL', 'MQY']``
        virt_flag: If true, it will use virtual correctors.
                   **Flags**: --virt_flag
                   **Action**: ``store_true``
        weights (float): Weights to apply to each measured quantity. Input in order of optics_params.
                         **Flags**: --weights
    """

    LOG.info("Starting Iterative Global Correction.")
    with logging_tools.DebugMode(active=opt.debug):
        not_implemented_params = [k for k in opt.optics_params
                                  if k not in _get_measurement_filters()]
        if any(not_implemented_params):
            raise NotImplementedError("Correct iterative is not equipped for parameters:"
                                      "'{:s}'".format(not_implemented_params))

        # ######### Preparations ######### #

        # check on opt
        opt = _check_opt(opt)
        meth_opt = _get_method_opt(opt)

        # get accelerator class
        accel_cls = manager.get_accel_class(accel_opt)
        accel_inst = accel_cls(model_dir=opt.model_dir)
        if opt.optics_file is not None:
            accel_inst.optics_file = opt.optics_file

        # convert numbers to dictionaries
        w_dict = dict(zip(opt.optics_params, opt.weights))
        mcut_dict = dict(zip(opt.optics_params, opt.modelcut))
        ecut_dict = dict(zip(opt.optics_params, opt.errorcut))

        # read data from files
        vars_list = _get_varlist(accel_cls, opt.variable_categories, opt.virt_flag)
        optics_params, meas_dict = _get_measurment_data(
            opt.optics_params,
            opt.meas_dir, opt.beta_file_name,
            w_dict,
        )

        if opt.fullresponse_path is not None:
            resp_dict = _load_fullresponse(opt.fullresponse_path, vars_list)
        else:
            resp_dict = _create_response(accel_inst, vars_list, optics_params)

        # the model in accel_inst is modified later, so save nominal model here to variables
        nominal_model = _maybe_add_coupling_to_model(accel_inst.get_model_tfs(), optics_params)

        # apply filters to data
        meas_dict = _filter_measurement(
            optics_params, meas_dict, nominal_model,
            opt.use_errorbars, w_dict, ecut_dict, mcut_dict
        )
        meas_dict = _append_model_to_measurement(nominal_model, meas_dict, optics_params)
        resp_dict = _filter_response_index(resp_dict, meas_dict, optics_params)
        resp_matrix = _join_responses(resp_dict, optics_params, vars_list)

        # _dump(os.path.join(opt.output_path, "measurement_dict.bin"), meas_dict)
        delta = tfs.TfsDataFrame(0, index=vars_list, columns=["DELTA"])
        change_params_path = os.path.join(opt.output_path, "changeparameters.madx")
        change_params_correct_path = os.path.join(opt.output_path, "changeparameters_correct.madx")

        # ######### Iteration Phase ######### #

        for iteration in range(opt.max_iter + 1):
            LOG.info("Correction Iteration {:d} of {:d}.".format(iteration, opt.max_iter))

            # ######### Update Model and Response ######### #
            if iteration > 0:
                LOG.debug("Updating model via MADX.")
                corr_model_path = os.path.join(opt.output_path, "twiss_" + str(iteration) + ".dat")

                _create_corrected_model(corr_model_path, change_params_path, accel_inst, opt.debug)

                corr_model_elements = tfs.read_tfs(corr_model_path, index="NAME")
                corr_model = corr_model_elements.loc[tfs.get_bpms(corr_model_elements), :]
                corr_model = _maybe_add_coupling_to_model(corr_model, optics_params)

                meas_dict = _append_model_to_measurement(corr_model, meas_dict, optics_params)
                if opt.update_response:
                    LOG.debug("Updating response.")
                    # please look away for the next two lines.
                    accel_inst._model = corr_model
                    accel_inst._elements = corr_model_elements
                    resp_dict = _create_response(accel_inst, vars_list, optics_params)
                    resp_dict = _filter_response_index(resp_dict, meas_dict, optics_params)
                    resp_matrix = _join_responses(resp_dict, optics_params, vars_list)

            # ######### Actual optimization ######### #
            delta += _calculate_delta(
                resp_matrix, meas_dict, optics_params, vars_list, opt.method, meth_opt)

            writeparams(change_params_path, delta)
            writeparams(change_params_correct_path, -delta)
            LOG.debug("Cumulative delta: {:.5e}".format(
                np.sum(np.abs(delta.loc[:, "DELTA"].values))))

        write_knob(os.path.join(opt.output_path, "changeparameters.knob"), delta)
    LOG.info("Finished Iterative Global Correction.")

# Main function helpers #######################################################


def _check_opt(opt):
    """ Check on options and put in missing values """
    # get unset paths from other paths
    if opt.output_path is None:
        opt.output_path = opt.meas_dir

    iotools.create_dirs(opt.output_path)

    bb_root = iotools.get_absolute_path_to_betabeat_root()
    opt.template_file_path = os.path.join(bb_root, "correction", "job.twiss_python.madx")

    # check cuts and weights:
    def_dict = _get_default_values()
    if opt.modelcut is None:
        opt.modelcut = [def_dict["modelcut"][p] for p in opt.optics_params]
    elif len(opt.optics_params) != len(opt.modelcut):
        raise ValueError("The length of modelcut is not the same as of the optical parameters!")

    if opt.errorcut is None:
        opt.errorcut = [def_dict["errorcut"][p] for p in opt.optics_params]
    elif len(opt.optics_params) != len(opt.errorcut):
        raise ValueError("The length of errorcut is not the same as of the optical parameters!")

    if opt.weights is None:
        opt.weights = [def_dict["weights"][p] for p in opt.optics_params]
    elif len(opt.optics_params) != len(opt.weights):
        raise ValueError("The length of the weights is not the same as of the optical parameters!")

    return opt


def _get_method_opt(opt):
    """ Slightly unnecessary function to separate method-options
    for easier debugging and readability """
    meth_opt = DotDict(
        svd_cut=opt.svd_cut,
    )
    return meth_opt


def _print_rms(meas, diff_w, r_delta_w):
    """ Prints current RMS status """
    f_str = "{:>20s} : {:.5e}"

    LOG.debug("RMS Model - Measure (before correction, w/o weigths):")
    for key in meas:
        LOG.debug(f_str.format(key, _rms(meas[key].loc[:, 'DIFF'].values)))

    LOG.info("RMS Model - Measure (before correction, w/ weigths):")
    for key in meas:
        LOG.info(f_str.format(
            key, _rms(meas[key].loc[:, 'DIFF'].values * meas[key].loc[:, 'WEIGHT'].values)))
    LOG.info(f_str.format("Model - Measure", _rms(diff_w)))
    LOG.debug(f_str.format("R * delta", _rms(r_delta_w)))
    LOG.debug("(Model - Measure) - (R * delta)   ")
    LOG.debug(f_str.format("", _rms(diff_w - r_delta_w)))


def _load_fullresponse(full_response_path, variables):
    """
    Full response is dictionary of optics-parameter gradients upon
    a change of a single quadrupole strength
    """
    LOG.debug("Starting loading Full Response optics")
    with open(full_response_path, "r") as full_response_file:
        full_response_data = pickle.load(full_response_file)

    for param in full_response_data:
        df = full_response_data[param]
        # fill with zeros
        not_found_vars = pd.Index(variables).difference(df.columns)
        if len(not_found_vars) > 0:
            LOG.debug(("Variables not in fullresponse {:s} " +
                       "(To be filled with zeros): {:s}").format(param, not_found_vars.values))
            df = df.assign(**dict.fromkeys(not_found_vars, 0.0))  # one-liner to add zero-columns
        # order variables
        full_response_data[param] = df.loc[:, variables]

    LOG.debug("Loading ended")
    return full_response_data


def _create_response(accel_inst, vars_list, optics_params):
    """ Create response via TwissResponse """
    varmap_path = check_varmap_file(accel_inst)
    tr = TwissResponse(varmap_path, accel_inst.get_elements_tfs(), vars_list)
    return tr.get_response_for(optics_params)


def _get_measurment_data(keys, meas_dir, beta_file_name, w_dict):
    """ Retruns a dictionary full of get_llm data """
    measurement = {}
    filtered_keys = [k for k in keys if w_dict[k] != 0]

    getllm_data = GetLlmMeasurement(meas_dir)
    for key in filtered_keys:
        if key == "MUX":
            measurement['MUX'] = getllm_data.phase_x
        elif key == 'MUY':
            measurement['MUY'] = getllm_data.phase_y
        elif key == "DX":
            measurement['DX'] = getllm_data.disp_x
        elif key == "DX":
            measurement['DY'] = getllm_data.disp_y
        elif key == "NDX":
            measurement['NDX'] = getllm_data.norm_disp
        elif key in ('F1001R', 'F1001I', 'F1010R', 'F1010I'):
            measurement[key] = getllm_data.coupling
        elif key == "Q":
            measurement["Q"] = pd.DataFrame({
                # Just fractional tunes:
                'VALUE': np.remainder([getllm_data.phase_x['Q1'],
                                       getllm_data.phase_y['Q2']], [1, 1]),
                # TODO measured errors not in the file
                'ERROR': np.array([0.001, 0.001])
            }, index=['Q1', 'Q2'])
        else:
            # a beta key
            if beta_file_name == "getbeta":
                if key in ("BBX", "BETX"):
                    measurement[key] = getllm_data.beta_x
                elif key in ("BBY", "BETY"):
                    measurement[key] = getllm_data.beta_y

            elif beta_file_name == "getampbeta":
                if key in ("BBX", "BETX"):
                    measurement[key] = getllm_data.amp_beta_x
                elif key in ("BBY", "BETY"):
                    measurement[key] = getllm_data.amp_beta_y

            elif beta_file_name == "getkmodbeta":
                if key in ("BBX", "BETX"):
                    measurement[key] = getllm_data.kmod_beta_x
                elif key in ("BBY", "BETY"):
                    measurement[key] = getllm_data.kmod_beta_y
    return filtered_keys, measurement


def _get_varlist(accel_cls, variables, virt_flag):  # TODO: Virtual?
    varlist = np.array(accel_cls.get_variables(classes=variables))
    if len(varlist) == 0:
        raise ValueError("No variables found! Make sure your categories are valid!")
    return varlist


def _maybe_add_coupling_to_model(model, keys):
    if any([key for key in keys if key.startswith("F1")]):
        tw_opt = TwissOptics(model)
        couple = tw_opt.get_coupling(method="cmatrix")
        model["F1001R"] = couple["F1001"].apply(np.real).astype(np.float64)
        model["F1001I"] = couple["F1001"].apply(np.imag).astype(np.float64)
        model["F1010R"] = couple["F1010"].apply(np.real).astype(np.float64)
        model["F1010I"] = couple["F1010"].apply(np.imag).astype(np.float64)
    return model


# Parameter filtering ########################################################


def _filter_measurement(keys, meas, model, errorbar, w_dict, e_dict, m_dict):
    """ Filters measurements and renames columns to VALUE, ERROR, WEIGHT"""
    filters = _get_measurement_filters()
    new = dict.fromkeys(keys)
    for key in keys:
        new[key] = filters[key](key, meas[key], model, errorbar, w_dict[key],
                                 modelcut=m_dict[key], errorcut=e_dict[key])
    return new


def _get_filtered_generic(key, meas, model, erwg, weight, modelcut, errorcut):
    common_bpms = meas.index.intersection(model.index)
    meas = meas.loc[common_bpms, :]

    # value
    new = tfs.TfsDataFrame(index=common_bpms)
    new.loc[:, "VALUE"] = meas[key]

    # errors
    if ("ERR" + key) in meas.columns.values:  # usually beta
        if ('STD' + key) in meas.columns.values:  # Old files or k-mod
            new['ERROR'] = np.sqrt(np.square(meas['ERR' + key].values) +
                                   np.square(meas['STD' + key].values))
        else:
            new['ERROR'] = meas['ERR' + key]

    else:
        key2num = {'1001': '1', '1010': '2'}
        if key[1:-1] in key2num:  # coupling
            new.loc[:, 'ERROR'] = meas['FWSTD' + key2num[key[1:-1]]]
        else:
            new.loc[:, 'ERROR'] = meas['STD'+key]

    # weights
    new.loc[:, 'WEIGHT'] = weight
    if erwg:
        new.loc[:, 'WEIGHT'] = new.loc[:, 'WEIGHT'].values / new.loc[:, 'ERROR'].values

    # filter cuts
    error_filter = new.loc[:, 'ERROR'].values < errorcut
    try:
        model_filter = np.abs(new.loc[:, 'VALUE'].values -
                              meas[key + 'MDL'].values) < modelcut
    except KeyError:
        # Why is there no standard for where "MDL" is attached to the name???
        model_filter = np.abs(new['VALUE'].values -
                              meas['MDL' + key].values) < modelcut

    good_bpms = error_filter & model_filter
    LOG.debug("Number of BPMs with {:s}: {:d}".format(key, np.sum(good_bpms)))
    return new.loc[good_bpms, :]


def _get_filtered_phases(key, meas, model, erwg, weight, modelcut, errorcut):
    common_bpms = meas.index.intersection(model.index)
    meas = meas.loc[common_bpms, :]

    col_val = "PHASE" + key[-1]
    col_err = "STDPH" + key[-1]
    col_mdl = "PH" + key[-1] + "MDL"

    # value
    new = tfs.TfsDataFrame(index=common_bpms)
    new.loc[:, "VALUE"] = meas[col_val]

    # errors
    new.loc[:, 'ERROR'] = meas[col_err]

    # weights
    new.loc[:, 'WEIGHT'] = weight
    if erwg:
        new.loc[:, 'WEIGHT'] = new['WEIGHT'] / new['ERROR']

    # filter cuts
    error_filter = new['ERROR'] < errorcut
    model_filter = np.abs(new['VALUE'] - meas[col_mdl]) < modelcut

    new.loc[:, 'NAME2'] = meas['NAME2']
    second_bpm_in = np.in1d(new['NAME2'].values, new.index.values)
    good_bpms = error_filter & model_filter & second_bpm_in
    good_bpms[-1] = False

    LOG.debug("Number of BPMs with {:s}: {:d}".format(key, np.sum(good_bpms)))
    return new.loc[good_bpms, :]


def _get_filtered_betabeat(key, meas, model, erwg, weight, modelcut, errorcut):
    # Beta-beating and its error RELATIVE as shown in GUI
    common_bpms = meas.index.intersection(model.index)
    meas = meas.loc[common_bpms, :]

    col_val = "BET" + key[-1]
    col_std = "STDBET" + key[-1]
    col_err = 'ERRBET' + key[-1]
    col_mdl = "BET" + key[-1] + "MDL"

    # value
    new = tfs.TfsDataFrame(index=common_bpms)
    new.loc[:, 'VALUE'] = meas[col_val]

    # errors
    if col_std in new.columns.values:  # Old files or k-mod
        new.loc[:, 'ERROR'] = np.sqrt(np.square(meas[col_err]) + np.square(meas[col_std]))
    else:
        new.loc[:, 'ERROR'] = meas[col_err]

    # weights
    new.loc[:, 'WEIGHT'] = weight
    if erwg:
        new.loc[:, 'WEIGHT'] = new['WEIGHT'] * meas[col_mdl] / new['ERROR']

    # filter cuts
    model_filter = np.abs(new['VALUE'] - meas[col_mdl]) / meas[col_mdl] < modelcut
    error_filter = new['ERROR'] / meas[col_mdl] < errorcut
    good_bpms = model_filter & error_filter

    LOG.debug("Number of BPMs with {:s}: {:d}".format(key, np.sum(good_bpms)))
    return new.loc[good_bpms, :]


def _get_tunes(key, meas, model, erwg, weight, modelcut=0.1, errorcut=0.027):
    meas.loc[:, 'WEIGHT'] = weight
    if erwg:
        meas.loc[:, 'WEIGHT'] = meas['WEIGHT'] / meas['ERROR']
    LOG.debug("Number of tune measurements: " + str(len(meas.index.values)))
    return meas


# Response filtering ##########################################################


def _filter_response_index(response, measurement, keys):
    not_in_response = [k for k in keys if k not in response]
    if len(not_in_response) > 0:
        raise KeyError("The following optical parameters are not present in current"
                       "response matrix: {:s}".format(not_in_response))

    filters = _get_response_filters()
    new_resp = {}
    for key in keys:
        new_resp[key] = filters[key](response[key], measurement[key])
    return new_resp


def _get_generic_response(resp, meas):
    return resp.loc[meas.index.values, :]


def _get_phase_response(resp, meas):
    phase1 = resp.loc[meas.index.values, :]
    phase2 = resp.loc[meas.loc[:, 'NAME2'].values, :]
    return -phase1.sub(phase2.values)  # phs2-phs1 but with idx of phs1


def _get_tune_response(resp, meas):
    return resp


# Model appending #############################################################


def _append_model_to_measurement(model, measurement, keys):
    appenders = _get_model_appenders()
    meas = {}
    for key in keys:
        meas[key] = appenders[key](model, measurement[key], key)
    return meas


def _get_model_generic(model, meas, key):
    with log_pandas_settings_with_copy(LOG.debug):
        meas.loc[:, 'MODEL'] = model.loc[meas.index.values, key].values
        meas.loc[:, 'DIFF'] = meas['VALUE'] - meas['MODEL']
    return meas


def _get_model_phases(model, meas, key):
    with log_pandas_settings_with_copy(LOG.debug):
        meas.loc[:, 'MODEL'] = (model.loc[meas['NAME2'].values, key].values -
                             model.loc[meas.index.values, key].values)
        meas.loc[:, 'DIFF'] = meas['VALUE'] - meas['MODEL']
    return meas


def _get_model_betabeat(model, meas, key):
    col = "BETX" if key == "BBX" else "BETY"
    with log_pandas_settings_with_copy(LOG.debug):
        meas.loc[:, 'MODEL'] = model.loc[meas.index.values, col].values
        meas.loc[:, 'DIFF'] = (meas['VALUE'] - meas['MODEL']) / meas['MODEL']
    return meas


def _get_model_norm_disp(model, meas, key):
    col = key[1:]
    beta = "BET" + key[-1]
    with log_pandas_settings_with_copy(LOG.debug):
        meas.loc[:, 'MODEL'] = (
            model.loc[meas.index.values, col].values /
            np.sqrt(model.loc[meas.index.values, beta].values)
        )
        meas.loc[:, 'DIFF'] = meas['VALUE'] - meas['MODEL']
    return meas


def _get_model_tunes(model, meas, key):
    # We want just fractional tunes
    with log_pandas_settings_with_copy(LOG.debug):
        meas.loc[:, 'MODEL'] = np.remainder([model['Q1'], model['Q2']], [1, 1])
        meas.loc[:, 'DIFF'] = meas['VALUE'] - meas['MODEL']
    return meas


# Main Calculation ################################################################


def _calculate_delta(resp_matrix, meas_dict, keys, vars_list, method, meth_opt):
    """ Get the deltas for the variables.

    Output is Dataframe with one column 'DELTA' and vars_list index. """
    weight_vector = _join_columns('WEIGHT', meas_dict, keys)
    diff_vector = _join_columns('DIFF', meas_dict, keys)

    resp_weighted = resp_matrix.mul(weight_vector, axis="index")
    diff_weighted = diff_vector * weight_vector

    delta = _get_method_fun(method)(resp_weighted, diff_weighted, meth_opt)
    delta = tfs.TfsDataFrame(delta, index=vars_list, columns=["DELTA"])

    update = np.dot(resp_weighted, delta["DELTA"])
    _print_rms(meas_dict, diff_weighted, update)
    return delta


def _get_method_fun(method):
    funcs = {
        "pinv": _pseudo_inverse
    }
    return funcs[method]


def _pseudo_inverse(response_mat, diff_vec, opt):
    """ Calculates the pseudo-inverse of the response via svd. (numpy) """
    return np.dot(np.linalg.pinv(response_mat, opt.svd_cut), diff_vec)


# MADX related ###############################################################


def _create_corrected_model(twiss_out, change_params, accel_inst, debug):
    """ Use the calculated deltas in changeparameters.madx to create a corrected model """
    # create script from template
    with open(accel_inst.get_update_correction_tmpl(), "r") as template:
        madx_template = template.read()

    replace_dict = {
        "LIB": accel_inst.MACROS_NAME,
        "MAIN_SEQ": accel_inst.load_main_seq_madx(),
        "OPTICS_PATH": accel_inst.optics_file,
        "CROSSING_ON": 0,  # TODO: Crossing
        "NUM_BEAM": accel_inst.get_beam(),
        "PATH_TWISS": twiss_out,
        "DPP": 0.0,
        "QMX": accel_inst.nat_tune_x,
        "QMY": accel_inst.nat_tune_y,
        "CORRECTIONS": change_params,
    }
    madx_script = madx_template % replace_dict

    # run madx
    if debug:
        with logging_tools.TempFile("correct_iter_madxout.tmp", LOG.debug) as log_file:
            madx_wrapper.resolve_and_run_string(
                madx_script,
                log_file=log_file,
            )
    else:
        madx_wrapper.resolve_and_run_string(
            madx_script,
            log_file=os.devnull,
        )


def write_knob(knob_path, delta):
    a = datetime.datetime.fromtimestamp(time.time())
    delta_out = - delta.loc[:, ["DELTA"]]
    delta_out.headers["PATH"] = os.path.dirname(knob_path)
    delta_out.headers["DATE"] = str(a.ctime())
    tfs.write_tfs(knob_path, delta_out, save_index="NAME")


def writeparams(path_to_file, delta):
    with open(path_to_file, "w") as madx_script:
        for var in delta.index.values:
            value = delta.loc[var, "DELTA"]
            madx_script.write("{var:s} = {var:s} {value:+e};\n".format(var=var, value=value))


# Small Helpers ################################################################


def _rms(a):
    return np.sqrt(np.mean(np.square(a)))


def _dump(path_to_dump, content):
    with open(path_to_dump, 'wb') as dump_file:
        cPickle.Pickler(dump_file, -1).dump(content)


def _join_responses(resp, keys, varslist):
    """ Returns matrix #BPMs * #Parameters x #variables """
    return pd.concat([resp[k] for k in keys],  # dataframes
                     axis="index",  # axis to join along
                     join_axes=[pd.Index(varslist)]
                     # other axes to use (pd Index obj required)
                     )


def _join_columns(col, meas, keys):
    """ Retuns vector: N= #BPMs * #Parameters (BBX, MUX etc.) """
    return np.concatenate([meas[key].loc[:, col].values for key in keys], axis=0)


# Related Public Methods #####################################################


def check_varmap_file(accel_inst):
    """ Checks on varmap file and creates it if not in model folder.
    THIS SHOULD BE REPLACED WITH A CALL TO JAIMES DATABASE, IF IT BECOMES AVAILABLE """
    varmapfile_name = accel_inst.NAME.lower() + "b" + str(accel_inst.get_beam())
    varmap_path = os.path.join(accel_inst.model_dir, varmapfile_name + "." + VARMAP_EXT)
    if not os.path.isfile(varmap_path):
        LOG.debug("Variable mapping not found. Creating it with madx/sequence_parser.")
        varmap_path = varmap_path.replace("." + VARMAP_EXT, ".seq")
        save_sequence_jobfile = os.path.join(accel_inst.model_dir, "job.save_sequence.madx")

        if os.path.isfile(save_sequence_jobfile):
            with logging_tools.TempFile("save_sequence_madxout.tmp", LOG.debug) as log_file:
                madx_wrapper.resolve_and_run_file(save_sequence_jobfile, log_file=log_file)

        else:
            LOG.warning("'job.save_sequence.madx' not found. Using standard 'main.seq'")
            iotools.copy_item(accel_inst.get_sequence_file(), varmap_path)

    return varmap_path


# Main invocation ############################################################


if __name__ == "__main__":
    global_correction()


