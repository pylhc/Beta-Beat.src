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
import shutil
import time

import numpy as np
import pandas as pd

import madx_wrapper
from model import manager
from segment_by_segment.segment_by_segment import GetLlmMeasurement
from utils import logging_tools
from utils import tfs_pandas as tfs, iotools
from utils.dict_tools import DotDict
from utils.entrypoint import entrypoint, EntryPointParameters
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
        name="meas_dir_path",
        required=True,
    )
    params.add_parameter(
        flags="--model_dir",
        help="Path to the model to use.",
        name="model_or_twiss_path",
        required=True,
    )
    params.add_parameter(
        flags="--fullresponse",
        help="Path to the fullresponse binary file.",
        name="fullresponse_path",
        required=True,
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
        name="weights_on_quantities",
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
        fullresponse_path: Path to the fullresponse binary file.
                           **Flags**: --fullresponse
        meas_dir_path: Path to the directory containing the measurement files.
                       **Flags**: --meas_dir
        model_or_twiss_path: Path to the model to use.
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
        weights_on_quantities (float): Weight to apply to each measured quantity.
                                       Input in order of optics_params.
                                       **Flags**: --weights
    """

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

        # convert numbers to dictionaries
        w_dict = dict(zip(opt.optics_params, opt.weights_on_quantities))
        m_dict = dict(zip(opt.optics_params, opt.modelcut))
        e_dict = dict(zip(opt.optics_params, opt.errorcut))

        # read data from files
        vars_list = _get_varlist(accel_cls, opt.variable_categories, opt.virt_flag)
        resp_dict = _load_fullresponse(opt.fullresponse_path, vars_list)

        optics_params, meas_dict = _get_measurment_data(
            opt.optics_params,
            opt.meas_dir_path, opt.beta_file_name,
            w_dict,
        )

        nominal_model = _load_model(opt.model_path, optics_params)

        # apply filters to data
        meas_dict = _filter_measurement(
            optics_params, meas_dict, nominal_model,
            opt.use_errorbars, w_dict, e_dict, m_dict
        )
        meas_dict = _append_model_to_measurement(nominal_model, meas_dict, optics_params)
        resp_dict = _filter_response_index(resp_dict, meas_dict, optics_params)

        _dump(os.path.join(opt.output_path, "measurement_dict.bin"), meas_dict)

        delta = tfs.TfsDataFrame(0, index=vars_list, columns=["DELTA"])

        # ######### Actual optimization algorithm ######### #

        # as long as response is not recalculated leave outside loop:
        resp_matrix = _join_responses(resp_dict, optics_params, vars_list)

        for iteration in range(opt.max_iter + 1):
            LOG.debug("Correction Iteration {:d} of {:d}.".format(iteration, opt.max_iter))

            if iteration > 0:
                LOG.debug("Updating model via MADX")
                corr_model_path = _create_corrected_model(iteration, accel_cls, nominal_model,
                                                          opt.optics_file, opt.template_file_path,
                                                          opt.output_path, opt.debug)

                corr_model = _load_model(corr_model_path, optics_params)
                meas_dict = _append_model_to_measurement(corr_model, meas_dict, optics_params)

            if opt.debug:
                _print_rms(meas_dict, optics_params)

            # new deltas
            delta += _calculate_delta(
                resp_matrix, meas_dict, optics_params, vars_list, opt.method, meth_opt)

            writeparams(opt.output_path, delta)
            LOG.debug("Cumulative delta: {:.5e}".format(
                np.sum(np.abs(delta.loc[:, "DELTA"].values))))

        write_knob(opt.output_path, delta)


# Main function helpers #######################################################


def _check_opt(opt):
    """ Check on options and put in missing values """
    # get unset paths from other paths
    if os.path.isdir(opt.model_or_twiss_path):
        opt.model_dir = opt.model_or_twiss_path
        opt.model_path = os.path.join(opt.model_dir, "twiss.dat")
    else:
        opt.model_dir = os.path.dirname(opt.model_or_twiss_path)
        opt.model_path = opt.model_or_twiss_path

    if opt.optics_file is None:
        opt.optics_file = os.path.join(opt.model_dir, "modifiers.madx")

    if opt.output_path is None:
        opt.output_path = opt.meas_dir_path

    iotools.create_dirs(opt.output_path)

    bb_root = iotools.get_absolute_path_to_betabeat_root()
    opt.template_file_path = os.path.join(bb_root, "correction", "job.twiss_python.madx")

    # check cuts and weights:
    def_dict = _get_default_values()
    if opt.modelcut is None:
        opt.modelcut = [def_dict["modelcut"][p] for p in opt.optics_params]

    if opt.errorcut is None:
        opt.errorcut = [def_dict["errorcut"][p] for p in opt.optics_params]

    if opt.weights_on_quantities is None:
        opt.weights_on_quantities = [def_dict["weights"][p] for p in opt.optics_params]
    return opt


def _get_method_opt(opt):
    """ Slightly unnecessary function to separate method-options
    for easier debugging and readability """
    meth_opt = DotDict(
        svd_cut=opt.svd_cut,
    )
    return meth_opt


def _print_rms(meas, keys):
    """ Prints current RMS status """
    for key in keys:
        LOG.debug("{:s} RMS: {:.5e}".format(key, _rms(meas[key].loc[:, 'DIFF'].values)))
        LOG.debug("{:s} Weighted RMS: {:.5e}".format(
            key, _rms(meas[key].loc[:, 'DIFF'].values * meas[key].loc[:, 'WEIGHT'].values)))


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


def _get_measurment_data(keys, meas_dir_path, beta_file_name, w_dict):
    """ Retruns a dictionary full of get_llm data """
    measurement = {}
    filtered_keys = [k for k in keys if w_dict[k] != 0]

    getllm_data = GetLlmMeasurement(meas_dir_path)
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
                'NAME': pd.Categorical(['Q1', 'Q2']),
                # Just fractional tunes:
                'VALUE': np.remainder([getllm_data.phase_x['Q1'],
                                       getllm_data.phase_y['Q2']], [1, 1]),
                # TODO measured errors not in the file
                'ERROR': np.array([0.001, 0.001])
            })
        else:
            # a beta key
            if beta_file_name == "getbeta":
                if key in ("BBX", "BETX"):
                    measurement[key] = getllm_data.beta_x
                elif key == ("BBY", "BETY"):
                    measurement[key] = getllm_data.beta_y

            elif beta_file_name == "getampbeta":
                if key in ("BBX", "BETX"):
                    measurement[key] = getllm_data.amp_beta_x
                elif key == ("BBY", "BETY"):
                    measurement[key] = getllm_data.amp_beta_y

            elif beta_file_name == "getkmodbeta":
                if key in ("BBX", "BETX"):
                    measurement[key] = getllm_data.kmod_beta_x
                elif key == ("BBY", "BETY"):
                    measurement[key] = getllm_data.kmod_beta_y
    return filtered_keys, measurement


def _get_varlist(accel_cls, variables, virt_flag):  # TODO: Virtual?
    return np.array(accel_cls.get_variables(classes=variables))


def _load_model(model_path, keys):
    model = tfs.read_tfs(model_path).set_index("NAME", drop=False)
    if any([key for key in keys if key.startswith("F1")]):
        model = _add_coupling_to_model(model)
    return model


def _add_coupling_to_model(model):
    tw_opt = TwissOptics(model)
    couple = tw_opt.get_coupling(method="cmatrix")
    model.loc[:, "F1001R"] = couple["F1001"].apply(np.real).astype(np.float64)
    model.loc[:, "F1001I"] = couple["F1001"].apply(np.imag).astype(np.float64)
    model.loc[:, "F1010R"] = couple["F1010"].apply(np.real).astype(np.float64)
    model.loc[:, "F1010I"] = couple["F1010"].apply(np.imag).astype(np.float64)
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

    # name and value
    new = meas.loc[common_bpms, ['NAME', key]]
    new.columns = ['NAME', 'VALUE']

    # errors
    if ("ERR" + key) in meas.columns.values:  # usually beta
        if ('STD' + key) in meas.columns.values:  # Old files or k-mod
            new['ERROR'] = np.sqrt(np.square(new.loc[common_bpms, 'ERR' + key].values) +
                                   np.square(new.loc[common_bpms, 'STD' + key].values))
        else:
            new['ERROR'] = meas.loc[common_bpms, 'ERR' + key]

    else:
        key2num = {'1001': '1', '1010': '2'}
        if key[1:-1] in key2num:  # coupling
            new['ERROR'] = meas.loc[common_bpms, 'FWSTD' + key2num[key[1:-1]]]
        else:
            new['ERROR'] = meas.loc[common_bpms, 'STD'+key]

    # weights
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = new.loc[:, 'WEIGHT'].values / new.loc[:, 'ERROR'].values

    # filter cuts
    error_filter = new.loc[:, 'ERROR'].values < errorcut
    try:
        model_filter = np.abs(new.loc[:, 'VALUE'].values -
                              meas.loc[common_bpms, key + 'MDL'].values) < modelcut
    except KeyError:
        # Why is there no standard for where "MDL" is attached to the name???
        model_filter = np.abs(new.loc[:, 'VALUE'].values -
                              meas.loc[common_bpms, 'MDL' + key].values) < modelcut

    good_bpms = error_filter & model_filter
    LOG.debug("Number of BPMs with {:s}: {:d}".format(key, np.sum(good_bpms)))
    return new.loc[good_bpms, :]


def _get_filtered_phases(key, meas, model, erwg, weight, modelcut, errorcut):
    common_bpms = meas.index.intersection(model.index)

    col_val = "PHASE" + key[-1]
    col_err = "STDPH" + key[-1]
    col_mdl = "PH" + key[-1] + "MDL"

    # name and value
    new = meas.loc[common_bpms, ['NAME', col_val]]
    new.columns = ['NAME', 'VALUE']

    # errors
    new['ERROR'] = meas.loc[common_bpms, col_err]

    # weights
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = new.loc[:, 'WEIGHT'].values / new.loc[:, 'ERROR'].values

    # filter cuts
    error_filter = new.loc[:, 'ERROR'].values < errorcut
    model_filter = np.abs(new.loc[:, 'VALUE'].values -
                          meas.loc[common_bpms, col_mdl].values) < modelcut

    new['NAME2'] = meas.loc[common_bpms, 'NAME2']
    second_bpm_in = np.in1d(new.loc[:, 'NAME2'].values,
                            new.loc[:, 'NAME'].values)
    good_bpms = error_filter & model_filter & second_bpm_in
    good_bpms[-1] = False

    LOG.debug("Number of BPMs with {:s}: {:d}".format(key, np.sum(good_bpms)))
    return new.loc[good_bpms, :]


def _get_filtered_betabeat(key, meas, model, erwg, weight, modelcut, errorcut):
    # Beta-beating and its error RELATIVE as shown in GUI
    # TODO: Rewrite with intersections instead of merge (like generic)
    new = pd.merge(meas, model, how='inner', on='NAME', suffixes=('', 'm'))
    if 'BETYMDL' in new.columns.values:
        plane = 'Y'
    else:
        plane = 'X'
    new['BETAMOD'] = new.loc[:, 'BET' + plane + 'MDL'].values
    new['VALUE'] = new.loc[:, 'BET' + plane].values
    if ('STDBET' + plane) in new.columns.values:  # Old files or k-mod
        new['ERROR'] = np.sqrt(np.square(new.loc[:, 'ERRBET' + plane].values) +
                               np.square(new.loc[:, 'STDBET' + plane].values))
    else:
        new['ERROR'] = new.loc[:, 'ERRBET' + plane].values
    model_close = (
        np.abs(new.loc[:, 'VALUE'].values - new.loc[:, 'BETAMOD'].values) /
        new.loc[:, 'BETAMOD'].values
    ) < modelcut
    error_low = (
        new.loc[:, 'ERROR'].values /
        new.loc[:, 'BETAMOD'].values
    ) < errorcut
    good_bpms = model_close & error_low
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = (
            new.loc[:, 'WEIGHT'].values * new.loc[:, 'BETAMOD'].values /
            new.loc[:, 'ERROR'].values
        )
    LOG.debug("Number of BPMs with beta in plane" +
              plane + ": " + str(np.sum(good_bpms)))
    return new.loc[good_bpms, ['NAME', 'VALUE', 'ERROR', 'WEIGHT']]


def _get_tunes(key, meas, model, erwg, weight, modelcut=0.1, errorcut=0.027):
    meas['WEIGHT'] = weight
    if erwg:
        meas['WEIGHT'] = (meas.loc[:, 'WEIGHT'].values /
                          meas.loc[:, 'ERROR'].values)
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
    meas['MODEL'] = model.loc[meas.loc[:, 'NAME'].values, key].values
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
    return meas


def _get_model_phases(model, meas, key):
    meas['MODEL'] = (model.loc[meas.loc[:, 'NAME2'].values, key].values -
                     model.loc[meas.loc[:, 'NAME'].values, key].values)
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
    return meas


def _get_model_betabeat(model, meas, key):
    col = "BETX" if key == "BBX" else "BETY"
    meas['MODEL'] = model.loc[meas.loc[:, 'NAME'].values, col].values
    meas['DIFF'] = (
        (meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values) /
        meas.loc[:, 'MODEL'].values
    )
    return meas


def _get_model_norm_disp(model, meas, key):
    col = key[1:]
    beta = "BET" + key[-1]
    meas['MODEL'] = (
        model.loc[meas.loc[:, 'NAME'].values, col].values /
        np.sqrt(model.loc[meas.loc[:, 'NAME'].values, beta].values)
    )
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
    return meas


def _get_model_tunes(model, meas, key):
    # We want just fractional tunes
    meas['MODEL'] = np.remainder([model['Q1'], model['Q2']], [1, 1])
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
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

    LOG.debug("Delta calculation: ")
    LOG.debug("RMS weightened, Model-Measure: {:.5e}".format(_rms(diff_weighted)))
    LOG.debug("RMS weightened, (Model-Measure) - R * delta: {:.5e}".format(
        _rms(diff_weighted - np.dot(resp_weighted, delta))))
    delta = tfs.TfsDataFrame(delta, index=vars_list, columns=["DELTA"])
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


def _create_corrected_model(iteration, accel_cls, nominal_model, optics_file,
                            template_file_path, output_path, debug):
    madx_script = _create_madx_script(accel_cls, nominal_model, optics_file,
                                      template_file_path, output_path)
    _call_madx(madx_script, debug)
    corr_model_path = os.path.join(output_path, "twiss_" + str(iteration) + ".dat")
    shutil.copy2(os.path.join(output_path, "twiss_corr.dat"), corr_model_path)
    return corr_model_path


def _create_madx_script(accel_cls, nominal_model, optics_file,
                        template_path, output_path):
    with open(template_path, "r") as template:
        madx_template = template.read()
        qx, qy = nominal_model.Q1, nominal_model.Q2
        beam = int(nominal_model.SEQUENCE[-1])
        replace_dict = {
            "LIB": accel_cls.MACROS_NAME,
            "MAIN_SEQ": accel_cls.load_main_seq_madx(),
            "OPTICS_PATH": optics_file,
            "CROSSING_ON": 0,  # TODO: Crossing
            "NUM_BEAM": beam,
            "PATH": output_path,
            "DPP": 0.0,
            "QMX": qx,
            "QMY": qy,
            "COR": "changeparameters.madx",
        }
        madx_script = madx_template % replace_dict
    return madx_script


def _call_madx(madx_script, debug):
    """ Call MADX and log into file """
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


def write_knob(path, delta):
    a = datetime.datetime.fromtimestamp(time.time())
    changeparameters_path = os.path.join(path, 'changeparameters.knob')
    delta_out = - delta.loc[:, ["DELTA"]]
    delta_out.headers["PATH"] = path
    delta_out.headers["DATE"] = str(a.ctime())
    tfs.write_tfs(changeparameters_path, delta_out, save_index="NAME")


def writeparams(path, delta):
    for correct in (True, False):
        filename = "changeparameters_correct.madx" if correct else "changeparameters.madx"
        with open(os.path.join(path, filename), "w") as madx_script:
            for var in delta.index.values:
                value = -delta.loc[var, "DELTA"] if correct else delta.loc[var, "DELTA"]
                madx_script.write(var + " = " + var
                                  + (" + " if value > 0 else " ")
                                  + str(value) + ";\n"
                                  )


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


# Main invocation ############################################################

if __name__ == "__main__":
    global_correction()


