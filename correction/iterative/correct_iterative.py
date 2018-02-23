"""
Iterative Correction Scheme





Possible problems and notes:
- do we need error cut, when we use error-based weights? probably not anymore
- error-based weights default? likely - but be carefull with low tune errors vs
svd cut in pseudoinverse
- manual creation of pd.DataFrame varslist, deltas? maybe
tunes in tfs_pandas single value or a column?
- Minimal strength removed
- Check the output files and when they are written
- There should be some summation/renaming for iterations
- For two beam correction
- The two beams can be treated separately until the calcultation of correction
- Values missing in the response (i.e. correctors of the other beam) shall be
treated as zeros
- Missing a part that treats the output from LSA
"""

# noinspection PyUnresolvedReferences
import __init__
import os
import shutil
import time
import datetime
import pickle
import cPickle
import numpy as np
import pandas as pd

from madx import madx_wrapper  # noqa
from utils import tfs_pandas as tfs, iotools  # noqa
from utils.contexts import timeit
from utils import logging_tools
from utils.dict_tools import DotDict
from utils.entrypoint import entrypoint, EntryPointParameters
from model import manager  # noqa
from segment_by_segment.segment_by_segment import GetLlmMeasurement

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
        flags="--meas",
        help="Path to the directory containing the measurement files.",
        name="meas_dir_path",
        required=True,
    )
    params.add_parameter(
        flags="--model",
        help="Path to the model to use.",
        name="model_twiss_path",
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
        nargs="*",
        default=DEFAULT_ARGS["optics_params"],
    )
    params.add_parameter(
        flags="--optics_file",
        help=("Path to the optics file to use, usually modifiers.madx. If "
              "not present will default to model_path/modifiers.madx"),
        name="optics_file",
    )
    params.add_parameter(
        flags="--output",
        help=("Path to the directory where to write the ouput files, will "
              "default to the --meas input path."),
        name="output_path",
        default=DEFAULT_ARGS["output_path"],
    )
    params.add_parameter(
        flags="--svd_cut",
        help=("Cutoff for small singular values of the pseudo inverse. (Method: 'pinv')"
              "Singular values smaller than `rcond`*largest_singular_value are set to zero"),
        name="svd_cut",
        type=float,
        default=DEFAULT_ARGS["svd_cut"],
    )
    params.add_parameter(
        flags="--model_cut",
        help=("Reject BPMs whose deviation to the model is higher than the "
              "correspoding input. Input in order of optics_params."),
        name="modelcut",
    )
    params.add_parameter(
        flags="--error_cut",
        help=("Reject BPMs whose error bar is higher than the "
              "correspoding input. Input in order of optics_params."),
        name="errorcut",
        nargs="*",
        type=float,
    )
    params.add_parameter(
        flags="--weights",
        help=("Weight to apply to each measured quatity. "
              "Input in order of optics_params."),
        name="weights_on_quantities",
        nargs="*",
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
        help="Comma separated names of the variables classes to use.",
        name="variable_categories",
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
        help="Maximum number of correction iterations to perform.",
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
    """ Do the global correction.

    Keyword Args:
        variable_categories: Comma separated names of the variables classes to use.
        errorcut: Reject BPMs whose error bar is higher than the correspoding input.
                  Input should be: Phase,Betabeat,NDx
        optics_file: Path to the optics file to use, usually modifiers.madx.
                     If not present will default to model_path/modifiers.madx
        model_twiss_path: Path to the model to use.
        num_reiteration: Number of correction iterations to perform.
        weights_on_quantities: Weight to apply to each measured quatity.
                               Input shoud be: PhaseX,PhaseY,BetaX,BetaY,NDx,Q
        beta_file_name: Prefix of the beta file to use. E.g.: getkmodbeta
        use_errorbars: If present, it will take into account the measured errorbars in the correction.
        modelcut: Reject BPMs whose deviation to the model is higher than the correspoding input. Input should be: Phase,Betabeat,NDx
        output_path: Path to the directory where to write the ouput files, will default to the --meas input path.
        virt_flag: If present, it will use virtual correctors.
        debug: Print debug information.
        meas_dir_path: Path to the directory containing the measurement files.
        singular_value_cut:
        fullresponse_path: Path to the fullresponse binary file.
    """

    if opt.debug:
        logging_tools.start_debug_mode()

    not_implemented_params = [k for k in opt.optics_params if k not in _get_measurement_filters()]
    if any(not_implemented_params):
        raise NotImplementedError("Correct iterative is not equipped for parameters:"
                                  "'{:s}'".format(not_implemented_params))

    with timeit(lambda t: LOG.debug("  Total time for Global Correction: {:f}s".format(t))):
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
        nominal_model = tfs.read_tfs(opt.model_twiss_path).set_index("NAME", drop=False)
        vars_list = _get_varlist(accel_cls, opt.variable_categories, opt.virt_flag)
        resp_dict = _load_fullresponse(opt.fullresponse_path, vars_list)

        optics_params, meas_dict = _get_measurment_data(
            opt.optics_params,
            opt.meas_dir_path, opt.beta_file_name,
            w_dict,
        )

        # apply filters to data
        meas_dict = _filter_measurement(
            optics_params, meas_dict, nominal_model,
            opt.use_errorbars, w_dict, e_dict, m_dict
        )
        meas_dict = _append_model_to_measurement(nominal_model, meas_dict, optics_params)
        resp_dict = _filter_response_index(resp_dict, meas_dict, optics_params)

        # as long as response is not recalculated leave outside loop:
        resp_matrix = _join_responses(resp_dict, optics_params, vars_list)

        if opt.debug:
            _print_rms(meas_dict, optics_params)

        _dump(os.path.join(opt.output_path, "measurement_dict.bin"), meas_dict)
        delta = _calculate_delta(resp_matrix, meas_dict, optics_params, vars_list, opt.method,
                                 meth_opt)
        writeparams(opt.output_path, delta)

        for _iter in range(opt.max_iter):
            LOG.debug("Running MADX, iteration {:d} of {:d}".format(_iter + 1, opt.max_iter))

            # update model
            madx_script = _create_madx_script(accel_cls, nominal_model, opt.optics_file,
                                              opt.template_file_path, opt.output_path)
            _call_madx(madx_script, opt.debug)
            new_model_path = os.path.join(opt.output_path, "twiss_" + str(_iter) + ".dat")
            shutil.copy2(os.path.join(opt.output_path, "twiss_corr.dat"), new_model_path)
            new_model = tfs.read_tfs(new_model_path, index="NAME")
            meas_dict = _append_model_to_measurement(new_model, meas_dict, optics_params)

            if opt.debug:
                _print_rms(meas_dict, optics_params)

            # get new deltas
            delta += _calculate_delta(
                resp_matrix, meas_dict, optics_params, vars_list, opt.method, meth_opt)

            writeparams(opt.output_path, delta, append=True)
            LOG.debug("Cumulative delta: {:.5e}".format(
                np.sum(np.abs(delta.loc[:, "DELTA"].values))))
        write_knob(opt.output_path, delta)


# Main function helpers #######################################################


def _check_opt(opt):
    """ Check on options and put in missing values """
    # get unset paths from other paths
    if opt.optics_file is None:
        opt.optics_file = os.path.join(os.path.dirname(opt.model_twiss_path),
                                       "modifiers.madx")
    if opt.output_path is None:
        opt.output_path = opt.meas_dir_path

    iotools.create_dirs(opt.output_path)

    opt.template_file_path = os.path.join(os.path.dirname(__file__), "job.twiss_python.madx")

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
    Full response is dictionary of MUX/Y, BBX/Y, NDX and Q gradients upon
    a change of a single quadrupole strength
    """
    LOG.debug("Starting loading Full Response optics")
    with open(full_response_path, "r") as full_response_file:
        full_response_data = pickle.load(full_response_file)

    for param in full_response_data:
        df = full_response_data[param]
        # fill with zeros
        not_found_vars = [var for var in variables if var not in df.columns.values]
        if len(not_found_vars) > 0:
            LOG.debug(("Variables not in fullresponse {:s} " +
                       "(To be filled with zeros): {:s}").format(param, not_found_vars))
            df = df.assign(**dict.fromkeys(not_found_vars, 0.0))  # one-liner to add zero-columns
        # order variables
        full_response_data[param] = df.loc[:, variables]

    LOG.debug("Loading ended")
    return full_response_data


def _get_measurment_data(keys, meas_dir_path, beta_file_name, w_dict):
    measurement = {}
    filtered_keys = [k for k in keys if w_dict[k] != 0]
    getllm_data = GetLlmMeasurement(meas_dir_path)

    if any(k in filtered_keys for k in ["MUX", "MUY", "Q"]):
        measurement['MUX'] = getllm_data.phase_x
        measurement['MUY'] = getllm_data.phase_y

        # tune
        measurement['Q'] = pd.DataFrame({
            'NAME': pd.Categorical(['Q1', 'Q2']),
            # Just fractional tunes:
            'VALUE': np.remainder([measurement['MUX']['Q1'],
                                   measurement['MUY']['Q2']], [1, 1]),
            # TODO measured errors not in the file
            'ERROR': np.array([0.001, 0.001])
        })

    if any(k in filtered_keys for k in ["BBX", "BBY"]):
        # it's the same data as BETX, BETY but after filtering it will be different
        if beta_file_name == "getgetbeta":
            measurement['BBX'] = getllm_data.beta_x
            measurement['BBY'] = getllm_data.beta_y
        elif beta_file_name == 'getampbeta':
            measurement['BBX'] = getllm_data.amp_beta_x
            measurement['BBY'] = getllm_data.amp_beta_y
        elif beta_file_name == "getkmodbeta":
            measurement['BBX'] = getllm_data.kmod_beta_x
            measurement['BBY'] = getllm_data.kmod_beta_y
        else:
            raise ValueError("Beta filename '{:s}' not recognized".format(beta_file_name))

    if any(k in filtered_keys for k in ["BETX", "BETY"]):
        if beta_file_name == "getgetbeta":
            measurement['BETX'] = getllm_data.beta_x
            measurement['BETY'] = getllm_data.beta_y
        elif beta_file_name == 'getampbeta':
            measurement['BETX'] = getllm_data.amp_beta_x
            measurement['BETY'] = getllm_data.amp_beta_y
        elif beta_file_name == "getkmodbeta":
            measurement['BETX'] = getllm_data.kmod_beta_x
            measurement['BETY'] = getllm_data.kmod_beta_y
        else:
            raise ValueError("Beta filename '{:s}' not recognized".format(beta_file_name))

    if any(k in filtered_keys for k in ["DX", "DY"]):
        measurement['DX'] = getllm_data.disp_x
        measurement['DY'] = getllm_data.disp_y

    if "NDX" in filtered_keys:
        try:
            measurement['NDX'] = getllm_data.norm_disp
        except IOError:
            LOG.warning("No good dispersion or inexistent file getNDx")
            LOG.warning("Correction will not take into account NDx")
            filtered_keys.remove('NDX')

    return filtered_keys, measurement


def _get_varlist(accel_cls, variables, virt_flag):  # TODO: Virtual?
    return np.array(accel_cls.get_variables(classes=variables))


def _read_tfs(path, file_name):
    return tfs.read_tfs(os.path.join(path, file_name))


# Parameter filtering ########################################################


def _filter_measurement(keys, meas, model, errorbar, w_dict, e_dict, m_dict):
    """ Filteres measurements and renames columns to VALUE, ERROR, WEIGHT"""
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
    model_filter = np.abs(new.loc[:, 'VALUE'].values -
                          meas.loc[common_bpms, key + 'MDL'].values) < modelcut

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
                     join_axes=[pd.Index(varslist)]  # other axes to use (pd Index obj required)
                     )


def _join_columns(col, meas, keys):
    """ Retuns vector: N= #BPMs * #Parameters (BBX, MUX etc.) """
    return np.concatenate([meas[key].loc[:, col].values for key in keys], axis=0)


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


def _create_madx_script(accel_cls, nominal_model, optics_file,
                        template_file_path, output_path):
    with open(template_file_path, "r") as textfile:
        madx_template = textfile.read()
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


def write_knob(path, delta):
    a = datetime.datetime.fromtimestamp(time.time())
    changeparameters_path = os.path.join(path, 'changeparameters.knob')
    delta_out = - delta.loc[:, ["DELTA"]]
    delta_out.headers["PATH"] = path
    delta_out.headers["DATE"] = str(a.ctime())
    tfs.write_tfs(changeparameters_path, delta_out, save_index="NAME")


def writeparams(path, delta, append=False):
    mode = 'a' if append else 'w'
    with open(os.path.join(path, "changeparameters.madx"), mode) as madx_script:
        for var in delta.index.values:
            value = delta.loc[var, "DELTA"]
            madx_script.write(var + " = " + var
                              + (" + " if value > 0 else " ")
                              + str(value) + ";\n"
                              )


# Main invocation ############################################################

if __name__ == "__main__":
    global_correction()


