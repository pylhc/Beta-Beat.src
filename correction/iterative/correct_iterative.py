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
from utils import tfs_pandas as tfs  # noqa
from utils.contexts import timeit
from utils import logging_tools
from utils.entrypoint import entrypoint, EntryPointParameters
from model import manager  # noqa

LOG = logging_tools.get_logger(__name__)

DEV_NULL = os.devnull


# Configuration ##################################################################


DEFAULTS = {
    "optics_file": None,
    "output_path": None,
    "singular_value_cut": 0.01,
    "optics_params": ['MUX', 'MUY', 'BBX', 'BBY', 'NDX', 'Q'],
    "modelcut": [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
    "errorcut": [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
    "weights_on_quantities": [1., 1., 0., 0., 0., 10.],
    "use_errorbars": False,
    "variables": ["MQM", "MQT", "MQTL", "MQY"],
    "beta_file_name": "getbeta",
    "virt_flag": False,
    "method": "newton",
    "max_iter": 3,
    "eps": None,
}


# Define functions here, to new optics params
def _get_measurement_filters():
    return {
        'MUX': _get_filtered_phases, 'MUY': _get_filtered_phases,
        'BBX': _get_filtered_betas, 'BBY': _get_filtered_betas,
        'NDX': _get_filtered_disp, 'Q': _get_tunes
        }


def _get_response_filters():
    return {
        'MUX': _get_phase_response, 'MUY': _get_phase_response,
        'BBX': _get_betabeat_response, 'BBY': _get_betabeat_response,
        'NDX': _get_disp_response, 'Q': _get_tune_response
        }


def _get_model_appenders():
    return {
        'MUX': _get_model_phases, 'MUY': _get_model_phases,
        'BBX': _get_model_betas, 'BBY': _get_model_betas,
        'NDX': _get_model_disp, 'Q': _get_model_tunes
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
        help="List of parameters to correct upon (e.g. [BETX, BETY])",
        name="optics_params",
        type=str,
        nargs="*",
        default=DEFAULTS["optics_params"],
    )
    params.add_parameter(
        flags="--optics_file",
        help=("Path to the optics file to use, usually modifiers.madx. If "
              "not present will default to model_path/modifiers.madx"),
        name="optics_file",
        default=DEFAULTS["optics_file"],
    )
    params.add_parameter(
        flags="--output",
        help=("Path to the directory where to write the ouput files, will "
              "default to the --meas input path."),
        name="output_path",
        default=DEFAULTS["output_path"],
    )
    params.add_parameter(
        flags="--svd_cut",
        help="",  # TODO
        name="singular_value_cut",
        type=float,
        default=DEFAULTS["singular_value_cut"],
    )
    params.add_parameter(
        flags="--model_cut",
        help=("Reject BPMs whose deviation to the model is higher than the "
              "correspoding input. Input should be: Phase,Betabeat,NDx"),
        name="modelcut",
        default=DEFAULTS["modelcut"],
    )
    params.add_parameter(
        flags="--error_cut",
        help=("Reject BPMs whose error bar is higher than the "
              "correspoding input. Input should be: Phase,Betabeat,NDx"),
        name="errorcut",
        nargs="*",
        type=float,
        default=DEFAULTS["errorcut"],
    )
    params.add_parameter(
        flags="--weights",
        help=("Weight to apply to each measured quatity. Input shoud be: "
              "PhaseX,PhaseY,BetaX,BetaY,NDx,Q"),
        name="weights_on_quantities",
        nargs="*",
        type=float,
        default=DEFAULTS["weights_on_quantities"],
    )
    params.add_parameter(
        flags="--use_errorbars",
        help=("If True, it will take into account the measured errorbars "
              "in the correction."),
        name="use_errorbars",
        action="store_" + str(not DEFAULTS["use_errorbars"]).lower(),
    )
    params.add_parameter(
        flags="--variables",
        help="Comma separated names of the variables classes to use.",
        name="variable_categories",
        default=DEFAULTS["variables"],
    )
    params.add_parameter(
        flags="--beta_file_name",
        help="Prefix of the beta file to use. E.g.: getkmodbeta",
        name="beta_file_name",
        default=DEFAULTS["beta_file_name"],
    )
    params.add_parameter(
        flags="--virt_flag",
        help="If true, it will use virtual correctors.",
        name="virt_flag",
        action="store_" + str(not DEFAULTS["virt_flag"]).lower(),
    )
    params.add_parameter(
        flags="--method",
        help="Optimization method to use.",
        name="method",
        type=str,
        default=DEFAULTS["method"],
        choices=["newton"]
    )
    params.add_parameter(
        flags="--max_iter",
        help="Maximum number of correction iterations to perform.",
        name="max_iter",
        type=int,
        default=DEFAULTS["max_iter"],
    )
    params.add_parameter(
        flags="--eps",
        help="Convergence criterion. If <|delta(PARAM * WEIGHT)|> < eps, stop iteration.",
        name="eps",
        type=float,
        default=DEFAULTS["eps"],
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
        _setup_debug()

    with timeit(lambda t: LOG.debug("  Total time for Global Correction: {:f}s".format(t))):
        # get unset paths from other paths
        if opt.optics_file is None:
            opt.optics_file = os.path.join(os.path.dirname(opt.model_twiss_path),
                                       "modifiers.madx")
        if opt.output_path is None:
            opt.output_path = opt.meas_dir_path

        template_file_path = os.path.join(os.path.dirname(__file__) , "job.twiss_python.madx")

        # get accelerator class
        accel_cls = manager.get_accel_class(accel_opt)

        # convert numbers to dictionaries
        w_dict = dict(zip(opt.optics_params, opt.weights_on_quantities))
        m_dict = dict(zip(opt.optics_params, opt.modelcut))
        e_dict = dict(zip(opt.optics_params, opt.errorcut))

        # read data from files
        nominal_model = tfs.read_tfs(opt.model_twiss_path)
        full_response = _load_fullresponse(opt.fullresponse_path)
        varslist = _get_varlist(accel_cls, opt.variable_categories, opt.virt_flag)

        # apply filters to data
        optics_params, meas_dict = _scan_meas_dir(
            opt.optics_params,
            opt.meas_dir_path, opt.beta_file_name,
            w_dict,
        )
        optics_params, meas_dict = _filter_measurement(
            optics_params, meas_dict, opt.nominal_model,
            opt.use_errorbars, w_dict, e_dict, m_dict
        )
        full_response = _filter_response_columns(full_response, meas_dict, optics_params)
        meas_dict = _append_model_to_measurement(nominal_model, meas_dict, optics_params)

        if opt.debug:
            _print_rms(meas_dict, optics_params)

        _dump(os.path.join(opt.output_path, "measurement_dict.bin"), meas_dict)
        deltas = _calculate_deltas(full_response, meas_dict, optics_params, varslist,
                                   cut=opt.singular_value_cut)
        writeparams(deltas, varslist, opt.output_path)

        for idx in range(opt.max_iter):
            LOG.debug("Running MADX, iteration {:d} of {:d}".format(idx + 1, opt.max_iter))

            madx_script = _create_madx_script(accel_cls, nominal_model, opt.optics_file,
                                              template_file_path, opt.output_path)
            _callMadx(madx_script)  # TODO
            new_model_path = os.path.join(opt.output_path, "twiss_" + str(idx) + ".dat")
            shutil.copy2(os.path.join(opt.output_path, "twiss_corr.dat"),
                         new_model_path)
            new_model = tfs.read_tfs(new_model_path)
            meas_dict = _append_model_to_measurement(new_model, meas_dict, optics_params)
            if opt.debug:
                _print_rms(meas_dict, optics_params)
            ideltas = _calculate_deltas(
                full_response, meas_dict, optics_params, varslist,
                cut=opt.singular_value_cut,
            )
            writeparams(ideltas, varslist, opt.output_path, append=True)
            deltas = deltas + ideltas
            LOG.debug("Cumulative deltas:" + str(np.sum(np.abs(deltas))))
        write_knob(deltas, varslist)


# Helper functions #############################################################


def _print_rms(meas, keys):
    """ Prints current RMS status """
    for key in keys:
        message = key + " RMS: " + str(np.std(meas[key].loc[:, 'DIFF'].values))
        message += "\n"
        message += key + " Weighted RMS: " + str(np.sqrt(
            np.average(np.square(meas[key].loc[:, 'DIFF'].values),
                       weights=np.square(meas[key].loc[:, 'WEIGHT'].values))
        ))
        LOG.debug(message)


def _load_fullresponse(full_response_path):
    """
    Full response is dictionary of MUX/Y, BBX/Y, NDX and Q gradients upon
    a change of a single quadrupole strength
    """
    LOG.debug("Starting loading Full Response optics")
    with open(full_response_path, "r") as full_response_file:
        full_response_data = pickle.load(full_response_file)
    LOG.debug("Loading ended")
    return full_response_data


def _scan_meas_dir(keys, meas_dir_path, beta_file_name, w_dict):
    measurement = {}
    filtered_keys = keys[:]  # Clone original keys list
    if w_dict['MUX'] or w_dict['MUY'] or w_dict['Q']:
        try:
            measurement['MUX'] = _read_tfs(meas_dir_path, "getphasex_free.out")
            measurement['MUY'] = _read_tfs(meas_dir_path, "getphasey_free.out")
        except IOError:
            measurement['MUX'] = _read_tfs(meas_dir_path, "getphasex.out")
            measurement['MUY'] = _read_tfs(meas_dir_path, "getphasey.out")
        measurement['Q'] = pd.DataFrame({
            'NAME': pd.Categorical(['Q1', 'Q2']),
            # Just fractional tunes:
            'VALUE': np.remainder([measurement['MUX']['Q1'],
                                   measurement['MUY']['Q2']], [1, 1]),
            # TODO measured errors not in the file
            'ERROR': np.array([0.001, 0.001])
        })
    else:
        filtered_keys.remove('MUX')
        filtered_keys.remove('MUY')
        filtered_keys.remove('Q')

    if w_dict['BBX'] or w_dict['BBY']:
        try:
            measurement['BBX'] = _read_tfs(meas_dir_path,
                                           beta_file_name + "x_free.out")
            measurement['BBY'] = _read_tfs(meas_dir_path,
                                           beta_file_name + "y_free.out")
        except IOError:
            measurement['MUX'] = _read_tfs(meas_dir_path,
                                           beta_file_name + "x.out")
            measurement['MUY'] = _read_tfs(meas_dir_path,
                                           beta_file_name + "y.out")
    else:
        filtered_keys.remove('BBX')
        filtered_keys.remove('BBY')

    if w_dict['NDX']:
        try:
            measurement['NDX'] = _read_tfs(meas_dir_path, "getNDx.out")
        except IOError:
            LOG.warning("No good dispersion or inexistent file getNDx")
            LOG.warning("Correction will not take into account NDx")
            filtered_keys.remove('NDX')
    else:
        filtered_keys.remove('NDX')

    return filtered_keys, measurement


def _get_varlist(accel_cls, variables, virt_flag):  # TODO: Virtual?
    return np.array(accel_cls.get_variables(classes=variables))


def _read_tfs(path, file_name):
    return tfs.read_tfs(os.path.join(path, file_name))


# Parameter filtering ########################################################

# TODO: Add error and model cuts:
def _filter_measurement(keys, measurement, model, w, w_dict,  e_dict, m_dict):
    filters = _get_measurement_filters()
    meas = {}
    new_keys = []
    for key in keys:
        if w_dict[key] > 0.0:
            meas[key] = filters[key](
                measurement[key], model, w_dict[key], w,
                modelcut=m_dict[key], errorcut=e_dict[key]
            )
            new_keys.append(key)
    return new_keys, meas


def _get_filtered_phases(meas, model, weight, erwg,
                         modelcut=0.05, errorcut=0.035):
    new = pd.merge(meas, model, how='inner', on='NAME', suffixes=('', 'm'))
    second_bpm_in = np.in1d(new.loc[:, 'NAME2'].values,
                            new.loc[:, 'NAME'].values)
    if 'PHYMDL' in new.columns.values:
        plane = 'Y'
    else:
        plane = 'X'
    new['PHMOD'] = new.loc[:, 'PH' + plane + 'MDL'].values
    new['VALUE'] = new.loc[:, 'PHASE' + plane].values
    new['ERROR'] = new.loc[:, 'STDPH' + plane].values
    good_phase = (
        (new.loc[:, 'ERROR'].values < errorcut) &
        (np.abs(new.loc[:, 'VALUE'].values - new.loc[:, 'PHMOD'].values) <
            modelcut)
    )
    good_bpms = good_phase & second_bpm_in
    # Removing the last pair as it goes beyond the end of lattice
    good_bpms[-1] = False
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = (new.loc[:, 'WEIGHT'].values /
                         new.loc[:, 'ERROR'].values)
    LOG.debug("Number of BPM pairs in plane " +
              plane + ": " + str(np.sum(good_bpms)))
    return new.loc[good_bpms, ['NAME', 'NAME2', 'VALUE', 'ERROR', 'WEIGHT']]


def _get_filtered_betas(meas, model, weight, erwg,
                        modelcut=0.2, errorcut=0.02):
    # Beta-beating and its error RELATIVE as shown in GUI
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


def _get_filtered_disp(meas, model, weight, erwg, modelcut=0.2, errorcut=0.02):
    new = pd.merge(meas, model, how='inner', on='NAME', suffixes=('', 'm'))
    new['VALUE'] = new.loc[:, 'NDX'].values
    new['ERROR'] = new.loc[:, 'STDNDX'].values
    good_bpms = (
        (new.loc[:, 'ERROR'].values < errorcut) &
        (np.abs(new.loc[:, 'VALUE'].values -
                new.loc[:, 'NDXMDL'].values) < modelcut)
    )
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = (new.loc[:, 'WEIGHT'].values /
                         new.loc[:, 'ERROR'].values)
    LOG.debug("Number of BPMs with NDx : " + str(np.sum(good_bpms)))
    return new.loc[good_bpms, ['NAME', 'VALUE', 'ERROR', 'WEIGHT']]


def _get_tunes(meas, mod, weight, erwg, modelcut=0.1, errorcut=0.027):
    meas['WEIGHT'] = weight
    if erwg:
        meas['WEIGHT'] = (meas.loc[:, 'WEIGHT'].values /
                          meas.loc[:, 'ERROR'].values)
    LOG.debug("Number of tune measurements: " + str(len(meas.index.values)))
    return meas


# Response filtering ##########################################################


def _filter_response_columns(response, measurement, keys):
    filters = _get_response_filters()
    new_resp = {}
    for key in keys:
        new_resp[key] = filters[key](response[key], measurement[key])
    return new_resp


def _get_phase_response(resp, meas):
    new = resp.loc[:, meas.loc[:, 'NAME'].values]
    new.sub(resp.as_matrix(columns=meas.loc[:, 'NAME2'].values), axis=0)
    return -new  # As we did subtraction name-name2


def _get_betabeat_response(resp, meas):
    return resp.loc[:, meas.loc[:, 'NAME'].values]


def _get_disp_response(resp, meas):
    return resp.loc[:, meas.loc[:, 'NAME'].values]


def _get_tune_response(resp, meas):
    return resp


# Model appending #############################################################


def _append_model_to_measurement(model, measurement, keys):
    appenders = _get_model_appenders()
    meas = {}
    for key in keys:
        meas[key] = appenders[key](model, measurement[key], key)
    return meas


def _get_model_phases(model, meas, key):
    new = model.set_index('NAME')
    meas['MODEL'] = (new.loc[meas.loc[:, 'NAME2'].values, key].values -
                     new.loc[meas.loc[:, 'NAME'].values, key].values)
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
    return meas


def _get_model_betas(model, meas, key):
    new = model.set_index('NAME')
    if key == 'BBX':
        meas['MODEL'] = new.loc[meas.loc[:, 'NAME'].values, 'BETX'].values
    else:
        meas['MODEL'] = new.loc[meas.loc[:, 'NAME'].values, 'BETY'].values
    meas['DIFF'] = (
        (meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values) /
        meas.loc[:, 'MODEL'].values
    )
    return meas


def _get_model_disp(model, meas, key):
    new = model.set_index('NAME')
    meas['MODEL'] = (
        new.loc[meas.loc[:, 'NAME'].values, 'DX'].values /
        np.sqrt(new.loc[meas.loc[:, 'NAME'].values, 'BETX'].values)
    )
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
    return meas


def _get_model_tunes(model, meas, key):
    # We want just fractional tunes
    meas['MODEL'] = np.remainder([model['Q1'], model['Q2']], [1, 1])
    meas['DIFF'] = meas.loc[:, 'VALUE'].values - meas.loc[:, 'MODEL'].values
    return meas


# Small Helpers ################################################################


def _dump(path_to_dump, content):
    with open(path_to_dump, 'wb') as dump_file:
        cPickle.Pickler(dump_file, -1).dump(content)


def _calculate_deltas(resp, meas, keys, varslist, cut=0.01, append=False):
    # TODO: think about output form
    # Return NumPy array: each row for magnet, columns for measurements
    response_matrix = _filter_response_rows_and_join(resp, keys, varslist)
    weight_vector = _join_weights(meas, keys)
    diff_vector = _join_diffs(meas, keys)
    delta = np.dot(
        np.linalg.pinv(np.transpose(response_matrix * weight_vector), cut),
        diff_vector * weight_vector
    )
    LOG.debug("Delta calculation: ")
    LOG.debug(np.std(np.abs(diff_vector * weight_vector)))
    LOG.debug(np.std(np.abs((diff_vector * weight_vector) -
                            np.dot(delta, response_matrix * weight_vector))))
    return delta


def _filter_response_rows_and_join(resp, keys, varslist):
    return pd.concat([resp[key] for key in keys],
                     axis=1,
                     join_axes=[pd.Index(varslist)])


def _join_diffs(meas, keys):
    return np.concatenate([meas[key].loc[:, 'DIFF'].values for key in keys])


def _join_weights(meas, keys):
    return np.concatenate([meas[key].loc[:, 'WEIGHT'].values for key in keys])


def _load_model(path_to_optics_files_dir, iteration):
    return _read_tfs(path_to_optics_files_dir,
                     "twiss_" + str(iteration) + ".dat")


def _callMadx(madx_script, debug):
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
            "CROSSING_ON": 0,  # TODO
            "NUM_BEAM": beam,
            "PATH": output_path,
            "DPP": 0.0,
            "QMX": qx,
            "QMY": qy,
            "COR": "changeparameters.madx",
        }
        madx_script = madx_template % replace_dict
    return madx_script


def write_knob(deltas, varslist, path="./"):
    # TODO: Use tfs_file_writer
    a = datetime.datetime.fromtimestamp(time.time())
    changeparameters_path = os.path.join(path, 'changeparameters.knob')
    with open(changeparameters_path, "w") as changeparameters:
        changeparameters.write("@ PATH %s" + path + "\n")
        changeparameters.write("@ DATE %s" + str(a.ctime()) + "\n")
        changeparameters.write("* NAME DELTA\n")
        changeparameters.write("$ %s %le\n")
        for i, var in enumerate(varslist):
            changeparameters.write(var + '   ' + str(-deltas[i]) + '\n')


def writeparams(deltas, variables, path, append=False):
    mode = 'a' if append else 'w'
    mad_script = open(os.path.join(path, "changeparameters.madx"), mode)
    for i, var in enumerate(variables):
        if deltas[i] > 0:
            mad_script.write(var + " = " + var + " + " +
                             str(deltas[i]) + ";\n")
        else:
            mad_script.write(var + " = " + var + " " +
                             str(deltas[i]) + ";\n")
    mad_script.close()


def _setup_debug():
    """ Setup Logger for debugging mode """
    import datetime
    LOG.setLevel(logging_tools.DEBUG)
    LOG.debug("Running in Debug-Mode.")

    now = str(datetime.datetime.now().isoformat())
    log_file = now + os.path.abspath(__file__).replace(".pyc", "").replace(".py", "") + ".log"

    LOG.debug("Writing log to file '{:s}'.".format(log_file))
    file_handler = logging_tools.file_handler(log_file)
    logging_tools.add_root_handler(file_handler)


# Main invocation ############################################################

if __name__ == "__main__":
    global_correction()


