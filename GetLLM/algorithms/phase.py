"""
.. module: tune

Created on 18/07/18

:author: Lukas Malina

It computes betatron phase advances and provides structures to store them.
"""
import sys
from os.path import join
import numpy as np
import pandas as pd
import compensate_excitation
from model.accelerators.accelerator import AccExcitationMode
from utils import logging_tools, stats, tfs_pandas

DEBUG = sys.flags.debug  # True with python option -d! ("python -d GetLLM.py...") (vimaier)
LOGGER = logging_tools.get_logger(__name__)
PLANES = ("X", "Y")


def calculate_phase(measure_input, input_files, tunes, header_dict):
    """
    Calculates phase advances and fills the following files:
        getphase(tot)(x/y)(_free).out

    Parameters:
        measure_input: the input object including settings and the accelerator class
        input_files: includes measurement tfs_pandas
        tunes: TunesDict object containing measured and model tunes
        header_dict: part of the header common for all output files

    Returns:
        an instance of PhaseDict filled with the results of get_phases
    """
    phase_d = PhaseDict()
    accelerator = measure_input.accelerator
    try:
        model_free = accelerator.get_best_knowledge_model_tfs()
    except AttributeError:
        model_free = accelerator.get_model_tfs()
    if accelerator.excitation == AccExcitationMode.FREE:
        model_of_measurement = model_free
    else:
        model_of_measurement = accelerator.get_driven_tfs()
    LOGGER.info('Calculating phase')
    for plane in PLANES:
        LOGGER.debug("tune of measurement files = {}".format(tunes[plane]["Q"]))
        phase_d[plane]["F"], output_dfs = get_phases(measure_input, input_files, model_of_measurement, plane)
        headers = _get_headers(header_dict, tunes, plane)
        for head, df in zip(headers, output_dfs):
            tfs_pandas.write_tfs(join(measure_input.outputdir, head['FILENAME']), df, head)
        if measure_input.accelerator.excitation != AccExcitationMode.FREE:
            phase_d[plane]["D"] = phase_d[plane]["F"]
            phase_d[plane]["ac2bpm"] = compensate_excitation.phase_ac2bpm(
                phase_d[plane]["F"]["MODEL"], tunes[plane]["Q"], tunes[plane]["QF"], plane,
                measure_input.accelerator)
            phase_d[plane]["F"], output_dfs = get_phases(measure_input, input_files, model_free,
                        plane, (tunes[plane]["Q"], tunes[plane]["QF"], phase_d[plane]["ac2bpm"]))
            headers = _get_headers(header_dict, tunes, plane, free=True)
            for head, df in zip(headers, output_dfs):
                tfs_pandas.write_tfs(join(measure_input.outputdir, head['FILENAME']), df, head)
            # phase_d[plane]["F2"]  = _get_free_phase(phase_d[plane]["F"], tune_d[plane]["Q"], tune_d[plane]["QF"], bpmsx, model_driven, model, plane)
    return phase_d


def get_phases(meas_input, input_files, model, plane, compensate=None, no_errors=True):
    """
    Computes phase advances among all BPMs.

    Args:
        meas_input: OpticsInput object
        input_files: InputFiles object
        model: model tfs_panda to be used
        plane: "X" or "Y"
        compensate: (driven_tune,free_tune,ac2bpm object)
        no_errors: if True measured errors shall not be propagated (only their spread)

    Returns:
        dictionary of DataFrames indexed (BPMi x BPMj) yielding phase advance phi_ij
            "MEAS" measured phase advances
            "ERRMEAS" errors of measured phase advances
            "MODEL" model phase advances

            +------++--------+--------+--------+--------+
            |      ||  BPM1  |  BPM2  |  BPM3  |  BPM4  |
            +======++========+========+========+========+
            | BPM1 ||   0    | phi_12 | phi_13 | phi_14 |
            +------++--------+--------+--------+--------+
            | BPM2 || phi_21 |    0   | phi_23 | phi_24 |
            +------++--------+--------+--------+--------+
            | BPM3 || phi_31 | phi_32 |   0    | phi_34 |
            +------++--------+--------+--------+--------+

            The phase advance between BPM_i and BPM_j can be obtained via:
                phase_advances["MEAS"].loc[BPMi,BPMj]
        list of output dataframes(for files)
    """
    phase_frame = pd.DataFrame(model).loc[:, ['S', 'MU' + plane]]
    how = 'outer' if meas_input.union else 'inner'
    phase_frame = pd.merge(phase_frame,
                           input_files.get_joined_frame(plane, ['MU' + plane, 'ERR_MU' + plane],
                                                        zero_dpp=True, how=how),
                           how='inner', left_index=True, right_index=True)
    phases_mdl = phase_frame.loc[:, 'MU' + plane].values
    phase_advances = {"MODEL": _get_square_data_frame((phases_mdl[np.newaxis, :] - phases_mdl[:, np.newaxis]) % 1.0, phase_frame.index)}
    phases_meas = input_files.get_data(phase_frame, 'MU' + plane) * meas_input.accelerator.get_beam_direction()
    phases_errors = input_files.get_data(phase_frame, 'ERR_MU' + plane)

    if compensate is not None:
        (driven_tune, free_tune, ac2bpmac) = compensate
        k_bpmac = ac2bpmac[2]
        phase_corr = ac2bpmac[1] - phases_meas[k_bpmac] + (0.5 * driven_tune)
        phases_meas = phases_meas + phase_corr[np.newaxis, :]
        r = np.sin(np.pi * (driven_tune - free_tune)) / np.sin(np.pi * ((driven_tune + free_tune) % 1.0))
        LOGGER.debug(plane + " compensation lambda = {}".format(r))
        LOGGER.debug(plane + " k_bpmac = {}".format(k_bpmac))
        LOGGER.debug(plane + " psid_ac2bpmac = {}".format(ac2bpmac[1]))
        LOGGER.debug(plane + " bpmac = {}".format(ac2bpmac[0]))
        phases_meas[k_bpmac:, :] = phases_meas[k_bpmac:, :] - driven_tune
        psi = (np.arctan((1 - r) / (1 + r) * np.tan(2 * np.pi * phases_meas)) / (2 * np.pi)) % 0.5
        phases_meas = np.where(phases_meas % 1.0 > 0.5, psi + .5, psi)
        phases_meas[k_bpmac:, :] = phases_meas[k_bpmac:, :] + free_tune

    if phases_meas.ndim < 2:
        phase_advances["MEAS"] = _get_square_data_frame((phases_meas[np.newaxis, :] - phases_meas[:, np.newaxis]) % 1.0, phase_frame.index)
        phase_advances["ERRMEAS"] = _get_square_data_frame(np.zeros((len(phases_meas), len(phases_meas))), phase_frame.index)
        return phase_advances
    if meas_input.union:
        mask = np.isnan(phases_meas)
        phases_meas[mask], phases_errors[mask] = 0.0, np.inf
        if no_errors:
            phases_errors[~mask] = 1e-10
    elif no_errors:
        phases_errors = None
    phases_3d = phases_meas[np.newaxis, :, :] - phases_meas[:, np.newaxis, :]
    if phases_errors is not None:
        errors_3d = phases_errors[np.newaxis, :, :] + phases_errors[:, np.newaxis, :]
    else:
        errors_3d = None
    phase_advances["MEAS"] = _get_square_data_frame(stats.circular_mean(phases_3d, period=1, errors=errors_3d,
                                                 axis=2) % 1.0, phase_frame.index)
    phase_advances["ERRMEAS"] = _get_square_data_frame(stats.circular_error(phases_3d, period=1, errors=errors_3d,
                                                     axis=2), phase_frame.index)
    return phase_advances, [create_output_df(phase_advances, phase_frame, plane), create_output_df(phase_advances, phase_frame, plane, tot=True)]


def create_output_df(phase_advances, model, plane, tot=False):
    meas = phase_advances["MEAS"]
    mod = phase_advances["MODEL"]
    err = phase_advances["ERRMEAS"]
    output_data = model.loc[:, ["S", "MU"+plane]].iloc[:-1, :]
    output_data.rename(columns={'MU' + plane: 'MU' + plane + 'MDL'}, inplace=True)
    if tot:
        output_data["NAME"] = model.index[1:].values
        output_data = output_data.assign(S2=model.at[model.index[0], "S"], NAME2=model.index[0])
        output_data["PHASE" + plane] = meas.values[0, 1:]
        output_data["STDPH" + plane] = err.values[0, 1:]
        output_data["PH{}MDL".format(plane)] = mod.values[0, 1:]
    else:
        output_data["NAME"] = output_data.index
        output_data = output_data.assign(S2=model.loc[:, "S"].values[1:], NAME2=model.index[1:].values)
        output_data["PHASE" + plane] = np.diag(meas.values, k=1)
        output_data["STDPH" + plane] = np.diag(err.values, k=1)
        output_data["PH{}MDL".format(plane)] = np.diag(mod.values, k=1)
    dif = (output_data.loc[:, "PHASE" + plane].values -
           output_data.loc[:, "PH{}MDL".format(plane)].values) % 1.0
    output_data["DELTAPHASE" + plane] = np.where(dif > 0.5, dif - 1.0, dif)
    return output_data


def _get_headers(header_dict, tunes, plane, free=False):
    header = header_dict.copy()
    header['Q1'] = tunes["X"]["QF"] if free else tunes["X"]["Q"]
    header['Q2'] = tunes["Y"]["QF"] if free else tunes["Y"]["Q"]
    header_tot = header.copy()
    header['FILENAME'] = "getphase" + plane.lower() + free * "_free" + ".out"
    header_tot['FILENAME'] = "getphasetot" + plane.lower() + free * "_free" + ".out"
    return [header, header_tot]


def _get_square_data_frame(data, index):
    return pd.DataFrame(data=data, index=index, columns=index)


class PhaseDict(dict):
    """
    Used as data structure to hold phase advances
    """
    def __init__(self):
        init_phases = {"ac2bpm": None, "D": None, "F": None, "F2": None}
        super(PhaseDict, self).__init__(zip(PLANES, (init_phases, init_phases)))


class _PhaseData(object):
    def __init__(self, phase_dict):
        self.ac2bpmac_x = phase_dict["X"]["ac2bpm"]
        self.ac2bpmac_y = phase_dict["Y"]["ac2bpm"]

        self.phase_advances_x = phase_dict["X"]["D"]
        self.phase_advances_free_x = phase_dict["X"]["F"]
        self.phase_advances_free2_x = phase_dict["X"]["F2"]
        self.phase_advances_y = phase_dict["Y"]["D"]
        self.phase_advances_free_y = phase_dict["Y"]["F"]
        self.phase_advances_free2_y = phase_dict["Y"]["F2"]


"""
def write_special_phase_file(plane, phase_advances, tune_x, tune_y, accel):
    plane_mu = "MU" + plane
    plane_tune = tune_x if plane == "X" else tune_y
    meas = phase_advances["MEAS"]
    err = phase_advances["ERRMEAS"]
    bd = accel.get_beam_direction()
    elements = accel.get_elements_tfs()
    for elem1, elem2 in accel.get_important_phase_advances():
        mus1 = elements.loc[elem1, plane_mu] - elements.loc[:, plane_mu]
        minmu1 = abs(mus1.loc[meas.index]).idxmin()
        mus2 = elements.loc[:, plane_mu] - elements.loc[elem2, plane_mu]
        minmu2 = abs(mus2.loc[meas.index]).idxmin()
        try:
            bpm_phase_advance = meas.loc[minmu1, minmu2]
            model_value = elements.loc[elem2, plane_mu] - elements.loc[elem1, plane_mu]
            if (elements.loc[elem1, "S"] - elements.loc[elem2, "S"]) * bd > 0.0:
                bpm_phase_advance += plane_tune
                model_value += plane_tune
            bpm_err = err.loc[minmu1, minmu2]
            phase_to_first = -mus1.loc[minmu1]
            phase_to_second = -mus2.loc[minmu2]
            ph_result = ((bpm_phase_advance + phase_to_first + phase_to_second) * bd)
            model_value = (model_value * bd)
            resultdeg = ph_result % .5 * 360
            if resultdeg > 90:
                resultdeg -= 180
            modeldeg = model_value % .5 * 360
            if modeldeg > 90:
                modeldeg -= 180
            model_desc = [elem1 + "__to__" + elem2 + "___MODL",
                          "{:8.4f}     {:6s} = {:6.2f} deg".format(model_value % 1, "",
                                                                   modeldeg)]
            result_desc = [elem1 + "__to__" + elem2 + "___MEAS",
                           "{:8.4f}  +- {:6.4f} = {:6.2f} +- {:3.2f} deg ({:8.4f} + {:8.4f} [{}, {}])".format(
                               ph_result % 1, bpm_err, resultdeg, bpm_err * 360,
                               bpm_phase_advance,
                               phase_to_first + phase_to_second,
                               minmu1, minmu2) ]
            tfs_file.add_string_descriptor(*model_desc)
            tfs_file.add_string_descriptor(*result_desc)
            LOGGER.debug("")
            LOGGER.debug("::" + " : ".join(model_desc))
            LOGGER.debug("::" + " : ".join(result_desc))
        except KeyError as e:
            LOGGER.error("Couldn't calculate the phase advance because " + e)
"""
