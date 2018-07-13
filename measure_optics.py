import os
from os.path import isfile, join
import sys
import traceback
import datetime
from time import time
import re
from collections import OrderedDict
import pandas as pd
import numpy as np
from GetLLM import algorithms, optics_input
from GetLLM.GetLLMError import GetLLMError, CriticalGetLLMError
from model.accelerators.accelerator import AccExcitationMode
from utils import tfs_pandas, logging_tools, iotools


VERSION = 'V3.0.0 Dev'
DEBUG = sys.flags.debug  # True with python option -d! ("python -d measure_optics.py...") (vimaier)
LOGGER = logging_tools.get_logger(__name__)
PLANES = ('X', 'Y')


def measure_optics(input_files, measure_input):
    LOGGER.info("Starting GetLLM " + VERSION)
    global __getllm_starttime
    __getllm_starttime = time()
    iotools.create_dirs(measure_input.outputdir)
    logging_tools.add_module_handler(logging_tools.file_handler(os.path.join(measure_input.outputdir, "getllm.log")))
    header_dict = _get_header()
    if sys.flags.debug:
        LOGGER.info("     DEBUG ON")

    tune_d = _TuneData()
    tune_d.q1mdl = measure_input.accelerator.get_model_tfs().headers["Q1"]
    tune_d.q2mdl = measure_input.accelerator.get_model_tfs().headers["Q2"]

    """
    Construct pseudo-double plane BPMs
    TODO This should be in accelerator class
    if (accelerator.__name__ == "SPS" or "RHIC" in accelerator.__name__) and twiss_d.has_zero_dpp_x() and twiss_d.has_zero_dpp_y():
        [pseudo_list_x, pseudo_list_y] = algorithms.helper.pseudo_double_plane_monitors(accelerator.get_model_tfs(), twiss_d.zero_dpp_x, twiss_d.zero_dpp_y, bpm_dictionary)
    else:
        Initialize variables otherwise calculate_coupling would raise an exception(vimaier)
        pseudo_list_x = None
        pseudo_list_y = None
    """
    print_time("BEFORE_PHASE", time() - __getllm_starttime)
    #-------- START Phase for beta calculation with best knowledge model in ac phase compensation
    try:
        phase_d_bk, tune_d = algorithms.phase.calculate_phase(measure_input, input_files)
    except:
        _tb_()
        # if phase crashed, none of the subsequent algorithms can run. Thus
        raise CriticalGetLLMError("get phase crashed. None of the following algorithms can work hence GetLLM will crash now. Good bye!")
    print_time("AFTER_PHASE", time() - __getllm_starttime)
    #-------- START coupling.
    try:
        tune_d = algorithms.coupling.calculate_coupling(measure_input, input_files, phase_d_bk, tune_d)
    except:
        _tb_()
    if measure_input.only_coupling:
        LOGGER.info("GetLLM was only calculating coupling. Skipping the rest and returning ...")
        return
    try:
        beta_d, beta_driven_x, beta_free_x = algorithms.beta.calculate_beta_from_phase(measure_input, input_files, tune_d, phase_d_bk)
    except:
        _tb_()
    if measure_input.three_bpm_method:
        print_time("AFTER_BETA_FROM_PHASE", time() - __getllm_starttime)
    else:
        print_time("AFTER_A_NBPM", time() - __getllm_starttime)
    # try:
    #     algorithms.lobster.get_local_observable( phase_d_bk, getllm_d.accelerator.get_model_tfs(), files_dict, tune_d.q1f)
    # except:
    #     _tb_()
    try:
        beta_d = algorithms.beta_from_amplitude.calculate_beta_from_amplitude(measure_input, input_files, tune_d, beta_d)
    except:
        _tb_()
    # in the following functions, nothing should change, so we choose the models now
    mad_twiss = measure_input.accelerator.get_model_tfs()
    mad_elements = measure_input.accelerator.get_elements_tfs()
    if measure_input.accelerator.excitation != AccExcitationMode.FREE:
        mad_ac = measure_input.accelerator.get_driven_tfs()
    else:
        mad_ac = mad_twiss
    try:
        algorithms.interaction_point.betastar_from_phase(measure_input.accelerator, phase_d_bk, mad_twiss)
    except:
        _tb_()
    try:
        algorithms.dispersion.calculate_orbit_and_dispersion(input_files, tune_d, mad_twiss, header_dict, measure_input.orbit_unit, measure_input.max_closed_orbit, beta_driven_x, measure_input.outputdir)
    except:
        _tb_()
    #------ Start get Q,JX,delta
    try:
        inv_x, inv_y = algorithms.kick.calculate_kick(mad_twiss, mad_ac, measure_input, input_files, beta_d, phase_d_bk, measure_input.outputdir, header_dict)
    except:
        _tb_()
    if measure_input.nonlinear:
        try:
            algorithms.resonant_driving_terms.calculate_RDTs(mad_twiss, measure_input, input_files, phase_d_bk, tune_d, inv_x, inv_y)
        except:
            _tb_()
        # TODO: what does this?
        # files_dict = algorithms._calculate_getsextupoles(twiss_d, phase_d_bk, accelerator.get_model_tfs(), files_dict, tune_d.q1f)
        # files_dict = algorithms.chi_terms.calculate_chiterms(getllm_d, twiss_d, accelerator.get_model_tfs(), files_dict)
    print_time("FINISH", time() - __getllm_starttime)


def _get_header():
    return OrderedDict([('GetLLMVersion', VERSION),
                        ('Command', sys.executable + " '" + "' '".join([] + sys.argv) + "'"),
                        ('CWD', os.getcwd()),
                        ('Date', datetime.datetime.today().strftime("%d. %B %Y, %H:%M:%S"))])
    # TODO add model directory


def _tb_():
    if sys.stdout.isatty():
        err_exc = re.sub(r"line\s([0-9]+)", "\33[1mline \33[38;2;80;160;255m\\1\33[0m\33[21m", traceback.format_exc())
        err_exc = re.sub("File\\s\"([^\"]+)\",", "File \33[38;2;0;255;100m\\1\33[0m", err_exc)
        err_excs = err_exc.split("\n")
        for line in err_excs:
            LOGGER.error(line)
    else:
        LOGGER.error(traceback.format_exc())


def print_time(index, t):
    LOGGER.debug(":::  GetLLM time  >>>>>>>>>> {:8.3f} s".format(t))


class _TuneData(object):
    """
    Used as data structure to hold tunes and phase advances.
    """
    def __init__(self):
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
            for plane in PLANES:
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
        return self[plane][np.logical_not(self.get_dpps(plane))]

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

    def calibrate(self, calibs):
        if calibs is None:
            pass
        for plane in PLANES:
            for i in range(len(self[plane])):
                data = pd.merge(self[plane][i].loc[:, "AMP" + plane], calibs[plane], how='left',
                                left_index=True, right_index=True).fillna(
                    value={"CALIBRATION": 1., "ERROR_CALIBRATION": 0.})
                self[plane][i]["AMP" + plane] = self[plane][i].loc[:, "AMP" + plane] * data.loc[:,"CALIBRATION"]
                self[plane][i]["ERRAMP" + plane] = data.loc[:, "ERROR_CALIBRATION"]  # TODO

    def use_average_tune(self, no_averaged_tune):
        if no_averaged_tune:
            pass
        for plane in PLANES:
            for i in range(len(self[plane])):
                self[plane][i].rename(columns={"AVG_MU" + plane: "MU" + plane, "MU" + plane: "OLD_MU" + plane})

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


def _copy_calibration_files(outputdir, calibrationdir):
    if calibrationdir is None:
        return None
    calibs = {}
    for plane in PLANES:
        cal_file = "calibration_{}.out".format(plane.lower())
        iotools.copy_item(join(calibrationdir, cal_file), join(outputdir, cal_file))
        calibs[plane] = tfs_pandas.read_tfs(join(outputdir, cal_file)).set_index("NAME")
    return calibs


if __name__ == "__main__":
    arguments = optics_input.parse_args()
    inputs = InputFiles(arguments.files)
    iotools.create_dirs(arguments.outputdir)
    calibrations = _copy_calibration_files(arguments.outputdir, arguments.calibrationdir)
    inputs.calibrate(calibrations)
    inputs.use_average_tune(arguments.no_averaged_tune)
    measure_optics(inputs, *optics_input.parse_args())
