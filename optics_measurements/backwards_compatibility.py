"""
.. module: backwards_compatibility

Created on 26/07/18

:author: Lukas Malina

It contains wrappers to new data structures such that they can be used in old functions.
"""
import warnings
from utils import logging_tools
from tfs_files import tfs_file_writer

LOGGER = logging_tools.get_logger(__name__)
warnings.simplefilter('always', DeprecationWarning)


def _get_output_tfs_files(header_dict, filename, outpath):
    tfs_file = tfs_file_writer.TfsFileWriter(filename, outputpath=outpath)
    for (key, value) in header_dict.items():
        tfs_file.add_string_descriptor(key, value)
    tfs_file.add_string_descriptor("FILENAME", filename)
    warnings.warn("Manual construction of tfs_file is deprecated, "
                  "consider using tfs_files.tfs_pandas",
                  DeprecationWarning, 2)
    return tfs_file


class _TwissData(object):
    warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                  DeprecationWarning, 2)

    def __init__(self, inputfiles):

        self.zero_dpp_x = inputfiles.zero_dpp_frames("X")
        self.zero_dpp_y = inputfiles.zero_dpp_frames("Y")
        self.non_zero_dpp_x = inputfiles["X"]
        self.non_zero_dpp_y = inputfiles["Y"]
        warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                      DeprecationWarning, 2)

    def has_zero_dpp_x(self):
        warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                      DeprecationWarning, 2)
        return 0 != len(self.zero_dpp_x)

    def has_non_zero_dpp_x(self):
        warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                      DeprecationWarning, 2)
        return 0 != len(self.non_zero_dpp_x)

    def has_zero_dpp_y(self):
        warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                      DeprecationWarning, 2)
        return 0 != len(self.zero_dpp_y)

    def has_non_zero_dpp_y(self):
        warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                      DeprecationWarning, 2)
        return 0 != len(self.non_zero_dpp_y)

    def has_no_input_files(self):
        warnings.warn("_TwissData is deprecated, consider using measure_optics.InputFiles",
                      DeprecationWarning, 2)
        return (not self.has_zero_dpp_x() and not self.has_zero_dpp_y()
                and not self.has_non_zero_dpp_x() and not self.has_non_zero_dpp_y())


class _TuneData(object):
    warnings.warn("_TuneData is deprecated, consider using optics_measurements.tune.TuneDict",
                  DeprecationWarning, 2)

    def __init__(self, tune_dict):
        self.q1 = tune_dict["X"]["Q"]  # Driven horizontal tune
        self.q2 = tune_dict["Y"]["Q"]  # Driven vertical tune
        # Free is from analytic equation
        self.q1f = tune_dict["X"]["QF"]  # Free horizontal tune
        self.q2f = tune_dict["Y"]["QF"]  # Free vertical tune
        self.q1mdl = tune_dict["X"]["QM"]
        self.q2mdl = tune_dict["Y"]["QM"]
        self.q1mdlf = tune_dict["X"]["QFM"]
        self.q2mdlf = tune_dict["Y"]["QFM"]
        warnings.warn("_TuneData is deprecated, consider using optics_measurements.tune.TuneDict",
                      DeprecationWarning, 2)


class _PhaseData(object):
    warnings.warn("_PhaseData is deprecated, consider using optics_measurements.phase.PhaseDict",
                  DeprecationWarning, 2)

    def __init__(self, phase_dict):
        self.ac2bpmac_x = phase_dict["X"]["ac2bpm"]
        self.ac2bpmac_y = phase_dict["Y"]["ac2bpm"]

        self.phase_advances_x = phase_dict["X"]["D"]
        self.phase_advances_free_x = phase_dict["X"]["F"]
        self.phase_advances_free2_x = phase_dict["X"]["F2"]
        self.phase_advances_y = phase_dict["Y"]["D"]
        self.phase_advances_free_y = phase_dict["Y"]["F"]
        self.phase_advances_free2_y = phase_dict["Y"]["F2"]
        warnings.warn("_PhaseData is deprecated, consider using optics_measurements.phase.PhaseDict",
                      DeprecationWarning, 2)
