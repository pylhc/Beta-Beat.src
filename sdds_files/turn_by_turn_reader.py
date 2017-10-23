from __future__ import print_function
import sys
import os
import logging
from datetime import datetime
import sdds_reader
import ascii_reader
import numpy as np
import pandas as pd

LOGGER = logging.getLogger(__name__)

TIMESTAMP_NAME = "acqStamp"
NUM_BUNCHES_NAME = "nbOfCapBunches"
NUM_TURNS_NAME = "nbOfCapTurns"
ALL_HOR_POSITIONS_NAME = "horPositionsConcentratedAndSorted"
ALL_VER_POSITIONS_NAME = "verPositionsConcentratedAndSorted"
BPM_NAMES_NAME = "bpmNames"
HOR_BUNCH_ID_NAME = "horBunchId"
HOR_FAILED_BUNCH_ID_NAME = "horBunchIdFailsInTurn"
VER_BUNCH_ID_NAME = "verBunchId"
VER_FAILED_BUNCH_ID_NAME = "verBunchIdFailsInTurn"

HOR, VER = 0, 1

# Number of decimal digits to print in the ASCII file
PRINT_PRECISION = 7
FORMAT_STRING = "{0:." + str(PRINT_PRECISION) + "f}"


# Public ###################

def read_tbt_file(file_path):
    if ascii_reader.is_ascii_file(file_path):
        return [TbtFile.create_from_matrices(
            *ascii_reader.read_ascii_file(file_path)
        )]  # If ASCII return only one TbtFile
    tbt_files = _TbtReader(file_path).read_file()
    return tbt_files


def transform_tbt_to_ascii(file_path, model_path, output_path):
    tbt_files = _TbtReader(file_path).read_file()
    _TbtAsciiWriter(tbt_files,
                    model_path,
                    output_path).transform_tbt_to_ascii()


def write_ascii_file(model_path, output_path,
                     bpm_names_x, matrix_x,
                     bpm_names_y, matrix_y, date,
                     headers_dict=None):
    if headers_dict is None:
        headers_dict = {}
    new_tbt_file = TbtFile.create_from_matrices(
        bpm_names_x, matrix_x,
        bpm_names_y, matrix_y, date
    )
    _TbtAsciiWriter([new_tbt_file],
                    model_path,
                    output_path,
                    headers_dict).transform_tbt_to_ascii()


class TbtFile(object):
    def __init__(self, date, num_bunches, num_turns, bunch_id):
        self.date = date
        self.num_bunches = num_bunches
        self.num_turns = num_turns
        self.bunch_id = bunch_id
        self._samples_matrix = {HOR: {}, VER: {}}

    @staticmethod
    def create_from_matrices(bpm_names_x, matrix_x,
                             bpm_names_y, matrix_y, date):
        if date is None:
            date = datetime.now()
        new_tbt_file = TbtFile(bpm_names_x, bpm_names_y,
                               date, 1, len(matrix_x[0]), 0)
        new_tbt_file._samples_matrix[HOR] = pd.DataFrame(
            index=bpm_names_x,
            data=matrix_x,
            dtype=float,
        )
        new_tbt_file._samples_matrix[VER] = pd.DataFrame(
            index=bpm_names_y,
            data=matrix_y,
            dtype=float,
        )
        return new_tbt_file

    def get_x_samples(self, bpm_name):
        """
        Returns the horizontal samples captured by bpm_name.
        """
        return self._get(bpm_name, HOR)

    def get_y_samples(self, bpm_name):
        """
        Returns the vertical samples captured by bpm_name.
        """
        return self._get(bpm_name, VER)

    def write_to_ascii(self, model_path, output_path):
        """
        Writes this tbt file to output_path in ascii format.
        """
        _TbtAsciiWriter([self],
                        model_path,
                        output_path).transform_tbt_to_ascii()

    @property
    def samples_matrix_x(self):
        """
        Returns the vertical samples matrix, a Pandas DataFrame that uses BPM
        names as index.
        E.g.: a.samples_matrix_y.loc["BPM12", 3] -> Sample 3 of monitor BPM12.
        """
        return self._samples_matrix[HOR]

    @property
    def samples_matrix_y(self):
        """
        Returns the vertical samples matrix, a Pandas DataFrame that uses BPM
        names as index.
        E.g.: a.samples_matrix_x.loc["BPM12", 3] -> Sample 3 of monitor BPM12.
        """
        return self._samples_matrix[VER]

    def _get(self, bpm_name, plane):
        return self._samples_matrix[plane].loc[bpm_name]


# Private ###################

class _TbtReader(object):
    def __init__(self, file_path):
        sdds_file = sdds_reader.read_sdds_file(file_path)
        parameters = sdds_file.get_parameters()
        arrays = sdds_file.get_arrays()
        self._timestamp = parameters[TIMESTAMP_NAME].value
        self._timestamp_to_date()
        self._num_bunches = parameters[NUM_BUNCHES_NAME].value
        self._num_turns = parameters[NUM_TURNS_NAME].value
        self._all_samples = {
            HOR: arrays[ALL_HOR_POSITIONS_NAME].values,
            VER: arrays[ALL_VER_POSITIONS_NAME].values
        }
        self._bpm_names = arrays[BPM_NAMES_NAME].values
        self._num_monitors = len(self._bpm_names)
        self._bunch_id = {
            HOR: arrays[HOR_BUNCH_ID_NAME].values,
            VER: arrays[VER_BUNCH_ID_NAME].values
        }
        self._failed_bunch_id = {
            HOR: arrays[HOR_FAILED_BUNCH_ID_NAME].values,
            VER: arrays[VER_FAILED_BUNCH_ID_NAME].values
        }
        samples_matrix_shape = (self._num_monitors,
                                self._num_bunches,
                                self._num_turns)
        self._all_samples[HOR].shape = samples_matrix_shape
        self._all_samples[VER].shape = samples_matrix_shape
        self._tbt_files = []
        for index in range(self._num_bunches):
            self._tbt_files.append(
                TbtFile(self._date, self._num_bunches, self._num_turns,
                        self._bunch_id[HOR][index])  # TODO: does plane matter?
            )

    def _timestamp_to_date(self):
        # The timestamp is in nanoseconds
        self._date = datetime.fromtimestamp(self._timestamp / 1e9)

    def read_file(self):
        self._timestamp_to_date()
        for plane in (HOR, VER):
            self._read_bpms(plane)
        return self._tbt_files

    def _read_bpms(self, plane):
        num_bunches = self._num_bunches
        all_samples = self._all_samples[plane]
        for bunch_index in range(num_bunches):
            tbt_file = self._tbt_files[bunch_index]
            tbt_file._samples_matrix[plane] = pd.DataFrame(
                index=self._bpm_names,
                data=all_samples[:, bunch_index, :],
                dtype=float,
            )


class _TbtAsciiWriter(object):
    def __init__(self, tbt_files, model_path, output_path, headers_dict=None):
        if headers_dict is None:
            headers_dict = {}
        self._tbt_files = tbt_files
        self._model_path = model_path
        self._output_path = output_path
        self._headers_dict = headers_dict

    def transform_tbt_to_ascii(self):
        np.set_printoptions(precision=PRINT_PRECISION)
        model_data = self._load_model()
        for tbt_file in self._tbt_files:
            suffix = ""
            if len(self._tbt_files) > 1:
                suffix = "_" + str(tbt_file.bunch_id)
            with open(self._output_path + suffix, "w") as output_file:
                self._write_header(tbt_file, output_file, model_data)
                self._write_tbt_data(tbt_file, output_file, model_data)

    def _load_model(self):
        _append_beta_beat_to_path()
        from Utilities import tfs_pandas
        return tfs_pandas.read_tfs(self._model_path)

    def _write_header(self, tbt_file, output_file, model_data):
        output_file.write("#SDDSASCIIFORMAT v1\n")
        try:
            output_file.write("#Beam: " + model_data.SEQUENCE + "\n")
        except AttributeError:
            pass
        output_file.write("#Created: " + datetime.now().strftime(
            "%Y-%m-%d at %H:%M:%S By: Python SDDS converter"
        ) + "\n")
        output_file.write("#Number of turns: " +
                          str(tbt_file.num_turns) + "\n")
        output_file.write("#Number of horizontal monitors: " +
                          str(tbt_file.num_monitors_x) + "\n")
        output_file.write("#Number of vertical monitors: " +
                          str(tbt_file.num_monitors_y) + "\n")
        output_file.write("#Acquisition date: " + tbt_file.date.strftime(
            "%Y-%m-%d at %H:%M:%S"
        ) + "\n")
        for name, value in self._headers_dict.iteritems():
            output_file.write("#" + name + ": " + str(value) + "\n")

    def _write_tbt_data(self, tbt_file, output_file, model_data):
        for bpm_index in range(len(model_data.NAME)):
            bpm_name = model_data.NAME[bpm_index]
            bpm_s = str(np.fromstring(model_data.S[bpm_index])[0])
            try:
                bpm_samples_x = tbt_file.get_x_samples(bpm_name)
                bpm_samples_y = tbt_file.get_y_samples(bpm_name)
            except KeyError:
                LOGGER.debug(bpm_name + " not found in measurement file")
                continue
            output_str_x = " ".join((str(HOR), bpm_name, str(bpm_s), " ",
                                     " ".join(map(str, bpm_samples_x)), "\n"))
            output_str_y = " ".join((str(VER), bpm_name, str(bpm_s), " ",
                                     " ".join(map(str, bpm_samples_y)), "\n"))
            output_file.write(output_str_x)
            output_file.write(output_str_y)


def _append_beta_beat_to_path(self):
    parent_path = os.path.abspath(os.path.join(
        os.path.dirname(__file__), ".."
    ))
    if os.path.basename(parent_path) == "Beta-Beat.src":
        sys.path.append(parent_path)
    else:
        LOGGER.warn("Not in Beta-Beat.src, using lintrack metaclass")
        if "win" in sys.platform and sys.platform != "darwin":
            afs_root = "\\\\AFS"
        else:
            afs_root = "/afs"
        sys.path.append(
            os.path.join(afs_root, "cern.ch", "eng",
                         "sl", "lintrack", "Beta-Beat.src")
        )


if __name__ == "__main__":
    _, _file_path, _model_path, _output_path = sys.argv
    transform_tbt_to_ascii(_file_path, _model_path, _output_path)
