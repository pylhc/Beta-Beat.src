from __future__ import print_function
import sys
import os
import sdds_reader
from datetime import datetime
import numpy as np


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
    tbt_files = _TbtReader(file_path).read_file()
    return tbt_files


def transform_tbt_to_ascii(file_path, model_path, output_path):
    tbt_files = _TbtReader(file_path).read_file()
    _TbtAsciiWriter(tbt_files,
                    model_path,
                    output_path).transform_tbt_to_ascii()


class TbtFile(object):
    def __init__(self, bpm_names, date,
                 num_bunches, num_turns, num_monitors, bunch_id):
        self.date = date
        self.num_bunches = num_bunches
        self.num_turns = num_turns
        self.num_monitors = num_monitors
        self.bunch_id = bunch_id
        self.bpm_names = bpm_names
        self._bpm_samples = {HOR: {}, VER: {}}
        self._samples_matrix = {HOR: {}, VER: {}}

    @property
    def bpm_samples_x(self):
        """
        Returns horizontal samples dictionary:
        bpm_name -> array_of_samples
        """
        return self._bpm_samples[HOR]

    @property
    def bpm_samples_y(self):
        """
        Returns vertical samples dictionary:
        bpm_name -> array_of_samples
        """
        return self._bpm_samples[VER]

    @property
    def samples_matrix_x(self):
        """
        Returns the horizontal samples matrix, guaranteeing the
        same monitor order than in self.bpm_names.
        E.g.: a.samples_matrix_x[5][3] -> Sample 3 of monitor 5
        """
        return self._samples_matrix[HOR]

    @property
    def samples_matrix_y(self):
        """
        Returns the vertical samples matrix, guaranteeing the
        same monitor order than in self.bpm_names.
        E.g.: a.samples_matrix_x[5][3] -> Sample 3 of monitor 5
        """
        return self._samples_matrix[VER]


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
                TbtFile(self._bpm_names, self._date, self._num_bunches,
                        self._num_turns, self._num_monitors,
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
        bpm_names = self._bpm_names
        num_bunches = self._num_bunches
        all_samples = self._all_samples[plane]
        for bunch_index in range(num_bunches):
            tbt_file = self._tbt_files[bunch_index]
            bpms_with_samples = tbt_file._bpm_samples[plane]
            for bpm_index in range(len(bpm_names)):
                bpm_name = bpm_names[bpm_index]
                bpms_with_samples[bpm_name] = all_samples[
                    bpm_index, bunch_index, :
                ]
            tbt_file._samples_matrix[plane] = all_samples[
                :, bunch_index, :
            ]
        return bpms_with_samples


class _TbtAsciiWriter(object):
    def __init__(self, tbt_files, model_path, output_path):
        self._tbt_files = tbt_files
        self._model_path = model_path
        self._output_path = output_path

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
        self._append_beta_beat_to_path()
        from Python_Classes4MAD import metaclass
        return metaclass.twiss(self._model_path)

    def _append_beta_beat_to_path(self):
        parent_path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..")
        )
        if os.path.basename(parent_path) == "Beta-Beat.src":
            sys.path.append(parent_path)
        else:
            print("Not in Beta-Beat.src, using lintrack metaclass")
            if "win" in sys.platform and not sys.platform == "darwin":
                afs_root = "\\\\AFS"
            else:
                afs_root = "/afs"
            sys.path.append(
                os.path.join(afs_root, "cern.ch", "eng",
                             "sl", "lintrack", "Beta-Beat.src")
            )

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
        output_file.write("#Number of monitors: " +
                          str(tbt_file.num_monitors) + "\n")
        output_file.write("#Acquisition date: " + tbt_file.date.strftime(
            "%Y-%m-%d at %H:%M:%S"
        ) + "\n")

    def _write_tbt_data(self, tbt_file, output_file, model_data):
        for bpm_index in range(len(model_data.NAME)):
            bpm_name = model_data.NAME[bpm_index]
            bpm_s = str(np.fromstring(model_data.S[bpm_index])[0])
            try:
                bpm_samples_x = tbt_file.bpm_samples_x[bpm_name]
                bpm_samples_y = tbt_file.bpm_samples_y[bpm_name]
            except KeyError:
                print(bpm_name + " not found in measurement file",
                      file=sys.stderr)
                continue
            output_str_x = " ".join((str(HOR), bpm_name, str(bpm_s), " ",
                                     " ".join(map(str, bpm_samples_x)), "\n"))
            output_str_y = " ".join((str(VER), bpm_name, str(bpm_s), " ",
                                     " ".join(map(str, bpm_samples_y)), "\n"))
            output_file.write(output_str_x)
            output_file.write(output_str_y)


if __name__ == "__main__":
    _, _file_path, _model_path, _output_path = sys.argv
    transform_tbt_to_ascii(_file_path, _model_path, _output_path)
