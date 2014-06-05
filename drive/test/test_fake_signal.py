import __init__  # @UnusedImport
import os
import sys
import subprocess
import unittest
import time
import datetime
import numpy as np
import Utilities.iotools
import random as rnd

CURRENT_PATH = os.path.dirname(__file__)
MAXIMUM_ABS_ERR = 1e-8
MAXIMUM_REL_ERR = 1e-8

# For drive seems to be really hard to find the phases and only gives this precision
MAXIMUM_PHASE_ABS_ERR = 1e-2
MAXIMUM_PHASE_REL_ERR = 1e-2

# This values must be the same as in the ./fake_signal/Drive.inp
NATURAL_TUNE = {"x": 0.27, "y": 0.32}
ISTUN = 0.005


class TestFakeData(unittest.TestCase):

    fake_data_file_path = os.path.join(CURRENT_PATH, "fake_signal", "test.sdds.cleaned")

    path_to_drive = os.path.join(CURRENT_PATH, "..", "Drive_God_lin")
    path_to_test = os.path.join(CURRENT_PATH, "fake_signal")
    path_to_driving_terms = os.path.join(CURRENT_PATH, "fake_signal", "DrivingTerms")
    path_to_sdds = os.path.join(CURRENT_PATH, "fake_signal", "test.sdds.cleaned")

    def setUp(self):
        self._generate_fake_data()
        self._run_drive()

    def tearDown(self):
        self._delete_drive_output()
        self._delete_fake_data()

    def testFakeData(self):
        self._compare_drive_output_with_fake_input()

    def _generate_fake_data(self):
        self.qx = 0.28
        self.qy = 0.31

        self.bpms_data = {
            'BPMX.1.BTEST': {'pos': 0.0, 'plane': 'x', 'phaseOffset': 0.0},
            'BPMY.1.BTEST': {'pos': 0.0, 'plane': 'y', 'phaseOffset': 0.0},
            'BPMX.2.BTEST': {'pos': 1.0, 'plane': 'x', 'phaseOffset': 0.0},
            'BPMY.2.BTEST': {'pos': 1.0, 'plane': 'y', 'phaseOffset': 0.0},
            'BPMX.3.BTEST': {'pos': 2.0, 'plane': 'x', 'phaseOffset': 0.0},
            'BPMY.3.BTEST': {'pos': 2.0, 'plane': 'y', 'phaseOffset': 0.0},
            'BPMX.4.BTEST': {'pos': 3.0, 'plane': 'x', 'phaseOffset': 0.0},
            'BPMY.4.BTEST': {'pos': 3.0, 'plane': 'y', 'phaseOffset': 0.0},
        }

        #  Resonances (random examples, may be adjusted)
        self.resonances_dict = {"x": {}, "y": {}}

        ResonanceLineData.add(self, "x", 1, 0, 100.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 0, 1, 10.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 0, 2, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 2, 0, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 1, 1, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 1, -2, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 1, 2, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", -2, 1, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 3, 0, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", 2, -2, 0.01 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "x", -1, 3, 0.01 * rnd.random(), np.pi * (2 * rnd.random() - 1))

        ResonanceLineData.add(self, "y", 0, 1, 100.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 1, 0, 10.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 0, 2, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 2, 0, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 1, 1, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 1, -1, 1.0 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 1, -2, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 2, 1, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", 0, 3, 0.1 * rnd.random(), np.pi * (2 * rnd.random() - 1))
        ResonanceLineData.add(self, "y", -1, 3, 0.01 * rnd.random(), np.pi * (2 * rnd.random() - 1))

        self.turns = 2000
        self.kick_turn = 100

        bpm_to_signal_dict = self._get_bpm_to_signal_dict(
                                                        self.bpms_data,
                                                        self.resonances_dict,
                                                        self.turns,
                                                        self.kick_turn
                                                        )

        write_result = self._write_fake_data_to_file(bpm_to_signal_dict, self.bpms_data)
        if write_result:
            print('Wrote data file to: %s' % self.fake_data_file_path)
        else:
            print('Something went wrong while writing the data to file!')

    def _get_bpm_to_signal_dict(self, bpms_data, resonances_dict, turns=2000, kick_turn=100):
        bpm_dict = {}

        for bpm_name in bpms_data:
            bpm_dict[bpm_name] = []

        for bpm_name in bpms_data:
            plane = bpms_data[bpm_name]['plane']
            max_signal = -500000.0
            min_signal = 500000.0
            co = 0.0
            corms = 0.0

            for turn in range(turns):
                bpm_signal = self._compute_signal_of_bpm(turn, plane, bpms_data[bpm_name]["phaseOffset"])
                bpm_dict[bpm_name].append(bpm_signal)

                co += bpm_signal
                corms += bpm_signal ** 2
                if(bpm_signal > max_signal):
                    max_signal = bpm_signal
                if(bpm_signal < min_signal):
                    min_signal = bpm_signal

            co = co / turns
            corms = np.ma.sqrt(corms / turns - co ** 2)

            bpms_data[bpm_name]['pk2pk'] = max_signal - min_signal
            bpms_data[bpm_name]['co'] = co
            bpms_data[bpm_name]['corms'] = corms
        return bpm_dict

    def _compute_signal_of_bpm(self, turn, plane, bpm_offset):
        pi2 = np.pi * 2

        signal = 0
        for resonance_data in self.resonances_dict[plane].values():
            signal += resonance_data.amplitude * np.cos(pi2 * resonance_data.tune * turn +
                                                        resonance_data.phase +
                                                        bpm_offset)
        return signal

    def _write_fake_data_to_file(self, bpmData, bpms):
        header = ("#SDDSASCIIFORMAT v1\n"
                    "#Beam: Test\n"
                    "#Created: %s By: Drive Fake Signal Generator\n"
                    "#bunchid :0\n"
                    "#number of turns :%d\n"
                    "#number of monitors :%d\n"
                    "#NTURNS calculated: %d\n") % (
                    datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d#%H-%M-%S'),
                    len(bpmData[bpmData.keys()[0]]),
                    len(bpms),
                    len(bpmData[bpmData.keys()[0]])
        )

        with open(self.fake_data_file_path, 'w') as out:
            out.write(header)
            for bpm in bpmData:
                if   bpms[bpm]['plane'] == 'x':
                    plane = 0
                elif bpms[bpm]['plane'] == 'y':
                    plane = 1
                else:
                    raise ValueError('Unknown plane: %s' % bpms[bpm]['plane'])

                out.write('%d %s %f %s\n' % (plane, bpm, bpms[bpm]['pos'], ' '.join(format(x.real, '.5f') for x in bpmData[bpm])))
        if os.path.isfile(self.fake_data_file_path):
            return True
        return False

    def _run_drive(self):
        file_driving_terms = open(self.path_to_driving_terms, "w")
        print >> file_driving_terms, self.path_to_sdds, "1", "2000"
        file_driving_terms.close()

        call_command = os.path.abspath(self.path_to_drive) + " " + os.path.abspath(self.path_to_test)
        process = subprocess.Popen(call_command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               shell=True)

        (out_stream, err_stream) = process.communicate()

        errcode = process.returncode

        print "Printing output:-------------------------"
        print out_stream
        if errcode != 0:
            print "Error running command:", call_command
            print >> sys.stderr, "Printing error output:-------------------"
        print >> sys.stderr, err_stream

    def _delete_drive_output(self):
        Utilities.iotools.delete_item(os.path.join(self.path_to_test, "BPM"))
        Utilities.iotools.delete_item(os.path.join(self.path_to_test, "sussix_v4.inp"))
        Utilities.iotools.delete_item(os.path.join(self.path_to_test, "test.sdds.cleaned_linx"))
        Utilities.iotools.delete_item(os.path.join(self.path_to_test, "test.sdds.cleaned_liny"))
        Utilities.iotools.delete_item(os.path.join(self.path_to_test, "test.sdds.cleaned_rejectedBpms_x"))
        Utilities.iotools.delete_item(os.path.join(self.path_to_test, "test.sdds.cleaned_rejectedBpms_y"))

    def _delete_fake_data(self):
        Utilities.iotools.delete_item(self.fake_data_file_path)
        Utilities.iotools.delete_item(self.path_to_driving_terms)

    def _compare_drive_output_with_fake_input(self):
        drive_output_content_dicts_linx = self._get_drive_data_dicts_for_file("test.sdds.cleaned_linx")
        drive_output_content_dicts_liny = self._get_drive_data_dicts_for_file("test.sdds.cleaned_liny")
        for values_dict in drive_output_content_dicts_linx:
            self._test_line_values(values_dict)
        for values_dict in drive_output_content_dicts_liny:
            self._test_line_values(values_dict)

    def _get_drive_data_dicts_for_file(self, file_name):
        dict_list = []
        header_list = []
        types_list = []
        with open(os.path.join(self.path_to_test, file_name), "r") as drive_file_data:
            drive_file_lines = drive_file_data.readlines()
            for drive_line in drive_file_lines:
                if drive_line.startswith("@"):
                    pass
                elif drive_line.startswith("*"):
                    header_list = self._clean_and_split_line(drive_line)
                    del header_list[0]
                elif drive_line.startswith("$"):
                    types_list = self._clean_and_split_line(drive_line)
                    del types_list[0]
                else:
                    line_dict = {}
                    line_values = self._clean_and_split_line(drive_line)
                    self.assertTrue(len(header_list) == len(types_list) == len(line_values),
                                    "Incorrect number of columns in file: " + file_name)
                    for col_index in range(len(line_values)):
                        if types_list[col_index] == "%s":
                            line_dict[header_list[col_index]] = line_values[col_index].strip('"')
                        elif types_list[col_index] == "%le":
                            line_dict[header_list[col_index]] = float(line_values[col_index])
                    dict_list.append(line_dict)
        return dict_list

    def _clean_and_split_line(self, line):
        return line.strip().split()

    def _test_line_values(self, values_dict):
        bpm_name = values_dict["NAME"]
        self.assertTrue(bpm_name in self.bpms_data, bpm_name + " not in the fake bpm list")
        bpm_data = self.bpms_data[bpm_name]
        for (value_name, value) in values_dict.iteritems():
            if value_name == "S":
                self.assertTrue(value == bpm_data["pos"],
                                self._get_diff_message(value_name, bpm_data["pos"], value))
            elif value_name == "BINDEX":
                pass  # Irrelevant, just the BPM number
            elif value_name == "SLABEL":
                pass  # Irrelevant, always 1
            elif value_name == "TUNEX":
                self.assertTrue(self._similar(value, self.qx),
                                self._get_diff_message(value_name, self.qx, value))
            elif value_name == "TUNEY":
                self.assertTrue(self._similar(value, self.qy),
                                self._get_diff_message(value_name, self.qy, value))
            elif value_name == "NOISE":
                pass  # Irrelevant, there is not noise in the fake data
            elif value_name == "PK2PK":
                self.assertTrue(self._similar(value, bpm_data["pk2pk"]),
                                self._get_diff_message(value_name, bpm_data["pk2pk"], value))
            elif value_name == "CO":
                self.assertTrue(self._similar(value, bpm_data["co"]),
                                self._get_diff_message(value_name, bpm_data["co"], value))
            elif value_name == "CORMS":
                self.assertTrue(self._similar(value, bpm_data["corms"]),
                                self._get_diff_message(value_name, bpm_data["corms"], value))
            elif value_name == "AMPX":
                self.assertTrue(self._similar(value, self.resonances_dict["x"][1, 0].amplitude / 2),
                                self._get_diff_message(value_name,
                                                       self.resonances_dict["x"][1, 0].amplitude / 2, value))
            elif value_name == "AMPY":
                self.assertTrue(self._similar(value, self.resonances_dict["y"][0, 1].amplitude / 2),
                                self._get_diff_message(value_name,
                                                       self.resonances_dict["y"][0, 1].amplitude / 2, value))
            elif value_name == "MUX":
                self.assertTrue(self._similar(abs(value), abs(self.resonances_dict["x"][1, 0].phase / (2 * np.pi))),
                                self._get_diff_message(value_name,
                                                        abs(self.resonances_dict["x"][1, 0].phase / (2 * np.pi)),
                                                        abs(value)))
            elif value_name == "MUY":
                self.assertTrue(self._similar(abs(value), abs(self.resonances_dict["y"][0, 1].phase / (2 * np.pi))),
                                self._get_diff_message(value_name,
                                                        abs(self.resonances_dict["y"][0, 1].phase / (2 * np.pi)),
                                                        abs(value)))
            elif value_name.startswith("AMP"):
                i, j = self._extract_index(value_name.replace("AMP", ""))
                if "AMPX" in values_dict:
                    fake_amp = values_dict["AMPX"]
                    plane = "x"
                elif "AMPY" in values_dict:
                    fake_amp = values_dict["AMPY"]
                    plane = "y"
                real_amp = 2 * fake_amp * value
                if (i, j) in self.resonances_dict[plane]:
                    drive_amp = self.resonances_dict[plane][i, j].amplitude
                elif (-i, -j) in self.resonances_dict[plane]:
                    drive_amp = self.resonances_dict[plane][-i, -j].amplitude
                else:
                    drive_amp = 0.0
                self.assertTrue(self._similar(real_amp, drive_amp),
                                self._get_diff_message(value_name, drive_amp, real_amp))
            elif value_name.startswith("PHASE"):
                i, j = self._extract_index(value_name.replace("PHASE", ""))
                self._test_resonance_phase(i, j, values_dict, value_name, value)
            elif value_name.startswith("NATTUNE"):
                plane = value_name.replace("NATTUNE", "").lower()
                natural_tune = self._get_natural_tune(bpm_name, plane)
                self.assertTrue(self._similar(natural_tune, value),
                                self._get_diff_message(value_name, natural_tune, value))
            elif value_name.startswith("NATAMP"):
                plane = value_name.replace("NATAMP", "").lower()
                natural_amplitude = self._get_natural_amplitude(bpm_name, plane)
                self.assertTrue(self._similar(natural_amplitude, value),
                                self._get_diff_message(value_name, natural_amplitude, value))

    def _get_diff_message(self, col_name, fake_value, drive_value):
        return ("Difference found in column " + col_name +
                " fake: " + str(fake_value) +
                " drive: " + str(drive_value))

    def _similar(self, value1, value2, abs_err=MAXIMUM_ABS_ERR, relative_err=MAXIMUM_PHASE_REL_ERR):
        min_val = min(abs(value1), abs(value2))
        if min_val == 0:
            min_val = 1
        return (abs(value1 - value2) < MAXIMUM_PHASE_ABS_ERR or
                (abs(value1 - value2) / min_val) < MAXIMUM_PHASE_REL_ERR)

    def _extract_index(self, raw_string):
        sign = 1
        is_set_i = False
        i = 0
        j = 0
        for char in raw_string:
            if char == "_":
                sign = -1
            elif not is_set_i:
                i = sign * int(char)
                sign = 1
                is_set_i = True
            else:
                j = sign * int(char)
        return i, j

    def _test_resonance_phase(self, i, j, values_dict, value_name, value):
        if "MUX" in values_dict:
            plane = "x"
        elif "MUY" in values_dict:
            plane = "y"
        amp_string = "AMP"
        if i >= 0:
            amp_string += str(i)
        else:
            amp_string += "_" + str(abs(i))
        if j >= 0:
            amp_string += str(j)
        else:
            amp_string += "_" + str(abs(j))
        drive_amp = values_dict[amp_string]
        #  Sometimes drive outputs a very small amplitude instead of 0.0, in that cases it also outputs a phase that
        #  in the general case won't be similar to 0.0, so we have to ignore it.
        if not self._similar(drive_amp, 0.0):
            if (i, j) in self.resonances_dict[plane]:
                fake_phase = abs(self.resonances_dict[plane][i, j].phase / (2 * np.pi))
            elif (-i, -j) in self.resonances_dict[plane]:
                fake_phase = abs(self.resonances_dict[plane][-i, -j].phase / (2 * np.pi))
            else:
                fake_phase = 0.0
            self.assertTrue(self._similar(fake_phase, abs(value), MAXIMUM_PHASE_ABS_ERR, MAXIMUM_PHASE_REL_ERR),
                            self._get_diff_message(value_name, abs(value), fake_phase))

    def _get_natural_tune(self, bpm_name, plane):
        maximum_params = self._get_maximum_amplitude_bpm_params(bpm_name, plane)
        return float(maximum_params[0])

    def _get_natural_amplitude(self, bpm_name, plane):
        maximum_params = self._get_maximum_amplitude_bpm_params(bpm_name, plane)
        return float(maximum_params[1])

    def _get_maximum_amplitude_bpm_params(self, bpm_name, plane):
        max_amplitude = -500000.0
        result_params = []
        path_to_bpm_spectrum = os.path.join(self.path_to_test, "BPM", bpm_name + "." + plane)
        with open(path_to_bpm_spectrum) as spectrum_file:
            for line in spectrum_file.readlines():
                if not (line.startswith("*") or line.startswith("$")):
                    bpm_params = self._clean_and_split_line(line)
                    amplitude = float(bpm_params[1])
                    tune = float(bpm_params[0])
                    if (amplitude > max_amplitude and
                        NATURAL_TUNE[plane] - ISTUN < tune and
                        tune < NATURAL_TUNE[plane] + ISTUN):
                            max_amplitude = amplitude
                            result_params = [tune, amplitude]
        if len(result_params) == 0:
            result_params[0] = -100.0
            result_params[1] = -100.0
        return result_params


class ResonanceLineData(object):

        def __init__(self, father, tune_i, tune_j, amplitude, phase):
            self._father = father
            self.amplitude = amplitude
            self.tune = self._get_tune(tune_i, tune_j)
            self.phase = phase

        def _get_tune(self, tune_i, tune_j):
            tune = tune_i * self._father.qx + tune_j * self._father.qy
            if not -0.5 < tune < 0.5:
                tune = 1 - tune
            return tune

        @staticmethod
        def add(father, plane, tune_i, tune_j, amplitude, phase):
            father.resonances_dict[plane][tune_i, tune_j] = ResonanceLineData(father, tune_i, tune_j, amplitude, phase)


if __name__ == '__main__':
    unittest.main()
