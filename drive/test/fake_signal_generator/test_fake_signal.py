import os
import sys
import subprocess
import unittest
import time
import datetime
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import Utilities.iotools

CURRENT_PATH = os.path.dirname(__file__)
MAXIMUM_ABS_ERR = 1e-3
MAXIMUM_REL_ERR = 1e-3

class TestFakeData(unittest.TestCase):

    plot_fake_data = False

    fake_data_file_path = os.path.join(CURRENT_PATH, "fake_signal", "test.sdds.cleaned")

    path_to_drive = os.path.join(CURRENT_PATH, "..", "..", "Drive_God_lin")
    path_to_test = os.path.join(CURRENT_PATH, "fake_signal")
    path_to_driving_terms = os.path.join(CURRENT_PATH, "fake_signal", "DrivingTerms")
    path_to_sdds = os.path.join(CURRENT_PATH, "fake_signal", "test.sdds.cleaned")
            
    def setUp(self):
        self._generate_fake_data()
        self._run_drive()

    def tearDown(self):
        # self._delete_drive_output()
        self._delete_fake_data()

    def testFakeData(self):
        self._compare_drive_output_with_fake_input()

    def _generate_fake_data(self):
        self.qx = 0.28
        self.qy = 0.31

        # BPMs
        self.bpms_data = {
            'BPMX.1.BTEST':{'pos':0.0, 'plane':'x', 'phaseOffset':0.0},
            'BPMY.1.BTEST':{'pos':0.0, 'plane':'y', 'phaseOffset':0.0},
            'BPMX.2.BTEST':{'pos':1.0, 'plane':'x', 'phaseOffset':0.5},
            'BPMY.2.BTEST':{'pos':1.0, 'plane':'y', 'phaseOffset':0.5},
            'BPMX.3.BTEST':{'pos':2.0, 'plane':'x', 'phaseOffset':0.0},
            'BPMY.3.BTEST':{'pos':2.0, 'plane':'y', 'phaseOffset':0.0},
            'BPMX.4.BTEST':{'pos':3.0, 'plane':'x', 'phaseOffset':-0.5},
            'BPMY.4.BTEST':{'pos':3.0, 'plane':'y', 'phaseOffset':-0.5},
        }
        
        # Resonances (random examples, may be adjusted)
        self.resonances_dict = {"x": {}, "y": {}}
        
        ResonanceLineData.add(self, "x", 1, 0, 100.0, 0.0)
        ResonanceLineData.add(self, "x", 0, 1, 10.0, 0.0)
        ResonanceLineData.add(self, "x", 1, 1, 1.0, 0.0)
        ResonanceLineData.add(self, "x", 2, 0, 0.1, 0.0)
        
        ResonanceLineData.add(self, "y", 0, 1, 100.0, 0.0)
        ResonanceLineData.add(self, "y", 1, 0, 10.0, 0.0)
        ResonanceLineData.add(self, "y", 1, 1, 1.0, 0.0)
        ResonanceLineData.add(self, "y", 1, -1, 0.1, 0.0)
    
        self.turns = 2000
        self.kick_turn = 100
    
        bpm_to_signal_dict = self._get_bpm_to_signal_dict(
                                                        self.bpms_data,
                                                        self.resonances_dict,
                                                        self.turns,
                                                        self.kick_turn
                                                        )
    
        if self.plot_fake_data:
            self._do_plot_fake_data(bpm_to_signal_dict, self.bpms_data, fft=True)
    
        write_result = self._write_fake_data_to_file(bpm_to_signal_dict, self.bpms_data)
        if write_result:
            print('Wrote data file to: %s' % self.fake_data_file_path)
        else:
            print('Something went wrong while writing the data to file!')
    
    def _do_plot_fake_data(self, bpmSignalDict, bpms, fft=True):
        ax = plt.subplot()[1]
        if fft:
            x = np.fft.fftfreq(2000, d=1.)
            y = []
            for bpm in bpms:
                y.append(abs(np.fft.fft(bpmSignalDict[bpm])))
                plt.plot(x, y[-1], alpha=.7, label=bpm)
            plt.xlim(-.5,.5)
            plt.xticks()
            plt.semilogy()
            majorLocator = MultipleLocator(0.1)
            minorLocator = MultipleLocator(0.05)
            ax.xaxis.set_major_locator(majorLocator)
            ax.xaxis.set_minor_locator(minorLocator)
            plt.xlabel('Tune [1/Turns]')
            plt.ylabel('Amplitude [a.u.]')
        else:
            for bpm in bpms:
                plt.plot(bpmSignalDict[bpm], alpha=.7, label=bpm)
            plt.xlabel('Turns [1]')
            plt.ylabel('Amplitude [a.u.]')
    
        plt.subplots_adjust(left=0.04, right=0.97, top=0.97, bottom=0.04)
        plt.grid(which='both')
        plt.legend()
        plt.show()

    def _get_bpm_to_signal_dict(self, bpms_data, resonances_dict, turns=2000, kick_turn=100):
        bpm_dict = {}
        # initialize dict
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
        pi2 = np.pi*2
        
        signal = 0
        for resonance_data in self.resonances_dict[plane].values():
            signal += resonance_data.amplitude * np.cos(pi2 * resonance_data.tune * turn + resonance_data.phase + bpm_offset)
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
                if   bpms[bpm]['plane'] == 'x': plane = 0
                elif bpms[bpm]['plane'] == 'y': plane = 1
                else: raise ValueError('Unknown plane: %s' % bpms[bpm]['plane'])
                
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
                pass  # Irrelevant, allways 1
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
                self.assertTrue(self._similar(value,self.resonances_dict["y"][0, 1].amplitude / 2), 
                                self._get_diff_message(value_name, 
                                                       self.resonances_dict["y"][0, 1].amplitude / 2, value))
            elif value_name == "MUX":
                pass  # TODO: How do we get this value?
            elif value_name == "MUY":
                pass  # TODO: How do we get this value?
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
                pass  # TODO: How do we get this value?
            elif value_name.startswith("NATTUNE"):
                pass  # TODO: How do we get this value?
            elif value_name.startswith("NATAMP"):
                pass  # TODO: How do we get this value?
            
    def _get_diff_message(self, col_name, fake_value, drive_value):
        return ("Difference found in column " + col_name + 
                " fake: " + str(fake_value) + 
                " drive: " + str(drive_value))
            
    def _similar(self, value1, value2):
        min_val = min(abs(value1), abs(value2))
        if min_val == 0:
            min_val = 1
        return (abs(value1 - value2) < MAXIMUM_ABS_ERR or 
                (abs(value1 - value2) / min_val) < MAXIMUM_REL_ERR)
        
    def _extract_index(self, raw_string):
        sign = 1;
        is_set_i = False;
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
