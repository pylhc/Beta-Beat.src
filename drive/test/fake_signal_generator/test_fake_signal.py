import os
import sys
import subprocess
import unittest
import time
import datetime
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

CURRENT_PATH = os.path.dirname(__file__)

class test_fake_data(unittest.TestCase):
	
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
		self._delete_drive_output()
		self._delete_fake_data()
		
	def testFakeData(self):
		pass
		
	def _generate_fake_data(self):
		qx = 0.28
		qy = 0.31
	
		# BPMs
		bpms_data = {
			'BPMX.1.BTEST':{'pos':0.0, 'plane':'x', 'phaseOffset':0.},
			'BPMY.1.BTEST':{'pos':0.0, 'plane':'y', 'phaseOffset':0.},
			'BPMX.2.BTEST':{'pos':1.0, 'plane':'x', 'phaseOffset':.5},
			'BPMY.2.BTEST':{'pos':1.0, 'plane':'y', 'phaseOffset':.5},
			'BPMX.3.BTEST':{'pos':2.0, 'plane':'x', 'phaseOffset':0.},
			'BPMY.3.BTEST':{'pos':2.0, 'plane':'y', 'phaseOffset':0.},
			'BPMX.4.BTEST':{'pos':3.0, 'plane':'x', 'phaseOffset':-.5},
			'BPMY.4.BTEST':{'pos':3.0, 'plane':'y', 'phaseOffset':-.5},
		}
		
		# Resonances (random examples, may be adjusted)
		amps   = {
			'x':[100., 10.,        1.,        .1,     1.],
			'y':[100., 10.,        1.,        .1]
		}
		tunes  = {
			'x':[  qx,  qy, 1-(qx+qy),      2*qx, qx-.01],
			'y':[  qy,  qx, 1-(qx+qy), 1-(qx-qy)]
		}
		phases = {
			'x':[  0.,  0.,        0.,        0.,     0.],
			'y':[  0.,  0.,        0.,        0.]
		}
	
		turns    = 2000
		kickTurn =  100
	
		bpm_to_signal_dict = self._get_bpm_to_signal_dict(bpms_data, amps, tunes, phases, turns=turns, kick_turn=kickTurn)
	
		if self.plot_fake_data:
			self._do_plot_fake_data(bpm_to_signal_dict, bpms_data, fft=True)
	
		wroteToFile = self._write_fake_data_to_file(bpm_to_signal_dict, bpms_data)
		if wroteToFile:
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
	
	def _get_bpm_to_signal_dict(self, bpms_data, amps, tunes, phases, turns=2000, kick_turn=100):
		bpm_dict = {}
		# initialize dict
		for bpm_name in bpms_data:
			bpm_dict[bpm_name] = []
		
		for bpm_name in bpms_data:
			plane = bpms_data[bpm_name]['plane']
			for turn in range(turns):
				ph = [phases[plane][i] + bpms_data[bpm_name]['phaseOffset'] for i in range(len(phases[plane]))]
				bpm_dict[bpm_name].append(self._compute_signal_of_bpm(turn, amps[plane], tunes[plane], ph))
			
		return bpm_dict
	
	def _compute_signal_of_bpm(self, turn, amps, tunes, phases):
		assert len(amps) == len(tunes) and len(amps) == len(phases)
		pi2 = np.pi*2
		
		signal = 0
		for i in range(len(amps)):
			signal += amps[i] * np.cos(pi2*tunes[i]*turn + phases[i])
		return signal
	
	def _write_fake_data_to_file(self, bpmData, bpms):
		header = '''#SDDSASCIIFORMAT v1
					#Beam: Test
					#Created: %s By: Drive Fake Signal Generator
					#bunchid :0
					#number of turns :%d
					#number of monitors :%d
					#NTURNS calculated: %d
				''' % (
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
	
		if 0 != errcode:
			print "Error running command:", call_command
			print >> sys.stderr, "Printing error output:-------------------"
			print >> sys.stderr, err_stream
		print "Printing output:-------------------------"
		print out_stream


if __name__=='__main__':
	unittest.main()
