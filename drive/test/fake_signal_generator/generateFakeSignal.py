# for creating and writing the fake BPM signals
from signalGenerator import signalGenerator
from fileWriter import writeBPMfile
import os
import sys
import subprocess

CURRENT_PATH = os.path.dirname(__file__)

def makePlot(bpmSignalDict, bpms, fft=True):
	import matplotlib.pyplot as plt
	from matplotlib.ticker import MultipleLocator
	import numpy as np
	
	ax = plt.subplots()[1]
	if fft:
		x = np.fft.fftfreq(2000, d=1.)
		y = []
		for bpm in bpms:
			y.append(abs(np.fft.fft(bpmSignalDict[bpm])))
			plt.plot(x, y[-1], alpha=.7, label=bpm)
		plt.xlim(-.5,.5)
		#plt.ylim(1.)
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
	
def main():
	outFile = os.path.join(CURRENT_PATH, "fake_signal", "test.sdds.cleaned")

	qx = 0.28
	qy = 0.31

	# BPMs
	bpms = {
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

	# create signals
	bpmSignalDict = signalGenerator(bpms, amps, tunes, phases, turns=turns, kickTurn=kickTurn)

	# plot
	plot = False
	if plot:
		makePlot(bpmSignalDict, bpms, fft=True)

	# write data to file
	wroteToFile = writeBPMfile(outFile, bpmSignalDict, bpms)
	if wroteToFile:
		print('Wrote data file to: %s' % outFile)
	else:
		print('Something went wrong while writing the data to file!')


def test():
	path_to_drive = os.path.join(CURRENT_PATH, "..", "..", "Drive_God_lin")
	path_to_test = os.path.join(CURRENT_PATH, "fake_signal")
	path_to_driving_terms = os.path.join(CURRENT_PATH, "fake_signal", "DrivingTerms")
	path_to_sdds = os.path.join(CURRENT_PATH, "fake_signal", "test.sdds.cleaned")
	
	file_driving_terms = open(path_to_driving_terms, "w")
	print >> file_driving_terms, path_to_sdds, "1", "2000"
	file_driving_terms.close()
	
	call_command = os.path.abspath(path_to_drive) + " " + os.path.abspath(path_to_test)
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
	main()
	test()
