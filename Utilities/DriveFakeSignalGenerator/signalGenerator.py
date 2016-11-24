#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import random

def noise(amp):
	return (random.random()*2.-1.)*amp

def signal(turn, amps, tunes, phases, noise_amp):
	assert len(amps) == len(tunes) and len(amps) == len(phases)
	pi2 = np.pi*2
	
	signal = 0
	for i in range(len(tunes)):
		print np.cos(pi2*tunes[i]*turn + phases[i])
		print amps[i]
		signal += amps[i] * np.cos(pi2*tunes[i]*turn + phases[i])
	signal += noise(noise_amp/2.)
	return signal

def signalGenerator(bpms, amps, tunes, phases, turns=2000, noise_amp=.03):
	bpmDict = {}
	# initialize dict
	for bpm in bpms:
		bpmDict[bpm] = []
	
	for bpm in bpms:
		plane = bpms[bpm]['plane']
		for turn in range(turns):
			ph = [phases[plane][i]+bpms[bpm]['phaseOffset'] for i in range(len(phases[plane]))]
			print plane, len(tunes[plane])
			bpmDict[bpm].append(signal(turn, amps[plane], tunes[plane], ph, noise_amp))
		
	return bpmDict