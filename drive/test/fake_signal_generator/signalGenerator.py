#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np

def signal(turn, amps, tunes, phases):
	assert len(amps) == len(tunes) and len(amps) == len(phases)
	pi2 = np.pi*2
	
	signal = 0
	for i in range(len(amps)):
		signal += amps[i] * np.cos(pi2*tunes[i]*turn + phases[i])
	return signal

def signalGenerator(bpms, amps, tunes, phases, turns=2000, kickTurn=100):
	bpmDict = {}
	# initialize dict
	for bpm in bpms:
		bpmDict[bpm] = []
	
	for bpm in bpms:
		plane = bpms[bpm]['plane']
		for turn in range(turns):
			ph = [phases[plane][i]+bpms[bpm]['phaseOffset'] for i in range(len(phases[plane]))]
			bpmDict[bpm].append(signal(turn, amps[plane], tunes[plane], ph))
		
	return bpmDict