import numpy as np
import sys
import os
from os.path import abspath, join, dirname
new_path = abspath(join(dirname(abspath(__file__)), os.pardir, os.pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

from utils.tfs_pandas import read_tfs


C = 299792458.0 
PI2 = np.pi * 2
#First all per BPM stuff is applied
def generate(infile, nturns=6600, plane ='X',       #single plane with no coupling
            max_dpoverp = 0.0002,                   #amplitude of dp over p modulation
            chroma = 5,                            #Chroma not taken from the file, since it is generally not matched, and may have crazy values
            rf_modulation_freqency = 5.0,                        #Hz frequency of RF modulation, i. e. change of beam energy
            deltaQ = -0.01,                         # difference between AC-Dipole driven tune and the natural tune, if 0.0 -> kicker case
            strength = 1.02,                        # of ac-dipole or kicker
            sigmaQ = 0.0001,                        # width of tune ... sigma
            damp_time=1500.0,                        # damping time of kicker induced betatron oscillations 
            bpm_noise=0.1                           #BPM noise in mm
            ):  
    twiss = read_tfs(infile)
    if plane == 'X':
        Q = np.remainder(twiss.headers['Q1'],1)
    else: 
        Q = np.remainder(twiss.headers['Q2'],1)
    nBPMs = len(twiss.index.values)
    dims = (nturns, nBPMs)
    phase0z = np.random.rand()            #initial phase of synchrotron oscilation
    phase0 = np.random.rand()             #initial phase of betatron oscilation at the beginning of the lattice 
    Qz = rf_modulation_freqency * twiss.headers['LENGTH'] / C     #synchrotron tune
    dpp = max_dpoverp * np.cos(PI2 * Qz * np.arange(nturns, dtype=float) + PI2 * phase0z) #dpoverp for all turns
    tune = dpp * chroma + Q                     #betatron tune for all turns
    dQ = 1/np.abs(dpp * chroma - deltaQ + _get_noise(PI2 * sigmaQ, d1=nturns)) 
        
    
    signal = 1000 * np.outer(twiss.loc[:,'D' + plane].values , dpp) #synchrotron motion nBPMs x nturns
    nat_tune = np.transpose(0.01 * np.sqrt(twiss.loc[:,'BET'+ plane].values) * np.cos(_get_tbt_phases(tune, twiss, plane, phase0, nturns, nBPMs, sigmaQ)))    #
    
    if not deltaQ:  # kicker case     
        nat_tune= nat_tune * 10 * np.exp(-1.0/damp_time * np.arange(nturns, dtype=float)) 
    else:       
        dQ = 1/np.abs(dpp * chroma - deltaQ + _get_noise(PI2 * sigmaQ, d1=nturns))
        dQ[np.where(dQ>10000.0)] = 10000.0        
        dr_tune = dQ * np.transpose(0.0006 * strength * np.sqrt(twiss.loc[:,'BET'+ plane].values) * np.cos(_get_tbt_phases_ac(Q + deltaQ, twiss, plane, phase0, nturns, nBPMs, sigmaQ)))
        signal += dr_tune
    signal = signal + nat_tune + _get_noise(bpm_noise, d1=nBPMs, d2=nturns)
    #plot_bpm_tbt('', signal, 40)
    return signal

def plot_bpm_tbt(path, signal, num):
    import matplotlib.pyplot as plt
    f,ax = plt.subplots()
    ax.plot(signal[num])
    f.savefig(os.path.join(path, str(num) +'.png'))

def save(file_name,tbt_samples):
    np.save(file_name,tbt_samples)

def load(file_name):
    return np.load(file_name)

#returns tbt phases (in rads) of shape nturns x nBPMs with phase advance between BPMs as in twiss assuming tune spread sigmaQ
def _get_tbt_phases(tune, twiss, plane, phase0, nturns, nBPMs, sigmaQ): 
    return np.cumsum(_get_noise(PI2 * sigmaQ, d1=nturns, d2=nBPMs) + _get_tbt_tunes(tune, nBPMs), axis=0) + _get_betetron_phases(twiss, plane, phase0)

def _get_tbt_phases_ac(q, twiss, plane, phase0, nturns, nBPMs, sigmaQ): 
    return np.cumsum(_get_noise(PI2 * sigmaQ, d1=nturns, d2=nBPMs) + _get_tbt_tunes_ac((nturns, nBPMs), q), axis=0) + _get_betetron_phases(twiss, plane, phase0)

def _get_tbt_tunes(tune,nBPM):
    return PI2 * np.outer(tune,np.ones(nBPM))

def _get_tbt_tunes_ac(dims,q):
    return np.full(dims,PI2 * q)

def _get_betetron_phases(twiss, plane, ph0): # 
    return PI2 * (twiss.loc[:,'MU' + plane].values + ph0)

def _get_noise(amp, d1=None, d2=None):
    if d1 == None:
        return amp * np.random.randn()
    elif d2 == None:
        return amp * np.random.randn(d1)
    else:
        return amp * np.random.randn(d1, d2)

if __name__ == '__main__':
    save('test.npy',generate('twiss.dat'))
