import os
import sys
import scipy.io
sys.path.append("/afs/cern.ch/work/g/garciata/public/beta_beat_src_good/Beta-Beat.src/")
from drive import drive_runner
from Utilities import tfs_pandas
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import *
from optparse import OptionParser
Darkred1 = '#8B0000'
Darkred2 = '#FF0000'
Brown1 = '#FAA460'
Brown2 = '#8B4513'
SkyBlue1 = '#87CEEB'
SkyBlue2 = '#00BFFF'
Gray1 = '#808080'
Gray2 = '#696969'
Gold1 = '#FFD700'
Gold2 = '#EEE8AA'
Violet1 = '#8A2BE2'
Violet2 = '#4B0082'
Olive1 = '#6B8E23'
Olive2 = '#556B2F'
Orange1 = '#FFA500'
Orange2 = '#FF8C00'
COLOR_LIST_LIGHT = [Darkred1,Orange1,Brown1,SkyBlue1,Gray1,Gold1,Violet1,Olive1,Orange1]
COLOR_LIST_DARK = [Darkred2,Orange2,Brown2,SkyBlue2,Gray2,Gold2,Violet2,Olive2,Orange2]
import Utilities.tfs_file_writer as tfs_writer
BPM_NAMES_STR = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"]
def _parse_args():
    parser = OptionParser()
    parser.add_option("-r", "--results",
                    help="results folder",
                    metavar="results",dest="results")
    parser.add_option("-p", "--plane",
                    help="plane to analyze",
                    metavar="PLANE",dest="plane_to_analyze")
    parser.add_option("-o", "--output_path",
                    help="output path",
                    metavar="OUTPUT",dest="output_path")

    parser.add_option("-f", "--free",type="int",
                    help="1: free.out or 2: free2.out",
                    metavar="freeNo",dest="freeNo")
	
    options, _ = parser.parse_args()
    
    print options.freeNo
    
    return options.results, options.plane_to_analyze,options.output_path, options.freeNo

def _configure_plots():
    matplotlib.rc('font', **{'family': 'sans-serif', 'serif': ['Computer Modern']})
    params = {'backend': 'pdf',
              'axes.labelsize': 12,
              'font.size': 14,
              'axes.titlesize': 12,
              'legend.fontsize': 16,
              'xtick.labelsize': 12,
              'ytick.labelsize': 12,
              'text.usetex': False,
              'axes.unicode_minus': True,
              'xtick.major.pad': 10,
              'ytick.major.pad': 10,
              'xtick.minor.pad': 10,
              'ytick.minor.pad': 10,
              }
    matplotlib.rcParams.update(params)



def plotting_automatic_beta_ratio_abs(file_input,file_output,plane,freeNo):
    name_output = os.path.join(file_output,"beta_beating_ratio_abs" + str(plane.lower()) + ".pdf")
    gs = gridspec.GridSpec(1, 1,height_ratios=[1])
    ax1 = subplot(gs[0])
    beta = "BET" + str(plane.upper())
    beta_error_phase = "ERRBET" + str(plane.upper())
    beta_error_amplitude = "BET" + str(plane.upper()) + "STD"
    beta_model = "BET" + str(plane.upper()) + "MDL"
    if plane == "x":
       ax1.set_ylabel(r'$\beta_A,x / \beta_{\Phi,x}$ [-]')
    elif plane == "y":
       ax1.set_ylabel(r'$\beta_A,y / \beta_{\Phi,y}$ [-]')
    ax1.set_xlabel(r'Longitudinal location [m]')
    #for i in range(len(results)):
    freeVer = ".out"
    if freeNo > 1:
        freeVer = "_free{}.out".format(freeNo)
    elif freeNo == 1:
        freeVer = "_free.out"
    
    
    print("freeNo = ",freeNo)
    print("freeVer = ",freeVer)
        
    beta_file_phase = tfs_pandas.read_tfs(os.path.join(file_input,"getbeta" + str(plane.lower()) + freeVer))
    beta_file_amplitude = tfs_pandas.read_tfs(os.path.join(file_input,"getampbeta" + str(plane.lower()) + freeVer))
    beta_error_phase_values = getattr(beta_file_phase,beta_error_phase)
    beta_error_amplitude_values = getattr(beta_file_amplitude,beta_error_amplitude)
    beta_error_total = ((beta_error_amplitude_values/getattr(beta_file_phase,beta))**2 + (beta_error_phase_values*getattr(beta_file_amplitude,beta)/getattr(beta_file_phase,beta)**2)**2)**0.5
    print "BETA RATIO ABS ######"
    print getattr(beta_file_amplitude,beta) / getattr(beta_file_phase,beta)
    print beta_error_total
    ax1.errorbar(beta_file_phase.S, getattr(beta_file_amplitude,beta) / getattr(beta_file_phase,beta) , yerr = (beta_error_total)**0.5 , fmt='o', color=COLOR_LIST_LIGHT[0], markersize=4, markeredgecolor=COLOR_LIST_DARK[0])
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc=4, frameon=True, numpoints = 1, ncol=1)
    tight_layout()
    savefig(name_output)


def main(file_results,plane_to_analyze,output_path,freeNo):
    plotting_automatic_beta_ratio_abs(file_results,output_path,plane_to_analyze,freeNo)
    
if __name__ == "__main__":
    (_file_results,_plane_to_analyze,_output_path,_freeNo) = _parse_args()
    main(_file_results,_plane_to_analyze,_output_path,_freeNo)
      
