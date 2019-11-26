import os
import sys
import scipy.io
from os.path import abspath, join, dirname, pardir

new_path = abspath(join(dirname(abspath(__file__)), pardir))
if new_path not in sys.path:
    sys.path.append(new_path)
from drive import drive_runner
from tfs_files import tfs_pandas
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import *
from optparse import OptionParser
import tfs_files.tfs_file_writer as tfs_writer
from scipy.optimize import curve_fit

#beam 1
BPM_IP1_B1 = ["BPMR.5L1.B1","BPMYA.4L1.B1","BPMWB.4L1.B1", "BPMSY.4L1.B1", "BPMS.2L1.B1", "BPMSW.1L1.B1", "BPMSW.1R1.B1", "BPMS.2R1.B1", "BPMSY.4R1.B1", "BPMWB.4R1.B1","BPMYA.4R1.B1"]
BPM_IP5_B1 = ["BPMYA.4L5.B1","BPMWB.4L5.B1", "BPMSY.4L5.B1", "BPMS.2L5.B1", "BPMSW.1L5.B1", "BPMSW.1R5.B1", "BPMS.2R5.B1", "BPMWB.4R5.B1","BPMYA.4R5.B1","BPM.5R5.B1"]
#beam 2
BPM_IP1_B2 = ["BPM.5L1.B2","BPMYA.4L1.B2","BPMWB.4L1.B2", "BPMSY.4L1.B2", "BPMS.2L1.B2", "BPMSW.1L1.B2", "BPMSW.1R1.B2", "BPMS.2R1.B2", "BPMSY.4R1.B2", "BPMWB.4R1.B2","BPMYA.4R1.B2"]
BPM_IP5_B2 = ["BPMYA.4L5.B2","BPMWB.4L5.B2", "BPMSY.4L5.B2", "BPMS.2L5.B2", "BPMSW.1L5.B2", "BPMSW.1R5.B2", "BPMS.2R5.B2", "BPMSY.4R5.B2", "BPMWB.4R5.B2","BPMYA.4R5.B2","BPMR.5R5.B2"]
BPMS_B1 = [BPM_IP1_B1,BPM_IP5_B1]
BPMS_B2 = [BPM_IP1_B2,BPM_IP5_B2]
BPM_dispersion_IP1_B1 = ["BPMSY.4L1.B1", "BPMS.2L1.B1", "BPMSW.1L1.B1", "BPMSW.1R1.B1", "BPMS.2R1.B1", "BPMSY.4R1.B1"]
BPM_dispersion_IP5_B1 = ["BPMSY.4L5.B1", "BPMS.2L5.B1", "BPMSW.1L5.B1", "BPMSW.1R5.B1", "BPMS.2R5.B1"]
BPM_dispersion_IP1_B2 = ["BPMSY.4L1.B2", "BPMS.2L1.B2", "BPMSW.1L1.B2", "BPMSW.1R1.B2", "BPMS.2R1.B2", "BPMSY.4R1.B2"]
BPM_dispersion_IP5_B2 = ["BPMSY.4L5.B2", "BPMS.2L5.B2", "BPMSW.1L5.B2", "BPMSW.1R5.B2", "BPMS.2R5.B2", "BPMSY.4R5.B2"]
BPM_dispersion_B1 = [BPM_dispersion_IP1_B1,BPM_dispersion_IP5_B1]
BPM_dispersion_B2 = [BPM_dispersion_IP1_B2,BPM_dispersion_IP5_B2]
COLUMN_NAMES = ["NAME", "S", "CALIBRATION", "ERROR_CALIBRATION","CAL_FIT","CAL_FIT_STD"]
label_size = 10
plt.rcParams.update({"font.size": 20})
Butter1 = '#FCE94F'
Butter2 = '#EDD400'
Butter3 = '#C4A000'
Orange1 = '#FCAF3E'
Orange2 = '#F57900'
Orange3 = '#CE5C00'
Chocolate1 = '#E9B96E'
Chocolate2 = '#C17D11'
Chocolate3 = '#8F5902'
Chameleon1 = '#8AE234'
Chameleon2 = '#73D216'
Chameleon3 = '#4E9A06'
SkyBlue1 = '#729FCF'
SkyBlue2 = '#3465A4'
SkyBlue3 = '#204A87'
Plum1 = '#AD7FA8'
Plum2 = '#75507B'
Plum3 = '#5C3566'
ScarletRed1 = '#EF2929'
ScarletRed2 = '#CC0000'
ScarletRed3 = '#A40000'
Aluminium1 = '#EEEEEC'
Aluminium2 = '#D3D7CF'
Aluminium3 = '#BABDB6'
Aluminium4 = '#888A85'
Aluminium5 = '#555753'
Aluminium6 = '#2E3436'

def _parse_args():
    parser = OptionParser()
    parser.add_option("-i", "--input",
                    help="Input measurement path",
                    metavar="INPUT",dest="input_path")
    parser.add_option("-m", "--model",
                    help="Input model path",
                    metavar="INPUT_MODEL", default=None, dest="input_path_model")
    parser.add_option("-o", "--output",
                    help="Output path",
                    metavar="OUTPUT", dest="output_path")
    options, _ = parser.parse_args()
    if options.output_path is None:
        options.output_path = options.input_path
    return options.input_path, options.input_path_model, options.output_path


def _get_beam_from_model(nominal_model):
    if nominal_model.SEQUENCE.endswith("B1"):
        return 1
    elif nominal_model.SEQUENCE.endswith("B2"):
        return 2
    else:
        raise ValueError("Wrong sequence in model.")

def _get_tfs_file(output_path,beam):
    name_file =  "calibration_dispersion_x.out"
    file_path = os.path.join(output_path, name_file)
    tfs_file_writer_calibration = tfs_writer.TfsFileWriter.open(file_path)
    tfs_file_writer_calibration.add_string_descriptor("PLANE", "X")
    tfs_file_writer_calibration.add_column_names(COLUMN_NAMES)
    tfs_file_writer_calibration.add_column_datatypes(["%s"] + ["%le"] * (len(COLUMN_NAMES) - 1))
    return tfs_file_writer_calibration


def print_files(bpm_set,position_phase,calibration_factor_dispersion,calibration_factor_dispersion_error,calibration_factors_from_fit,calibration_factors_from_fit_error,tfs_file):
    for index in range(len(bpm_set)):
        list_row_entries = [bpm_set[index],position_phase[index],calibration_factor_dispersion[index],calibration_factor_dispersion_error[index],calibration_factors_from_fit[index],calibration_factors_from_fit_error[index]]
        tfs_file.add_table_row(list_row_entries)

def get_dispersion_files(input_path):
    dispersion_file = os.path.join(input_path,"getDx.out")
    normalized_dispersion_file = os.path.join(input_path,"getNDx.out")
    beta_amp_file = os.path.join(input_path,"getampbetax.out")
    beta_phase_file = os.path.join(input_path,"getbetax.out")
    return dispersion_file,normalized_dispersion_file,beta_amp_file,beta_phase_file

def get_beta_from_phase_ir(beta_phase_file,bpm_set):
    beta_phase_measurement = tfs_pandas.read_tfs(beta_phase_file)
    position_phase = beta_phase_measurement.set_index("NAME").loc[bpm_set,"S"]
    beta_phase = beta_phase_measurement.set_index("NAME").loc[bpm_set,"BETX"]
    beta_phase_error = beta_phase_measurement.set_index("NAME").loc[bpm_set,"ERRBETX"]
    return position_phase,beta_phase,beta_phase_error 
    
    
def get_beta_from_amplitude_ir(beta_amp_file,bpm_set):
    beta_amp_measurement = tfs_pandas.read_tfs(beta_amp_file)
    position_amp = beta_amp_measurement.set_index("NAME").loc[bpm_set,"S"]
    beta_amp  = beta_amp_measurement.set_index("NAME").loc[bpm_set,"BETX"]
    beta_amp_error = beta_amp_measurement.set_index("NAME").loc[bpm_set,"BETXSTD"]
    return position_amp, beta_amp,beta_amp_error

def get_dispersion_normalized_dispersion_ir(dispersion_file,normalized_dispersion_file,bpm_set):
    dispersion_measurement = tfs_pandas.read_tfs(dispersion_file)
    normalized_dispersion_measurement = tfs_pandas.read_tfs(normalized_dispersion_file)
    dispersion = dispersion_measurement.set_index("NAME").loc[bpm_set,"DX"]
    dispersion_error = dispersion_measurement.set_index("NAME").loc[bpm_set,"STDDX"] 
    normalized_dispersion = normalized_dispersion_measurement.set_index("NAME").loc[bpm_set,"NDX"]
    normalized_dispersion_error = normalized_dispersion_measurement.set_index("NAME").loc[bpm_set,"STDNDX"]
    return dispersion,dispersion_error,normalized_dispersion,normalized_dispersion_error
    
def get_dispersion_using_phase(normalized_dispersion,normalized_dispersion_error,beta_phase,beta_phase_error):
    dispersion_from_phase = normalized_dispersion*(beta_phase)**0.5
    dispersion_from_phase_error = ((normalized_dispersion_error*((beta_phase)**0.5))**2 + (0.5*normalized_dispersion*beta_phase_error*(beta_phase)**(-0.5))**2)**0.5
    return dispersion_from_phase, dispersion_from_phase_error
    
def get_dispersion_fit(position,normalized_dispersion,normalized_dispersion_error):
    print normalized_dispersion
    #popt, pcov = curve_fit(func,position,normalized_dispersion, sigma=1/(np.array(normalized_dispersion_error)**2))
    popt, pcov = curve_fit(func,position,normalized_dispersion)
    perr = np.sqrt(np.diag(pcov))
    # prepare confidence level curves
    nstd = 1.0 # to draw 5-sigma intervals
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    return popt,perr,popt_up,popt_dw
    
    
def get_calibration_factors_from_dispersion(dispersion_from_phase,dispersion_from_phase_error,dispersion,dispersion_error):
    calibration_factor_dispersion =  dispersion_from_phase/dispersion
    calibration_factor_dispersion_error = ((dispersion_from_phase_error/dispersion)**2 + (dispersion_error*dispersion_from_phase/(dispersion)**2)**2)**0.5
    return calibration_factor_dispersion,calibration_factor_dispersion_error
    
def get_normalized_dispersion_from_fit(position_phase,popt,popt_up,popt_dw):
   normalized_dispersion_from_fit = func(position_phase, *popt)
   normalized_dispersion_from_fit_positive_error = func(position_phase, *popt_up)
   normalized_dispersion_from_fit_negative_error = func(position_phase, *popt_dw)
   normalized_dispersion_from_fit_error = normalized_dispersion_from_fit_positive_error - normalized_dispersion_from_fit
   return normalized_dispersion_from_fit, normalized_dispersion_from_fit_error,normalized_dispersion_from_fit_positive_error,normalized_dispersion_from_fit_negative_error
   
def get_calibration_from_normalized_dispersion_fit(dispersion_from_fit,dispersion_from_fit_error,dispersion_from_fit_phase,dispersion_from_fit_phase_error):    
   calibration_factors_from_fit = dispersion_from_fit_phase/dispersion_from_fit
   calibration_factors_from_fit_error = ((dispersion_from_fit_phase_error/dispersion_from_fit)**2 + (dispersion_from_fit_error*dispersion_from_fit_phase/(dispersion_from_fit)**2)**2)**0.5
   return calibration_factors_from_fit, calibration_factors_from_fit_error 
       
def func(x, a, b):
    return a * x + b

    
def plotting_dispersion(position_phase,position_phase_fit,normalized_dispersion,normalized_dispersion_error,normalized_dispersion_from_phase,normalized_dispersion_from_phase_error,normalized_dispersion_from_fit,normalized_dispersion_from_fit_error,normalized_dispersion_from_fit_phase, normalized_dispersion_from_fit_error_phase,normalized_dispersion_from_fit_positive_error,normalized_dispersion_from_fit_negative_error,normalized_dispersion_from_fit_positive_error_phase,normalized_dispersion_from_fit_negative_error_phase,output_path,bpm_set):
    file_name = "calibration_dispersion_beam_no_normalized_" + bpm_set[0][-1] + "_ir_" + bpm_set[0][-4] + "_with_fit.pdf"
    file_name_path = os.path.join(output_path,file_name)
    fig, ax = plt.subplots(2)
    gs = matplotlib.gridspec.GridSpec(1, 1, height_ratios=[1])
    ax = plt.subplot(gs[0])
    ax.axes.get_xaxis().set_ticks([])
    ax.set_xlim(np.min(position_phase)-10,np.max(position_phase)+10)
    plt.grid(False)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.errorbar(position_phase,normalized_dispersion_from_phase, yerr=normalized_dispersion_from_phase_error, xerr=0,fillstyle="none",markersize=6, elinewidth=2,capsize=4,capthick=2, color=SkyBlue3,markeredgecolor=SkyBlue3, fmt='o',label=r'Dx (calibration independent)')
    ax.errorbar(position_phase,normalized_dispersion, yerr=normalized_dispersion_error, xerr=0, color=ScarletRed3,fillstyle="none",markersize=6, elinewidth=2,capsize=4,capthick=2,markeredgecolor=ScarletRed3, fmt='o', label=r'Dx (calibration dependent)')
    ax.plot(position_phase_fit, normalized_dispersion_from_fit_phase, SkyBlue1, lw=2, label='_nolegend_')
    xlabel('Position [m]', fontsize=18)
    ylabel('Dx [m]', fontsize=18)
    xticks(np.arange(np.min(position_phase)+10,np.max(position_phase)-10, step=100))
    legend(loc='best',fontsize=16, numpoints=1)
    matplotlib.pyplot.savefig(file_name_path, bbox_inches='tight')   
    
def main(input_path,input_path_model,output_path):
   nominal_model = tfs_pandas.read_tfs(os.path.join(input_path_model,"twiss.dat"))
   beam = _get_beam_from_model(nominal_model)
   [dispersion_file,normalized_dispersion_file,beta_amp_file,beta_phase_file] = get_dispersion_files(input_path)
   if beam == 1:
      bpm_sets = BPMS_B1
      bpm_dispersion = BPM_dispersion_B1
      tfs_file = _get_tfs_file(output_path, beam)
   elif beam == 2:
      bpm_sets = BPMS_B2
      bpm_dispersion = BPM_dispersion_B2
      tfs_file = _get_tfs_file(output_path, beam)
   for i in range(len(bpm_sets)):
       [position_phase,beta_phase,beta_phase_error] = get_beta_from_phase_ir(beta_phase_file,bpm_sets[i])
       [position_phase_dispersion,beta_phase_dispersion,beta_phase_error_dispersion] = get_beta_from_phase_ir(beta_phase_file,bpm_sets[i])
       [position_phase_fit,beta_phase_fit,beta_phase_error_fit] = get_beta_from_phase_ir(beta_phase_file,bpm_dispersion[i])
       [position_amp,beta_amp,beta_amp_error] = get_beta_from_amplitude_ir(beta_amp_file,bpm_sets[i])
       [dispersion,dispersion_error,normalized_dispersion,normalized_dispersion_error] = get_dispersion_normalized_dispersion_ir(dispersion_file,normalized_dispersion_file,bpm_sets[i])
       [dispersion_fit,dispersion_error_fit,normalized_dispersion_fit,normalized_dispersion_error_fit] = get_dispersion_normalized_dispersion_ir(dispersion_file,normalized_dispersion_file,bpm_dispersion[i])
       [dispersion_from_phase,dispersion_from_phase_error] = get_dispersion_using_phase(normalized_dispersion,normalized_dispersion_error,beta_phase,beta_phase_error)
       [dispersion_from_phase_fit,dispersion_from_phase_error_fit] = get_dispersion_using_phase(normalized_dispersion_fit,normalized_dispersion_error_fit,beta_phase_fit,beta_phase_error_fit)
       [calibration_factor_dispersion,calibration_factor_dispersion_error] = get_calibration_factors_from_dispersion(dispersion_from_phase,dispersion_from_phase_error,dispersion,dispersion_error)
       [popt,perr,popt_up,popt_dw] = get_dispersion_fit(position_phase_fit,dispersion_fit,dispersion_error_fit)
       [popt_phase,perr_phase,popt_up_phase,popt_dw_phase] = get_dispersion_fit(position_phase_fit,dispersion_from_phase_fit,dispersion_from_phase_error_fit)
       [dispersion_from_fit,dispersion_from_fit_error,dispersion_from_fit_positive_error,dispersion_from_fit_negative_error] = get_normalized_dispersion_from_fit(position_phase_fit,popt,popt_up,popt_dw)
       [dispersion_from_fit_phase,dispersion_from_fit_error_phase,dispersion_from_fit_positive_error_phase,dispersion_from_fit_negative_error_phase] = get_normalized_dispersion_from_fit(position_phase_fit,popt_phase,popt_up_phase,popt_dw_phase)
       [calibration_factors_from_fit, calibration_factors_from_fit_error] = get_calibration_from_normalized_dispersion_fit(dispersion,dispersion_error,dispersion_from_fit_phase,dispersion_from_fit_error_phase)
       print_files(bpm_sets[i],position_phase,calibration_factor_dispersion,calibration_factor_dispersion_error,calibration_factors_from_fit, calibration_factors_from_fit_error,tfs_file)
       plotting_dispersion(position_phase,position_phase_fit,dispersion,dispersion_error,dispersion_from_phase,dispersion_from_phase_error,dispersion_from_fit,dispersion_from_fit_error,dispersion_from_fit_phase, dispersion_from_fit_error_phase,dispersion_from_fit_positive_error,dispersion_from_fit_negative_error,dispersion_from_fit_positive_error_phase,dispersion_from_fit_negative_error_phase,output_path,bpm_sets[i])
   tfs_file.write_to_file()  
    
    

if __name__ == "__main__":
    (_input_path, _input_path_model, _output_path)=_parse_args()
    main(_input_path, _input_path_model, _output_path) 

