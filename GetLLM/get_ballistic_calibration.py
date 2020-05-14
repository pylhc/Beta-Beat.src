# -*- coding: utf8 -*-

import sys
import os
from os.path import abspath, join, dirname, pardir
new_path = abspath(join(dirname(abspath(__file__)), pardir))
if new_path not in sys.path:
    sys.path.append(new_path)
from tfs_files import tfs_pandas
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
#import tfs_files.tfs_file_writer as tfs_writer
import utils
import pandas as pd
#from Python_Classes4MAD import metaclass
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from optparse import OptionParser
from scipy.interpolate import interp1d


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

OUTPUT_FILE_PREFIX_PLOT = "plot_beta_ip"
OUTPUT_FILE_PREFIX_CALIBRATION_FILE = "calibration"

#BPM_PREFIXES = "BPMR.5L{ip}.B{beam},BPMYA.4L{ip}.B{beam},BPMWB.4L{ip}.B{beam},BPMSY.4L{ip}.B{beam},BPMS.2L{ip}.B{beam},BPMSW.1L{ip}.B{beam},BPMSW.1R{ip}.B{beam},BPMS.2R{ip}.B{beam},BPMSY.4R{ip}.B{beam},BPMWB.4R{ip}.B{beam},BPMYA.4R{ip}.B{beam}".format(ip=ip, beam=beam).split(",")
#COLUMN_NAMES = ["NAME", "S", "CAL_AMP", "CAL_AMP_STD", "CALIBRATION", "ERROR_CALIBRATION", "CAL_BETA", "CAL_BETA_STD"]
LABELS = ["S","CALIBRATION","ERROR_CALIBRATION","CALIBRATION_PHASE_FIT","ERROR_CALIBRATION_PHASE_FIT"]
LABELS_BETA = ["CALIBRATION_BETA","ERROR_CALIBRATION_BETA"]
IPS = [1, 4, 5]
PLANES = ["X", "Y"]

INITIAL_BETA_STAR_ESTIMATION = 200
LEGEND_POSITION = 1


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


def main(input_path, input_path_model, output_path) :
    utils.iotools.create_dirs(output_path)
    (beta_from_phase_x, beta_from_amp_x) = _get_twiss_for_one_of(input_path, "getbetax_free.out", "getbetax.out")
    (beta_from_phase_y, beta_from_amp_y) = _get_twiss_for_one_of(input_path, "getbetay_free.out", "getbetay.out")
    nominal_model = tfs_pandas.read_tfs(input_path_model)
    _configure_plots()
    files_phase = (beta_from_phase_x, beta_from_phase_y)
    files_amplitude = (beta_from_amp_x, beta_from_amp_y)
    beam = _get_beam_from_model(nominal_model)
    for i in range(len(PLANES)):
        for ip in IPS:
            (summary_amplitude,summary_beta) = _compute_calibration_for_ip_and_plane(ip, PLANES[i], files_phase[i], files_amplitude[i], nominal_model, beam, output_path)
            name_file = OUTPUT_FILE_PREFIX_CALIBRATION_FILE + "_" + str(PLANES[i].lower()) + ".out"
            calibration_file_path = os.path.join(output_path,name_file)
            try:
                summary_amplitude_tot = tfs_pandas.read_tfs(calibration_file_path)
                summary_amplitude_tot = summary_amplitude_tot.append(summary_amplitude)
                tfs_pandas.write_tfs(calibration_file_path,summary_amplitude_tot,save_index=False)
            except:
                tfs_pandas.write_tfs(calibration_file_path,summary_amplitude,save_index=False)


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


def _get_beam_from_model(nominal_model):
    if nominal_model.SEQUENCE.endswith("B1"):
        return 1
    elif nominal_model.SEQUENCE.endswith("B2"):
        return 2
    else:
        raise ValueError("Wrong sequence in model.")


def func_phase(x, A, B):
    return A + ((x - B) ** 2) / A


def print_files(names_range, positions_range, amplitude_ratio_fit, error_amplitude_ratio_fit, amplitude_ratio_measured, error_amplitude_ratio_measured, beta_ratio, error_beta_ratio, beta_ratio_phasefit, error_beta_ratio_phasefit, tfs_file):
    for index in range(len(names_range)):
        list_row_entries = [names_range[index], positions_range[index], amplitude_ratio_fit[index], error_amplitude_ratio_fit[index], amplitude_ratio_measured[index], error_amplitude_ratio_measured[index], beta_ratio[index], error_beta_ratio[index]]
        tfs_file.add_table_row(list_row_entries)


def _compute_calibration_for_ip_and_plane(
    ip, plane, file_phase, file_amplitude, nominal_model, beam, output_path
):

    # Define BPMs to be used for a combination of IP and Beam
    bpms = {1: {1: ['BPMR.5L1.B1',
                    'BPMYA.4L1.B1',
                    'BPMWB.4L1.B1',
                    'BPMSY.4L1.B1',
                    'BPMS.2L1.B1',
                    'BPMSW.1L1.B1',
                    'BPMSW.1R1.B1',
                    'BPMS.2R1.B1',
                    'BPMSY.4R1.B1',
                    'BPMWB.4R1.B1',
                    'BPMYA.4R1.B1'],
                2: ['BPM.5L1.B2', 
                    'BPMYA.4L1.B2', 
                    'BPMWB.4L1.B2', 
                    'BPMSY.4L1.B2', 
                    'BPMS.2L1.B2', 
                    'BPMSW.1L1.B2', 
                    'BPMSW.1R1.B2', 
                    'BPMS.2R1.B2', 
                    'BPMSY.4R1.B2', 
                    'BPMWB.4R1.B2', 
                    'BPMYA.4R1.B2']
                },
            4: {1: [
                    #'BPMYA.5L4.B1',
                    'BPMWI.A5L4.B1',
                    'BPMWA.B5L4.B1',
                    'BPMWA.A5L4.B1',
                    'BPMWA.A5R4.B1',
                    'BPMWA.B5R4.B1',
                    #'BPMYB.5R4.B1',
                    #'BPMYA.6R4.B1',
                    ],
                2: [ 
                    #'BPMYB.5L4.B2',
                    'BPMWA.B5L4.B2',
                    'BPMWA.A5L4.B2',
                    'BPMWA.A5R4.B2',
                    'BPMWA.B5R4.B2',
                    'BPMWI.A5R4.B2',
                    #'BPMYA.5R4.B2',
                    #'BPMYB.6R4.B2'
                    ]
                },
            5: {1: ['BPMYA.4L5.B1', 
                    'BPMWB.4L5.B1', 
                    'BPMSY.4L5.B1', 
                    'BPMS.2L5.B1', 
                    'BPMSW.1L5.B1', 
                    'BPMSW.1R5.B1', 
                    'BPMS.2R5.B1', 
                    'BPMSY.4R5.B1', 
                    'BPMWB.4R5.B1', 
                    'BPMYA.4R5.B1', 
                    'BPM.5R5.B1'],
                2: ['BPMYA.4L5.B2', 
                    'BPMWB.4L5.B2', 
                    'BPMSY.4L5.B2', 
                    'BPMS.2L5.B2', 
                    'BPMSW.1L5.B2', 
                    'BPMSW.1R5.B2', 
                    'BPMS.2R5.B2', 
                    'BPMSY.4R5.B2', 
                    'BPMWB.4R5.B2', 
                    'BPMYA.4R5.B2', 
                    'BPMR.5R5.B2']
                }
            }

    print('Current data: IP{} and Beam{} / Plane {}'.format(ip, beam, plane))
    
    # Get the BPM for a combination of IP and Beam
    if ip in bpms.keys() and beam in bpms[ip].keys():
        bpm_names = bpms[ip][beam]
    else:
        raise Exception("No BPMs defined for IP or beam combination: {ip}, "
                        "{beam}")

    # Filter the list of BPMs we got in our dataframe
    bpm_names_filter = file_phase.set_index("NAME").reindex(bpm_names)\
                                 .dropna().index

    IR_positions_common = nominal_model.set_index("NAME")\
                                       .loc[bpm_names_filter, "S"].dropna()
    IR_positions_phase = file_phase.set_index("NAME")\
                                       .loc[bpm_names_filter, "S"].dropna()

    # Get the maximum and minimum positions of the elements in the IR
    IR_minimum = nominal_model.set_index("NAME").loc[bpm_names_filter[0], "S"]
    IR_maximum = nominal_model.set_index("NAME").loc[bpm_names_filter[-1], "S"]
    IP_position = (IR_maximum - IR_minimum) / 2

    print(nominal_model)
    print(nominal_model.columns)

    beta = "BET" + plane

    beta_phase_error = "ERRBET" + str(plane.upper())
    beta_amplitude_error = "BET" + str(plane.upper()) + "STD"
    beta_range = file_phase.set_index("NAME").loc[bpm_names_filter, beta].dropna()
    beta_range_err = file_phase.set_index("NAME").loc[bpm_names_filter, beta_phase_error].dropna()
    beta_range_amp = file_amplitude.set_index("NAME").loc[bpm_names_filter, beta].dropna()
    beta_range_amp_err = file_amplitude.set_index("NAME").loc[bpm_names_filter, beta_amplitude_error].dropna()

    initial_values = [INITIAL_BETA_STAR_ESTIMATION, IP_position]
    beta_phasefit_curve, beta_phasefit_curve_err = curve_fit(func_phase, IR_positions_phase, beta_range, p0=initial_values, sigma=beta_range_err, maxfev=1000000)
    beta_phasefit = func_phase(IR_positions_common, beta_phasefit_curve[0], beta_phasefit_curve[1])
    beta_phasefit_max = func_phase(IR_positions_common, beta_phasefit_curve[0] + beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] + beta_phasefit_curve_err[1, 1] ** 0.5)
    beta_phasefit_min = func_phase(IR_positions_common, beta_phasefit_curve[0] - beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] - beta_phasefit_curve_err[1, 1] ** 0.5)
    beta_phasefit_err = (beta_phasefit_max - beta_phasefit_min) / 2
    amplitude_ratio_phasefit = pd.DataFrame((beta_phasefit / beta_range_amp) ** 0.5)
    amplitude_ratio_measured = ((beta_range / beta_range_amp) ** 0.5)
    error_amplitude_ratio_phasefit = (((beta_phasefit_err) ** 2 * 1 / (beta_range_amp * beta_phasefit * 4) + beta_phasefit_err ** 2 * beta_phasefit / (beta_range_amp ** 3 * 4)) ** 0.5)
    error_amplitude_ratio_measured = (((beta_range_err) ** 2 * 1 / (beta_range_amp * beta_range * 4) + beta_range_amp_err**2 * beta_range / (beta_range_amp ** 3 * 4)) ** 0.5)
    beta_ratio = beta_phasefit / beta_range_amp
    error_beta_ratio = ((beta_range_amp_err * beta_phasefit / beta_range_amp ** 2) ** 2 + (beta_phasefit_err / beta_range_amp) ** 2) ** 0.5
    summary_calibration_factors_amplitude = pd.concat([IR_positions_common,amplitude_ratio_measured.dropna(), error_amplitude_ratio_measured.dropna(),amplitude_ratio_phasefit.dropna(),error_amplitude_ratio_phasefit.dropna()], axis=1) 
    summary_calibration_factors_amplitude.columns = LABELS
    summary_calibration_factors_amplitude = summary_calibration_factors_amplitude.reset_index()
    summary_calibration_factors_beta = pd.concat([beta_ratio,error_beta_ratio], axis=1)
    summary_calibration_factors_beta.columns = LABELS_BETA
    summary_calibration_factors_beta = summary_calibration_factors_beta.reset_index()


    _plot_calibration_fit(output_path, beam, plane, ip, beta_phasefit_curve, beta_phasefit_curve_err, IR_positions_common, IR_positions_common, beta_range, beta_range_err, beta_range_amp, beta_range_amp_err, nominal_model)

    return (summary_calibration_factors_amplitude, summary_calibration_factors_beta)


def _plot_calibration_fit(output_path, beam, plane, ip, beta_phasefit_curve, beta_phasefit_curve_err, position_fit, IR_positions_common, beta_range, beta_range_err, beta_range_amp, beta_range_amp_err, nominal_model):
    name_file_pdf = OUTPUT_FILE_PREFIX_PLOT + str(ip) + "_" + "B" + str(beam) + "_" + str(plane) + ".png"
    file_path_pdf = os.path.join(output_path, name_file_pdf)
    gs = matplotlib.gridspec.GridSpec(1, 1, height_ratios=[1])
    ax = plt.subplot(gs[0])

    # Add labels to the plot
    label_phase = r'$\beta^{\phi}$'
    label_phase_fit = r'$\beta^{fit}$'
    ax.set_ylabel(r'$\beta_{{{plane}}}$ [m]'.format(plane=plane), fontsize=19)

    beta_phasefit_allpositions = func_phase(position_fit, beta_phasefit_curve[0], beta_phasefit_curve[1])
    beta_phasefit_max_allpositions = func_phase(position_fit, beta_phasefit_curve[0] + beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] + beta_phasefit_curve_err[1, 1] ** 0.5)
    beta_phasefit_min_allpositions = func_phase(position_fit, beta_phasefit_curve[0] - beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] - beta_phasefit_curve_err[1, 1] ** 0.5)
    beta_phasefit_err_allpositions = (beta_phasefit_max_allpositions - beta_phasefit_min_allpositions) / 2
    xfine = np.linspace(position_fit[0], position_fit[len(position_fit) - 1], 2000)
    beta_fit = (func_phase(xfine, beta_phasefit_curve[0], beta_phasefit_curve[1]))
    plt.grid(False)
    ax.set_xlabel('position [m]')
    
    ax.errorbar(IR_positions_common, beta_range, yerr=beta_range_err, fmt='o', color=SkyBlue1, markersize=6, markeredgecolor=SkyBlue3, label= label_phase)


    # ------------------
    # Filter the nominal model and display it
    filter_s_min = nominal_model['S'] >= position_fit[0]
    filter_s_max  = nominal_model['S'] <= position_fit[-1]
    nominal_filtered = (nominal_model[filter_s_min])[filter_s_max]

    # Get the X and Y axes
    beta = "BET" + plane
    x_nominal = nominal_filtered.set_index('S').index
    beta_nominal_model = nominal_filtered.set_index(beta).index

    # Interpolate the model because we only have β values at the BPMs
    func = interp1d(x_nominal, beta_nominal_model, kind='cubic')
    beta_model_interp = func(xfine)
    
    # And plot
    plt.plot(xfine, beta_model_interp, label=r'$\beta^{mdl}$')
    print('Length of data points: {}'.format(len(beta_nominal_model)))
    # ------------------

    # ------------------
    # Compute the χ²
    chi = chisquare(f_obs=beta_fit, f_exp=beta_model_interp)
    ax.set_title(u'IP{} B{} / χ²: {}'.format(ip, beam, chi[0] / len(xfine)))

    # ------------------
    ax.errorbar(xfine,beta_fit, fmt='r', color=ScarletRed3, markersize=4, markeredgecolor=ScarletRed3,label=label_phase_fit)
    ax.legend(numpoints=1, ncol=1, loc="best", fontsize=14)
    matplotlib.pyplot.savefig(file_path_pdf, bbox_inches='tight', dpi=150)


def _get_twiss_for_one_of(directory, *file_names):
    for file_name in file_names:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            file_name_amplitude = "getampbeta" + file_name[7:]
            file_path_amplitude = os.path.join(directory, file_name_amplitude)
            return (tfs_pandas.read_tfs(file_path), tfs_pandas.read_tfs(file_path_amplitude))
    raise IOError("None of the files exist:\n\t" + "\n\t".join(file_names))


if __name__ == "__main__":
    (_input_path, _input_path_model, _output_path)=_parse_args()
    main(_input_path, _input_path_model, _output_path)
