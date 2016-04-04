
import os
import __init__  # @UnusedImport
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import Utilities.tfs_file_writer as tfs_writer
import Utilities.bpm

from Python_Classes4MAD import metaclass
from scipy.optimize import curve_fit
from optparse import OptionParser


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

BPM_PREFIXES = ["BPMWB.4L", "BPMSY.4L", "BPMS.2L", "BPMSW.1L", "BPMSW.1L", "BPMSW.1R", "BPMS.2R", "BPMSY.4R", "BPMWB.4R"]
COLUMN_NAMES = ["NAME", "S", "CAL_AMP", "CAL_AMP_STD", "CALIBRATION", "ERROR_CALIBRATION", "CAL_BETA", "CAL_BETA_STD"]

IPS = [1, 5]
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
    Utilities.iotools.create_dirs(output_path)
    (beta_from_phase_x, beta_from_amp_x) = _get_twiss_for_one_of(input_path, "getbetax_free.out", "getbetax.out")
    (beta_from_phase_y, beta_from_amp_y) = _get_twiss_for_one_of(input_path, "getbetay_free.out", "getbetay.out")
    nominal_model = metaclass.twiss(input_path_model)
    _configure_plots()
    files_phase = (beta_from_phase_x, beta_from_phase_y)
    files_amplitude = (beta_from_amp_x, beta_from_amp_y)
    beam = _get_beam_from_model(nominal_model)
    for i in range(len(PLANES)):
        tfs_file = _get_tfs_file(output_path, PLANES[i], beam)
        for ip in IPS:
            (names_range, positions_range, amplitude_ratio_phasefit, error_amplitude_ratio_phasefit, amplitude_ratio_measured, error_amplitude_ratio_measured, beta_ratio_phasefit, error_beta_ratio_phasefit) = _compute_calibration_for_ip_and_plane(ip, PLANES[i], files_phase[i], files_amplitude[i], nominal_model, beam, output_path)
            print_files(names_range, positions_range, amplitude_ratio_phasefit, error_amplitude_ratio_phasefit, amplitude_ratio_measured, error_amplitude_ratio_measured, beta_ratio_phasefit, error_beta_ratio_phasefit, beta_ratio_phasefit, error_beta_ratio_phasefit, tfs_file)
        tfs_file.write_to_file()


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

def _get_tfs_file(output_path, plane, beam):

    name_file = OUTPUT_FILE_PREFIX_CALIBRATION_FILE + "_" + str(plane.lower()) + ".out"
    file_path = os.path.join(output_path, name_file)
    tfs_file_writer_calibration = tfs_writer.TfsFileWriter.open(file_path)
    tfs_file_writer_calibration.add_string_descriptor("PLANE", plane)
    tfs_file_writer_calibration.add_column_names(COLUMN_NAMES)
    tfs_file_writer_calibration.add_column_datatypes(["%s"] + ["%le"] * (len(COLUMN_NAMES) - 1))
    return tfs_file_writer_calibration


def func_phase(x, A, B):
    return A + ((x - B) ** 2) / A


def print_files(names_range, positions_range, amplitude_ratio_fit, error_amplitude_ratio_fit, amplitude_ratio_measured, error_amplitude_ratio_measured, beta_ratio, error_beta_ratio, beta_ratio_phasefit, error_beta_ratio_phasefit, tfs_file):
    for index in range(len(names_range)):
        list_row_entries = [names_range[index], positions_range[index], amplitude_ratio_fit[index], error_amplitude_ratio_fit[index], amplitude_ratio_measured[index], error_amplitude_ratio_measured[index], beta_ratio[index], error_beta_ratio[index]]
        tfs_file.add_table_row(list_row_entries)


def _compute_calibration_for_ip_and_plane(
    ip, plane, file_phase, file_amplitude, nominal_model, beam, output_path
):
    IR_positions_common = []
    beta_range = []
    beta_range_amp = []
    beta_range_err = []
    beta_range_amp_err = []
    names_range_IP = []
    position_fit = []
    BPM_names = []
    # define if its beam 1 or beam 2
    for bpm_prefix in BPM_PREFIXES:
        bpm_prefix = bpm_prefix + str(ip) + ".B" + str(beam)
        BPM_names.append(bpm_prefix)
        position_fit.append(nominal_model.S[nominal_model.NAME.index(bpm_prefix)])
    IR_minimum = getattr(nominal_model, "S")[nominal_model.NAME.index(BPM_names[0])]
    IR_maximum = nominal_model.S[nominal_model.NAME.index(BPM_names[len(BPM_names) - 1])]
    IP_position = (IR_maximum - IR_minimum) / 2

    # just in case a BPM is missing in one file but not in the other ( BPM from ampl but not from phase )
    commonbpms = Utilities.bpm.intersect([file_amplitude, file_phase])

    for i in range(0, len(commonbpms)):
        common_names_amplitude_phase = str.upper(commonbpms[i][1])

        if nominal_model.S[nominal_model.NAME.index(common_names_amplitude_phase)] >= IR_minimum and nominal_model.S[nominal_model.NAME.index(common_names_amplitude_phase)] <= IR_maximum:
            IR_positions_common.append(nominal_model.S[nominal_model.NAME.index(common_names_amplitude_phase)])
            beta_range.append(getattr(file_phase, "BET" + plane)[i])
            beta_range_err.append(((getattr(file_phase, "STDBET" + plane)[i] ** 2) + (getattr(file_phase, "ERRBET" + plane)[i] ** 2)) ** 0.5)
            beta_range_amp.append(getattr(file_amplitude, "BET" + plane)[i])
            beta_range_amp_err.append(getattr(file_amplitude, "BET" + plane + "STD")[i])
            names_range_IP.append(file_phase.NAME[i])

    initial_values = [INITIAL_BETA_STAR_ESTIMATION, IP_position]
    beta_phasefit_curve, beta_phasefit_curve_err = curve_fit(func_phase, IR_positions_common, beta_range, p0=initial_values, sigma=beta_range_err, maxfev=1000000)
    beta_phasefit = []
    beta_phasefit_max = []
    beta_phasefit_min = []
    beta_phasefit_err = []
    amplitude_ratio_phasefit = []
    amplitude_ratio_measured = []
    error_amplitude_ratio_phasefit = []
    error_amplitude_ratio_measured = []
    beta_ratio = []
    error_beta_ratio = []
    for i in range(len(IR_positions_common)):
        beta_phasefit.append(func_phase(IR_positions_common[i], beta_phasefit_curve[0], beta_phasefit_curve[1]))
        beta_phasefit_max.append(func_phase(IR_positions_common[i], beta_phasefit_curve[0] + beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] + beta_phasefit_curve_err[1, 1] ** 0.5))
        beta_phasefit_min.append(func_phase(IR_positions_common[i], beta_phasefit_curve[0] - beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] - beta_phasefit_curve_err[1, 1] ** 0.5))
        beta_phasefit_err.append((beta_phasefit_max[i] - beta_phasefit_min[i]) / 2)
        amplitude_ratio_phasefit.append((beta_phasefit[i] / beta_range_amp[i]) ** 0.5)
        amplitude_ratio_measured.append((beta_range[i] / beta_range_amp[i]) ** 0.5)
        error_amplitude_ratio_phasefit.append(((beta_phasefit_err[i]) ** 2 * 1 / (beta_range_amp[i] * beta_phasefit[i] * 4) + beta_phasefit_err[i] ** 2 * beta_phasefit[i] / (beta_range_amp[i] ** 3 * 4)) ** 0.5)
        error_amplitude_ratio_measured.append(((beta_range_err[i]) ** 2 * 1 / (beta_range_amp[i] * beta_range[i] * 4) + beta_range_amp_err[i]**2 * beta_range[i] / (beta_range_amp[i] ** 3 * 4)) ** 0.5)
        beta_ratio.append(beta_phasefit[i] / beta_range_amp[i])
        error_beta_ratio.append(((beta_range_amp_err[i] * beta_phasefit[i] / beta_range_amp[i] ** 2) ** 2 + (beta_phasefit_err[i] / beta_range_amp[i]) ** 2) ** 0.5)
   
    _plot_calibration_fit(output_path, beam, plane, ip, beta_phasefit_curve, beta_phasefit_curve_err, position_fit, IR_positions_common, beta_range, beta_range_err, beta_range_amp, beta_range_amp_err)
    return(names_range_IP, IR_positions_common, amplitude_ratio_phasefit, error_amplitude_ratio_phasefit, amplitude_ratio_measured, error_amplitude_ratio_measured, beta_ratio, error_beta_ratio)


def _plot_calibration_fit(output_path, beam, plane, ip, beta_phasefit_curve, beta_phasefit_curve_err, position_fit, IR_positions_common, beta_range, beta_range_err, beta_range_amp, beta_range_amp_err):
    # beta from fit for ALL BPMS in the IR
    name_file_pdf = OUTPUT_FILE_PREFIX_PLOT + str(ip) + "_" + "B" + str(beam) + "_" + str(plane) + ".pdf"
    file_path_pdf = os.path.join(output_path, name_file_pdf)
    beta_phasefit_allpositions = []
    beta_phasefit_max_allpositions = []
    beta_phasefit_min_allpositions = []
    beta_phasefit_err_allpositions = []
    beta_mdl = []
    j = 0
    if plane == "X":
        label_amp = r'$\beta$ from amplitude (x) '
        label_phase = r'$\beta$ from phase (x) '
        label_phase_fit = r'$\beta$ fit (x) '
    elif plane == "Y":
        label_amp = r'$\beta$ from amplitude (y) '
        label_phase = r'$\beta$ from phase (y) '
        label_phase_fit = r'$\beta$ fit (y) '
    for i in position_fit:
        beta_phasefit_allpositions.append(func_phase(i, beta_phasefit_curve[0], beta_phasefit_curve[1]))
        beta_phasefit_max_allpositions.append(func_phase(i, beta_phasefit_curve[0] + beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] + beta_phasefit_curve_err[1, 1] ** 0.5))
        beta_phasefit_min_allpositions.append(func_phase(i, beta_phasefit_curve[0] - beta_phasefit_curve_err[0, 0] ** 0.5, beta_phasefit_curve[1] - beta_phasefit_curve_err[1, 1] ** 0.5))
        beta_phasefit_err_allpositions.append((beta_phasefit_max_allpositions[j] - beta_phasefit_min_allpositions[j]) / 2)
        j = j + 1
    xfine = np.linspace(position_fit[0], position_fit[len(position_fit) - 1], 2000)
    for i in xfine:
        beta_mdl.append(func_phase(i, beta_phasefit_curve[0], beta_phasefit_curve[1]))
    gs = matplotlib.gridspec.GridSpec(1, 1, height_ratios=[1])
    ax2 = plt.subplot(gs[0])
    plt.grid(False)

    ax2.set_ylabel(r'$\beta [m]$', fontsize=19)
    ax2.set_xlabel('position [m]')
    ax2.errorbar(IR_positions_common, beta_range, yerr=beta_range_err, fmt='o', color=SkyBlue1, markersize=4, markeredgecolor=SkyBlue3, label= label_phase)
    ax2.errorbar(IR_positions_common, beta_range_amp, yerr=beta_range_amp_err, fmt='o', color=ScarletRed1, markersize=3, markeredgecolor=ScarletRed3, label= label_amp)
    ax2.errorbar(position_fit, beta_phasefit_allpositions, yerr=beta_phasefit_err_allpositions, fmt='o', color=Orange1, markersize=4, markeredgecolor=Orange3, label=label_phase_fit )
    ax2.errorbar(xfine, beta_mdl, fmt='r', color=Orange1, markersize=4, markeredgecolor=Orange3)
    ax2.legend(numpoints=1, ncol=1, loc=LEGEND_POSITION, fontsize=14)
    matplotlib.pyplot.savefig(file_path_pdf, bbox_inches='tight')


def _get_twiss_for_one_of(directory, *file_names):
    for file_name in file_names:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            file_name_amplitude = "getampbeta" + file_name[7:]
            file_path_amplitude = os.path.join(directory, file_name_amplitude)
            return (metaclass.twiss(file_path), metaclass.twiss(file_path_amplitude))
    raise IOError("None of the files exist:\n\t" + "\n\t".join(file_names))


if __name__ == "__main__":
    (_input_path, _input_path_model, _output_path)=_parse_args()
    main(_input_path, _input_path_model, _output_path)
