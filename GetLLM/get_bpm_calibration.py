import os
import sys
from os.path import abspath, join, dirname

new_path = abspath(join(dirname(abspath(__file__)), os.pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

import matplotlib
import matplotlib.pyplot as plt
import logging
import tfs_files.tfs_file_writer as tfs_writer
import utils.bpm

from Python_Classes4MAD import metaclass
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


COLUMN_NAMES = ["NAME", "S", "CALIBRATION", "ERROR_CALIBRATION"]

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
    
    if options.input_path is None:
        options.input_path = os.path.dirname(os.path.realpath(__file__))
        
    if options.output_path is None:
        options.output_path = options.input_path
    return options.input_path, options.input_path_model, options.output_path


def main(input_path, input_path_model, output_path) :
    print("output_path = ", output_path)
    print("input_path_model = ", input_path_model)
    utils.iotools.create_dirs(output_path)
    (beta_from_phase_x, beta_from_amp_x) = _get_twiss_for_one_of(input_path, "getbetax_free.out", "getbetax.out")
    (beta_from_phase_y, beta_from_amp_y) = _get_twiss_for_one_of(input_path, "getbetay_free.out", "getbetay.out")
    
    #nominal_model = metaclass.twiss(input_path_model)
    
    _configure_plots()
    
    files_phase = (beta_from_phase_x, beta_from_phase_y)
    files_amplitude = (beta_from_amp_x, beta_from_amp_y)
    
    calibs = [0]*2
    calibs_err = [0]*2
    names = [0]*2
    spos = [0]*2
    
    for i in range(len(PLANES)):
        tfs_file = _get_tfs_file(output_path, PLANES[i])

        (names[i], spos[i],\
              calibs[i], \
              calibs_err[i])   = _compute_calibration(PLANES[i], files_phase[i], files_amplitude[i], output_path)
            
        print_files(names[i], spos[i],  calibs[i], calibs_err[i], tfs_file)
        tfs_file.write_to_file()
        

    _doplots(names, spos,  calibs, calibs_err)
    
def _doplots(names, spos, cals, cals_err):
    
    tmp = names[0]
    nn = tmp[:]
    i=0
    for name in names[0]:
        # remove PR.
        efflen = len(name)-3
        newlen = (efflen)*2   
        vertname = [""]*newlen
        jj=0
        for j in range(3,len(name)):
            vertname[jj] = name[j]
            jj=jj+1
            vertname[jj] = "\n"
            jj=jj+1
            
        nn[i] = vertname
        #print(vertname)
        i=i+1
        
        
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.errorbar(spos[0], cals[0], yerr=cals_err[0], color='dodgerblue')
    
    ax = fig.add_subplot(212)
    ax.errorbar(spos[1], cals[1], yerr=cals_err[1], color='dodgerblue')
    #plt.xticks(spos[0], nn)
    plt.xticks(spos[0], names[0],rotation='vertical')
    
    print(names[0])
    
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    
    plt.show()


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



def _get_tfs_file(output_path, plane):

    name_file = OUTPUT_FILE_PREFIX_CALIBRATION_FILE + "_" + str(plane.lower()) + ".out"
    file_path = os.path.join(output_path, name_file)
    tfs_file_writer_calibration = tfs_writer.TfsFileWriter.open(file_path)
    tfs_file_writer_calibration.add_string_descriptor("PLANE", plane)
    tfs_file_writer_calibration.add_column_names(COLUMN_NAMES)
    tfs_file_writer_calibration.add_column_datatypes(["%s"] + ["%le"] * (len(COLUMN_NAMES) - 1))
    return tfs_file_writer_calibration


def print_files(names_range, positions_range, beta_ratio, error_beta_ratio, tfs_file):
    for index in range(len(names_range)):
        list_row_entries = [names_range[index], positions_range[index], beta_ratio[index], error_beta_ratio[index]]
        tfs_file.add_table_row(list_row_entries)

#_________________________________________
''' Calculates the calibration   '''
def _compute_calibration(plane, file_phase, file_amplitude, output_path):
     
    commonbpms = utils.bpm.intersect([file_amplitude, file_phase])
    nbpms = len(commonbpms)     
    
    names_range =      [""] * nbpms
    positions_common = [0] * nbpms
    beta_ratio       = [0] * nbpms
    error_beta_ratio = [0] * nbpms

    
    #print("file_phase = ", file_phase.NAME)
    #print("file_amplitude = ", file_amplitude.NAME)
    
    #sysbet  = getattr(file_phase, "SYSBET" + plane)
    #statbet = getattr(file_phase, "STATBET" + plane)
    betp_err = getattr(file_phase, "ERRBET" + plane)
    beta_err = getattr(file_amplitude, "BET" + plane + "STD")
    
    betp_arr = getattr(file_phase, "BET" + plane)
    beta_arr = getattr(file_amplitude, "BET" + plane )
    
    for i in range(0, nbpms):
        bpmname = str.upper(commonbpms[i][1])
        print(bpmname)
        names_range[i] = bpmname
        idxa = file_amplitude.NAME.index(bpmname);
        idxp = file_phase.NAME.index(bpmname);

        positions_common[i] =file_amplitude.S[i]  

        #print(i, bpmname, idxa, idxp, file_phase.NAME[idxp] )
        
        beta = beta_arr[idxa];
        betp = betp_arr[idxp];
        beta_ratio[i] = beta / betp;
        
        # print(i, " bet amp = ",beta, betp, beta_ratio[i])
        
        #betp_err = ((sysbet[i] ** 2) + (statbet[i] ** 2)) ** 0.5
        
        error_beta_ratio[i] = (( beta_err[i] * betp / beta**2 )**2   +  (betp_err[i] / beta)**2)**0.5
            
    #bphase = metaclass.twiss(file_phase)
    #bamp = metaclass.twiss(file_amplitude)
    #print(bphase.S)
    #print(bamp.S)
    
    return (names_range, positions_common, \
            beta_ratio, error_beta_ratio)




def _get_twiss_for_one_of(directory, *file_names):
    for file_name in file_names:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            file_name_amplitude = "getampbeta" + file_name[7:]
            file_path_amplitude = os.path.join(directory, file_name_amplitude)
            return (metaclass.twiss(file_path), metaclass.twiss(file_path_amplitude))
    raise IOError("None of the files :\n\t\t" + "\n\t\t".join(file_names) + "\n\t exist in \n\t\t" + directory )


if __name__ == "__main__":
    logging.basicConfig()
    (_input_path, _input_path_model, _output_path)=_parse_args()
    main(_input_path, _input_path_model, _output_path)
