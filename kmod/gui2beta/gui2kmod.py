#!/afs/cern.ch/work/o/omc/anaconda/bin/python

import sys
import os
from os.path import abspath, join, dirname, pardir
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

new_path = abspath(join(dirname(abspath(__file__)), pardir, pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

from kmod.gui2beta.read_Timber_output import merge_data
from kmod.gui2beta import Magnet_definitions, KModUtilities
from Python_Classes4MAD import metaclass
from kmod.gui2beta.make_fit_plots import plot_fitting
from tfs_files import tfs_file_writer
from utils import outliers

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

# TODO: Short term: Think about the accelerator class here for positions and Ks
# TODO: Immediately: Use a logger for logging
# TODO: Immediately: get rid of repetive code and use loops and functions
# TODO: Immediately: Use tfs_pandas instead of metaclass

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--BetastarAndWaist',
                        help='Estimated beta star of measurements and waist shift',
                        action='store', type=str, dest='betastar')
    parser.add_argument('--working_directory',
                        help='path to working directory with stored KMOD measurement files',
                        action='store', type=str, dest='work_dir')
    parser.add_argument('--cminus',
                        help='C Minus',
                        action='store', type=float, dest='cminus', default=argparse.SUPPRESS)
    parser.add_argument('--misalignment',
                        help='misalignment of the modulated quadrupoles in m',
                        action='store', type=float, dest='misalign', default=argparse.SUPPRESS)
    parser.add_argument('--errorK',
                        help='error in K of the modulated quadrupoles, unit m^-2',
                        action='store', type=float, dest='ek', default=argparse.SUPPRESS)
    parser.add_argument('--Tuneuncertainty',
                        help='tune measurement uncertainty',
                        action='store', type=float, dest='tunemeasuncertainty', default=2.5e-5)
    parser.add_argument('--beam',
                        help='define beam used: b1 or b2',
                        action='store', type=str, dest='beam', choices=['b1', 'b2', 'B1', 'B2'], required=True)
    parser.add_argument('--instruments',
                        help='define instruments (use keywords from twiss) at which beta should be calculated , separated by comma, e.g. MONITOR,RBEND,INSTRUMENT,TKICKER',
                        action='store', type=str, dest='instruments', default='MONITOR,SBEND,TKICKER,INSTRUMENT')
    parser.add_argument('--log',
                        help='flag for creating a log file',
                        action='store_true', dest='log')
    parser.add_argument('--noautoclean',
                        help='flag for manually cleaning data',
                        action='store_true', dest='a_clean')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument( '--circuit',
                       help='circuit names of the modulated quadrupoles',
                       action='store', type=str, dest='magnets')
    group.add_argument( '--interaction_point',
                       help='define interaction point',
                       action='store', type=str, dest='ip', choices=['ip1', 'ip2', 'ip5', 'ip8', 'IP1', 'IP2', 'IP5', 'IP8'])

    options = parser.parse_args()

    return options


class clicker_class(object):
    def __init__(self, ax, data, pix_err=1):
        self.canvas = ax.get_figure().canvas
        self.cid = None
        self.data = data
        self.pt_lst = []
        self.pt_plot = ax.plot([], [], marker='o',
                               linestyle='-', zorder=5)[0]
        self.cl_plot = ax.plot([], [], color='r', marker='o',
                               linestyle='', zorder=5)[0]
        self.tr_plot = ax.plot([], [], color='g', marker='o',
                               linestyle='-', zorder=5)[0]
        self.pix_err = pix_err
        self.connect_sf()

    def set_visible(self, visible):
        '''sets if the curves are visible '''
        self.pt_plot.set_visible(visible)

    def clear(self):
        '''Clears the points'''
        self.pt_lst = []
        x, y = [], []
        self.pt_plot.set_xdata(x)
        self.pt_plot.set_ydata(y)
        self.cl_plot.set_xdata(x)
        self.cl_plot.set_ydata(y)
        self.tr_plot.set_xdata(x)
        self.tr_plot.set_ydata(y)

        self.canvas.draw()

    def connect_sf(self):
        if self.cid is None:
            self.cid = self.canvas.mpl_connect('button_press_event',
                                               self.click_event)
            self.cid = self.canvas.mpl_connect('key_press_event',
                                               self.key_event)

    def disconnect_sf(self):
        if self.cid is not None:
            self.canvas.mpl_disconnect(self.cid)
            # print self.data
            # print self.cleaned_data
            self.cid = None

    def key_event(self, event):
        ''' Extracts locations from the user'''
        if event.key == 'c':
            self.cleaned_data = self.data[0]
            self.disconnect_sf()
            plt.close()
            return

    def click_event(self, event):
        ''' Extracts locations from the user'''
        if event.key == 'shift':
            self.pt_lst = []
            self.redraw()
            return
        if event.xdata is None or event.ydata is None:
            return
        if event.button == 1:
            self.pt_lst.append((event.xdata, event.ydata))
            if len(self.pt_lst) > 4:
                self.disconnect_sf()
                plt.close()
                return
        elif event.button == 3:
            self.clear()
        self.redraw()
        if len(self.pt_lst) > 3:
            self.start_clean()

    def start_clean(self):
        self.cleaned_data = clean(self.data, self.pt_lst)
        self.cl_plot.set_xdata(self.cleaned_data[:, 0])
        self.cl_plot.set_ydata(self.cleaned_data[:, 1])

        pt_list = self.pt_lst
        pt_list.append(pt_list[0])
        ptdata = zip(*pt_list)

        self.tr_plot.set_xdata(ptdata[0])
        self.tr_plot.set_ydata(ptdata[1])
        self.canvas.draw()

    def remove_pt(self, loc):
        if len(self.pt_lst) > 0:
            self.pt_lst.pop(np.argmin(map(lambda x:
                                          np.sqrt((x[0] - loc[0]) ** 2 +
                                                  (x[1] - loc[1]) ** 2),
                                          self.pt_lst)))

    def redraw(self):
        if len(self.pt_lst) > 0:
            x, y = zip(*self.pt_lst)
        else:
            x, y = [], []
        self.pt_plot.set_xdata(x)
        self.pt_plot.set_ydata(y)

        self.canvas.draw()

    def return_clean_data(self):
        '''Returns the clicked points in the format the rest of the
        code expects'''
        return self.cleaned_data


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p) >= 0


def clean(data, trapezium):
    mask = in_hull([data[0, :, 0:2]], trapezium)
    cleaned_data = data[mask]
    return cleaned_data


def start_cleaning_data(k, tune_data, tune_data_err):
    data = np.dstack((k, tune_data, tune_data_err))
    plt.figure(figsize=(15, 15))
    plt.xlabel('K')
    plt.ylabel('Tune')
    plt.title('Left click: Select corners,  Right click: Cancel selection,  c: Skip')
    plt.errorbar(k, tune_data, yerr=tune_data_err, fmt='o')
    ax = plt.gca()
    cc = clicker_class(ax, data)
    plt.show()
    return cc.return_clean_data()


def automatic_cleaning_data(k,tune_data, tune_data_err, limit=1e-5):
    data = np.dstack((k, tune_data, tune_data_err))
    mask = outliers.get_filter_mask(tune_data, x_data=k, limit=limit)
    return data[0, mask, :]


def run_analysis_simplex(path, beam, magnet1, magnet2, hor_bstar, vert_bstar, waist, working_directory, instruments, ek, misalign,
                         cminus, twiss, log, logfile, auto_clean):

    fitx_2, fitx_1, fity_2, fity_1, errx_1, erry_1, errx_2, erry_2, K1, K2, dK, Qx1, Qy1, Qx2, Qy2 = lin_fit_data(path, beam,
                                                                                                      working_directory,
                                                                                                      magnet1, magnet2, log, logfile, auto_clean)

    fitx_2 = fitx_2 * dK
    fitx_1 = fitx_1 * dK
    fity_2 = fity_2 * dK
    fity_1 = fity_1 * dK

    if (Magnet_definitions.MagnetPolarity(magnet1, beam, twiss) == 1. and
        Magnet_definitions.MagnetPolarity(magnet2, beam, twiss) == -1.):

        if log == True:
            logfile.write('Focussing magnet: %s  \n' % (magnet1))
            logfile.write('\n')

        fitx_foc = fitx_1
        fitx_def = fitx_2

        fity_foc = fity_2
        fity_def = fity_1

        errx_foc = errx_1
        errx_def = errx_2

        erry_foc = erry_2
        erry_def = erry_1

        Qx_foc = Qx1
        Qy_foc = Qy1

        Qx_def = Qx2
        Qy_def = Qy2

        K_foc = K1
        K_def = K2

        l_foc = Magnet_definitions.MagnetLength(magnet1, beam, twiss)
        l_def = Magnet_definitions.MagnetLength(magnet2, beam, twiss)

    elif (Magnet_definitions.MagnetPolarity(magnet1, beam, twiss) == -1. and
          Magnet_definitions.MagnetPolarity(magnet2, beam, twiss) == 1.):
        if log == True:
            logfile.write('Focussing magnet: %s  \n' % (magnet2))
            logfile.write('\n')

        fitx_foc = fitx_2
        fitx_def = fitx_1

        fity_foc = fity_1
        fity_def = fity_2

        errx_foc = errx_2
        errx_def = errx_1

        erry_foc = erry_1
        erry_def = erry_2

        Qx_foc = Qx2
        Qy_foc = Qy2

        Qx_def = Qx1
        Qy_def = Qy1        

        K_foc = K2
        K_def = K1

        l_foc = Magnet_definitions.MagnetLength(magnet2, beam, twiss)
        l_def = Magnet_definitions.MagnetLength(magnet1, beam, twiss)

    L_star = Magnet_definitions.Lstar(magnet1, magnet2, beam, twiss)

    resultsx = KModUtilities.analysis(Qx_foc, Qx_def, Qy_foc, L_star, misalign, K_foc, dK, l_foc, K_def, dK, l_def, fitx_foc, errx_foc,
                                      fitx_def, errx_def, ek, ek, cminus, hor_bstar, waist,
                                      (magnet1 + '-' + magnet2 + '.' + beam) + '.X', log, logfile)
    resultsy = KModUtilities.analysis(Qy_foc, Qy_def, Qx_foc, L_star, misalign, K_foc, dK, l_foc, K_def, dK, l_def, fity_foc, erry_foc,
                                      fity_def, erry_def, ek, ek, cminus, vert_bstar, waist,
                                      (magnet1 + '-' + magnet2 + '.' + beam) + '.Y', log, logfile)

    results = tfs_file_writer.TfsFileWriter.open(os.path.join(path, get_results_filename()))
    results.set_column_width(15)
    results.add_column_names(
        ['LABEL', 'BETAWAIST', 'BETAWAIST_ERR', 'WAIST', 'WAIST_ERR', 'BETA_AV_FOC', 'BETA_AV_FOC_ERR', 'BETA_AV_DEF',
         'BETA_AV_DEF_ERR'])
    results.add_column_datatypes(['%s', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])

    results.add_table_row(resultsx)
    results.add_table_row(resultsy)

    results.write_to_file()

    calc_beta_star(path, magnet1, magnet2, beam, L_star, twiss)

    for instr in instruments:
        calc_beta_instr(path, magnet1, magnet2, beam, instr, log, logfile, twiss)


def calc_beta_instr(path, magnet1, magnet2, beam, instr, log, logfile, twiss):
    if instr == 'MONITOR':
        name = 'BPM'
    else:
        name = instr

    if Magnet_definitions.FindKeywordBetweenMagnets(magnet1, magnet2, instr, beam, twiss):

        if log:
            logfile.write('%s found, calculating Betas\n' % (name))

        names, positions = Magnet_definitions.ReturnDataofBPMinBetweenMagnets(magnet1, magnet2, instr, beam, twiss)

        L_star_position = Magnet_definitions.LstarPosition(magnet1, magnet2, beam, twiss)

        result_waist = metaclass.twiss(os.path.join(path, get_results_filename()))
        beta_waist = result_waist.BETAWAIST
        beta_waist_err = result_waist.BETAWAIST_ERR
        waist = result_waist.WAIST
        waist_err = result_waist.WAIST_ERR

        waist = waist * (1, -1)

        Magnet1Pos = Magnet_definitions.MagnetPosition(magnet1, beam, twiss)
        Magnet2Pos = Magnet_definitions.MagnetPosition(magnet2, beam, twiss)

        if Magnet1Pos > Magnet2Pos and Magnet_definitions.MagnetPolarity(magnet1, beam, twiss) == -1:
            waist = - waist

        elif Magnet2Pos > Magnet1Pos and Magnet_definitions.MagnetPolarity(magnet2, beam, twiss) == -1:
            waist = - waist

        Waist_pos = L_star_position + waist

        beta_bpm_x = beta_waist[0] + abs(Waist_pos[0] - positions) ** 2 / beta_waist[0]

        beta_waist_err_x = np.linspace(-beta_waist_err[0], beta_waist_err[0], 2) + beta_waist[0]
        waist_err_x = np.linspace(-waist_err[0], waist_err[0], 2) + waist[0]
        Waist_pos_err_x = L_star_position + waist_err_x

        beta_err_x = np.zeros((4, len(positions)))
        n = 0
        for i in range(2):
            for j in range(2):
                beta_err_x[n] = beta_waist_err_x[i] + abs(Waist_pos_err_x[j] - positions) ** 2 / beta_waist_err_x[i]
                n += 1

        err_x = (abs(np.nanmax(beta_err_x, axis=0) - np.nanmin(beta_err_x, axis=0))) / 2.

        beta_bpm_y = beta_waist[1] + abs(Waist_pos[1] - positions) ** 2 / beta_waist[1]

        beta_waist_err_y = np.linspace(-beta_waist_err[1], beta_waist_err[1], 2) + beta_waist[1]
        waist_err_y = np.linspace(-waist_err[1], waist_err[1], 2) + waist[1]
        Waist_pos_err_y = L_star_position + waist_err_y

        beta_err_y = np.zeros((4, len(positions)))
        n = 0
        for i in range(2):
            for j in range(2):
                beta_err_y[n] = beta_waist_err_y[i] + abs(Waist_pos_err_y[j] - positions) ** 2 / beta_waist_err_y[i]
                n += 1

        err_y = (abs(np.nanmax(beta_err_y, axis=0) - np.nanmin(beta_err_y, axis=0))) / 2.

        if name == 'BPM':
            xdata = tfs_file_writer.TfsFileWriter.open(os.path.join(path, get_beta_filename("x")))
        else:
            xdata = tfs_file_writer.TfsFileWriter.open(os.path.join(path, 'Beta_%s_X.out' % name))
        xdata.set_column_width(20)
        xdata.add_column_names(['NAME', 'S', 'COUNT', 'BETX', 'BETXSTD', 'BETXMDL', 'MUXMDL', 'BETXRES', 'BETXSTDRES'])
        xdata.add_column_datatypes(['%s', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])

        if name == 'BPM':
            ydata = tfs_file_writer.TfsFileWriter.open(os.path.join(path, get_beta_filename("y")))
        else:
            ydata = tfs_file_writer.TfsFileWriter.open(os.path.join(path, 'Beta_%s_Y.out' % name))
        ydata.set_column_width(20)
        ydata.add_column_names(['NAME', 'S', 'COUNT', 'BETY', 'BETYSTD', 'BETYMDL', 'MUYMDL', 'BETYRES', 'BETYSTDRES'])
        ydata.add_column_datatypes(['%s', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])

        for i in range(len(names)):
            xdata.add_table_row([names[i], 0, 0, beta_bpm_x[i], err_x[i], 0, 0, 0, 0])
            ydata.add_table_row([names[i], 0, 0, beta_bpm_y[i], err_y[i], 0, 0, 0, 0])
        xdata.write_to_file()
        ydata.write_to_file()

    else:
        if log:
            logfile.write('No %s found in between magnets\n' % (name))


def calc_beta_star(path, magnet1, magnet2, beam, lstar, twiss):
    if Magnet_definitions.FindParentBetweenMagnets(magnet1, magnet2, 'OMK', beam, twiss):

        results_write = tfs_file_writer.TfsFileWriter.open(
            os.path.join(path, get_beta_star_filename()))
        results_write.set_column_width(20)
        results_write.add_column_names(
            ['LABEL', 'BETASTAR', 'BETASTAR_ERR', 'WAIST', 'WAIST_ERR', 'BETAWAIST', 'BETAWAIST_ERR', 'PHASEADV', 'ERRPHASEADV'])
        results_write.add_column_datatypes(['%s', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])

        results = metaclass.twiss(os.path.join(path, get_results_filename()))
        beta_w = results.BETAWAIST

        for i, b_w in enumerate(beta_w):
            label = results.LABEL[i]

            waist = results.WAIST[i]

            beta_w_err = results.BETAWAIST_ERR[i]
            waist_err = results.WAIST_ERR[i]

            beta_star = b_w + waist ** 2 / b_w

            beta_star_err = np.zeros(4)

            beta_star_err[0] = b_w + beta_w_err + waist ** 2 / (b_w + beta_w_err)
            beta_star_err[1] = b_w - beta_w_err + waist ** 2 / (b_w - beta_w_err)
            beta_star_err[2] = b_w + (waist + waist_err) ** 2 / b_w
            beta_star_err[3] = b_w + (waist - waist_err) ** 2 / b_w

            std = math.sqrt(max(abs(beta_star_err[0] - beta_star), abs(beta_star_err[1] - beta_star)) ** 2 +
                            max(abs(beta_star_err[2] - beta_star), abs(beta_star_err[3] - beta_star)) ** 2)

            phadv, ephadv = KModUtilities.phase_adv_from_kmod(lstar, beta_star, std, waist, waist_err)

            results_write.add_table_row([label, beta_star, std, waist, waist_err, b_w, beta_w_err, phadv, ephadv])

        results_write.write_to_file()


def lin_fit_data(path, beam, working_directory, magnet1, magnet2, log, logfile, auto_clean):
    file_path_1 = working_directory + '/' + magnet1 + '.' + beam + '.dat'
    file_path_2 = working_directory + '/' + magnet2 + '.' + beam + '.dat'

    right_data = metaclass.twiss(file_path_1)
    left_data = metaclass.twiss(file_path_2)

    if auto_clean == True:
        cleaned_x1 = start_cleaning_data(right_data.K, right_data.TUNEX, right_data.TUNEX_ERR)
        cleaned_y1 = start_cleaning_data(right_data.K, right_data.TUNEY, right_data.TUNEY_ERR)
        cleaned_x2 = start_cleaning_data(left_data.K, left_data.TUNEX, left_data.TUNEX_ERR)
        cleaned_y2 = start_cleaning_data(left_data.K, left_data.TUNEY, left_data.TUNEY_ERR)
    else:
        cleaned_x1 = automatic_cleaning_data(right_data.K, right_data.TUNEX, right_data.TUNEX_ERR)
        cleaned_y1 = automatic_cleaning_data(right_data.K, right_data.TUNEY, right_data.TUNEY_ERR)
        cleaned_x2 = automatic_cleaning_data(left_data.K, left_data.TUNEX, left_data.TUNEX_ERR)
        cleaned_y2 = automatic_cleaning_data(left_data.K, left_data.TUNEY, left_data.TUNEY_ERR)

    fitx_1, covx_1 = np.polyfit(cleaned_x1[:, 0], cleaned_x1[:, 1], 1, cov=True, w=1 / cleaned_x1[:, 2] ** 2)
    fity_1, covy_1 = np.polyfit(cleaned_y1[:, 0], cleaned_y1[:, 1], 1, cov=True, w=1 / cleaned_y1[:, 2] ** 2)
    fitx_2, covx_2 = np.polyfit(cleaned_x2[:, 0], cleaned_x2[:, 1], 1, cov=True, w=1 / cleaned_x2[:, 2] ** 2)
    fity_2, covy_2 = np.polyfit(cleaned_y2[:, 0], cleaned_y2[:, 1], 1, cov=True, w=1 / cleaned_y2[:, 2] ** 2)

    plot_fitting(fitx_2, fitx_1, fity_2, fity_1, left_data, right_data, path)

    dK = 1.0e-5
    K2 = np.average(cleaned_x2[:, 0])
    K1 = np.average(cleaned_x1[:, 0])

    Qx1 = np.average(cleaned_x1[:, 1])
    Qy1 = np.average(cleaned_y1[:, 1])

    Qx2 = np.average(cleaned_x2[:, 1])
    Qy2 = np.average(cleaned_y2[:, 1])


    errx_1 = np.sqrt(np.diag(covx_1)[0]) * dK
    erry_1 = np.sqrt(np.diag(covy_1)[0]) * dK
    errx_2 = np.sqrt(np.diag(covx_2)[0]) * dK
    erry_2 = np.sqrt(np.diag(covy_2)[0]) * dK


    return fitx_2[0], fitx_1[0], fity_2[0], fity_1[
        0], errx_1, erry_1, errx_2, erry_2, K1, K2, dK, Qx1, Qy1, Qx2, Qy2  # kmod_data  # Array with all dQ's (slopes of fit scaled with dK) and the dK spread. [xR, xL, yR, yL, dK ]


def returnmagnetname(circuit, beam, twiss):
    circuit = circuit.split('.')

    number = circuit[0][-1]
    side = circuit[1][0]
    ip = circuit[1][1]

    searchstring = '.'+str(number)+str(side)+str(ip)

    magnet = Magnet_definitions.findQuadrupoleType(searchstring, beam, twiss)
    return magnet


def returncircuitname(magnet, beam):
    magnet = magnet.split('.')
    number = magnet[1][0]
    side = magnet[1][1]
    ip = magnet[1][2]
    name = 'RQ' + str(number) + '.' + str(side) + str(ip)
    if magnet[0] != 'MQXA':
        name += beam.upper()

    return name


def _main():
    options = parse_args()

    IP_default_err = {'cminus': 1e-3, 'misalign': 0.006, 'ek': 0.001}
    Circuit_default_err = {'cminus': 1e-3, 'misalign': 0.001, 'ek': 0.001}

    if "cminus" not in options:
        if options.ip is not None:
            cminus = IP_default_err['cminus']
        else:
            cminus = Circuit_default_err['cminus']
    else:
        cminus = options.cminus

    if "ek" not in options:
        if options.ip is not None:
            ek = IP_default_err['ek']
        else:
            ek = Circuit_default_err['ek']
    else:
        ek = options.ek

    if "misalign" not in options:
        if options.ip is not None:
            misalign = IP_default_err['misalign']
        else:
            misalign = Circuit_default_err['misalign']
    else:
        misalign = options.misalign

    working_directory = options.work_dir
    beam = options.beam.upper()

    instruments = options.instruments.split(',')
    instruments = [x.upper() for x in instruments]

    bs = options.betastar
    bs = bs.split(",")
    
    if len(bs)==2:
        hor_bstar = bs[0]
        vert_bstar = bs[0]
        waist = bs[1]
    if len(bs)==3:
        hor_bstar = bs[0]
        vert_bstar = bs[1]
        waist = bs[2]
    
    
    auto_clean = options.a_clean
    command = open(working_directory + '/command.run', 'a')
    command.write(str(' '.join(sys.argv)))
    command.write('\n')
    command.close()

    if beam == 'B1':
        twissfile = os.path.join(CURRENT_PATH, "sequences", "twiss_lhcb1.dat")
    else:
        twissfile = os.path.join(CURRENT_PATH, "sequences", "twiss_lhcb2.dat")
    twiss = metaclass.twiss(twissfile)

    if options.ip is not None:
        if options.ip == 'ip1' or options.ip == 'IP1':
            magnet1, magnet2 = 'MQXA.1L1', 'MQXA.1R1'
        elif options.ip == 'ip5' or options.ip == 'IP5':
            magnet1, magnet2 = 'MQXA.1L5', 'MQXA.1R5'
        elif options.ip == 'ip8' or options.ip == 'IP8':
            magnet1, magnet2 = 'MQXA.1L8', 'MQXA.1R8'
        elif options.ip == 'ip2' or options.ip == 'IP2':
            magnet1, magnet2 = 'MQXA.1L2', 'MQXA.1R2'

    else:
        circuits = options.magnets
        circuits = circuits.split(",")
        circuit1, circuit2 = circuits
        magnet1 = returnmagnetname(circuit1, beam, twiss)
        magnet2 = returnmagnetname(circuit2, beam, twiss)

    path = os.path.join(working_directory, magnet1 + '.' + magnet2 + '.' + beam)

    if not os.path.exists(path):
        os.makedirs(path)
    if options.log == True:
        logdata = open(path + '/data.log', 'w')

    merge_data(working_directory, magnet1, returncircuitname(magnet1, beam), magnet2, returncircuitname(magnet2, beam),
               beam, options.ip, options.tunemeasuncertainty)

    run_analysis_simplex(path, beam, magnet1, magnet2, hor_bstar, vert_bstar, waist, working_directory, instruments, ek,
                         misalign, cminus, twiss, options.log, logdata, auto_clean)

    logdata.close()


def get_beta_filename(plane):
    return "getkmodbeta{:s}.out".format(plane)


def get_beta_star_filename():
    return "beta_star.out"


def get_results_filename():
    return "results.out"



if __name__ == '__main__':
    _main()

