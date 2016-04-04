#!/afs/cern.ch/work/o/omc/anaconda/bin/python

# import __init__
import sys
import os
new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../"))
if new_path not in sys.path:
    sys.path.append(new_path)
import numpy as np
import matplotlib

matplotlib.use('qt4agg')  # THIS BACKEND IS NEEDED FOR THE CLEANING !

import matplotlib.pyplot as plt
from Python_Classes4MAD import metaclass
from scipy import stats
from make_fit_plots import plot_fitting
from scipy.spatial import Delaunay
from optparse import OptionParser
from Utilities import tfs_file_writer
from read_Timber_output import merge_data
from IR_planes import IR_definitions

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

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
        self.cl_plot.set_xdata(self.cleaned_data[:,0])
        self.cl_plot.set_ydata(self.cleaned_data[:,1])

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
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0


def clean(data, trapezium):
    mask =  in_hull([data[0,:,0:2]],trapezium)
    cleaned_data = data[mask]
    return cleaned_data


def start_cleaning_data(k,tune_data,tune_data_err):
    data = np.dstack((k,tune_data,tune_data_err))
    plt.figure(figsize=(15,15))
    plt.xlabel('K')
    plt.ylabel('Tune')
    plt.title('Left click: Select corners,  Right click: Cancel selection,  c: Skip')
    plt.errorbar(k,tune_data, yerr=tune_data_err, fmt='o')
    ax = plt.gca()
    cc = clicker_class(ax, data)
    plt.show()
    return cc.return_clean_data()


def run_analysis_simplex(path, beam, ip, bs, working_directory):
    fitx_L, fitx_R, fity_L, fity_R, errx_r, erry_r, errx_l, erry_l, K, dK, Q1, Q2 = lin_fit_data(path)
    
    fitx_L = fitx_L*dK
    fitx_R = fitx_R*dK
    fity_L = fity_L*dK
    fity_R = fity_R*dK
    
    pathx = 'python '+ os.path.join(CURRENT_PATH,'..','GetBetaStarFromKmod.py') + ' -Q 0.28 -L 22.965 -l 6.37 -K %s -D ' %str(abs(K))
    commandx = pathx + str(dK) +','+ str(dK)+ ' -d ' + str(abs(fitx_L))+','+str(abs(fitx_R)) + ' -e ' + str(abs(errx_l))+','+str(abs(errx_r)) + ' -b ' + bs + ' -t ' + (ip+beam)+'.X' 
    pathy = 'python '+ os.path.join(CURRENT_PATH,'..','GetBetaStarFromKmod.py') + ' -Q 0.31 -L 22.965 -l 6.37 -K %s -D ' %str(abs(K))
    commandy = pathy + str(dK) +','+ str(dK)+ ' -d ' + str(abs(fity_L))+','+str(abs(fity_R)) + ' -e ' + str(abs(erry_l))+','+str(abs(erry_r)) + ' -b ' + bs + ' -t ' + (ip+beam)+'.Y' 

    os.system(commandx)
    os.system(commandy)

    
    with open(os.path.join(working_directory,'command.run'), 'a') as commands_write:
        commands_write.write(commandx +' \n')
        commands_write.write(commandy)

    with open('BetaStarResults.dat', 'r') as results_read:
        content = results_read.read()
    
    with open(os.path.join(path,'%s.results' %(ip+beam)), 'w') as results_write:
        results_write.write('* LABEL   BETASTAR   BETASTAR_ERR   WAIST       WAIST_ERR   BETAWAIST   BETAWAIST_ERR   BBEAT_FOC   BBEAT_FOC_ERR   BBEAT_DEF   BBEAT_DEF_ERR \n')
        results_write.write('$ %s      %le        %le            %le         %le         %le         %le             %le         %le             %le         %le \n')
        results_write.write(content)

    os.remove('BetaStarResults.dat')
    calc_BPM_beta(path, ip, beam)


def lin_fit_data(path):
    file_path_R = path+'R.dat'
    file_path_L = path+'L.dat'
    right_data = metaclass.twiss(file_path_R)
    left_data  = metaclass.twiss(file_path_L)

    cleaned_xR = start_cleaning_data(right_data.K, right_data.TUNEX,right_data.TUNEX_ERR)
    cleaned_yR = start_cleaning_data(right_data.K, right_data.TUNEY,right_data.TUNEY_ERR)
    cleaned_xL = start_cleaning_data(left_data.K, left_data.TUNEX,left_data.TUNEX_ERR)
    cleaned_yL = start_cleaning_data(left_data.K, left_data.TUNEY,left_data.TUNEY_ERR)

    fitx_R, covx_r = np.polyfit(cleaned_xR[:,0], cleaned_xR[:,1], 1, cov=True, w = 1/cleaned_xR[:,2]**2)
    fity_R, covy_r = np.polyfit(cleaned_yR[:,0], cleaned_yR[:,1], 1, cov=True, w = 1/cleaned_yR[:,2]**2)
    fitx_L, covx_l = np.polyfit(cleaned_xL[:,0], cleaned_xL[:,1], 1, cov=True, w = 1/cleaned_xL[:,2]**2)
    fity_L, covy_l = np.polyfit(cleaned_yL[:,0], cleaned_yL[:,1], 1, cov=True, w = 1/cleaned_yL[:,2]**2)

    plot_fitting(fitx_L,fitx_R,fity_L,fity_R,left_data,right_data,path)

    dK = 1.0e-5
    K = np.average(left_data.K)
    Q1 = np.average(right_data.TUNEX)
    Q2 = np.average(right_data.TUNEY)

    errx_r = np.sqrt(np.diag(covx_r)[0]) * dK
    erry_r = np.sqrt(np.diag(covy_r)[0]) * dK
    errx_l = np.sqrt(np.diag(covx_l)[0]) * dK
    erry_l = np.sqrt(np.diag(covy_l)[0]) * dK

    return fitx_L[0], fitx_R[0], fity_L[0], fity_R[0], errx_r, erry_r, errx_l, erry_l, K, dK, Q1, Q2 #kmod_data  # Array with all dQ's (slopes of fit scaled with dK) and the dK spread. [xR, xL, yR, yL, dK ]


def which_bpms(ips):
    bpm_dict = {'ip1b1': ['"BPMSW.1L1.B1"', '"BPMSW.1R1.B1"'],
                'ip2b1': ['"BPMSW.1L2.B1"', '"BPMSW.1R2.B1"'],
                'ip5b1': ['"BPMSW.1L5.B1"', '"BPMSW.1R5.B1"'],
                'ip8b1': ['"BPMSW.1L8.B1"', '"BPMSW.1R8.B1"'],
                'ip1b2': ['"BPMSW.1L1.B2"', '"BPMSW.1R1.B2"'],
                'ip2b2': ['"BPMSW.1L2.B2"', '"BPMSW.1R2.B2"'],
                'ip5b2': ['"BPMSW.1L5.B2"', '"BPMSW.1R5.B2"'],
                'ip8b2': ['"BPMSW.1L8.B2"', '"BPMSW.1R8.B2"'],
                }
    return bpm_dict[ips]


def calc_BPM_beta(path, ip, beam):
    if ip == 'ip1' or ip == 'ip5':
        L_ip2bpm = 21.564  # Length between IP and first BPM:  BPMSW.1L1.B1, BPMSW.1R1.B1, BPMSW.1L5.B1, BPMSW.1R5.B1
    elif ip == 'ip8' or ip == 'ip2':
        L_ip2bpm = 21.595    
    ip_data = metaclass.twiss(os.path.join(path,'%s.results' %(ip+beam)))
    bw     = ip_data.BETAWAIST
    bw_err = ip_data.BETAWAIST_ERR
    w      = ip_data.WAIST
    w_err  = ip_data.WAIST_ERR
    label  = ip_data.LABEL
    
    IR = ip + beam 

    if IR_definitions[IR+'.X'][0] == 'foc':
        w[0] = -w[0] 
    elif IR_definitions[IR+'.Y'][0] == 'foc':
        w[1] = -w[1]


    b_bpmR = bw + (L_ip2bpm - w)**2/bw
    b_bpmL = bw + (L_ip2bpm + w)**2/bw

    b_bpmR_err = []
    b_bpmL_err = []
    count = 0

    for l in range(len(label)):
        count += 1
        rerr = []
        lerr = []
        we  = np.linspace(-w_err[l] , w_err[l] ,2) + w[l]
        bwe = np.linspace(-bw_err[l], bw_err[l],2) + bw[l]
        for i in range(2):
            for j in range(2):
                rerr.append(bwe[i] + (L_ip2bpm - we[j])**2/bwe[i])
                lerr.append(bwe[i] + (L_ip2bpm + we[j])**2/bwe[i])
        b_bpmR_err.append((max(rerr)-min(rerr))/2.)
        b_bpmL_err.append((max(lerr)-min(lerr))/2.)

    beta_bpm = np.transpose(np.vstack((b_bpmL,b_bpmR)))
    beta_bpm_err = np.transpose(np.vstack((b_bpmL_err,b_bpmR_err)))
    bpms = which_bpms(ip+beam)

    xdata = tfs_file_writer.TfsFileWriter.open(os.path.join(path, 'getkmodbetax_%s.out' %ip))
    xdata.set_column_width(20)
    xdata.add_column_names(['NAME', 'S'  , 'COUNT',   'BETX',    'BETXSTD',       'BETXMDL'    ,        'MUXMDL'     ,      'BETXRES'    ,    'BETXSTDRES' ])
    xdata.add_column_datatypes(['%s', '%le','%le','%le', '%le', '%le', '%le', '%le', '%le'])

    ydata = tfs_file_writer.TfsFileWriter.open(os.path.join(path, 'getkmodbetay_%s.out' %ip))
    ydata.set_column_width(20)
    ydata.add_column_names(['NAME', 'S'  , 'COUNT',    'BETY',    'BETYSTD',   'BETYMDL'      ,      'MUYMDL'    ,       'BETYRES'    ,    'BETYSTDRES'])
    ydata.add_column_datatypes(['%s', '%le','%le','%le', '%le', '%le', '%le', '%le', '%le'])

    for i in range(len(bpms)):
        xdata.add_table_row([bpms[i], 0, 0, beta_bpm[0][i], beta_bpm_err[0][i], 0, 0, 0, 0 ])
        ydata.add_table_row([bpms[i], 0, 0, beta_bpm[1][i], beta_bpm_err[1][i], 0, 0, 0, 0 ])
    xdata.write_to_file()  
    ydata.write_to_file()  


def check_files(path):
    IR_files = [path+'L.dat', path+'R.dat']
    IR_check = True
    if os.path.exists(IR_files[0])==False or os.path.exists(IR_files[1])==False :
        IR_check = False
        print 'File missing for %s. Ommitting IR in analysis..' %path
    return IR_check


def parse_args():
    usage = 'Usage: %prog -w WORKING_DIR -o OUT_FILE_PATH [options]'
    parser = OptionParser(usage=usage)
    parser.add_option('-b', '--betastar', 
                            help='Estimated beta star of measurements', 
                            action='store', type='string', dest='betastar')
    parser.add_option('-w', '--working_directory', 
                            help='path to working directory with stored KMOD measurement files', 
                            action='store', type='string', dest='work_dir')
    parser.add_option('-i', '--interaction_point', 
                            help='define interaction point: ip1 or ip5', 
                            action='store', type='string', dest='ip')
    parser.add_option('-e', '--beam', 
                            help='define beam used: b1 or b2', 
                            action='store', type='string', dest='beam')
    (options, _) = parser.parse_args()
    return options


if __name__ == '__main__':
    options = parse_args()
    working_directory = options.work_dir
    beam = options.beam
    ip = options.ip
    bs = options.betastar

    IR = ip + beam

    merge_data(working_directory, ip ,beam)

    path = os.path.join(working_directory, IR)
    if check_files(path) == True:
        if not os.path.exists(path):
            os.makedirs(path)
            print '%s folder created..' % IR
        run_analysis_simplex(path, beam, ip, bs,working_directory)
