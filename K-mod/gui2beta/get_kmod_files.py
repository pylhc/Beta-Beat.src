#!/afs/cern.ch/work/o/omc/anaconda/bin/python

import os
import sys

new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../"))
if new_path not in sys.path:
    sys.path.append(new_path)

import numpy as np
from Python_Classes4MAD import metaclass
from optparse import OptionParser
from Utilities import tfs_file_writer


def write_global_files(beam, kmod_dir, res_dir, mod_path):
    bpms = np.array([])
    betx = np.array([])
    bety = np.array([])
    betx_std = np.array([])
    bety_std = np.array([])

    IPs = ['ip1', 'ip5', 'ip8']
    for ip in IPs:
        pathx, pathy = get_paths(kmod_dir, beam, ip)
        try:
            datax = metaclass.twiss(pathx)
            datay = metaclass.twiss(pathy)
        except IOError:
            print "Cannot find " + ip + " kmod data, won't copy those files."
            continue
        bpms = np.concatenate((bpms, datax.NAME))
        betx     = np.concatenate((betx, datax.BETX))
        bety     = np.concatenate((bety, datay.BETY))
        betx_std = np.concatenate((betx_std, datax.BETXSTD))
        bety_std = np.concatenate((bety_std, datay.BETYSTD))

    betx_mod, bety_mod = get_model_beta(bpms, mod_path)

    xdata = tfs_file_writer.TfsFileWriter.open(os.path.join(res_dir, 'getkmodbetax.out'))
    xdata.set_column_width(20)
    xdata.add_column_names(['NAME', 'S', 'COUNT', 'BETX', 'STDBETX', 'BETXMDL', 'ERRBETX', 'BETXRES', 'BETXSTDRES'])
    xdata.add_column_datatypes(['%s', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])

    ydata = tfs_file_writer.TfsFileWriter.open(os.path.join(res_dir, 'getkmodbetay.out'))
    ydata.set_column_width(20)
    ydata.add_column_names(['NAME', 'S', 'COUNT', 'BETY', 'STDBETY', 'BETYMDL', 'ERRBETY', 'BETYRES', 'BETYSTDRES'])
    ydata.add_column_datatypes(['%s', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])

    for i in range(len(bpms)):
        xdata.add_table_row([bpms[i], 0, 0, betx[i], betx_std[i], betx_mod[i], 0, 0, 0 ])
        ydata.add_table_row([bpms[i], 0, 0, bety[i], bety_std[i], bety_mod[i], 0, 0, 0 ])
    xdata.write_to_file()
    ydata.write_to_file()


def get_model_beta(bpms, mod_path):
    model_twiss = metaclass.twiss(mod_path)
    names = np.array(model_twiss.NAME)
    betx =  model_twiss.BETX
    bety =  model_twiss.BETY
    betasx = []
    betasy = []
    for bpm in bpms:
        index = np.where(names==bpm)[0][0]
        betasx.append(betx[index])
        betasy.append(bety[index])
    return np.array(betasx), np.array(betasy)


def get_paths(path, beam, ip):
    pathx = os.path.join(path, ip+beam, 'getkmodbetax_%s.out' %ip)
    pathy = os.path.join(path, ip+beam, 'getkmodbetay_%s.out' %ip)
    return pathx, pathy


def parse_args():
    usage = 'Usage: %prog -w WORKING_DIR -o OUT_FILE_PATH [options]'
    parser = OptionParser(usage=usage)
    parser.add_option('-k', '--kmod_directory',
                            help='path to kmod directory with stored KMOD measurement files', 
                            action='store', type='string', dest='kmod_dir')
    parser.add_option('-r', '--results_dir',
                            help='Specify results directory of optics analysis', 
                            action='store', type='string', dest='res_dir')
    parser.add_option('-m', '--model_path',
                            help='Specify path to current model',
                            action='store', type='string', dest='mod')
    parser.add_option('-b', '--beam',
                            help='define beam used: b1 or b2',
                            action='store', type='string', dest='beam')
    (options, _) = parser.parse_args()
    return options


if __name__ == '__main__':
    options = parse_args()
    beam = options.beam
    kmod_dir = options.kmod_dir
    res_dir  = options.res_dir
    mod_path = options.mod

    write_global_files(beam, kmod_dir, res_dir, mod_path)
