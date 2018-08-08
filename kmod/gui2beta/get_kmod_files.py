#!/afs/cern.ch/work/o/omc/anaconda/bin/python

import os
import sys

new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../"))
if new_path not in sys.path:
    sys.path.append(new_path)

from Python_Classes4MAD import metaclass
from optparse import OptionParser
from tfs_files import tfs_file_writer
import re


def write_global_files(beam, kmod_dir, res_dir, mod_path):    
    pattern = re.compile(".*R[0-9]\." + beam.upper())
    model_twiss = metaclass.twiss(mod_path)
    
    ip_dir_names = []
    for path,dirnames,files in os.walk(kmod_dir):
        for dirname in dirnames:
            if pattern.match(dirname):
                ip_dir_names.append(dirname)

    for plane in ("x", "y"):
        uplane = plane.upper()
        results_writer = tfs_file_writer.TfsFileWriter.open(os.path.join(res_dir, 'getkmodbeta' + plane + '.out'))
        write_headers(results_writer, uplane)
        bpm_names, bpm_s, betas, betas_std, betas_mdl = [], [], [], [], []
        for ip_dir_name in ip_dir_names:
            path = os.path.join(kmod_dir, ip_dir_name, 'getkmodbeta' + plane + '.out')
            try:
                data = metaclass.twiss(path)
            except IOError:
                print "Cannot find kmod data in " + path + ", won't copy those files."
                continue
            print("Kmod data found in: {}".format(path))
            for bpm_name in data.NAME:
                if bpm_name in model_twiss.NAME:
                    index = data.indx[bpm_name]
                    bpm_names.append(bpm_name)
                    betas.append(getattr(data, "BET" + uplane)[index])
                    betas_std.append(getattr(data, "BET" + uplane + "STD")[index])
                    model_index = model_twiss.indx[bpm_name]
                    betas_mdl.append(getattr(model_twiss, "BET" + uplane)[model_index])
                    bpm_s.append(getattr(model_twiss, "S")[model_index])

        for i in range(len(bpm_names)):
            results_writer.add_table_row([bpm_names[i], bpm_s[i], 0, betas[i],
                                          betas_std[i], betas_mdl[i], 0, 0, 0 ])
        results_writer.write_to_file()


def write_headers(resuts_writers, uplane):
    resuts_writers.set_column_width(20)
    resuts_writers.add_column_names(
        ['NAME', 'S', 'COUNT', 'BET' + uplane, 'STDBET' + uplane,
         'BET' + uplane + 'MDL', 'ERRBET' + uplane, 'BET' + uplane + 'RES',
         'BET' + uplane + 'STDRES']
    )
    resuts_writers.add_column_datatypes(['%s', '%le', '%le', '%le', '%le',
                                         '%le', '%le', '%le', '%le'])


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
