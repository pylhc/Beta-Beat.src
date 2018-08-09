"""
 :module: twiss_to_lin

 Created on 18/02/18

 :author: Lukas Malina


This module generates the .linx/y (both on-momentum and off-momentum) from two twiss files,
for free motion and for driven motion. The twisses should contain the chromatic functions as well.

"""

from collections import OrderedDict
import numpy as np
import pandas as pd
import sys
from os.path import abspath, join, dirname, pardir
new_path = abspath(join(dirname(abspath(__file__)), pardir, pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

from tfs_files.tfs_pandas import write_tfs, read_tfs


def optics_measurement_test_files(modeldir, outpath):
    generate_lin_files(modeldir, outfile=join(outpath, 'onmom1'), dpp=0.0)
    generate_lin_files(modeldir, outfile=join(outpath, 'onmom2'), dpp=0.0)
    generate_lin_files(modeldir, outfile=join(outpath, 'onmom3'), dpp=0.0)
    generate_lin_files(modeldir, outfile=join(outpath, 'offmom1'), dpp=-4e-4)
    generate_lin_files(modeldir, outfile=join(outpath, 'offmom2'), dpp=4e-4)
    return join(outpath, 'onmom1') + ',' + join(outpath, 'onmom2') + ',' + join(outpath, 'onmom3') \
           + ',' + join(outpath, 'offmom1') + ',' + join(outpath, 'offmom2')


def generate_lin_files(model_dir, outfile='test', dpp=0.0):
    free = read_tfs(join(model_dir, 'twiss.dat'))
    driven = read_tfs(join(model_dir, 'twiss_ac.dat'))
    nattune = {"X": np.remainder(free.headers['Q1'], 1), "Y": np.remainder(free.headers['Q2'], 1)}
    tune = {"X": np.remainder(driven.headers['Q1'], 1), "Y": np.remainder(driven.headers['Q2'], 1)}
    model = pd.merge(free, driven, how='inner', on='NAME', suffixes=('_f', '_d'))
    model['S'] = model.loc[:, 'S_f']
    nbpms = len(model.index.values)
    coup = 0.01

    for plane in ['X', 'Y']:
        lin = model.loc[:, ['NAME', 'S']]
        lin['NOISE'] = np.abs(np.random.randn(nbpms) * 0.0002 + 0.0002)
        lin['AVG_NOISE'] = lin.loc[:, 'NOISE']
        lin['CO'] = dpp * 1000 * model.loc[:, 'D' + plane + '_d']  # meters to millimeters
        lin['CORMS'] = np.abs(np.random.randn(nbpms) * 0.003 + 0.003)
        lin['PK2PK'] = 0.07 * np.sqrt(model.loc[:, 'BET' + plane + '_d'])
        lin['TUNE' + plane] = tune[plane] + 1e-7 * np.random.randn(nbpms)
        lin['MU' + plane] = np.remainder(model.loc[:, 'MU' + plane + '_d']
                                         + dpp * model.loc[:, 'DMU' + plane + '_d']
                                         + 0.0001 *np.random.randn(nbpms) + np.random.rand(), 1)
        lin['AMP' + plane] = 0.03 * np.sqrt(model.loc[:, 'BET' + plane + '_d'] *
                              (1 + dpp * np.sin(2 * np.pi * model.loc[:, 'PHI' + plane + '_d'])
                              * model.loc[:, 'W' + plane + '_d']))\
                              + 0.001 *np.random.randn(nbpms)
        lin['NATTUNE' + plane] = nattune[plane] + 1e-6 * np.random.randn(nbpms)
        lin['NATMU' + plane] = np.remainder(model.loc[:, 'MU' + plane + '_f']
                                            + 0.001 * np.random.randn(nbpms) + np.random.rand(), 1)
        lin['NATAMP' + plane] = 0.0003 * np.sqrt(model.loc[:, 'BET' + plane + '_f'])\
                                + 0.00001 *np.random.randn(nbpms)

        plane_number = {"X": "1", "Y": "2"}[plane]
        header = OrderedDict()
        header["Q" + plane_number] = tune[plane]
        header["Q" + plane_number + "RMS"] = 1e-7
        header["DPP"] = dpp
        header["NATQ" + plane_number] = nattune[plane]
        header["NATQ" + plane_number + "RMS"] = 1e-6
        if plane == 'X':
            lin['PHASE01'] = np.remainder(model.loc[:, 'MUY_d'] + dpp * model.loc[:, 'DMUY_d']
                            + 0.0001 *np.random.randn(nbpms) + np.random.rand(), 1)
            lin['AMP01'] = coup * np.sqrt(model.loc[:, 'BETY_d']
                            * (1 + dpp * np.sin(model.loc[:, 'PHIY_d']) * model.loc[:, 'WY_d'])) \
                            + 0.0001 * np.random.randn(nbpms)
            write_tfs(outfile + '.linx', lin, header)
        else:
            lin['PHASE10'] = np.remainder(model.loc[:, 'MUX_d']
                            + 0.0001 *np.random.randn(nbpms) + np.random.rand(), 1)
            lin['AMP10'] = coup * np.sqrt(model.loc[:, 'BETX_d'] \
                            * (1 + dpp * np.sin(model.loc[:, 'PHIX_d']) * model.loc[:, 'WX_d'])) \
                            + 0.0001 * np.random.randn(nbpms)
            write_tfs(outfile + '.liny', lin, header)


if __name__ == '__main__':
    optics_measurement_test_files('./', './')
