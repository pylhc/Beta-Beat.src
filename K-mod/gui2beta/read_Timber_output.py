import __init__
import os
import numpy as np
from Timber import *
from Python_Classes4MAD import metaclass
from Utilities import tfs_file_writer

TUNE_MEAS_PRECISION = 2.5e-5


def merge_data(working_directory, IR):
    planes = ['X','Y']
    sides  = ['L','R']
    for side in sides:
        result = {'X':[],'Y':[]}
        tdatax = metaclass.twiss(os.path.join(working_directory,IR+side+'X.tfs'))
        tdatay = metaclass.twiss(os.path.join(working_directory,IR+side+'Y.tfs'))
        kdata= metaclass.twiss(os.path.join(working_directory,side+'K.tfs'))
        K, Qx, Qxrms, Qy, Qyrms = pair(tdatax,tdatay,kdata)
        
        write_tfs_files(K, Qx, Qxrms, Qy, Qyrms, working_directory, IR, side)
        
        
def write_tfs_files(K, Qx, Qxrms, Qy, Qyrms, working_directory, IR, side):
    result = tfs_file_writer.TfsFileWriter.open(os.path.join(working_directory, IR+side+'.dat'))
    result.set_column_width(20)
    result.add_column_names(['K',    'TUNEX',     'TUNEX_ERR',     'TUNEY',     'TUNEY_ERR'])
    result.add_column_datatypes(['%le', '%le', '%le', '%le', '%le'])

    for i in range(len(K)):
        result.add_table_row([K[i], Qx[i], Qxrms[i], Qy[i], Qyrms[i] ])
    result.write_to_file()  
    
        
def pair(tdatax,tdatay,kdata):
    Qx    = []
    Qxrms = []
    Qy    = []
    Qyrms = []
    K = []
    
    if len(tdatax.TIME) > len(kdata.TIME):
        step = 300
        for i in range(len(kdata.TIME)):
            if kdata.TIME[i] > tdatax.TIME[0] and kdata.TIME[i] < tdatax.TIME[len(tdatax.TIME)-1]:
                new_timex = tdatax.TIME - kdata.TIME[i]
                maskx  = (new_timex**2<step**2)
                tunex_mask= tdatax.TUNE[maskx]
                aveQx = np.average(tunex_mask)
                Qrmsx = np.std(tunex_mask)
                
                new_timey = tdatay.TIME - kdata.TIME[i]
                masky  = (new_timey**2<step**2)
                tuney_mask= tdatay.TUNE[masky]
                aveQy = np.average(tuney_mask)
                Qrmsy = np.std(tuney_mask)

                if len(tunex_mask)>0 and len(tuney_mask)>0:
                    Qx.append(aveQx)
                    Qxrms.append(Qrmsx)
                    Qy.append(aveQy)
                    Qyrms.append(Qrmsy)
                    K.append(kdata.K[i])
                
    
    
    Qxrms = np.sqrt(np.array(Qxrms) ** 2 + TUNE_MEAS_PRECISION ** 2);
    Qyrms = np.sqrt(np.array(Qyrms) ** 2 + TUNE_MEAS_PRECISION ** 2);
    return K, Qx, Qxrms, Qy, Qyrms



if __name__=='__main__':
    merge_data()
