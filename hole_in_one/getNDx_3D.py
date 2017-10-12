r'''
.. module: getNDx_3D

Created on 31/07/17

:author: Lukas Malina  

It computes horizontal normalised dispersion from 3D kicks, 
it uses the model with AC-Dipole element in.
'''
import sys
import os
from optparse import OptionParser
import time
import numpy as np
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Utilities import tfs_pandas as tfs
from model import manager

PI2I = 2 * np.pi * complex(0, 1)

def _parse_args():
    parser = OptionParser()
    parser.add_option("--files",
                    help="List of comma-separated sdds files from analysis",
                    metavar="FILES", dest="files")
    parser.add_option("--model",
                    help="Model with ac-dipole in",
                    metavar="MODEL", default="", dest="model")
    parser.add_option("--output",
                    help="Output folder",
                    metavar="OUTPUT", default="", dest="output")
    options, _ = parser.parse_args()

    return options.files, options.model, options.output


# Important this assumes the spectral line amplitudes scaled by the main line amplitude  
def getNDX(files,model,output):
    #getting the list of BPMs and arcBPMs
    file_list = [(file_name.strip() + ".linx") for file_name in files.strip("\"").split(",")]
    file_dict = {}
    model_tfs=tfs.read_tfs(model)
    bpms=model_tfs.loc[:,"NAME"].values
    for file_name in file_list:
        filetfs = tfs.read_tfs(file_name)
        bpms = intersect(bpms,filetfs.loc[:,"NAME"].values)
    bpms=np.array(bpms)
    arc_bpms=get_arc_bpms(model_tfs, bpms)
    
    #merging the Dataframes
    model_panda=load_panda(model)
    for i,file_name in enumerate(file_list):
        file_panda = load_panda(file_name)
        model_panda = pd.merge(model_panda, file_panda, how='inner', on='NAME', suffixes=('', str(i+1)))
    model_panda['NDXMDL'] = (model_panda.loc[:, 'DX'].values / np.sqrt(model_panda.loc[:, 'BETX'].values))
    columns=['NAME','S', 'NDXMDL']
    for c in model_panda.columns.values:
        if c.startswith('AMPZ') or c.startswith('MUZ'):
            columns.append(c)
    results = model_panda.loc[:,columns]
    results.set_index("NAME", inplace=True)
    columns=['S', 'NDXMDL']
    cols=[]
    # scaling to the model, and getting the synchrotron phase in the arcs
    for c in results.columns.values:
        if c.startswith('MUZ'):
            results['sc'+c.replace('MU','AMP')]=results.loc[:,c.replace('MU','AMP')].values * np.sum(results.loc[arc_bpms,'NDXMDL']) / np.sum(results.loc[arc_bpms,c.replace('MU','AMP')])
            results['s'+c]=np.angle(np.exp(PI2I*results.loc[:,c].values))/(2*np.pi)
            field = results.loc[:,'s'+c].values- np.angle(np.sum(np.exp(PI2I * results.loc[arc_bpms,'s'+c])))/(2*np.pi)
            results['sc'+c]=np.abs(np.where(np.abs(field)>0.5,field-np.sign(field),field))
            d='sc'+c
            cols.append(d)
    #resolving the sign of dispersion
    for c in cols:
        results[c.replace('scMUZ','fNDX')]=results.loc[:,c.replace('MU','AMP')] * np.sign(0.25 - np.abs(np.angle(np.sum(np.exp(PI2I * results.loc[:,cols]),axis=1)))/(2*np.pi))
        columns.append(c.replace('scMUZ','fNDX'))
    forfile=results.loc[:,columns]
    f=[]
    #averaging over files and error calculation
    for c in forfile.columns.values:
        if c.startswith('fNDX'):
            f.append(c)
    if len(f)>1:
        forfile['STDNDX']=np.std(forfile.loc[:,f],axis=1)*t_value_correction(len(f))
    else:
        forfile['STDNDX']=0.0
    forfile['NDX']=np.mean(forfile.loc[:,f],axis=1)
    forfile['DNDX']=forfile.loc[:,'NDX']-forfile.loc[:,'NDXMDL']
    print np.mean(forfile.loc[:,'STDNDX'])
    tfs.write_tfs(forfile,{},os.path.join(output + "getNDx.out"))
    return 

def load_panda(model):
    return pd.DataFrame(tfs.read_tfs(model))
    


def get_arc_bpms(model_twiss, bpm_names):  # twiss_ac, intersected BPM names 
    model_twiss.set_index("NAME", inplace=True)
    sequence = model_twiss.headers["SEQUENCE"].lower().replace("b1", "").replace("b2", "")
    AccelClass = manager.get_accel_class(sequence)
    arc_bpms_mask = AccelClass.get_arc_bpms_mask(bpm_names)
    arc_bpm_names = bpm_names[arc_bpms_mask]
    return arc_bpm_names


def t_value_correction(num):
    
    correction_dict = {2:1.8394733927562799, 3:1.3224035682262103, 4:1.1978046912864673, 
                       5:1.1424650980932523, 6:1.1112993008590089, 7:1.0913332519214189, 
                       8:1.0774580800762166, 9:1.0672589736833817, 10:1.0594474783177483,
                       11:1.053273802733051, 12:1.0482721313740653, 13:1.0441378866779087,
                       14:1.0406635564353071, 15:1.0377028976401199, 16:1.0351498875115406,
                       17:1.0329257912610941, 18:1.0309709166064416, 19:1.029239186837585, 
                       20:1.0276944692596461}
    if num > 1 and num <=20:
        t_factor = correction_dict[num]
    else:
        t_factor = 1
    return t_factor


def intersect(a, b):
    return list(set(a) & set(b))

if __name__ == "__main__":
    timeStartGlobal = time.time()
    _files, _model, _output = _parse_args()
    getNDX(_files, _model, _output)
    timeGlobal = time.time() - timeStartGlobal
    print "Duration:", timeGlobal, "s"
