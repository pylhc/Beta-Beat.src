import sys
import os
import shutil
import optparse
import json
import time
import datetime
import logging
import pickle
import cPickle
import numpy as np
import pandas as pd
import iotools
import madx_wrapper
sys.path.append("/afs/cern.ch/work/l/lmalina/Beta-Beat.src/")
#import __init__  # @UnusedImport init will include paths
import Utilities.tfs_pandas as tfs

#LOGGER = logging.getLogger(__name__)



DEV_NULL = os.devnull 
  
#=======================================================================================================================
#Possible problems and notes:
#                   do we need error cut, when we use error-based weights? probably not anymore 
#                   error-based weights default? likely - but be carefull with low tune errors vs svd cut in pseudoinverse
#                   manual creation of pd.DataFrame varslist, deltas? maybe
#                   tunes in tfs_pandas single value or a column?
#                   Minimal strength removed
#                   Check the output files and when they are written
#                       There should be some summation/renaming for iterations
#                   For two beam correction 
#                       The two beams can be treated separately until the calcultation of correction 
#                   Values missing in the response (i.e. correctors of the other beam) shall be treated as zeros 
#                   Missing a part that treats the output from LSA                      
#=======================================================================================================================


def global_correction(
         output_path = "./",
         accel="LHCB1",
         singular_value_cut=0.01,
         errorcut="0.02,0.05,0.02", #phase, betabeat, ndx
         modelcut="0.2,0.2,0.2",    #phase, betabeat, ndx
         beta_beat_root='/afs/cern.ch/work/l/lmalina/online_model/',
         weights_on_corrections="1,1,0,0,0,10",
         error_weights=0,
         path_to_optics_files_dir="./",
         variables=["MQM","MQT","MQTL","MQY"],#["VQ1"],#["MQM","MQT","MQTL","MQY"],
         beta_file_name = "getbeta",
         virt_flag=False,
         num_reiteration=3
         ):
    keys = ['MUX', 'MUY', 'BBX', 'BBY', 'NDX', 'Q']
    w_dict = _get_weights_dictionary(weights_on_corrections)
    e_dict = _get_error_cut_dictionary(errorcut)
    m_dict = _get_model_cut_dictionary(modelcut)
    print w_dict
    keys, response, measurement, model, varslist = _load_input_files(keys, w_dict, path_to_optics_files_dir, output_path, beta_beat_root, beta_file_name, variables, virt_flag)
    keys, measurement = _filter_measurement(measurement, model, keys, w_dict, error_weights, e_dict, m_dict)
    response = _filter_response_columns(response, measurement, keys)

    measurement = _append_model_to_measurement(model, measurement, keys)
    _print_rms(measurement,keys)    
    _dump(os.path.join(output_path, "used_measurement"), measurement)
    init_meas = measurement
    deltas = _calculate_deltas(response, measurement, keys, varslist, cut=singular_value_cut, app=0, path=output_path)
    writeparams(deltas, varslist, app=0)      
    
    for i in range(num_iteration):
        print "\n\nRunning MADX:", i+1
        _callMadx(os.path.join(output_path,"job.LHCB1.twisscor.madx")) #TODO
        shutil.copy2(os.path.join(output_path,"twiss_corr.dat"),os.path.join(output_path,"twiss_"+str(i)+".dat"))
        model = _load_model(path_to_optics_files_dir, i)
        measurement = _append_model_to_measurement(model, measurement, keys)
        _print_rms(measurement,keys)
        ideltas = _calculate_deltas(response, measurement, keys, varslist, cut=singular_value_cut, app=0, path=output_path)
        writeparams(ideltas, varslist, app=1)    
        print np.sum(np.abs(ideltas))
        deltas= deltas + ideltas
        print "Cumulative deltas:", np.sum(np.abs(deltas))
    write_knob(deltas, varslist)   
    # _handle_data_for_accel(accel)
    

#=======================================================================================================================
# helper functions
#=======================================================================================================================
#

def _print_rms(meas,keys):
    for key in keys:
        print key, " RMS: ",  np.std(meas[key].loc[:,'DIFF'].values)
        print key, " Weighted RMS: ", np.sqrt(np.average(np.square(meas[key].loc[:,'DIFF'].values), weights=np.square(meas[key].loc[:,'WEIGHT'].values)) )

def _get_weights_dictionary(weights_on_corrections):
    w = weights_on_corrections.split(',')    
    return {'MUX': float(w[0]), 'MUY': float(w[1]), 'BBX': float(w[2]), 'BBY': float(w[3]), 'NDX': float(w[4]), 'Q' : float(w[5])}

def _get_errorcut_dictionary(errorcut):
    w = errorcut.split(',')    
    return {'MUX': float(w[0]), 'MUY': float(w[0]), 'BBX': float(w[1]), 'BBY': float(w[1]), 'NDX': float(w[2]), 'Q' : float(w[0])}

def _get_modelcut_dictionary(modelcut):
    w = modelcut.split(',')    
    return {'MUX': float(w[0]), 'MUY': float(w[0]), 'BBX': float(w[1]), 'BBY': float(w[1]), 'NDX': float(w[2]), 'Q' : float(w[0])}


def _load_input_files(keys, w_dict, path_to_optics_files_dir, output_path, beta_beat_root, beta_file_name, variables, virt_flag, flag_lsa):
    print "Starting loading Full Response optics"
    #Full response is dictionary of MUX/Y, BBX/Y, NDX and Q gradients upon a change of a single quadrupole strength
    full_response = pickle.load(open(os.path.join(path_to_optics_files_dir, 'FullResponse'), 'rb'))
    print "Loading ended"
    
    measurement = {}
    if (w_dict['MUX'] or w_dict['MUY'] or w_dict['Q']):
        try:
            measurement['MUX'] = _read_tfs(output_path, "getphasex_free.out")
            measurement['MUY'] = _read_tfs(output_path, "getphasey_free.out")
        except IOError:
            measurement['MUX'] = _read_tfs(output_path, "getphasex.out")
            measurement['MUY'] = _read_tfs(output_path, "getphasey.out")
        measurement['Q'] = pd.DataFrame({ 'NAME' : pd.Categorical(['Q1','Q2']),
                                          'VALUE' : np.remainder([measurement['MUX']['Q1'], measurement['MUY']['Q2']], [1, 1]), #just fractional tunes
                                          'ERROR' : np.array([0.001,0.001])})                                      #TODO measured errors not in the file
    else:
        keys.remove('MUX')
        keys.remove('MUY')
        keys.remove('Q')

    if w_dict['BBX'] or w_dict['BBY']:
        try:
            measurement['BBX'] = _read_tfs(output_path, beta_file_name + "x_free.out")
            measurement['BBY'] = _read_tfs(output_path, beta_file_name + "y_free.out")
        except IOError:
            measurement['MUX'] = _read_tfs(output_path, beta_file_name + "x.out")
            measurement['MUY'] = _read_tfs(output_path, beta_file_name + "y.out")
    else:
        keys.remove('BBX')
        keys.remove('BBY')
    
    if w_dict['NDX']:
        try:
            measurement['NDX'] = _read_tfs(output_path, "getNDx.out")      
        except IOError:
            print "WARNING: No good dispersion or inexistent file getNDx"
            print "WARNING: Correction will not take into account NDx"
            keys.remove('NDX')
    else:
        keys.remove('NDX')
    new_keys = []
    
    for key in keys:
        if w_dict[key]:
            new_keys.append(key)

    model = _read_tfs(path_to_optics_files_dir, "twiss.dat")
    if virt_flag:
        json_file="vqs.json"
    else:
        json_file="qs1.json"
    knobsdict = json.load(file(os.path.join(beta_beat_root, json_file), 'r')) 
    varslist = []
    for var in variables:
        try:
            varslist += knobsdict[var]
        except KeyError:
            print >> sys.stderr, "Given variable({0}) not in dictionary".format(var)
    print len(varslist)
    return new_keys, full_response, measurement, model, np.array(varslist)

#def _load_measurements():

def _read_tfs(path, file_name):
    return tfs.read_tfs(os.path.join(path, file_name))


def _filter_measurement(measurement, model, keys, w_dict, w, e_dict, m_dict):    # add error and model cuts
    meas = {}
    new_keys=[]   
    for key in keys:   
        if w_dict[key]>0.0:
            meas[key] = MEASUREMENT_FILTERS[key](measurement[key], model, w_dict[key], w, modelcut=m_dict[key],errorcut=e_dict[key])
            new_keys.append(key)
    return new_keys, meas


def _get_filtered_phases(meas, model, weight, erwg, modelcut=0.05, errorcut=0.035):
    new = pd.merge(meas, model, how='inner',on='NAME',suffixes=('', 'm'))    
    second_bpm_in = np.in1d(new.loc[:,'NAME2'].values, new.loc[:,'NAME'].values)
    if 'PHYMDL' in new.columns.values:
        plane = 'Y'
    else:
        plane = 'X'
    new['PHMOD'] = new.loc[:,'PH' + plane + 'MDL'].values
    new['VALUE'] = new.loc[:,'PHASE' + plane].values
    new['ERROR'] = new.loc[:,'STDPH' + plane].values      
    good_phase = ((new.loc[:,'ERROR'].values < errorcut) & (np.abs(new.loc[:,'VALUE'].values - new.loc[:,'PHMOD'].values) < modelcut)) 
    good_bpms = good_phase & second_bpm_in
    good_bpms[-1] = False   #removing the last pair as it goes beyond the end of lattice
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = new.loc[:,'WEIGHT'].values / new.loc[:,'ERROR'].values
    #print new['WEIGHT']
    #number_of_removed_bpms = len(x.NAME) - np.sum(good_bpms)
    #if num_of_removed_bpms > 0:
    #    print "Warning: ", num_removed_bpms, "BPM pairs removed from data for not beeing in the model or having too large error deviations: ", "PHMDL", modelcut, "STDPH", errorcut, "LEN", np.sum(good_bpms)
    print "Number of BPM pairs in plane", plane, ": ", np.sum(good_bpms)    
    return new.loc[good_bpms,['NAME','NAME2','VALUE','ERROR', 'WEIGHT']]


def _get_filtered_betas(meas, model, weight, erwg, modelcut=0.2, errorcut=0.02):      #Beta-beating and its error RELATIVE as shown in GUI
    new = pd.merge(meas, model, how='inner',on='NAME',suffixes=('', 'm'))
    if 'BETYMDL' in new.columns.values:
        plane = 'Y'
    else:
        plane = 'X'
    new['BETAMOD'] = new.loc[:,'BET' + plane + 'MDL'].values
    new['VALUE'] = new.loc[:,'BET' + plane].values
    if ('STDBET'+ plane) in new.columns.values: # old files or k-mod
        new['ERROR'] = np.sqrt(np.square(new.loc[:,'ERRBET' + plane].values) + np.square(new.loc[:,'STDBET' + plane].values) ) 
    else:
        new['ERROR'] = new.loc[:,'ERRBET' + plane].values
    model_close = (np.abs(new.loc[:,'VALUE'].values - new.loc[:,'BETAMOD'].values) / new.loc[:,'BETAMOD'].values) < modelcut
    error_low = ((new.loc[:,'ERROR'].values / new.loc[:,'BETAMOD'].values) < errorcut)         
    good_bpms = (model_close & error_low)
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = new.loc[:,'WEIGHT'].values * new.loc[:,'BETAMOD'].values / new.loc[:,'ERROR'].values
    #number_of_removed_bpms = len(x.NAME) - np.sum(good_bpms)
    #if num_of_removed_bpms > 0:
    #    print "Warning: ", num_of_removed_bpms, "BPMs removed from data for having too large error"
    print "Number of BPMs with beta in plane", plane, ": ", np.sum(good_bpms)
    return new.loc[good_bpms,['NAME','VALUE','ERROR', 'WEIGHT']]


def _get_filtered_disp(meas, model, weight, erwg, modelcut=0.2, errorcut=0.02):
    new = pd.merge(meas, model, how='inner',on='NAME',suffixes=('', 'm'))
    new['VALUE'] = new.loc[:,'NDX'].values
    new['ERROR'] = new.loc[:,'STDNDX'].values   
    good_bpms = ((new.loc[:,'ERROR'].values < errorcut) & (np.abs(new.loc[:,'VALUE'].values - new.loc[:,'NDXMDL'].values) < modelcut))
    new['WEIGHT'] = weight
    if erwg:
        new['WEIGHT'] = new.loc[:,'WEIGHT'].values / new.loc[:,'ERROR'].values
    #print "Number of x BPMs", len(x.NAME)
    #number_of_removed_bpms = len(x.NAME) - np.sum(good_bpms)
    #if num_of_removed_bpms > 0:
    #    print "Warning: ", num_of_removed_bpms, "BPMs removed from data for having too large error"
    print "Number of BPMs with NDx : ", np.sum(good_bpms)
    return new.loc[good_bpms,['NAME','VALUE','ERROR', 'WEIGHT']]


def _get_tunes(meas, mod, weight, erwg, modelcut=0.1, errorcut=0.027): 
    meas['WEIGHT'] = weight
    if erwg:
        meas['WEIGHT'] = meas.loc[:,'WEIGHT'].values / meas.loc[:,'ERROR'].values 
    print "Number of tune measurements: ", len(meas.index.values)
    return meas


def _filter_response_columns(response, measurement, keys):
    resp = {}    
    for key in keys:
        resp[key] = RESPONSE_FILTERS[key](response[key],measurement[key])
    return resp

def _get_phase_response(resp,meas):
    new = resp.loc[:,meas.loc[:,'NAME'].values]
    new.sub(resp.as_matrix(columns=meas.loc[:,'NAME2'].values),axis=0)
    return -new # as we did subtraction name-name2

def _get_betabeat_response(resp,meas):
    return resp.loc[:,meas.loc[:,'NAME'].values]

def _get_disp_response(resp,meas):
    return resp.loc[:,meas.loc[:,'NAME'].values]

def _get_tune_response(resp,meas):
    return resp

def _append_model_to_measurement(model, measurement, keys):
    meas = {}    
    for key in keys:
        meas[key] = MODEL_APPENDERS[key](model, measurement[key], key)
    return meas

def _get_model_phases(model, meas, key):
    new = model.set_index('NAME')
    meas['MODEL'] = new.loc[meas.loc[:,'NAME2'].values, key].values - new.loc[meas.loc[:,'NAME'].values, key].values
    meas['DIFF'] = meas.loc[:,'VALUE'].values - meas.loc[:,'MODEL'].values
    return meas

def _get_model_betas(model, meas, key):
    new = model.set_index('NAME')
    if key == 'BBX':
        meas['MODEL'] = new.loc[meas.loc[:,'NAME'].values, 'BETX'].values
    else:
        meas['MODEL'] = new.loc[meas.loc[:,'NAME'].values, 'BETY'].values
    meas['DIFF'] = (meas.loc[:,'VALUE'].values - meas.loc[:,'MODEL'].values) / meas.loc[:,'MODEL'].values
    
    return meas

def _get_model_disp(model, meas, key):
    new = model.set_index('NAME')
    meas['MODEL'] = new.loc[meas.loc[:,'NAME'].values, 'DX'].values / np.sqrt(new.loc[meas.loc[:,'NAME'].values, 'BETX'].values)
    meas['DIFF'] = meas.loc[:,'VALUE'].values - meas.loc[:,'MODEL'].values
    return meas


def _get_model_tunes(model, meas, key): 
    meas['MODEL'] = np.remainder([model['Q1'], model['Q2']], [1,1])    #we want just fractional tunes
    meas['DIFF'] = meas.loc[:,'VALUE'].values - meas.loc[:,'MODEL'].values
    return meas

def _dump(pathToDump, content):
    dumpFile = open(pathToDump, 'wb')
    cPickle.Pickler(dumpFile, -1).dump(content)
    dumpFile.close()

'''   
def _generate_changeparameters():   #TODO think about cross-talks with iterations using madx twiss  maybe we don't need this
    [deltas, varslist] = _calculate_deltas(response, measurement, varslist, cut=_InputData.singular_value_cut, app=0, path=_InputData.output_path)
    while np.sum(np.abs(deltas) < _InputData.min_strength):
        mask = (np.abs(deltas) > _InputData.min_strength)
        varslist = varslist[mask]        
        if len(varslist) == 0:
            print >> sys.stderr, "You want to correct with too high cut on the corrector strength"
            sys.exit(1)
        [deltas, varslist] = _calculate_deltas(response, measurement, varslist, cut=_InputData.singular_value_cut, app=0, path=_InputData.output_path)
        print "Initial correctors:", len(mask), ". Current: ", len(varslist), ". Removed for being lower than:", _InputData.min_strength
    if PRINT_DEBUG:
        print deltas
'''

def _calculate_deltas(resp, meas, keys, varslist, cut=0.01, app=0, path="./"):  #TODO think about output form
    response_matrix = _filter_response_rows_and_join(resp, keys, varslist)  #return NumPy array: each row for magnet, columns for measurements
    weight_vector = _join_weights(meas, keys)
    diff_vector = _join_diffs(meas, keys)
    delta = np.dot(np.linalg.pinv(np.transpose(response_matrix * weight_vector), cut), diff_vector * weight_vector)
    print "Delta calculation: "    
    print np.std(np.abs(diff_vector * weight_vector))
    print np.std(np.abs((diff_vector * weight_vector) - np.dot(delta,response_matrix * weight_vector)))
    #a = (np.std(np.abs(diff_vector * weight_vector - np.dot(delta,response_matrix * weight_vector)))/np.std(np.abs(diff_vector * weight_vector)))**4 * np.sum(np.abs(delta))
    #print a    
    return delta


def _filter_response_rows_and_join(resp, keys, varslist):
    return pd.concat([resp[key] for key in keys], axis=1, join_axes=[pd.Index(varslist)])


def _join_diffs(meas, keys):
    return np.concatenate([meas[key].loc[:,'DIFF'].values for key in keys])
    
def _join_weights(meas, keys):
    return np.concatenate([meas[key].loc[:,'WEIGHT'].values for key in keys])

def _load_model(path_to_optics_files_dir, iteration):
    return _read_tfs(path_to_optics_files_dir, "twiss_" + str(iteration) + ".dat")


def _callMadx(pathToTemplateFile): #TODO construct the file
    return madx_wrapper.resolve_and_run_file(pathToTemplateFile, log_file=DEV_NULL) #_create_madx_script(pathToTemplateFile), log_file=DEV_NULL)
'''
def _create_madx_script(pathToTemplateFile)
with open(template_file) as textfile:
            madx_template = textfile.read()
        iqx, iqy = cls._get_full_tunes(lhc_instance)
        replace_dict = {
            "RUN": lhc_instance.MACROS_NAME,
            "OPTICS_PATH": lhc_instance.optics_file,
            "NUM_BEAM": lhc_instance.get_beam(),
            "PATH": output_path,
            "DPP": lhc_instance.dpp,
            "QMX": iqx,
            "QMY": iqy,
            "COR": "changeparameters.madx",
        }
        madx_script = madx_template % replace_dict
        return madx_script

  
def _handle_data_for_accel(accel):
    if "LHC" in accel:  # .knob should always exist to be sent to LSA!
        src = os.path.join(os.path.join(_InputData.output_path, "changeparameters.tfs"))
        dst = os.path.join(os.path.join(_InputData.output_path, "changeparameters.knob"))
        iotools.copy_item(src, dst)  # madx table
        b = metaclass.twiss(os.path.join(_InputData.output_path, "changeparameters.tfs"))
        mad_script = open(os.path.join(_InputData.output_path, "beta_match.madx"), "w")
        names = getattr(b, "NAME", [])
        delta = getattr(b, "DELTA", [])
        for i in range(len(names)):
            if cmp(delta[i], 0) == 1:
                mad_script.write(names[i] + " = " + names[i] + " + " + str(delta[i]) + ";\n")
            else:
                mad_script.write(names[i] + " = " + names[i] + " " + str(delta[i]) + ";\n")
        mad_script.close()

'''
def write_knob(deltas, varslist, path="./"):
    a = datetime.datetime.fromtimestamp(time.time())
    f = open (os.path.join(path,'changeparameters.knob'),"w")
    print >> f, "@", "PATH", "%s", path
    print >> f, "@", "DATE", "%s", a.ctime()
    print >> f, "*", "NAME", "DELTA"
    print >> f, "$", "%s", "%le"
    for i, var in enumerate(varslist):
        f.write(var+'   '+str(-deltas[i])+'\n')
    f.close()

def writeparams(deltas, variables, app=0, path="./"):
    if (app == 0):
        mode = 'w'
    if (app == 1):
        mode = 'a'
    mad_script = open(os.path.join(path, "changeparameters.madx"), mode)    
    for i, var in enumerate(variables):
        if cmp(deltas[i], 0) == 1:
            mad_script.write(var + " = " + var + " + " + str(deltas[i]) + ";\n")
        else:
            mad_script.write(var + " = " + var + " " + str(deltas[i]) + ";\n")
    mad_script.close()


MEASUREMENT_FILTERS = {'MUX': _get_filtered_phases, 'MUY': _get_filtered_phases, 'BBX': _get_filtered_betas, 'BBY': _get_filtered_betas, 'NDX': _get_filtered_disp, 'Q' : _get_tunes}
RESPONSE_FILTERS = {'MUX': _get_phase_response, 'MUY': _get_phase_response, 'BBX': _get_betabeat_response, 'BBY': _get_betabeat_response, 'NDX': _get_disp_response, 'Q' : _get_tune_response}
MODEL_APPENDERS = {'MUX': _get_model_phases, 'MUY': _get_model_phases, 'BBX': _get_model_betas, 'BBY': _get_model_betas, 'NDX': _get_model_disp, 'Q' : _get_model_tunes}

#=======================================================================================================================
# main invocation
#=======================================================================================================================
def _start():

    timeStartGlobal = time.time()

    #options = _parse_args()
    global_correction()
    '''
         output_path=options.path,
         accel=options.ACCEL,
         singular_value_cut=options.cut,
         errorcut=options.errorcut,
         modelcut=options.modelcut,
         beta_beat_root=options.rpath,
         weights_on_corrections=options.WGT,
         error_weights=options.ErrWeight,
         path_to_optics_files_dir=options.opt,
         variables=options.var,
         index_of_num_of_beams_in_gui=options.JustOneBeam,
         beta_file_name=options.betaFile
         )
    '''
    timeGlobal = time.time() - timeStartGlobal
    print "Duration:", timeGlobal, "s"

if __name__ == "__main__":
    _start()
