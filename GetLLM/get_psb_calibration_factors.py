import sys
import os
from optparse import OptionParser
from collections import OrderedDict
import numpy as np
import pandas as pd

new_path = abspath(join(dirname(abspath(__file__)), pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

import tfs_files.tfs_file_writer as tfs_writer
from tfs_files import tfs_pandas


def _parse_args():
    parser = OptionParser()
    parser.add_option("-r", "--results",
                    help="results folder",
                    metavar="results",dest="results")
    parser.add_option("-o", "--output_path",
                    help="output path",
                    metavar="OUTPUT",dest="output_path")
    options, _ = parser.parse_args()
    return options.results,options.output_path

    

def get_beta_from_amplitude(beta_amplitude,plane):
    beta_file_amplitude = tfs_pandas.read_tfs(os.path.join(beta_amplitude,"getampbeta" + str(plane.lower()) + "_free.out"))
    return beta_file_amplitude


def get_beta_from_phase(beta_phase,plane):
    beta_file_phase = tfs_pandas.read_tfs(os.path.join(beta_phase,"getbeta" + str(plane.lower()) + "_free.out"))
    return beta_file_phase


def get_betas_values_with_errors(beta_file_phase,beta_file_amplitude,plane):
    beta_error_phase = "ERRBET" + str(plane.upper())
    beta_error_amplitude = "BET" + str(plane.upper()) + "STD"
    beta_model = "BET" + str(plane.upper()) + "MDL"
    beta = "BET" + str(plane.upper())
    beta_phase = getattr(beta_file_phase,beta)
    beta_model = getattr(beta_file_phase,beta_model)
    beta_amplitude = getattr(beta_file_amplitude,beta)
    beta_error_phase_values = getattr(beta_file_phase,beta_error_phase)
    beta_error_amplitude_values = getattr(beta_file_amplitude,beta_error_amplitude)
    names = beta_file_phase.NAME
    position = beta_file_phase.S
    return beta_phase,beta_error_phase_values,beta_amplitude,beta_error_amplitude_values,beta_model,names,position


def get_calibration_factors(beta_phase,beta_amplitude,beta_error_phase,beta_error_amplitude):
    beta_calibration = (beta_phase/beta_amplitude)**0.5
    beta_calibration_error = 0.25*(beta_phase/beta_amplitude**3)*beta_error_amplitude**2 +  0.25*(1/(beta_phase*beta_amplitude))*beta_error_phase**2
    return beta_calibration,beta_calibration_error


def writting_calibration_factors_to_file(name,position,beta_calibration,beta_calibration_error,file_output):
    tfs_calibration_pandas = file_output
    all_summary = pd.DataFrame(OrderedDict([('NAME', name),('S', position),('CALIBRATION',beta_calibration),('ERROR_CALIBRATION',beta_calibration_error)]))
    tfs_pandas.write_tfs(tfs_calibration_pandas,all_summary,save_index=False)


def main(results_directory,output_directory):
    planes = ["x","y"] 
    for plane in planes:
        file_output_name = "calibration_" + str(plane.lower()) + ".out"
        file_output_results = os.path.join(output_directory,file_output_name)
        beta_from_phase_file = get_beta_from_phase(results_directory,str(plane.lower())) 
        beta_from_amplitude_file = get_beta_from_amplitude(results_directory,str(plane.lower()))
        [beta_phase, beta_error_phase, beta_amplitude, beta_error_amplitude, beta_model, name, position] = get_betas_values_with_errors(beta_from_phase_file,beta_from_amplitude_file,plane)
        [beta_calibration, beta_calibration_error] = get_calibration_factors(beta_phase,beta_amplitude,beta_error_phase,beta_error_amplitude)
        writting_calibration_factors_to_file(name, position, beta_calibration, beta_calibration_error, file_output_results)
    
    
if __name__ == "__main__":
    (_file_results,_output_path) = _parse_args()
    main(_file_results,_output_path)
