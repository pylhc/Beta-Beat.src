from __future__ import print_function
import sys
import os
from datetime import datetime
import argparse

sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    ".."
)))

from sdds_files import turn_by_turn_reader
from harpy import harpy
import clean
import numpy as np
import time

RAW_SUFFIX = ".raw"
CLEAN_SUFFIX = ".clean"
DO_WRITE_RAW = False
DO_WRITE_CLEAN = False 
DO_WRITE_LIN = False


def run_all(bin_file, model_path, output_dir, tune_x = 0.28, tune_y = 0.31, 
         tolerance=0.01, nat_tune_x=None, nat_tune_y=None,
         tune_z=0.0, do_write_raw_sdds=DO_WRITE_RAW, do_write_clean_sdds=DO_WRITE_CLEAN, do_write_lin=DO_WRITE_LIN):
    start_time = time.time()
    if not ((do_write_raw_sdds) or (do_write_clean_sdds) or (do_write_lin)):
        print("No output selected, no output produced!")
        return
    tbt_files = turn_by_turn_reader.read_tbt_file(bin_file)
    for tbt_file in tbt_files:
        if do_write_raw_sdds:
            print("Writing raw sdds")
            tbt_file.write_to_ascii(model_path, get_outpath_with_suffix(bin_file, output_dir, RAW_SUFFIX)) 
        if ((do_write_clean_sdds) or (do_write_lin)):        
            dict_to_write_clean_sdds={}
            for plane in ["x", "y"]:
                bpm_names = np.array(getattr(tbt_file, "bpm_names_" + plane))
                bpm_data = getattr(tbt_file, "samples_matrix_" + plane)
                bpm_names, bpm_data, bad_bpms = clean.clean(
                    bpm_names, bpm_data,
                    min_peak_to_peak=0.0001, max_peak_cut=10.0, do_resync=resync_from_date(tbt_file.date)
                )
                bpm_names, bpm_data, bpm_res, bad_bpms_svd = clean.svd_clean(
                    bpm_names, bpm_data,
                    singular_values_amount_to_keep=12,
                    single_svd_bpm_threshold=0.925
                )
                if do_write_clean_sdds:
                    dict_to_write_clean_sdds["BPMs_" + plane] = bpm_names
                    dict_to_write_clean_sdds["SAMPLES_" + plane] = bpm_data
                    print("Writing clean sdds")
                bad_bpms_fft=[]
                if do_write_lin:
                    allowed = get_allowed_length(rang=[4000, bpm_data.shape[1]])[-1]
                    bpm_data = bpm_data[:, :allowed]
                    bad_bpms_fft = harmonic_analysis(bin_file, bpm_names, bpm_data, bpm_res, model_path,
                        tune_x, tune_y, tune_z, plane, output_dir, tolerance, nat_tune_x, nat_tune_y)
                write_bad_bpms_into_file(bin_file, bad_bpms + bad_bpms_svd + bad_bpms_fft, output_dir, plane)
            if do_write_clean_sdds:
                turn_by_turn_reader.write_ascii_file(
                     model_path, get_outpath_with_suffix(bin_file, output_dir, CLEAN_SUFFIX),
                     list(dict_to_write_clean_sdds["BPMs_x"]), dict_to_write_clean_sdds["SAMPLES_x"],
                     list(dict_to_write_clean_sdds["BPMs_y"]), dict_to_write_clean_sdds["SAMPLES_y"],
                     tbt_file.date)

    print(">> Total time for file: {0}s".format(time.time() - start_time))


def harmonic_analysis(bin_file, bpm_names, bpm_data, bpm_res, model_path, tune_x, tune_y, tune_z, plane,
                      output_dir, tolerance, nat_tune_x, nat_tune_y):
    time_start = time.time()
    tunes = (float(tune_x), float(tune_y), float(tune_z))
    nattunes = (nat_tune_x, nat_tune_y)
    if nat_tune_x is None or nat_tune_y is None:
        nattunes = None
    output_file = get_outpath_with_suffix(bin_file, output_dir, ".lin" + plane)
    drivemat = harpy.init_from_matrix(
        bpm_names, bpm_data, tunes, plane.upper(),
        output_file, model_path, nattunes=nattunes,
        tolerance=tolerance,
        start_turn=0, end_turn=None, sequential=False
    )
    bpms_after_fft = []
    for i in range(len(drivemat.bpm_results)):
        drivemat.bpm_results[i].bpm_resolution = bpm_res[np.where(bpm_names==drivemat.bpm_results[i].name)[0]][0]
        bpms_after_fft.append(drivemat.bpm_results[i].name)
    bad_bpms_fft = []    
    for bpm_name in bpm_names:
        if bpm_name not in bpms_after_fft:
            bad_bpms_fft.append(bpm_name + " Could not find the main resonance")
    drivemat.write_full_results()
    print(">> Time for harmonic_analysis: {0}s".format(
        time.time() - time_start)
    )
    return bad_bpms_fft


def get_allowed_length(rang=[300, 10000], p2max=14, p3max=9, p5max=6):
    ind = np.indices((p2max, p3max, p5max))
    nums = (np.power(2, ind[0]) *
            np.power(3, ind[1]) *
            np.power(5, ind[2])).reshape(p2max * p3max * p5max)
    nums = nums[(nums > rang[0]) & (nums <= rang[1])]
    return np.sort(nums)


def resync_from_date(acqdate):
    print("Resynchronizes BPMs")
    return acqdate > datetime(2016,4,1)


def write_bad_bpms_into_file(bin_path, bad_bpms_with_reasons, output_dir, plane):
    with open(get_outpath_with_suffix(bin_path, output_dir, ".bad_bpms_" + plane), 'w') as f:
        for line in bad_bpms_with_reasons:
            f.write(line + '\n')

    
def get_outpath_with_suffix(path, output_dir, suffix):
    return os.path.join(output_dir,os.path.basename(path)+suffix)
    

def _parse_args():
    """ Parses the arguments and returns args """
    parser = argparse.ArgumentParser(description="""Takes a binary SDDS and cleans it.\n """)

    parser.add_argument("-f", "--file", help="Binary File to clean", required=True, dest="binfile")
    parser.add_argument("-m", "--model", help="Model for BPM locations", dest="model")
    parser.add_argument("-o", "--outputdir", help="Output directory", dest="outputdir")
    parser.add_argument("-r", "--write_raw", help="Writing of raw sdds file", dest="do_write_raw", action="store_true")
    parser.add_argument("-c", "--write_clean", help="Writing of clean sdds file", dest="do_write_clean", action="store_true")
    parser.add_argument("-l", "--write_lin", help="Writing of lin files", dest="do_write_lin", action="store_true")
    options = parser.parse_args()
    return options


def _run_from_console(options):
    run_all(options.binfile, options.model, options.outputdir,
            do_write_raw_sdds=options.do_write_raw,
            do_write_clean_sdds=options.do_write_clean,
            do_write_lin=options.do_write_lin)


if __name__ == "__main__":
    _options = _parse_args()
    _run_from_console(_options)
