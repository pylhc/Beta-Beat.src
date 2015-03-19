'''
This script calculates the bet_deviations.npy binary file to improve the betas calculation.

@author: alangner, jcoellod
'''

import __init__  # @UnusedImport
import os
import multiprocessing
import time
import numpy as np
from optparse import OptionParser
from Utilities import iotools
from Python_Classes4MAD import madxrunner, metaclass
import math
import sys
import json

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
LHC_UNCERTAINTIES_FILE = os.path.abspath(os.path.join(CURRENT_PATH,
                                                      "..", "..", "MODEL", "LHCB_II", "model", "uncertainties.json"))

B2_ERRORS_TEMPLATE = \
'!!! %(NAME)s\n\
select, flag=error, clear; select, flag=error, %(SELECT)s;\n\
B2r = %(UNCERTAINTY)s;\n\
exec SetEfcomp_Q;\n'

TRANS_ALING_ERROR_TEMPLATE = \
"!!! %(NAME)s\n\
select, flag=error, clear; SELECT,FLAG=ERROR, %(SELECT)s; \n\
EALIGN, DX:=%(UNCERTAINTY)s*TGAUSS(GCUTR), DY:=%(UNCERTAINTY)s*TGAUSS(GCUTR);\n"

LONG_ALING_ERROR_TEMPLATE = \
"!!! %(NAME)s\n\
select, flag=error, clear; SELECT,FLAG=ERROR, %(SELECT)s; \n\
EALIGN, DS:=%(UNCERTAINTY)s*TGAUSS(GCUTR);\n"

NUM_SIMULATIONS = 1000  # Default number of simulations, used if no input argument is found
NUM_PROCESSES = multiprocessing.cpu_count()  # Default number of processes to use in simulations

LHC_RUNI_BASE_SEQ = \
              'call, file = "db5/as-built/V6.5.seq";\n' \
            + 'call, file = "db5/install_additional_elements.madx";'

BB_PATH = os.path.abspath(os.path.join(CURRENT_PATH, "..", "..", ))
LHC_RUNII_BASE_SEQ = \
              'call, file = "' + os.path.join(BB_PATH, "MODEL", "LHCB_II", "model", "base_sequence.madx") + '";'


def _parse_args():
    parser = OptionParser()
    parser.add_option("-m", "--model",
                    help="Model twiss file to use, the modifiers.madx file is assumed to be in the same directory.",
                    metavar="model", dest="model_twiss")
    parser.add_option("-o", "--output",
                    help="Output directory for the results.",
                    metavar="output", default="", dest="output_dir")
    parser.add_option("--error-tables",
                      help="Where to find the error tables, if not specified, will try to use the templates.",
                      metavar="ERROR", dest="errors_path")
    parser.add_option("-n", "--numsim",
                    help="Number of simulations to run",
                    metavar="numsim", default=NUM_SIMULATIONS, dest="num_simulations")
    parser.add_option("-p", "--processes",
                    help="Number of parallel processes to use in the simulation.",
                    metavar="numproc", default=NUM_PROCESSES, dest="num_processes")
    parser.add_option("-a", "--accelerator",
                    help="Accelerator to use: either LHCB1, LHCB2 or ALBA.",
                    metavar="ACCEL", dest="accelerator")
    parser.add_option("-x", "--tunex",
                    help="Horizontal tune.",
                    metavar="TUNEX", dest="tunex")
    parser.add_option("-y", "--tuney",
                    help="Vertical tune.",
                    metavar="TUNEY", dest="tuney")
    parser.add_option("-e", "--energy",
                    help="The energy of the beam.",
                    metavar="ENERGY", dest="energy")
    options, _ = parser.parse_args()

    if options.output_dir == "":
        options.output_dir = os.path.dirname(options.model_twiss)
    if options.accelerator is None or options.accelerator.upper() not in ["LHCB1", "LHCB2", "LHCB1RUNII", "LHCB2RUNII", "ALBA"]:
        print >> sys.stderr, "Accelerator sequence must be defined, it must be LHCB1, LHCB2 or ALBA."
        sys.exit(-1)
    if options.accelerator.startswith("LHCB"):
        if options.tunex is None or options.tuney is None or options.energy is None:
            print >> sys.stderr, "Vertical and horizontal tunes and energy must be defined."
            sys.exit(-1)
    return options.model_twiss, int(options.num_simulations), int(options.num_processes), options.output_dir, options.errors_path, options.accelerator, options.energy, options.tunex, options.tuney


def get_systematic_errors(model_twiss, num_simulations, num_processes, output_dir, errors_path, accelerator, energy, tunex, tuney):
    print "Started systematic error script, using " + str(num_processes) + " processes for " + str(num_simulations) + " simulations"
    model_dir_path = os.path.dirname(model_twiss)
    run_data_path = os.path.join(output_dir, "RUN_DATA")
    if accelerator.upper().startswith("LHCB"):
        if errors_path is None:
            beam = accelerator.upper().replace("LHC", "").replace("RUNII", "")
            errors_path = os.path.join(CURRENT_PATH, "..", "..", "MODEL", "LHC" + beam, "dipole_b2_errors")
    elif accelerator.upper() == "ALBA":
        if not errors_path is None:
            print >> sys.stdout, "ALBA doesn't need error tables, ignoring input"
    else:
        print >> sys.stderr, "No error table templates available for", accelerator, "specify an error tables path (--error-tables option)"
        sys.exit(-1)
    iotools.create_dirs(run_data_path)

    pool = multiprocessing.Pool(processes=num_processes)

    print "Running simulations..."
    start_time = time.time()
    times = _run_parallel_simulations(run_data_path, model_dir_path, num_simulations, errors_path, accelerator, energy, tunex, tuney, pool)
    _show_time_statistics(times)
    end_time = time.time()
    print "Done (" + str(end_time - start_time) + " seconds)\n"

    print "Calculating systematic error bars..."
    start_time = time.time()
    _parallel_get_systematic_errors_binary_file(model_twiss, run_data_path, output_dir, pool, num_simulations, num_processes)
    end_time = time.time()
    print "Done (" + str(end_time - start_time) + " seconds)\n"

    print "Cleaning output directory..."
    _create_summary(model_twiss, num_simulations, num_processes, output_dir, errors_path, accelerator, energy, tunex, tuney)
    iotools.delete_item(run_data_path)
    print "All done."


def _run_parallel_simulations(run_data_path, model_dir_path, num_simulations, errors_path, accelerator, energy, tunex, tuney, pool):
    times = []
    uncertainties = None
    if accelerator.upper() in ["LHCB1", "LHCB2", "LHCB1RUNII", "LHCB2RUNII"]:
        print "Using uncertainties file: ", LHC_UNCERTAINTIES_FILE
        uncertainties_json = json.load(open(LHC_UNCERTAINTIES_FILE, "r"))
        uncertainties = _extract_uncertainties(uncertainties_json, energy)
        simulation_function = _run_single_lhc_madx_simulation
    elif accelerator.upper() == "ALBA":
        simulation_function = _run_single_alba_madx_simulation
    else:
        print >> sys.stderr, "Accelerator", accelerator, "not yet implemented"
        sys.exit(-1)
    args = [(seed, run_data_path, model_dir_path, errors_path, accelerator, energy, tunex, tuney, uncertainties) for seed in range(1, num_simulations + 1)]
    tasks = pool.map_async(simulation_function, args, callback=times.append)
    tasks.wait()
    return times


def _extract_uncertainties(uncertainties_json, energy):
    if not energy in uncertainties_json:
        print >> sys.stderr, "The selected energy isn't specified in the uncertainties file."
    uncertainties_json = uncertainties_json[energy]
    uncertainties = ""
    for uncertainty_type_name, uncertainty_type in uncertainties_json.iteritems():
        print "Adding", uncertainty_type_name, "..."
        uncertainties += "!!! " + uncertainty_type_name + " !!!\n\n"
        if uncertainty_type_name == "B2_errors":
            template = B2_ERRORS_TEMPLATE
        elif uncertainty_type_name == "Alignment_longitudinal":
            template = LONG_ALING_ERROR_TEMPLATE
        elif uncertainty_type_name == "Alignment_transverse":
            template = TRANS_ALING_ERROR_TEMPLATE
        for name, error in uncertainty_type.iteritems():
            select = error["select"]
            uncertainty = error["uncertainty"]
            uncertainties += template % dict(NAME=name, SELECT=select, UNCERTAINTY=uncertainty) + "\n"
    print "Done adding uncertainties."
    return uncertainties


def _run_single_lhc_madx_simulation(seed_path_tuple):
    seed = seed_path_tuple[0]
    run_data_path = seed_path_tuple[1]
    model_dir_path = seed_path_tuple[2]
    errors_path = seed_path_tuple[3]
    accelerator = seed_path_tuple[4]
    energy = seed_path_tuple[5]
    tunex = seed_path_tuple[6]
    tuney = seed_path_tuple[7]
    uncertainties = seed_path_tuple[8]
    dict_for_replacing = dict(
            BB_PATH=BB_PATH,
            BASE_SEQ=LHC_RUNII_BASE_SEQ if accelerator.upper().endswith("RUNII") else LHC_RUNI_BASE_SEQ,
            ERR_NUM=str((seed % 60) + 1).zfill(4),
            UNCERTAINTIES=uncertainties,
            SEED=str(seed),
            RUN_DATA_PATH=run_data_path,
            ERRORS_PATH=errors_path,
            ACCEL=accelerator.upper().replace("RUNII", ""),
            BEAM=accelerator.upper().replace("LHC", "").replace("RUNII", ""),

            PATH=model_dir_path,
            QMX=tunex,
            QMY=tuney,
            ENERGY=energy
            )
    raw_madx_mask = iotools.read_all_lines_in_textfile(os.path.join(CURRENT_PATH, 'job.systematic.LHC.mask'))
    madx_job = raw_madx_mask % dict_for_replacing

    start_time = time.time()
    madxrunner.runForInputString(madx_job, stdout=open(os.devnull, "w"))
    end_time = time.time()
    return end_time - start_time


def _run_single_alba_madx_simulation(seed_path_tuple):
    seed = seed_path_tuple[0]
    run_data_path = seed_path_tuple[1]
    dict_for_replacing = dict(
            ALBA_MODEL=os.path.join(CURRENT_PATH, "..", "..", "MODEL", "ALBA"),
            SEED=str(seed),
            SEED1=str(1 + seed * 5),
            SEED2=str(2 + seed * 5),
            SEED3=str(3 + seed * 5),
            SEED4=str(4 + seed * 5),
            SEED5=str(5 + seed * 5),
            RUN_DATA_PATH=run_data_path)
    raw_madx_mask = iotools.read_all_lines_in_textfile(os.path.join(CURRENT_PATH, 'job.systematic.ALBA.mask'))
    madx_job = raw_madx_mask % dict_for_replacing

    start_time = time.time()
    madxrunner.runForInputString(madx_job, stdout=open(os.devnull, "w"))
    end_time = time.time()
    return end_time - start_time


def _show_time_statistics(times):
    times = times[0]
    average_time = 0
    max_time = 0
    min_time = sys.float_info.max
    for time in times:
        average_time += time
        if time > max_time:
            max_time = time
        if time < min_time:
            min_time = time
    average_time = average_time / len(times)
    print "Average simulation time: " + str(average_time) + " seconds"
    print "Max simulation time: " + str(max_time) + " seconds"
    print "Min simulation time: " + str(min_time) + " seconds"


def _parallel_get_systematic_errors_binary_file(model_twiss_path, run_data_path, output_path, pool, num_simulations, num_processes):
    model_twiss = metaclass.twiss(model_twiss_path)

    list_of_bpm = []
    for i in model_twiss.NAME:
        if "BPM" in i and i not in list_of_bpm:
            list_of_bpm.append(i)

    final_list = [{}, {}, 0]  # This list is used to pass the parameters "by reference" to the callback method

    chunksize = _compute_chunksize(num_simulations, num_processes)
    print "Using chunk size:", chunksize
    args = [(run_data_path, model_twiss, list_of_bpm, range(1, num_simulations + 1)[x:x + chunksize]) for x in xrange(0, num_simulations, chunksize)]
    for run_data_path, model_twiss, list_of_bpm, sim_nums in args:
        pool.apply_async(_get_error_bar_for_single_simulation, (run_data_path, model_twiss, list_of_bpm, sim_nums),
                         callback=lambda result: _merge_betas_dict(final_list, result))
    pool.close()
    pool.join()

    sys.stdout.write("\n")

    final_beta_hor, final_beta_ver, final_num_valid_data = final_list

    for key, value in final_beta_hor.iteritems():
        final_beta_hor[key] = value / final_num_valid_data
    for key, value in final_beta_ver.iteritems():
        final_beta_ver[key] = value / final_num_valid_data

    np.save(os.path.join(output_path, 'bet_deviations'), [final_beta_hor, final_beta_ver])


def _compute_chunksize(num_simulations, num_processes):
    if num_simulations <= num_processes:
        return 1
    else:
        return int(num_simulations / num_processes)


def _get_error_bar_for_single_simulation(run_data_path, model_twiss, list_of_bpm, sim_nums):
    final_list = [{}, {}, 0]
    for current_sim_index in range(len(sim_nums)):
        sim_num = sim_nums[current_sim_index]
        beta_hor = {}
        beta_ver = {}
        num_valid_data = 0
        BPM_RANGE = 13
        BPM_EACH_SIDE = int((BPM_RANGE - 1) / 2.)
        left = range(-1, -1 * (BPM_EACH_SIDE + 1), -1)
        right = range(1, BPM_EACH_SIDE + 1)
        left_comb = [(x, y) for x in left for y in left if x < y]
        right_comb = [(x, y) for x in right for y in right if x < y]
        mid_comb = [(x, y) for x in left for y in right]
        all_comb = left_comb + mid_comb + right_comb
        try:
            error_twiss = metaclass.twiss(os.path.join(run_data_path, 'twiss' + str(sim_num) + '.dat'))
            for probed_bpm in range(len(list_of_bpm)):
                err_betx = error_twiss.BETX[error_twiss.indx[list_of_bpm[probed_bpm]]]
                err_bety = error_twiss.BETY[error_twiss.indx[list_of_bpm[probed_bpm]]]
                probed_bpm_name = list_of_bpm[(probed_bpm) % len(list_of_bpm)]
                terms = {}
                for term in all_comb:
                    if term in left_comb:
                        s1_hor = (BetaFromPhase_BPM_right(probed_bpm_name,
                                                               list_of_bpm[(probed_bpm + term[0]) % len(list_of_bpm)],
                                                               list_of_bpm[(probed_bpm + term[1]) % len(list_of_bpm)],
                                                               model_twiss, error_twiss, 'H') - err_betx) / err_betx
                        s1_ver = (BetaFromPhase_BPM_right(probed_bpm_name,
                                                               list_of_bpm[(probed_bpm + term[0]) % len(list_of_bpm)],
                                                               list_of_bpm[(probed_bpm + term[1]) % len(list_of_bpm)],
                                                               model_twiss, error_twiss, 'V') - err_bety) / err_bety
                    elif term in mid_comb:
                        s1_hor = (BetaFromPhase_BPM_mid(probed_bpm_name,
                                                             list_of_bpm[(probed_bpm + term[0]) % len(list_of_bpm)],
                                                             list_of_bpm[(probed_bpm + term[1]) % len(list_of_bpm)],
                                                             model_twiss, error_twiss, 'H') - err_betx) / err_betx
                        s1_ver = (BetaFromPhase_BPM_mid(probed_bpm_name,
                                                             list_of_bpm[(probed_bpm + term[0]) % len(list_of_bpm)],
                                                             list_of_bpm[(probed_bpm + term[1]) % len(list_of_bpm)],
                                                             model_twiss, error_twiss, 'V') - err_bety) / err_bety
                    elif term in right_comb:
                        s1_hor = (BetaFromPhase_BPM_left(probed_bpm_name,
                                                              list_of_bpm[(probed_bpm + term[0]) % len(list_of_bpm)],
                                                              list_of_bpm[(probed_bpm + term[1]) % len(list_of_bpm)],
                                                              model_twiss, error_twiss, 'H') - err_betx) / err_betx
                        s1_ver = (BetaFromPhase_BPM_left(probed_bpm_name,
                                                              list_of_bpm[(probed_bpm + term[0]) % len(list_of_bpm)],
                                                              list_of_bpm[(probed_bpm + term[1]) % len(list_of_bpm)],
                                                              model_twiss, error_twiss, 'V') - err_bety) / err_bety
                    terms[term] = (s1_hor, s1_ver)
                for term1 in all_comb:
                    t1_index = all_comb.index(term1)
                    for term2 in all_comb[t1_index:]:
                        s1_hor, s1_ver = terms[term1]
                        s2_hor, s2_ver = terms[term2]
                        beta_hor[probed_bpm_name +
                                 list_of_bpm[(probed_bpm + term1[0]) % len(list_of_bpm)] +
                                 list_of_bpm[(probed_bpm + term1[1]) % len(list_of_bpm)] +
                                 list_of_bpm[(probed_bpm + term2[0]) % len(list_of_bpm)] +
                                 list_of_bpm[(probed_bpm + term2[1]) % len(list_of_bpm)]] = s1_hor * s2_hor
                        beta_ver[probed_bpm_name +
                                 list_of_bpm[(probed_bpm + term1[0]) % len(list_of_bpm)] +
                                 list_of_bpm[(probed_bpm + term1[1]) % len(list_of_bpm)] +
                                 list_of_bpm[(probed_bpm + term2[0]) % len(list_of_bpm)] +
                                 list_of_bpm[(probed_bpm + term2[1]) % len(list_of_bpm)]] = s1_ver * s2_ver
            num_valid_data = num_valid_data + 1
        except:
            pass
        if len(sim_nums) == 1:
            return beta_hor, beta_ver, num_valid_data
        else:
            _merge_betas_dict(final_list, (beta_hor, beta_ver, num_valid_data), False)
    return final_list


def _merge_betas_dict(final_list, result, print_num=True):
    beta_hor, beta_ver, num_valid_data = result
    for key, value in beta_hor.iteritems():
        if not key in final_list[0]:
            final_list[0][key] = value
        else:
            final_list[0][key] += value
    for key, value in beta_ver.iteritems():
        if not key in final_list[1]:
            final_list[1][key] = value
        else:
            final_list[1][key] += value
    final_list[2] += num_valid_data
    if print_num:
        sys.stdout.write("Done for: " + str(final_list[2]) + " files.\r")
        sys.stdout.flush()


def _TEST_compare_with_template(dict_list):
    template = np.load(os.path.join(CURRENT_PATH, "test_bet_deviations_good.npy"))
    t_beta_hor = template[0]
    t_beta_ver = template[1]
    beta_hor = dict_list[0]
    beta_ver = dict_list[1]
    n = 0
    for key, value in beta_hor.iteritems():
        if abs(t_beta_hor[key] - value) > 0.000000000001:
            n += 1
            print "Error for hor values: "
            print repr(t_beta_hor[key])
            print repr(value)
            print abs(t_beta_hor[key] - value)
    for key, value in beta_ver.iteritems():
        if abs(t_beta_ver[key] - value) > 0.000000000001:
            print "Error for ver values: "
            print repr(t_beta_hor[key])
            print repr(value)
            print abs(t_beta_ver[key] - value)
    print n


def BetaFromPhase_BPM_left(bn1, bn3, bn2, MADTwiss, ERRTwiss, plane):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is left of the other two BPMs
    :Parameters:
        'bn1':string
            Name of probed BPM
        'bn2':string
            Name of first BPM right of the probed BPM
        'bn3':string
            Name of second BPM right of the probed BPM
        'MADTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'betstd':float
            calculated error on beta function at probed BPM
        'alf':float
            calculated alfa function at probed BPM
        'alfstd':float
            calculated error on alfa function at probed BPM
    '''
    bn1_err_index = ERRTwiss.indx[bn1]
    bn2_err_index = ERRTwiss.indx[bn2]
    bn3_err_index = ERRTwiss.indx[bn3]
    bn1_mdl_index = MADTwiss.indx[bn1]
    bn2_mdl_index = MADTwiss.indx[bn2]
    bn3_mdl_index = MADTwiss.indx[bn3]
    if plane == 'H':
        ph2pi12 = 2. * np.pi * (ERRTwiss.MUX[bn2_err_index] - ERRTwiss.MUX[bn1_err_index]) % ERRTwiss.Q1
        ph2pi13 = 2. * np.pi * (ERRTwiss.MUX[bn3_err_index] - ERRTwiss.MUX[bn1_err_index]) % ERRTwiss.Q1
        phmdl12 = 2. * np.pi * (MADTwiss.MUX[bn2_mdl_index] - MADTwiss.MUX[bn1_mdl_index]) % MADTwiss.Q1
        phmdl13 = 2. * np.pi * (MADTwiss.MUX[bn3_mdl_index] - MADTwiss.MUX[bn1_mdl_index]) % MADTwiss.Q1

        betmdl1 = MADTwiss.BETX[bn1_mdl_index]
        betmdl2 = MADTwiss.BETX[bn2_mdl_index]
        betmdl3 = MADTwiss.BETX[bn3_mdl_index]
        alpmdl1 = MADTwiss.ALFX[bn1_mdl_index]
    elif plane == 'V':
        ph2pi12 = 2. * np.pi * (ERRTwiss.MUY[bn2_err_index] - ERRTwiss.MUY[bn1_err_index]) % ERRTwiss.Q2
        ph2pi13 = 2. * np.pi * (ERRTwiss.MUY[bn3_err_index] - ERRTwiss.MUY[bn1_err_index]) % ERRTwiss.Q2
        phmdl12 = 2. * np.pi * (MADTwiss.MUY[bn2_mdl_index] - MADTwiss.MUY[bn1_mdl_index]) % MADTwiss.Q2
        phmdl13 = 2. * np.pi * (MADTwiss.MUY[bn3_mdl_index] - MADTwiss.MUY[bn1_mdl_index]) % MADTwiss.Q2

        betmdl1 = MADTwiss.BETY[bn1_mdl_index]
        betmdl2 = MADTwiss.BETY[bn2_mdl_index]
        betmdl3 = MADTwiss.BETY[bn3_mdl_index]
        alpmdl1 = MADTwiss.ALFY[bn1_mdl_index]
    if betmdl3 < 0 or betmdl2 < 0 or betmdl1 < 0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)

    # Find beta1 and alpha1 from phases assuming model transfer matrix
    # Matrix M: BPM1-> BPM2
    # Matrix N: BPM1-> BPM3
    M11 = math.sqrt(betmdl2 / betmdl1) * (np.cos(phmdl12) + alpmdl1 * np.sin(phmdl12))
    M12 = math.sqrt(betmdl1 * betmdl2) * np.sin(phmdl12)
    N11 = math.sqrt(betmdl3 / betmdl1) * (np.cos(phmdl13) + alpmdl1 * np.sin(phmdl13))
    N12 = math.sqrt(betmdl1 * betmdl3) * np.sin(phmdl13)

    denom = M11 / M12 - N11 / N12 + 1e-16
    numer = 1 / np.tan(ph2pi12) - 1 / np.tan(ph2pi13)
    bet = numer / denom

    return bet


def BetaFromPhase_BPM_mid(bn2, bn1, bn3, MADTwiss, ERRTwiss, plane):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is between the other two BPMs
    :Parameters:
        'bn1':string
            Name of BPM left of the probed BPM
        'bn2':string
            Name of probed BPM
        'bn3':string
            Name of BPM right of the probed BPM
        'MADTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'betstd':float
            calculated error on beta function at probed BPM
        'alf':float
            calculated alfa function at probed BPM
        'alfstd':float
            calculated error on alfa function at probed BPM
    '''
    bn1_err_index = ERRTwiss.indx[bn1]
    bn2_err_index = ERRTwiss.indx[bn2]
    bn3_err_index = ERRTwiss.indx[bn3]
    bn1_mdl_index = MADTwiss.indx[bn1]
    bn2_mdl_index = MADTwiss.indx[bn2]
    bn3_mdl_index = MADTwiss.indx[bn3]
    if plane == 'H':
        ph2pi12 = 2. * np.pi * (ERRTwiss.MUX[bn2_err_index] - ERRTwiss.MUX[bn1_err_index]) % ERRTwiss.Q1
        ph2pi23 = 2. * np.pi * (ERRTwiss.MUX[bn3_err_index] - ERRTwiss.MUX[bn2_err_index]) % ERRTwiss.Q1
        phmdl12 = 2. * np.pi * (MADTwiss.MUX[bn2_mdl_index] - MADTwiss.MUX[bn1_mdl_index]) % MADTwiss.Q1
        phmdl23 = 2. * np.pi * (MADTwiss.MUX[bn3_mdl_index] - MADTwiss.MUX[bn2_mdl_index]) % MADTwiss.Q1

        betmdl1 = MADTwiss.BETX[bn1_mdl_index]
        betmdl2 = MADTwiss.BETX[bn2_mdl_index]
        betmdl3 = MADTwiss.BETX[bn3_mdl_index]
        alpmdl2 = MADTwiss.ALFX[bn2_mdl_index]
    elif plane == 'V':
        ph2pi12 = 2. * np.pi * (ERRTwiss.MUY[bn2_err_index] - ERRTwiss.MUY[bn1_err_index]) % ERRTwiss.Q2
        ph2pi23 = 2. * np.pi * (ERRTwiss.MUY[bn3_err_index] - ERRTwiss.MUY[bn2_err_index]) % ERRTwiss.Q2
        phmdl12 = 2. * np.pi * (MADTwiss.MUY[bn2_mdl_index] - MADTwiss.MUY[bn1_mdl_index]) % MADTwiss.Q2
        phmdl23 = 2. * np.pi * (MADTwiss.MUY[bn3_mdl_index] - MADTwiss.MUY[bn2_mdl_index]) % MADTwiss.Q2

        betmdl1 = MADTwiss.BETY[bn1_mdl_index]
        betmdl2 = MADTwiss.BETY[bn2_mdl_index]
        betmdl3 = MADTwiss.BETY[bn3_mdl_index]
        alpmdl2 = MADTwiss.ALFY[bn2_mdl_index]
    if betmdl3 < 0 or betmdl2 < 0 or betmdl1 < 0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)

    # Find beta2 and alpha2 from phases assuming model transfer matrix
    # Matrix M: BPM1-> BPM2
    # Matrix N: BPM2-> BPM3
    M22 = math.sqrt(betmdl1 / betmdl2) * (np.cos(phmdl12) - alpmdl2 * np.sin(phmdl12))
    M12 = math.sqrt(betmdl1 * betmdl2) * np.sin(phmdl12)
    N11 = math.sqrt(betmdl3 / betmdl2) * (np.cos(phmdl23) + alpmdl2 * np.sin(phmdl23))
    N12 = math.sqrt(betmdl2 * betmdl3) * np.sin(phmdl23)

    denom = M22 / M12 + N11 / N12 + 1e-16
    numer = 1 / np.tan(ph2pi12) + 1 / np.tan(ph2pi23)
    bet = numer / denom

    return bet


def BetaFromPhase_BPM_right(bn3, bn1, bn2, MADTwiss, ERRTwiss, plane):
    '''
    Calculates the beta/alfa function and their errors using the
    phase advance between three BPMs for the case that the probed BPM is right the other two BPMs
    :Parameters:
        'bn1':string
            Name of second BPM left of the probed BPM
        'bn2':string
            Name of first BPM left of the probed BPM
        'bn3':string
            Name of probed BPM
        'MADTwiss':twiss
            model twiss file
        'phase':dict
            measured phase advances
        'plane':string
            'H' or 'V'
    :Return:tupel (bet,betstd,alf,alfstd)
        'bet':float
            calculated beta function at probed BPM
        'betstd':float
            calculated error on beta function at probed BPM
        'alf':float
            calculated alfa function at probed BPM
        'alfstd':float
            calculated error on alfa function at probed BPM
    '''
    bn1_err_index = ERRTwiss.indx[bn1]
    bn2_err_index = ERRTwiss.indx[bn2]
    bn3_err_index = ERRTwiss.indx[bn3]
    bn1_mdl_index = MADTwiss.indx[bn1]
    bn2_mdl_index = MADTwiss.indx[bn2]
    bn3_mdl_index = MADTwiss.indx[bn3]
    if plane == 'H':
        ph2pi23 = 2. * np.pi * (ERRTwiss.MUX[bn3_err_index] - ERRTwiss.MUX[bn2_err_index]) % ERRTwiss.Q1
        ph2pi13 = 2. * np.pi * (ERRTwiss.MUX[bn3_err_index] - ERRTwiss.MUX[bn1_err_index]) % ERRTwiss.Q1
        phmdl23 = 2. * np.pi * (MADTwiss.MUX[bn3_mdl_index] - MADTwiss.MUX[bn2_mdl_index]) % MADTwiss.Q1
        phmdl13 = 2. * np.pi * (MADTwiss.MUX[bn3_mdl_index] - MADTwiss.MUX[bn1_mdl_index]) % MADTwiss.Q1

        betmdl1 = MADTwiss.BETX[bn1_mdl_index]
        betmdl2 = MADTwiss.BETX[bn2_mdl_index]
        betmdl3 = MADTwiss.BETX[bn3_mdl_index]
        alpmdl3 = MADTwiss.ALFX[bn3_mdl_index]
    elif plane == 'V':
        ph2pi23 = 2. * np.pi * (ERRTwiss.MUY[bn3_err_index] - ERRTwiss.MUY[bn2_err_index]) % ERRTwiss.Q2
        ph2pi13 = 2. * np.pi * (ERRTwiss.MUY[bn3_err_index] - ERRTwiss.MUY[bn1_err_index]) % ERRTwiss.Q2
        phmdl23 = 2. * np.pi * (MADTwiss.MUY[bn3_mdl_index] - MADTwiss.MUY[bn2_mdl_index]) % MADTwiss.Q2
        phmdl13 = 2. * np.pi * (MADTwiss.MUY[bn3_mdl_index] - MADTwiss.MUY[bn1_mdl_index]) % MADTwiss.Q2

        betmdl1 = MADTwiss.BETY[bn1_mdl_index]
        betmdl2 = MADTwiss.BETY[bn2_mdl_index]
        betmdl3 = MADTwiss.BETY[bn3_mdl_index]
        alpmdl3 = MADTwiss.ALFY[bn3_mdl_index]
    if betmdl3 < 0 or betmdl2 < 0 or betmdl1 < 0:
        print >> sys.stderr, "Some of the off-momentum betas are negative, change the dpp unit"
        sys.exit(1)

    # Find beta3 and alpha3 from phases assuming model transfer matrix
    # Matrix M: BPM2-> BPM3
    # Matrix N: BPM1-> BPM3
    M22 = math.sqrt(betmdl2 / betmdl3) * (np.cos(phmdl23) - alpmdl3 * np.sin(phmdl23))
    M12 = math.sqrt(betmdl2 * betmdl3) * np.sin(phmdl23)
    N22 = math.sqrt(betmdl1 / betmdl3) * (np.cos(phmdl13) - alpmdl3 * np.sin(phmdl13))
    N12 = math.sqrt(betmdl1 * betmdl3) * np.sin(phmdl13)

    denom = M22 / M12 - N22 / N12 + 1e-16
    numer = 1 / np.tan(ph2pi23) - 1 / np.tan(ph2pi13)
    bet = numer / denom

    return bet


def _create_summary(model_twiss, num_simulations, num_processes, output_dir, errors_path, accelerator, energy, tunex, tuney):
    with open(os.path.join(output_dir, "bet_deviation_summary.txt"), "w") as summary_file:
        summary_file.write("*** Systematic errors computed at: " + time.strftime("%Y-%m-%d %H:%M:%S") + " ***\n")
        summary_file.write("\n")
        summary_file.write("Accelerator: " + accelerator + "\n")
        if accelerator.lower().startswith("lhc"):
            summary_file.write("Energy: " + energy + "\n")
            summary_file.write("Tune X: " + tunex + "\n")
            summary_file.write("Tune Y: " + tuney + "\n")
            if errors_path is None:
                summary_file.write("Error tables from templates.\n")
            else:
                summary_file.write("Error tables from: \"" + os.path.abspath(errors_path) + "\"\n")
            modifiers_file_path = os.path.join(os.path.dirname(model_twiss), "modifiers.madx")
            if os.path.exists(modifiers_file_path):
                summary_file.write("Optics: \n")
                with open(modifiers_file_path, "r") as modifiers_file:
                    for modifier_line in modifiers_file:
                        summary_file.write("    " + modifier_line)
        summary_file.write("\n")
        summary_file.write("Computed using " + str(num_simulations) + " simulations." + "\n")


if __name__ == "__main__":
    _model_twiss, _num_simulations, _num_processes, _output_dir, _errors_path, _accelerator, _energy, _tunex, _tuney = _parse_args()
    get_systematic_errors(_model_twiss, _num_simulations, _num_processes, _output_dir, _errors_path, _accelerator, _energy, _tunex, _tuney)
