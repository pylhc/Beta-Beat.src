import cPickle as pickle
import json
import multiprocessing
import os
import random
from collections import OrderedDict
from Correction.iterative.response_pandas import _generate_fullresponse as generate_madx_resp

import numpy as np
import pandas as pd

from Utilities import logging_tools as logtool, iotools
from Utilities import tfs_pandas as tfs
from Utilities.dict_tools import DotDict
from madx import madx_wrapper as madx
from test_helpers import plot_df_comparison
from test_helpers import plot_single_magnet
from test_helpers import rms
from test_helpers import error_meanabs as error_fun
from twiss_optics.response_class import TwissResponse
from twiss_optics.response_class import get_delta

LOG = logtool.get_logger(__name__)

EXPONENT_RANGE = np.arange(-6, 0, .5)
N_RUNS_SINGLE = 25
N_RUNS_MULTI = 25
N_STATS = 5
N_MAGNETS = 10

DATA_DIR = os.path.abspath('test_data')
DATA_RESP = os.path.join(DATA_DIR, "response")
DATA_RESULT = os.path.join(DATA_DIR, "results")

"""
=============================   Tests   =============================
"""


def test_k1_change():
    LOG.info("Running K1 changing test.")

    # File Definitions
    madxfile_path = os.path.join(DATA_DIR, 'test_k1_change.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_k1_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    # Variable Changes
    dk1_start = pd.Series({
        "kq4.r1b1": 5e-5,
        "kq5.r1b1": 1e-4,
        "kq6.r1b1": -2e-5,
        "kq7.r1b1": 2e-5,
        "kq8.r1b1": 5e-5,
        "kq9.r1b1": -1e-4,
    })


    dk1 = dk1_start * 1e-1


    # Calculation
    twiss_delta, madx_delta = _get_deltas(madxfile_path,
                                          seqfile_path, seq_name,
                                          variables_path, dk1)

    # Results
    plot_df_comparison(twiss_delta, madx_delta,
                       title="Beta Response",
                       planes=["BETX", "BETY"],
                       ylabel="beta",
                       accel=madx_delta.SEQUENCE,
                       data_labels=["TwissResponse", "MAD-X"])

    plot_df_comparison(twiss_delta, madx_delta,
                       title="Phase Response",
                       planes=["MUX", "MUY"],
                       ylabel="phasetot",
                       accel=madx_delta.SEQUENCE,
                       data_labels=["TwissResponse", "MAD-X"])

    plot_df_comparison(twiss_delta, madx_delta,
                       title="Dispersion Response",
                       planes=["DX", "DY"],
                       ylabel="dispersion",
                       accel=madx_delta.SEQUENCE,
                       data_labels=["TwissResponse", "MAD-X"])

    for tune in ["QX", "QY"]:
        LOG.info("  {tune:s} absolute Error: {err:g}".format(
            tune=tune,
            err=(twiss_delta.headers[tune] - madx_delta.headers[tune])))

    LOG.info("  Mean absolute K1L change: {:g}".format(np.mean(np.abs(dk1))))


def test_dispersion():
    LOG.info("Running Dispersion Response Test.")
    LOG.info("Running K1 changing test.")
    random.seed(2011)

    # File Definitions
    madxfile_path = os.path.join(DATA_DIR, 'test_dispersion_response.madx')  #tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    # Variable Changes
    # dk_start = pd.Series({
    #     # K0
    #     "ad2.l2": 0,
    #     "ad2.l5": 5e-5,
    #     "ad2.l8": 5e-5,
    #     # K1
    #     "kq4.r1b1": 5e-5,
    #     "kq5.r1b1": 1e-4,
    #     "kq6.r1b1": -2e-5,
    #     "kq7.r1b1": 2e-5,
    #     "kq8.r1b1": 5e-5,
    #     "kq9.r1b1": -1e-4,
    #     # K1S
    #     "kqs.a23b1": 5e-5,
    #     "kqs.r7b1": 5e-5,
    #     "kqsx3.l1": 5e-5,
    #     "kqsx3.l2": 5e-5,
    # })


    # dk_start = pd.Series({
    #     # K0
    #     "ad2.l2": 0,
    #     "ad2.l5": 0,
    #     "ad2.l8": 0,
    #     # K1
    #     "kq4.r1b1": 0,
    #     "kq5.r1b1": 0,
    #     "kq6.r1b1": 0,
    #     "kq7.r1b1": 0,
    #     "kq8.r1b1": 0,
    #     "kq9.r1b1": 0,
    #     # K1S
    #     "kqs.a23b1": 1e-4,
    #     "kqs.r7b1": -1e-5,
    #     "kqsx3.l1": -1e-4,
    #     "kqsx3.l2": 1e-4,
    # })


    # dk_start = pd.Series({
    #     # K0
    #     "ad2.l2": 0,
    #     "ad2.l5": 0,
    #     "ad2.l8": 0,
    #     # K1
    #     "kq4.r1b1": 5e-5,
    #     "kq5.r1b1": 1e-4,
    #     "kq6.r1b1": -2e-5,
    #     "kq7.r1b1": 2e-5,
    #     "kq8.r1b1": 5e-5,
    #     "kq9.r1b1": -1e-4,
    #     # K1S
    #     "kqs.a23b1": 0,
    #     "kqs.r7b1": 0,
    #     "kqsx3.l1": 0,
    #     "kqsx3.l2": 0,
    # })


    # dk_start = pd.Series({
    #     "kq4.r1b1": 5e-5,
    #     "kq6.r1b1": -2e-5,
    #     "kq5.r1b1": 1e-5,
    # })


    # dk_start = pd.Series({
    #     # K0
    #     "ad2.l2": 5e-4,
    #     "ad2.l5": -1e-5,
    #     "ad2.l8": 5e-5,
    #     # K1
    #     "kq4.r1b1": 0,
    #     "kq5.r1b1": 0,
    #     "kq6.r1b1": 0,
    #     "kq7.r1b1": 0,
    #     "kq8.r1b1": 0,
    #     "kq9.r1b1": 0,
    #     # K1S
    #     "kqs.a23b1": 0,
    #     "kqs.r7b1": 0,
    #     "kqsx3.l1": 0,
    #     "kqsx3.l2": 0,
    # })


    dk_start = pd.Series({
        # K0
        "ad2.l2": 5e-4,
        "ad2.l5": -1e-5,
        "ad2.l8": 5e-5,
        # K1
        "kq4.r1b1": 0,
        "kq5.r1b1": 0,
        "kq6.r1b1": 0,
        "kq7.r1b1": 0,
        "kq8.r1b1": 0,
        "kq9.r1b1": 0,
        # K1S
        "kqs.a23b1": 1e-4,
        "kqs.r7b1": -1e-5,
        "kqsx3.l1": -1e-4,
        "kqsx3.l2": 1e-4,
    })

    dk = dk_start * 1e-1

    # Calculation
    twiss_delta, madx_delta = _get_deltas(madxfile_path,
                                          seqfile_path, seq_name,
                                          variables_path, dk, [])

    plot_df_comparison(twiss_delta, madx_delta,
                       title="Dispersion Response",
                       planes=["DX", "DY"],
                       ylabel="dispersion",
                       accel=madx_delta.SEQUENCE,
                       data_labels=["TwissResponse", "MAD-X"])

    LOG.info("  Mean absolute K1L change: {:g}".format(np.mean(np.abs(dk))))


def test_K0_single():
    LOG.info("Running K0 Single test.")
    random.seed(2012)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_k0_single.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_k0_single.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["K0L"]
    exponents = EXPONENT_RANGE
    n_magnets = 1
    n_magnet_runs = N_RUNS_SINGLE
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "$\delta$ K0L Single")


def test_K0_multi():
    LOG.info("Running K0 Multi test.")
    random.seed(2013)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_k0_multi.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_k0_multi.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["K0L"]
    exponents = EXPONENT_RANGE
    n_magnets = N_MAGNETS
    n_magnet_runs = N_RUNS_MULTI
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "$\delta$ K0L Multi")


def test_K1_single():
    LOG.info("Running K1 Single test.")
    random.seed(2014)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_k1_single.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_k1_single.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["K1L"]
    exponents = EXPONENT_RANGE
    n_magnets = 1
    n_magnet_runs = N_RUNS_SINGLE
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "$\delta$ K1L Single")


def test_K1_multi():
    LOG.info("Running K1 Multi test.")
    random.seed(2015)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_k1_multi.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_k1_multi.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["K1L"]
    exponents = EXPONENT_RANGE
    n_magnets = N_MAGNETS
    n_magnet_runs = N_RUNS_MULTI
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "$\delta$ K1L Multi")


def test_K1S_single():
    LOG.info("Running K1S Single test.")
    random.seed(2016)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_k1s_single.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_k1s_single.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["K1SL"]
    exponents = EXPONENT_RANGE
    n_magnets = 1
    n_magnet_runs = N_RUNS_SINGLE
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "$\delta$ K1SL Single")


def test_K1S_multi():
    LOG.info("Running K1S Multi test.")
    random.seed(2017)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_k1s_multi.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_k1s_multi.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["K1SL"]
    exponents = EXPONENT_RANGE
    n_magnets = N_MAGNETS
    n_magnet_runs = N_RUNS_MULTI
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "$\delta$ K1SL Multi")


def test_ALL_multi():
    LOG.info("Running All Multi test.")
    random.seed(2018)

    # File Definitions
    resultsfile_path = os.path.join(DATA_RESULT, 'results_test_all_multi.dat')
    madxfile_path = os.path.join(DATA_RESP, 'test_all_multi.madx')  # tbc
    variables_path = os.path.join(DATA_DIR, 'lhcb1_all_list.json')  # exists
    seqfile_path = os.path.join(DATA_DIR, 'lhcb1.seq')  # exists
    seq_name = "lhcb1"

    k_list = _load_varnames(variables_path)["all"]
    exponents = EXPONENT_RANGE
    n_magnets = 2 * N_MAGNETS
    n_magnet_runs = N_RUNS_MULTI
    n_stats = N_STATS

    results = _get_stats(exponents, k_list, n_magnets, n_magnet_runs, n_stats,
                         resultsfile_path, madxfile_path, seqfile_path, seq_name, variables_path)

    plot_single_magnet(results, "All $\delta$ K Multi")


"""
=============================   Subfunctions General  =============================
"""


def _get_stats(exponents, k_list, n_vars, n_runs_per_var, n_stats, resultsfile_path, madxfile_path,
               seqfile_path, seq_name, variables_path):
    """ Gathers results over variables, strengths and statistic-runs and saves them to a file """
    if os.path.exists(resultsfile_path):
        with open(resultsfile_path, "rb") as f:
            results = pickle.load(f)
    else:
        results = {}
        variables = [None] * n_runs_per_var
        # optimized to keep only the errors in memory and free space for next set of variables
        for idx_runs in range(n_runs_per_var):
            # use the same variables for every strength
            variables[idx_runs] = take_at_random(k_list, n_vars)

        for idx_exp, exp in enumerate(exponents):
            results_temp = {}
            for idx_runs in range(n_runs_per_var):
                results_temp[idx_runs] = {}
                for idx_stat in range(n_stats):
                    dk_list = pd.Series()
                    for var in variables[idx_runs]:
                        dk_list[var] = _get_random_strength(exp)

                    LOG.info("Strength-Scale {:d} of {:d}".format(idx_exp + 1, len(exponents)))
                    LOG.info("  {:d} variable(s):  {:d} of {:d}".format(
                        n_vars, idx_runs + 1, n_runs_per_var))
                    LOG.info(
                        "    Stats {:d} of {:d}".format(idx_stat+1, n_stats))
                    LOG.info("    Current dk \n {:}".format(dk_list))

                    # Calculation
                    twiss_delta, madx_delta, madx_resp_delta = _get_deltas(madxfile_path,
                                                                           seqfile_path, seq_name,
                                                                           variables_path, dk_list)

                    results_temp[idx_runs][idx_stat] = {
                        "twiss": twiss_delta,
                        "madx": madx_delta,
                        "madx_resp": madx_resp_delta,
                    }
            error_av = calc_error_averages(results_temp)

            # regroup
            for mad_x in ["madx_mean", "madx_rms", "madx_resp_mean", "madx_resp_rms"]:
                if mad_x not in results:
                    results[mad_x] = tfs.TfsDataFrame(None, columns=error_av.columns)
                results[mad_x].loc[exp, :] = error_av.loc[mad_x, :]

        iotools.create_dirs(os.path.dirname(resultsfile_path))
        with open(resultsfile_path, "wb") as f:
            pickle.dump(results, f, -1)

    return results


def _get_deltas(madxfile_path, seqfile_path, seq_name, variables_path, delta_k, exclude=None):
    """ Get the differences calculated by TwissResponse and by MAD-X """
    model_before, model_after = _calc_madx_models(madxfile_path, seqfile_path, seq_name, delta_k)
    madx_delta = _get_madx_delta(model_before, model_after)
    twiss_delta = _get_twiss_delta(model_before, seqfile_path, variables_path, delta_k, exclude)
    madx_resp_delta = _get_madx_response_delta(model_before, seqfile_path, seq_name,
                                               variables_path, delta_k, exclude)

    elmts = madx_delta.index.intersection(twiss_delta.index)  # as TR filters (M|BPM), MADX [MB]
    madx_delta = madx_delta.loc[elmts, :]
    twiss_delta = twiss_delta.loc[elmts, :]
    madx_resp_delta = madx_resp_delta.loc[elmts, :]
    return twiss_delta, madx_delta, madx_resp_delta


def _load_varnames(variables_path):
    """ Load variable Names from json file """
    with open(variables_path, "r") as varfile:
        varnames = json.load(varfile)
    return varnames


"""
=============================   Subfunctions Twiss  =============================
"""


def _get_twiss_delta(model_path, seqfile_path, variables_path, delta_k, exclude):
    """ Returns the differences as calculated by means of TwissResponse """
    # Load TwissResponse or calculate it
    folder = os.path.dirname(model_path)
    filename = "TwissResp_{seq:s}_{var:s}.dat".format(
        seq=_shorten(seqfile_path),
        var=_shorten(variables_path),
    )

    tr_path = os.path.join(folder, filename)
    if os.path.exists(tr_path):
        with open(tr_path, "rb") as tr_file:
            tr = pickle.load(tr_file)
        LOG.info("Loaded TwissResponse '{:s}'".format(tr_path))

    else:
        tr = TwissResponse(seqfile_path, model_path, variables_path,
                           exclude_categories=exclude, at_elements="all").get_fullresponse()

        iotools.create_dirs(os.path.dirname(tr_path))
        with open(tr_path, "wb") as tr_file:
            pickle.dump(tr, tr_file, -1)
        LOG.info("Saved TwissResponse '{:s}'".format(tr_path))

    # use response_class.get_delta() for the dot-product
    delta = get_delta(tr, delta_k)
    return delta


"""
=============================   Subfunctions MAD-X Calculation  =============================
"""


def _calc_madx_models(madxfile_path, seqfile_path, seq_name, delta_k):
    """ Calls madx to create the two models before and after changing the varaiable values """
    LOG.debug("Creating MADX models for '{:s}'".format(madxfile_path))
    mod_before_path = madxfile_path.replace(".madx", "_before.dat")
    mod_after_path = madxfile_path.replace(".madx", "_after.dat")
    logfile_path = madxfile_path.replace(".madx", ".log")

    madx_string = """
    call, file='{seq_file:s}';
    beam, sequence={seq_name:s}, particle=proton, energy=450, kbunch=1, npart=1.15E11, bv={bv:d};
    use, sequence={seq_name:s};
    select, flag=twiss, clear;
    select, flag=twiss, pattern="^[MB]", column=NAME, S, X, Y, BETX, BETY, DX, DY, MUX, MUY;
    twiss, sequence={seq_name:s}, file='{model_before:s}'; 
    """

    for k in delta_k.index:
        madx_string += "{var:s} = {var:s} + ({val:f});\n".format(var=k, val=delta_k[k])

    madx_string += "twiss, sequence={seq_name:s}, file='{model_after:s}';"

    iotools.create_dirs(os.path.dirname(madxfile_path))
    with open(madxfile_path, "w") as madxfile:
        madxfile.write(
            madx_string.format(
                seq_file=seqfile_path,
                seq_name=seq_name,
                model_before=mod_before_path,
                model_after=mod_after_path,
                bv=1 if seq_name[-1] == "1" else -1,
            )
        )
    madx.resolve_and_run_file(madxfile_path, log_file=logfile_path)
    return mod_before_path, mod_after_path


def _get_madx_delta(twiss_before, twiss_after):
    """ Returns the differences as calculated by means of MAD-X """
    df_before = tfs.read_tfs(twiss_before, index="NAME")
    df_after = tfs.read_tfs(twiss_after, index="NAME")

    columns = ["BETX", "BETY", "MUX", "MUY", "DX", "DY"]
    delta = df_after.loc[:, columns] - df_before.loc[:, columns]

    delta.headers["QX"] = df_after.Q1 - df_before.Q1
    delta.headers["QY"] = df_after.Q2 - df_before.Q2

    delta["S"] = df_before["S"]
    delta.headers["SEQUENCE"] = df_before.SEQUENCE
    return delta


"""
=============================   Subfunctions MAD-X Response  =============================
"""


def _get_madx_response_delta(model_path, seqfile_path, seq_name, variables_path, delta_k, exclude):
    """ Returns the differences as calculated by means of TwissResponse """
    # Load MadxResponse or calculate it
    filename = "MadxResp_{seq:s}_{var:s}.dat".format(
        seq=_shorten(seqfile_path),
        var=_shorten(variables_path),
    )

    mr_path = os.path.join(DATA_RESP, filename)
    if not os.path.exists(mr_path):
        temp_dir = os.path.join(DATA_RESP, "madx_temp")
        iotools.create_dirs(temp_dir)

        jobfile_path = os.path.join(temp_dir, "job.iterate.madx")

        madx_string = """
            call, file='{seq_file:s}';
            beam, sequence={seq_name:s}, particle=proton, energy=450, kbunch=1, npart=1.15E11, bv={bv:d};
            use, sequence={seq_name:s};
            select, flag=twiss, clear;
            select, flag=twiss, pattern="^[MB]", column=NAME, S, X, Y, BETX, BETY, DX, DY, MUX, MUY;
            
            %CALL_ITER_FILE%
        """
        iotools.create_dirs(os.path.dirname(jobfile_path))
        with open(jobfile_path, "w") as madxfile:
            madxfile.write(
                madx_string.format(
                    seq_file=seqfile_path,
                    seq_name=seq_name,
                    bv=1 if seq_name[-1] == "1" else -1,
                )
            )

        options = DotDict({
            "outfile_path": mr_path,
            "temp_dir": temp_dir,
            "varfile_path": variables_path,
            "original_jobfile_path": jobfile_path,
            "pattern": "%CALL_ITER_FILE%",
            "num_proc": multiprocessing.cpu_count(),
            "delta_k": 2e-6,
            "exclude_categories": [] if exclude is None else exclude,
        })
        generate_madx_resp(options)

    with open(mr_path, "rb") as tr_file:
        madx_resp = pickle.load(tr_file)
    LOG.info("Loaded MadxResponse '{:s}'".format(mr_path))

    madx_resp = madx_to_twiss_resp(madx_resp, model_path)

    # use response_class.get_delta() for the dot-product
    delta = get_delta(madx_resp, delta_k)
    return delta


def madx_to_twiss_resp(madx_resp, model_path):
    """ Converts the pandas output to twiss compatible output """
    twiss = tfs.read_tfs(model_path, index="NAME")
    return {
        "BETX": madx_resp["BBX"].mul(twiss["BETX"], axis="columns").transpose(),
        "BETY": madx_resp["BBY"].mul(twiss["BETY"], axis="columns").transpose(),
        "MUX": madx_resp["MUX"].transpose(),
        "MUY": madx_resp["MUY"].transpose(),
        "DX": madx_resp["NDX"].mul(np.sqrt(twiss["BETX"]), axis="columns").transpose(),
        "DY": madx_resp["NDY"].mul(np.sqrt(twiss["BETY"]), axis="columns").transpose(),
        # "DY": tfs.TfsDataFrame(None, index=twiss["BETX"].index),
        "Q":  madx_resp["Q"].transpose(),
    }


"""
=============================   Helpers   =============================
"""


def _get_random_strength(exp):
    """ Returns a random strength "gaussian" around 10**exp)

     "gaussian" means that mu=10**exp and sigma=.2*10**exp,
      but the tails are cut of at 3 sigma.
    """
    mu = 10 ** exp
    sigma = .2 * mu
    s_min = mu - 3 * sigma
    s_max = mu + 3 * sigma
    while True:
        strength = (np.sign(random.random() - .5) *
                    random.gauss(mu, sigma))
        if s_min <= np.abs(strength) <= s_max:
            return strength


def _shorten(string):
    """ Get only filename (w/o extension) with removed all vowels """
    string = os.path.basename(string)
    string = os.path.splitext(string)[0]
    string = string.replace("a", "").replace("e", "").replace(
        "i", "").replace("o", "").replace("u", "")
    return string


def take_at_random(items, n):
    """ Take n elements at random out of items

    (implementation tested in time against random.randrange() + pop())
    """
    items = list(items)
    random.shuffle(items)  # shuffles in place
    return items[:n]


def calc_error_averages(results):
    """ Create DataFrame for MADX and MADX_response vs TwissResponse errors """
    parameters = ["BETX", "BETY", "MUX", "MUY", "DX", "DY", "QX", "QY"]
    df = tfs.TfsDataFrame(None, index=["madx_mean", "madx_rms",
                                       "madx_resp_mean", "madx_resp_rms"])

    for mad_x in ["madx", "madx_resp"]:
        for param in parameters:
            rms_values = []
            mean_values = []
            for run in results:
                for stat in results[run]:
                    error_values = error_fun(
                        results[run][stat]["twiss"][param],
                        results[run][stat][mad_x][param]
                    )
                    rms_values.append(rms(error_values))
                    mean_values.append(np.mean(error_values))

            df.loc[mad_x + "_mean", param] = np.mean(mean_values)
            df.loc[mad_x + "_mean", param + '_min'] = np.min(mean_values)
            df.loc[mad_x + "_mean", param + '_max'] = np.max(mean_values)

            df.loc[mad_x + "_rms", param] = np.mean(rms_values)
            df.loc[mad_x + "_rms", param + '_min'] = np.min(rms_values)
            df.loc[mad_x + "_rms", param + '_max'] = np.max(rms_values)

    return df


"""
=============================   Main   =============================
"""


if __name__ == "__main__":
    test_K0_single()
    test_K0_multi()
    test_K1_single()
    test_K1_multi()
    test_K1S_single()
    test_K1S_multi()
    test_ALL_multi()
    # test_dispersion()
    # test_k1_change()
    # plt.show()
