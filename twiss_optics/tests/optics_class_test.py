""" Test cases for optics_class.py

Look in test_helpers.py to change plotting options.
"""

import os
import numpy as np
from twiss_optics.optics_class import TwissOptics
from utils.contexts import timeit
from utils import logging_tools
from utils import tfs_pandas as tfs
from test_helpers import plot_df_comparison

LOG = logging_tools.get_logger(__name__)

"""
====================== Test Cases ======================
"""


def test_rdt_two_octu(rdterms=None):
    LOG.info("Comparison TwissOptics with PTC Octupoles:")
    twissfile = os.path.join('test_data', 'twiss_2oct_elements.dat')
    ptc_file = os.path.join('test_data', 'ptc_normal_2octupoles.dat')
    _test_rdt_ptc(twissfile, ptc_file, rdterms)
    LOG.info('Octupole Test Finished\n\n')


def test_rdt_two_sextu(rdterms=None):
    LOG.info("Comparison TwissOptics with PTC Sextupoles:")
    twissfile = os.path.join('test_data', 'twiss_2sext_elements.dat')
    ptc_file = os.path.join('test_data', 'ptc_normal_2sextupoles.dat')
    _test_rdt_ptc(twissfile, ptc_file, rdterms)
    LOG.info('Sextupole Test Finished\n\n')


def test_linear_dispersion():
    LOG.info("Comparison TwissOptics with MADX Linear Dispersion:")
    twissfile = os.path.join('test_data', 'twiss_dispersion.dat')
    _test_linear_dispersion(twissfile)
    LOG.info("Linear Dispersion Test finished:")


def test_eight_coupling():
    LOG.info("Comparison TwissOptics with MADX Coupling:")
    twissfile = os.path.join('test_data', 'twiss_coupling.dat')
    _test_coupling(twissfile)
    LOG.info('Coupling Test Finished\n\n')


def test_eight_coupling_cmatrix():
    LOG.info("Comparison TwissOptrics RDT-Formula vs. CMatrix Coupling:")
    twissfile = os.path.join('test_data', 'twiss_coupling.dat')
    _test_coupling_cmatrix(twissfile)
    LOG.info('Coupling Test Finished\n\n')


def test_chromatic_beating_general():
    LOG.info("Comparison TwissOptics with MADX Chromatic Beating:")
    twissfile = os.path.join('test_data', 'twiss_chromatic_beating.dat')
    _test_chromatic_beating(twissfile)
    LOG.info('Coupling Test Finished\n\n')


"""
====================== Tests ======================
"""


def _test_rdt_ptc(twissfile, ptcfile, rdterms):
    LOG.info("RDT Test started.")
    LOG.info("  Twiss file: {:s}".format(twissfile))
    LOG.info("  PTC file: {:s}".format(ptcfile))

    ptc = tfs.read_tfs(ptcfile).set_index("NAME")

    if not rdterms:
        rdterms = [col.replace('_AMP', '').lower() for col in ptc if col.startswith('f')]

    with timeit(lambda t: LOG.info("  Time for import in TwissOptics: {:f}s".format(t))):
        twop = TwissOptics(twissfile)

    with timeit(lambda t: LOG.info("  Time to calculate RDTs in TwissOptics: {:f}s".format(t))):
        twop.calc_rdts(rdterms)

    for rdt in rdterms:
        tw_rdt = twop.get_rdts(rdt).dropna()
        intersect_bpms = ptc.index.intersection(tw_rdt.index)

        ptc_rdt = ptc.loc[intersect_bpms, [rdt.lower() + '_AMP']]
        ptc_rdt.columns = [rdt]
        tw_rdt = np.abs(tw_rdt.loc[intersect_bpms, ["S", rdt.upper()]])
        tw_rdt.columns = ["S", rdt]

        plot_df_comparison(tw_rdt, ptc_rdt,
                    "Resonance Driving Term", [rdt], 'absolute',
                           accel=twop.mad_twiss.SEQUENCE,
                           plane_labels=[rdt])


def _test_linear_dispersion(twissfile, ptcfile=None):
    LOG.info("Linear Dispersion Test started.")
    LOG.info("  Twiss file: {:s}".format(twissfile))

    with timeit(lambda t: LOG.info("  Time for import in TwissOptics: {:f}s".format(t))):
        twop = TwissOptics(twissfile)

    with timeit(lambda t:
                LOG.info("  Time to calculate linear dispersion in TwissOptics: {:f}s".format(t))):
        twop.calc_linear_dispersion(bpms_only=False)
        tw_ldisp = twop.get_linear_dispersion().dropna()

    if ptcfile is None:
        mad_twiss = twop.mad_twiss
        bpms = tw_ldisp.index
    else:
        mad_twiss = tfs.read_tfs(ptcfile, index='NAME')
        bpms = mad_twiss.index.intersection(tw_ldisp.index)
        tw_ldisp = tw_ldisp.loc[bpms, :]

    mad_ldisp = mad_twiss.loc[bpms, ['DX', 'DY']]

    plot_df_comparison(tw_ldisp, mad_ldisp,
                "Linear Dispersion", ['DX', 'DY'], 'dispersion',
                       twop.mad_twiss.SEQUENCE)


def _test_coupling(twissfile, ptcfile=None):
    LOG.info("Coupling Test started.")
    LOG.info("  Twiss file: {:s}".format(twissfile))
    couple_rdts = ['F1001', 'F1010']

    with timeit(lambda t: LOG.info("  Time for import in TwissOptics: {:f}s".format(t))):
        twop = TwissOptics(twissfile)

    with timeit(lambda t:
                LOG.info("  Time to calculate coupling in TwissOptics: {:f}s".format(t))):
        twop.calc_rdts(couple_rdts)
        tw_coupl = twop.get_rdts(couple_rdts).dropna()

    if ptcfile is None:
        mad_twiss = twop.mad_twiss
        bpms = tw_coupl.index
    else:
        mad_twiss = tfs.read_tfs(ptcfile, index='NAME')
        bpms = mad_twiss.index.intersection(tw_coupl.index)
        tw_coupl = tw_coupl.loc[bpms, :]

    tfs.add_coupling(mad_twiss)
    mad_coupl = mad_twiss.loc[bpms, [r.lower() for r in couple_rdts]]
    mad_coupl.columns = couple_rdts

    plot_df_comparison(np.abs(tw_coupl), np.abs(mad_coupl),
                "Coupling Term", couple_rdts, 'absolute',
                       accel=twop.mad_twiss.SEQUENCE,
                       plane_labels=couple_rdts)


def _test_coupling_cmatrix(twissfile):
    LOG.info("CMatrix Coupling Test started.")
    LOG.info("  Twiss file: {:s}".format(twissfile))
    couple_rdts = ['F1001', 'F1010']

    with timeit(lambda t: LOG.info("  Time for import in TwissOptics: {:f}s".format(t))):
        twop = TwissOptics(twissfile)

    with timeit(lambda t:
                LOG.info("  Time to calculate coupling via RDT-formulas: {:f}s".format(t))):
        twop.calc_rdts(couple_rdts)

    with timeit(lambda t:
                LOG.info("  Time to calculate coupling via C-Matrix: {:f}s".format(t))):
        twop.calc_cmatrix()

    rdt_coupl = twop.get_coupling(method='rdt').dropna()
    cm_couple = twop.get_coupling(method='cmatrix').loc[rdt_coupl.index, :]


    plot_df_comparison(np.abs(rdt_coupl), np.abs(cm_couple),
                "Coupling Term", couple_rdts, 'absolute',
                       accel=twop.mad_twiss.SEQUENCE,
                       plane_labels=couple_rdts,
                       data_labels=['rdt', 'cmatrix'])


def _test_chromatic_beating(twissfile, ptcfile=None):
    """
    Hint: "Notice also that in MAD-X, PT substitutes DELTAP as longitudinal variable. Dispersive and
    chromatic functions are hence derivatives with respect to PT. Since PT=BETA*DELTAP,
    where BETA is the relativistic Lorentz factor, those functions given by MAD-X must be multiplied
    by BETA a number of times equal to the order of the derivative to find the functions
    given in the literature" - MADX User Guide
    """

    LOG.info("Chromatic Beating Test started.")
    LOG.info("  Twiss file: {:s}".format(twissfile))

    with timeit(lambda t: LOG.info("  Time for import in TwissOptics: {:f}s".format(t))):
        twop = TwissOptics(twissfile)

    with timeit(lambda t:
                LOG.info("  Time to calculate chromatic beating in TwissOptics: {:f}s".format(t))):
        twop.calc_chromatic_beating()
        tw_beat = twop.get_chromatic_beating().dropna()

    if ptcfile is None:
        mad_twiss = twop.mad_twiss
        bpms = tw_beat.index
    else:
        mad_twiss = tfs.read_tfs(ptcfile, index='NAME')
        bpms = mad_twiss.index.intersection(tw_beat.index)
        tw_beat = tw_beat.loc[bpms, :]

    mad_beat = mad_twiss.loc[bpms, ['WX', 'WY']].mul(
        mad_twiss.loc[bpms, ['PHIX', 'PHIY']].apply(lambda x: np.cos(2*np.pi*x)).values,
        axis='index') * (1 - 1/(mad_twiss.GAMMA**2))**.5
    mad_beat.columns = ['DBEATX', 'DBEATY']

    plot_df_comparison(tw_beat, mad_beat,
                "Chromatic Beating", ['DBEATX', 'DBEATY'], 'dbetabeat',
                       twop.mad_twiss.SEQUENCE)


"""
====================== Main ======================
"""


if __name__ == '__main__':
    test_rdt_two_octu()
    test_rdt_two_sextu()
    test_linear_dispersion()
    test_eight_coupling()
    test_eight_coupling_cmatrix()
    test_chromatic_beating_general()


