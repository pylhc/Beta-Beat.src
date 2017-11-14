"""
    Provides Classes to calculate optics from twiss parameters.
    The calculation is based on formulas in [1,2].

    Resonance Driving Terms: Eq.A8 in [2]
    Linear Dispersion: Eq. 24 in [2]

    [1]
         First simultaneous measurement of sextupolar and octupolar resonance driving terms in a
         circular accelerator from turn-by-turn beam position monitor data,
         Phys. Rev. ST Accel. Beams 17, 074001 (2014).
         https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.17.074001

    [2]  A. Franchi et al.,
         Analytic formulas for the rapid evaluation of the orbit response matrix and chromatic
         functions from lattice parameters in circular accelerators
         NOT YET PUBLISHED

"""


from __future__ import print_function
import __init__
import os, sys
import numpy as np
import pandas as pd
import itertools
from math import factorial
from Utilities import tfs_pandas as tfs
from Utilities.contexts import timeit
import matplotlib.pyplot as plt
from Utilities.plotting import plot_style as pstyle
from Utilities import logging_tools as logtool

import re


LOG = logtool.get_logger(__name__, level_console=0)


################################
#           CLASSES
################################

class TwissOptics(object):

    ################################
    #       init Functions
    ################################

    def __init__(self, model_path):
        LOG.debug("Creating TwissOptics from '{:s}'".format(model_path))
        self.mad_twiss = tfs.read_tfs(model_path).set_index('NAME')
        self._add_ip_pos()  # self._ip_pos
        self._remove_nonnecessaries()
        self._add_result_dataframe()  # self._results_df
        self._add_element_types()  # self._elements, self._elements_mapped
        self._add_phase_advances()  # self._phase_advances

    def _add_ip_pos(self):
        tw = self.mad_twiss
        self._ip_pos = tw.loc[regex_in("\AIP\d$", tw.index), 'S']

    def _remove_nonnecessaries(self):
        LOG.debug("Removing non-necessaries:")
        LOG.debug("  Entries total: {:d}".format(self.mad_twiss.shape[0]))
        with timeit(lambda t:
                    LOG.debug("  Removed in {:f}s".format(t))):
            self.mad_twiss = self.mad_twiss.loc[regex_in("\A[MB]", self.mad_twiss.index), :]
        LOG.debug("  Entries left: {:d}".format(self.mad_twiss.shape[0]))

    def _add_result_dataframe(self):
        LOG.debug("Adding results dataframe.")
        self._results_df = tfs.TfsDataFrame(
            data=self.mad_twiss['S'],
            index=self.mad_twiss.index,
        )

    def _add_element_types(self):
        LOG.debug("Sorting Elements into types:")

        tw = self.mad_twiss
        with timeit(lambda t:
                    LOG.debug("  Sorted in {:f}s".format(t))):
            self._elements = {
                'BPM': regex_in('\ABPM', tw.index),
                'MB': regex_in('\AMB[^S]', tw.index),
                'MBS': regex_in('\AMBS', tw.index),
                'MQ': regex_in('\AMQ[^ST]', tw.index),
                'MQS': regex_in('\AMQS', tw.index),
                'MQT': regex_in('\AMQT', tw.index),
                'MS': regex_in('\AMS[^S]', tw.index),
                'MSS': regex_in('\AMSS', tw.index),
                'MO': regex_in('\AMO[^S]', tw.index),
                'MOS': regex_in('\AMOS', tw.index),
                'MCB': regex_in('\AMCB', tw.index),
                'MCS': regex_in('\AMCS', tw.index),
                'MCO': regex_in('\AMCO', tw.index),
                'MCD': regex_in('\AMCD', tw.index),
            }

            el = self._elements
            self._elements_mapped = {
                'K0L': el['BPM'] | el['MB'],
                'K0SL': el['BPM'] | el['MBS'],
                'K1L': el['BPM'] | el['MQ'],
                'K1SL': el['BPM'] | el['MQS'],
                'K2L': el['BPM'] | el['MS'],
                'K2SL': el['BPM'] | el['MSS'],
                'K3L': el['BPM'] | el['MO'],
                'K3SL': el['BPM'] | el['MOS'],
            }

    def _add_phase_advances(self):
        """
        Calculate phase advances between all elements
        :return: Matrices similar to DPhi(i,j) = Phi(j) - Phi(i)
        """
        LOG.debug("Adding Phase Advances:")
        self._phase_advances = dict.fromkeys(['X', 'Y'])
        with timeit(lambda t:
                    LOG.debug("  Phase Advances calculated in {:f}s".format(t))):
            for plane in ['X', 'Y']:
                colmn_phase = "MU" + plane

                phases_mdl = self.mad_twiss.loc[self.mad_twiss.index, colmn_phase]
                phase_advances = pd.DataFrame((phases_mdl[:, None] - phases_mdl[None, :]))
                phase_advances.set_index(self.mad_twiss.index, inplace=True)
                phase_advances.columns = self.mad_twiss.index
                # Do not calculate dphi and tau here.
                # only slices of phase_advances as otherwise super slow
                self._phase_advances[plane] = phase_advances

    ################################
    #         Properties
    ################################

    @property
    def betas(self):
        """ Betas Function """
        return self._results_df.betas

    @betas.setter
    def betas(self, new):
        self._results_df.betas = new

    @property
    def dispersion(self):
        """ Dispersion Function """
        return self._results_df.dispersion

    @dispersion.setter
    def dispersion(self, new):
        self._results_df.dispersion = new

    @property
    def coupling(self):
        """ Coupling Function """
        return self._results_df.coupling

    @coupling.setter
    def coupling(self, new):
        self._results_df.coupling = new

    ################################
    #   Resonance Driving Terms
    ################################

    def calc_rdts(self, order):
        """ Calculate RDTs """
        if isinstance(order, int):
            order = get_all_rdts(order)
        elif not isinstance(order, list):
            order = [order]

        LOG.debug("Calculating RDTs: {:s}.".format(str(order)[1:-1]))
        with timeit(lambda t:
                    LOG.debug("  RDTs calculated in {:f}s".format(t))):

            i2pi = 2j * np.pi
            tw = self.mad_twiss
            phs_adv = self._phase_advances
            res_df = self._results_df

            for rdt in order:
                assertion(len(rdt) == 5 and rdt[0].upper() == 'F',
                          ValueError("'{:s}' does not seem to be a valid RDT name.".format(rdt)))

                conj_rdt = ''.join(['F', rdt[2], rdt[1], rdt[4], rdt[3]])

                if conj_rdt in self._results_df:
                    res_df[rdt.upper()] = np.conjugate(self._results_df[conj_rdt])
                else:
                    j, k, l, m = int(rdt[1]), int(rdt[2]), int(rdt[3]), int(rdt[4])
                    n = j + k + l + m

                    assertion(n >= 2, ValueError(
                        "The RDT-order has to be >1 but was {:d} for {:s}".format(n, rdt)))

                    denom = 1./(factorial(j) * factorial(k) * factorial(l) * factorial(m) *
                                  2**n * (1. - np.exp(i2pi * ((j-k) * tw.Q1 + (l-m) * tw.Q2))))

                    if (l + m) % 2 == 0:
                        src = 'K' + str(n-1) + 'L'
                        sign = -(1j ** (l+m))
                    else:
                        src = 'K' + str(n-1) + 'SL'
                        sign = -(1j ** (l+m+1))

                    k_mask = tw[src] != 0
                    el_mask = self._elements_mapped[src] | k_mask

                    if sum(k_mask) > 0:
                        # the next three lines determine the main order of speed, hence
                        # - mask as much as possible
                        # - additions are faster than multiplications (-> applymap last)
                        phx = self._dphi(phs_adv['X'].loc[k_mask, el_mask], tw.Q1)
                        phy = self._dphi(phs_adv['Y'].loc[k_mask, el_mask], tw.Q2)
                        phase_term = ((j-k) * phx + (l-m) * phy).applymap(lambda p: np.exp(i2pi*p))

                        beta_term = tw.loc[k_mask, src] * \
                                    tw.loc[k_mask, 'BETX'] ** ((j+k) / 2.) * \
                                    tw.loc[k_mask, 'BETY'] ** ((l+m) / 2.)

                        res_df.loc[el_mask, rdt.upper()] = sign * phase_term.multiply(
                            beta_term, axis="index").sum(axis=0).transpose() * denom

                        LOG.debug("  Average RDT amplitude |{:s}|: {:g}".format(rdt, np.mean(
                            np.abs(res_df.loc[el_mask, rdt.upper()]))))
                    else:
                        LOG.debug("  All {:s} == 0. RDT '{:s}' will be zero.".format(src, rdt))
                        res_df.loc[el_mask, rdt.upper()] = 0

    def get_rdts(self, rdt_names=None):
        """ Return dataframe with rdts """
        if rdt_names:
            if not isinstance(rdt_names, list):
                rdt_names = [rdt_names]
            return self._results_df.loc[:, ["S"] + [rdt.upper() for rdt in rdt_names]]
        else:
            return self._results_df.loc[:, regex_in('\A(S|F\d{4})$', self._results_df.columns)]

    def plot_rdts(self, rdt_names=None, apply_fun=np.abs, combined=True):
        rdts = self.get_rdts(rdt_names)
        is_s = regex_in('\AS$', rdts.columns)
        rdts = rdts.dropna()
        rdts.loc[:, ~is_s] = rdts.loc[:, ~is_s].applymap(apply_fun)
        pstyle.set_style('standard',
                         {u'lines.linestyle': '-',
                          u'lines.marker': ''})
        if combined:
            ax = rdts.plot(x='S')
            ax.set_title('Resonance Driving Terms')
            pstyle.set_yaxis_label(apply_fun.__name__, 'F_{{jklm}}', ax)
            self._nice_axes(ax)
        else:
            for rdt in rdts.loc[:, ~is_s]:
                ax = rdts.plot(x='S', y=rdt)
                ax.set_title('Resonance Driving Term ' + rdt)
                pstyle.set_yaxis_label(apply_fun.__name__, rdt, ax)
                self._nice_axes(ax)

    ################################
    #      Linear Dispersion
    ################################

    def calc_lin_dispersion(self):
        # helper functions
        def coeff_fun(beta, q):
            return np.sqrt(beta) / (2 * np.sin(np.pi * q))

        def sum_fun(k, j, d, beta, tau):
            # k, j, d , beta = columns -> convert to Series -> broadcasted
            # tau = Matrix as Frame
            calc_column = (k + j * d) * np.sqrt(beta)
            if isinstance(calc_column, pd.DataFrame):
                calc_column = calc_column.squeeze()
            return np.cos(2*np.pi * tau).mul(calc_column, axis='index').sum(axis='index')

        tw = self.mad_twiss
        phs_adv = self._phase_advances
        res_df = self._results_df
        tau = self._tau

        # Calculate
        LOG.debug("Calculate Linear Dispersion")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            # sources
            k0_mask = tw['K0L'] != 0
            k0s_mask = tw['K0SL'] != 0
            k1s_mask = tw['K1SL'] != 0

            # results
            mx_mask = k0_mask | k1s_mask  # magnets contributing to Dx,j (-> Dy,m)
            my_mask = k0s_mask | k1s_mask  # magnets contributing to Dy,j (-> Dx,m)

            if not any(mx_mask | my_mask):
                LOG.warning("  No linear dispersion contributions found. Values will be zero.")
                res_df['DX'] = 0
                res_df['DY'] = 0
                return

            # create temporary DataFrame for magnets with coefficients already in place
            df = tfs.TfsDataFrame(index=tw.index).join(
                coeff_fun(tw.loc[:, 'BETX'], tw.Q1)).join(
                coeff_fun(tw.loc[:, 'BETY'], tw.Q2))
            df.columns = ['COEFFX', 'COEFFY']

            LOG.debug("  Calculate uncoupled linear dispersion")
            df.loc[my_mask, 'DX'] = df.loc[my_mask, 'COEFFX'] * \
                                    sum_fun(tw.loc[mx_mask, 'K0L'],
                                            0,
                                            0,
                                            tw.loc[mx_mask, 'BETX'],
                                            tau(phs_adv['X'].loc[mx_mask, my_mask], tw.Q1)
                                            ).transpose()
            df.loc[mx_mask, 'DY'] = df.loc[mx_mask, 'COEFFY'] * \
                                    sum_fun(-tw.loc[my_mask, 'K0SL'],  # MINUS!
                                            0,
                                            0,
                                            tw.loc[my_mask, 'BETY'],
                                            tau(phs_adv['Y'].loc[my_mask, mx_mask], tw.Q2)
                                            ).transpose()

            LOG.debug("  Calculate full linear dispersion values")
            res_df.loc[:, 'DX'] = df.loc[:, 'COEFFX'] * \
                                sum_fun(tw.loc[mx_mask, 'K0L'],
                                        tw.loc[mx_mask, 'K1SL'],
                                        df.loc[mx_mask, 'DY'],
                                        tw.loc[mx_mask, 'BETX'],
                                        tau(phs_adv['X'].loc[mx_mask, :], tw.Q1)
                                        ).transpose()
            res_df.loc[:, 'DY'] = df.loc[:, 'COEFFY'] * \
                                sum_fun(-tw.loc[my_mask, 'K0SL'],  # MINUS!
                                        tw.loc[my_mask, 'K1SL'],
                                        df.loc[my_mask, 'DX'],
                                        tw.loc[my_mask, 'BETY'],
                                        tau(phs_adv['Y'].loc[my_mask, :], tw.Q2)
                                        ).transpose()

        LOG.debug("  Average linear dispersion Dx: {:g}".format(
                                                    np.mean(res_df.loc[:, 'DX'])))
        LOG.debug("  Average linear dispersion Dy: {:g}".format(
                                                    np.mean(res_df.loc[:, 'DY'])))

    def get_lin_dispersion(self):
        """ Return dataframe with rdts """
        return self._results_df.loc[:, ["S", "DX", "DY"]]

    def plot_lin_dispersion(self, combined=True, madx=False, error=False):
        tw = self.mad_twiss
        lin_disp = self.get_lin_dispersion().dropna()
        pstyle.set_style('standard',
                         {u'lines.linestyle': '-',
                          u'lines.marker': ''})

        if combined:
            ax_dx = lin_disp.plot(x='S')
            ax_dx.set_title('Linear Dispersion')
            pstyle.set_yaxis_label('dispersion', 'x,y', ax_dx)
            ax_dy = ax_dx
        else:
            ax_dx = lin_disp.plot(x='S', y='DX')
            ax_dx.set_title('Linear Dispersion')
            pstyle.set_yaxis_label('dispersion', 'x', ax_dx)

            ax_dy = lin_disp.plot(x='S', y='DY')
            ax_dy.set_title('Linear Dispersion')
            pstyle.set_yaxis_label('dispersion', 'y', ax_dy)

        if madx:
            ax_dx.plot(tw['S'], tw['DX'], label='DX MADX', linestyle='--')
            ax_dy.plot(tw['S'], tw['DY'], label='DY MADX', linestyle='--')

        for ax in (ax_dx, ax_dy):
            self._nice_axes(ax)
            ax.legend()

        if error:
            _, ax_errx = plt.subplots()
            ax_errx.set_title('Dispersion Estimation Error -  Analytical to MADX')
            if combined:
                ax_erry = ax_errx
                pstyle.set_yaxis_label('dispersion', 'x,y', ax_errx, delta=True)
            else:
                _, ax_erry = plt.subplots()
                pstyle.set_yaxis_label('dispersion', 'x', ax_errx, delta=True)
                pstyle.set_yaxis_label('dispersion', 'y', ax_erry, delta=True)
                ax_erry.set_title('Dispersion Estimation Error -  Analytical to MADX')

            ax_errx.plot(lin_disp['S'], (tw.loc[lin_disp.index, 'DX'] - lin_disp['DX']),
                         label='Error DX')
            ax_erry.plot(lin_disp['S'], (tw.loc[lin_disp.index, 'DY'] - lin_disp['DY']),
                         label='Error DY')

            for ax in (ax_errx, ax_erry):
                self._nice_axes(ax)
                ax.legend()

    ################################
    #      Phase Adv Shifts
    ################################

    def calc_shift_phase(self):
        pass

    def get_shift_phase(self):
        pass

    def plot_shift_phase(self):
        pass

    ################################
    #     Linear Chromaticity
    ################################

    def calc_lin_chromaticity(self):
        pass

    def get_lin_chromaticity(self):
        pass

    def plot_lin_chromaticity(self):
        pass

    ################################
    #     Chromatic Beating
    ################################

    def calc_chrom_beat(self):
        pass

    def get_chrom_beat(self):
        pass

    def plot_chrom_beat(self):
        pass

    ################################
    #  Chromatic Phase Adv Shifts
    ################################

    def calc_chrom_shift_phase(self):
        pass

    def get_chrom_shift_phase(self):
        pass

    def plot_chrom_shift_phase(self):
        pass

    ################################
    #  Second Order Dispersion
    ################################

    def calc_sec_dispersion(self):
        pass

    def get_sec_dispersion(self):
        pass

    def plot_sec_dispersion(self):
        pass

    ################################
    #  Chromatic Coupling
    ################################

    def calc_chrom_coupling(self):
        pass

    def get_chrom_coupling(self):
        pass

    def plot_chrom_coupling(self):
        pass

    ################################
    #       Class Helpers
    ################################

    @staticmethod
    def _dphi(data, q):
        return data + np.where(data < 0, q, 0)

    @staticmethod
    def _tau(data, q):
        return data + np.where(data < 0, q / 2, -q / 2)

    def _nice_axes(self, ax):
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 3))
        pstyle.set_xaxis_label(ax)
        try:
            pstyle.set_xLimits(self.mad_twiss.SEQUENCE, ax)
        except pstyle.ArgumentError:
            pass
        if self._ip_pos is not None and len(self._ip_pos) > 0:
            pstyle.show_ir(self._ip_pos, ax)


################################
#           Helpers
################################


def assertion(condition, exception):
    """ Raise Exception if condition is not fulfilled """
    if not condition:
        raise exception


def regex_in(regex, lst):
    """ Return boolean array of length lst, determining if that element starts with regex """
    return np.array([re.search(regex, element) is not None for element in lst])


def get_all_rdts(n):
    """ Returns list of all valid RDTs of order 2 to n """
    assertion(n > 1, ValueError("'n' must be greater 1 for resonance driving terms."))
    permut = [x for x in itertools.product(range(n + 1), repeat=4)
           if 1 < sum(x) <= n and not (x[0] == x[1] and x[2] == x[3])]
    return ['f{:d}{:d}{:d}{:d}'.format(j, k, l, m) for j, k, l, m in sorted(permut, key=sum)]


################################
#           MAIN
################################

def do_main():
    ''' Main function launcher '''
    model_2oct = os.path.abspath('tests/twiss_2oct_elements.dat')
    model_test = "/media/jdilly/Storage/Projects/2017_11_Analytical_Apporox_of_Beamparam/30_madx_tests/b1.twissBlubb"
    test_optics = TwissOptics(model_test)
    # test_optics.calc_rdts(['f1001', 'f1010'])
    test_optics.calc_lin_dispersion()
    test_optics.plot_lin_dispersion(combined=False, madx=True, error=True)
    plt.draw_all()
    plt.show()


if __name__ == '__main__':
    do_main()
