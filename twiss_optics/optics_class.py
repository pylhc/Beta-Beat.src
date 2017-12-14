"""
    Provides Classes to calculate optics from twiss parameters.
    The calculation is based on formulas in [1,2].

    Only works properly for on-orbit twiss files.

    Resonance Driving Terms: Eq. A8 in [2]
    Linear Dispersion: Eq. 24 in [2]
    Linear Chromaticity: Eq. 31 in [2]
    Chromatic Beating: Eq. 36 in [2]


    [1]  A. Franchi et al.,
         First simultaneous measurement of sextupolar and octupolar resonance driving terms in a
         circular accelerator from turn-by-turn beam position monitor data,
         Phys. Rev. ST Accel. Beams 17, 074001 (2014).
         https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.17.074001

    [2]  A. Franchi et al.,
         Analytic formulas for the rapid evaluation of the orbit response matrix and chromatic
         functions from lattice parameters in circular accelerators
         NOT YET PUBLISHED

"""

import os
import numpy as np
import pandas as pd
from math import factorial
import matplotlib.pyplot as plt
from Utilities.plotting import plot_style as pstyle
from Utilities import logging_tools as logtool
from Utilities import tfs_pandas as tfs
from Utilities.contexts import timeit
from twiss_optics.twiss_functions import get_phase_advances, tau, dphi
from twiss_optics.twiss_functions import assertion, regex_in, get_all_rdts

LOG = logtool.get_logger(__name__)

################################
#        TwissOptics
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

        self._add_result_dataframes()  # self._results_df
        self._add_element_types()  # self._elements, self._elements_mapped

        self._phase_advance = get_phase_advances(self.mad_twiss)
        self._define_plot_style()

    def _add_ip_pos(self):
        tw = self.mad_twiss
        self._ip_pos = tw.loc[regex_in(r"\AIP\d$", tw.index), 'S']

    def _remove_nonnecessaries(self):
        LOG.debug("Removing non-necessaries:")
        LOG.debug("  Entries total: {:d}".format(self.mad_twiss.shape[0]))
        with timeit(lambda t:
                    LOG.debug("  Removed in {:f}s".format(t))):
            self.mad_twiss = self.mad_twiss.loc[regex_in(r"\A[MB]", self.mad_twiss.index), :]
        LOG.debug("  Entries left: {:d}".format(self.mad_twiss.shape[0]))

    def _add_result_dataframes(self):
        LOG.debug("Adding results dataframes.")
        self._results_df = tfs.TfsDataFrame(index=self.mad_twiss.index)
        self._results_df["S"] = self.mad_twiss["S"]

    def _add_element_types(self):
        LOG.debug("Sorting Elements into types:")

        tw = self.mad_twiss
        with timeit(lambda t:
                    LOG.debug("  Sorted in {:f}s".format(t))):
            self._elements = {
                'BPM': regex_in(r'\ABPM', tw.index),
                'MB': regex_in(r'\AMB[^S]', tw.index),
                'MBS': regex_in(r'\AMBS', tw.index),
                'MQ': regex_in(r'\AMQ[^ST]', tw.index),
                'MQS': regex_in(r'\AMQS', tw.index),
                'MQT': regex_in(r'\AMQT', tw.index),
                'MS': regex_in(r'\AMS[^S]', tw.index),
                'MSS': regex_in(r'\AMSS', tw.index),
                'MO': regex_in(r'\AMO[^S]', tw.index),
                'MOS': regex_in(r'\AMOS', tw.index),
                'MCB': regex_in(r'\AMCB', tw.index),
                'MCS': regex_in(r'\AMCS', tw.index),
                'MCO': regex_in(r'\AMCO', tw.index),
                'MCD': regex_in(r'\AMCD', tw.index),
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

    def _define_plot_style(self):
        self._plot_style = 'standard'
        self._plot_manual = {u'lines.linestyle': '-',
                             u'lines.marker': '',
                             }

    ################################
    #         Properties
    ################################

    def get_coupling(self, method='rdt'):
        if method == 'rdt':
            # available after calc_rdts
            return self._results_df[['S', 'F1001', 'F1010']]
        elif method == 'cmatrix':
            # available after calc_cmatrix
            return self._results_df[['S', 'F1001_C', 'F1010_C']].rename(
                columns=lambda x: x.replace("_C", ""))
        else:
            raise ValueError("method '{:s}' not recognized.".format(method))

    def get_K(self, order):
        return self._results_df[['S', 'K{0}L'.format(order)]]

    def set_K(self, order, value):
        self._results_df[['S', 'K{0}L'.format(order)]] = value

    def get_J(self, order):
        return self._results_df[['S', 'K{0}SL'.format(order)]]

    def set_J(self, order, value):
        self._results_df[['S', 'K{0}SL'.format(order)]] = value

    def get_S(self):
        return self._results_df[["S"]]

    ################################
    #          C Matrix
    ################################

    def calc_cmatrix(self):
        """
        Calculates C matrix and Coupling and Gamma from it
        :return:
        """
        tw = self.mad_twiss
        res = self._results_df

        LOG.info("Calculating CMatrix.")
        with timeit(lambda t:
                    LOG.debug("  CMatrix calculated in {:f}s".format(t))):

            j = np.array([[0., 1.],
                          [-1., 0.]])
            rs = np.reshape(tw.as_matrix(columns=["R11", "R12",
                                                  "R21", "R22"]),
                            (len(tw), 2, 2))
            cs = np.einsum("ij,kjn,no->kio",
                           -j, np.transpose(rs, axes=(0, 2, 1)), j)
            cs = np.einsum("k,kij->kij", (1 / np.sqrt(1 + np.linalg.det(rs))), cs)

            g11a = 1 / np.sqrt(tw.loc[:, "BETX"])
            g12a = np.zeros(len(tw))
            g21a = tw.loc[:, "ALFX"] / np.sqrt(tw.loc[:, "BETX"])
            g22a = np.sqrt(tw.loc[:, "BETX"])
            gas = np.reshape(np.array([g11a, g12a,
                                       g21a, g22a]).T,
                             (len(tw), 2, 2))

            ig11b = np.sqrt(tw.loc[:, "BETY"])
            ig12b = np.zeros(len(tw))
            ig21b = -tw.loc[:, "ALFY"] / np.sqrt(tw.loc[:, "BETY"])
            ig22b = 1. / np.sqrt(tw.loc[:, "BETY"])
            igbs = np.reshape(np.array([ig11b, ig12b,
                                        ig21b, ig22b]).T,
                              (len(tw), 2, 2))
            cs = np.einsum("kij,kjl,kln->kin", gas, cs, igbs)
            gammas = np.sqrt(1 - np.linalg.det(cs))

            res.loc[:, "GAMMA_C"] = gammas
            res.loc[:, "F1001_C"] = ((cs[:, 0, 0] + cs[:, 1, 1]) * 1j +
                                     (cs[:, 0, 1] - cs[:, 1, 0])) / 4 / gammas
            res.loc[:, "F1010_C"] = ((cs[:, 0, 0] - cs[:, 1, 1]) * 1j +
                                     (-cs[:, 0, 1]) - cs[:, 1, 0]) / 4 / gammas

            res.loc[:, "C11"] = cs[:, 0, 0]
            res.loc[:, "C12"] = cs[:, 0, 1]
            res.loc[:, "C21"] = cs[:, 1, 0]
            res.loc[:, "C22"] = cs[:, 1, 1]

            LOG.debug("  Average coupling amplitude |F1001|: {:g}".format(np.mean(
                np.abs(res.loc[:, "F1001_C"]))))
            LOG.debug("  Average coupling amplitude |F1010|: {:g}".format(np.mean(
                np.abs(res.loc[:, "F1010_C"]))))
            LOG.debug("  Average gamma: {:g}".format(np.mean(
                np.abs(res.loc[:, "GAMMA_C"]))))

        self._log_added('GAMMA_C', 'F1001_C', 'F1010_C', 'C11', 'C12', 'C21', 'C22')

    def get_gamma(self):
        """ available after calc_cmatrix """
        return self._results_df[['S', 'GAMMA_C']].rename(columns=lambda x: x.replace("_C", ""))

    def get_cmatrix(self):
        """ available after calc_cmatrix """
        return self._results_df[['S', 'C11', 'C12', 'C21', 'C22']]

    ################################
    #   Resonance Driving Terms
    ################################

    def calc_rdts(self, order):
        """ Eq. A8 in [2] """
        if isinstance(order, int):
            order = get_all_rdts(order)
        elif not isinstance(order, list):
            order = [order]

        LOG.info("Calculating RDTs: {:s}.".format(str(order)[1:-1]))
        with timeit(lambda t:
                    LOG.debug("  RDTs calculated in {:f}s".format(t))):

            i2pi = 2j * np.pi
            tw = self.mad_twiss
            phs_adv = self._phase_advance
            res = self._results_df

            for rdt in order:
                assertion(len(rdt) == 5 and rdt[0].upper() == 'F',
                          ValueError("'{:s}' does not seem to be a valid RDT name.".format(rdt)))

                conj_rdt = ''.join(['F', rdt[2], rdt[1], rdt[4], rdt[3]])

                if conj_rdt in self._results_df:
                    res[rdt.upper()] = np.conjugate(self._results_df[conj_rdt])
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
                        phx = dphi(phs_adv['X'].loc[k_mask, el_mask], tw.Q1)
                        phy = dphi(phs_adv['Y'].loc[k_mask, el_mask], tw.Q2)
                        phase_term = ((j-k) * phx + (l-m) * phy).applymap(lambda p: np.exp(i2pi*p))

                        beta_term = tw.loc[k_mask, src] * \
                                    tw.loc[k_mask, 'BETX'] ** ((j+k) / 2.) * \
                                    tw.loc[k_mask, 'BETY'] ** ((l+m) / 2.)

                        res.loc[el_mask, rdt.upper()] = sign * phase_term.multiply(
                            beta_term, axis="index").sum(axis=0).transpose() * denom

                        LOG.debug("  Average RDT amplitude |{:s}|: {:g}".format(rdt, np.mean(
                            np.abs(res.loc[el_mask, rdt.upper()]))))
                    else:
                        LOG.debug("  All {:s} == 0. RDT '{:s}' will be zero.".format(src, rdt))
                        res.loc[el_mask, rdt.upper()] = 0

        self._log_added(*order)

    def get_rdts(self, rdt_names=None):
        """ available after calc_rdts """
        if rdt_names:
            if not isinstance(rdt_names, list):
                rdt_names = [rdt_names]
            return self._results_df[["S"] + [rdt.upper() for rdt in rdt_names]]
        else:
            return self._results_df.loc[:, regex_in(r'\A(S|F\d{4})$', self._results_df.columns)]

    def plot_rdts(self, rdt_names=None, apply_fun=np.abs, combined=True):
        """ available after calc_rdts """
        LOG.info("Plotting Resonance Driving Terms")
        rdts = self.get_rdts(rdt_names)
        is_s = regex_in(r'\AS$', rdts.columns)
        rdts = rdts.dropna()
        rdts.loc[:, ~is_s] = rdts.loc[:, ~is_s].applymap(apply_fun)
        pstyle.set_style(self._plot_style, self._plot_manual)

        if combined:
            ax = rdts.plot(x='S')
            ax.set_title('Resonance Driving Terms')
            pstyle.small_title(ax)
            pstyle.set_name('Resonance Driving Terms', ax)
            pstyle.set_yaxis_label(apply_fun.__name__, 'F_{{jklm}}', ax)
            self._nice_axes(ax)
        else:
            for rdt in rdts.loc[:, ~is_s]:
                ax = rdts.plot(x='S', y=rdt)
                ax.set_title('Resonance Driving Term ' + rdt)
                pstyle.small_title(ax)
                pstyle.set_name('Resonance Driving Term ' + rdt, ax)
                pstyle.set_yaxis_label(apply_fun.__name__, rdt, ax)
                self._nice_axes(ax)

    ################################
    #      Linear Dispersion
    ################################

    def calc_linear_dispersion(self, bpms_only=False):
        """ Eq. 24 in [2] """
        tw = self.mad_twiss
        phs_adv = self._phase_advance
        res = self._results_df
        coeff_fun = self._linear_dispersion_coeff
        sum_fun = self._linear_dispersion_sum

        # Calculate
        LOG.info("Calculate Linear Dispersion")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            # sources
            k0_mask = tw['K0L'] != 0
            k0s_mask = tw['K0SL'] != 0
            k1s_mask = tw['K1SL'] != 0

            mx_mask = k0_mask | k1s_mask  # magnets contributing to Dx,j (-> Dy,m)
            my_mask = k0s_mask | k1s_mask  # magnets contributing to Dy,j (-> Dx,m)

            if not any(mx_mask | my_mask):
                LOG.warning("  No linear dispersion contributions found. Values will be NaN.")
                res['DX'] = np.nan
                res['DY'] = np.nan
                return

            # results
            if bpms_only:
                res_mask = self._elements["BPM"]
            else:
                res_mask = np.ones(res.shape[0], dtype=bool)

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
            res.loc[res_mask, 'DX'] = df.loc[res_mask, 'COEFFX'] * \
                                sum_fun(tw.loc[mx_mask, 'K0L'],
                                        tw.loc[mx_mask, 'K1SL'],
                                        df.loc[mx_mask, 'DY'],
                                        tw.loc[mx_mask, 'BETX'],
                                        tau(phs_adv['X'].loc[mx_mask, res_mask], tw.Q1)
                                        ).transpose()
            res.loc[res_mask, 'DY'] = df.loc[res_mask, 'COEFFY'] * \
                                sum_fun(-tw.loc[my_mask, 'K0SL'],  # MINUS!
                                        tw.loc[my_mask, 'K1SL'],
                                        df.loc[my_mask, 'DX'],
                                        tw.loc[my_mask, 'BETY'],
                                        tau(phs_adv['Y'].loc[my_mask, res_mask], tw.Q2)
                                        ).transpose()

        LOG.debug("  Average linear dispersion Dx: {:g}".format(
                                                    np.mean(res['DX'])))
        LOG.debug("  Average linear dispersion Dy: {:g}".format(
                                                    np.mean(res['DY'])))
        self._log_added('DX', 'DY')

    def get_linear_dispersion(self):
        """ available after calc_linear_dispersion """
        return self._results_df.loc[:, ["S", "DX", "DY"]]

    def plot_linear_dispersion(self, combined=True):
        """ available after calc_linear_dispersion """
        LOG.info("Plotting Linear Dispersion")
        lin_disp = self.get_linear_dispersion().dropna()
        title = 'Linear Dispersion'
        pstyle.set_style(self._plot_style, self._plot_manual)

        if combined:
            ax_dx = lin_disp.plot(x='S')
            ax_dx.set_title(title)
            pstyle.set_yaxis_label('dispersion', 'x,y', ax_dx)
            ax_dy = ax_dx
        else:
            ax_dx = lin_disp.plot(x='S', y='DX')
            ax_dx.set_title(title)
            pstyle.set_yaxis_label('dispersion', 'x', ax_dx)

            ax_dy = lin_disp.plot(x='S', y='DY')
            ax_dy.set_title(title)
            pstyle.set_yaxis_label('dispersion', 'y', ax_dy)

        for ax in (ax_dx, ax_dy):
            self._nice_axes(ax)
            ax.legend()

    # helpers
    @staticmethod
    def _linear_dispersion_coeff(beta, q):
        return np.sqrt(beta) / (2 * np.sin(np.pi * q))

    @staticmethod
    def _linear_dispersion_sum(k, j, d, beta, tau):
        # k, j, d , beta = columns -> convert to Series -> broadcasted
        # tau = Matrix as Frame
        calc_column = (k + j * d) * np.sqrt(beta)
        if isinstance(calc_column, pd.DataFrame):
            calc_column = calc_column.squeeze()
        return np.cos(2 * np.pi * tau).mul(calc_column, axis='index').sum(axis='index')

    ################################
    #      Phase Adv Shifts
    ################################

    def plot_phase_advance(self, combined=True):
        """ Plots the phase advance, but only between two consecutive elements """
        LOG.info("Plotting Phase Advance")
        tw = self.mad_twiss
        pa = self._phase_advance
        dpa = self._dphase_advance
        phase_advx = np.append(pa['X'].iloc[0, -1] + tw.Q1, pa['X'].values.diagonal(offset=-1))
        dphase_advx = np.append(dpa['X'].iloc[0, -1], dpa['X'].values.diagonal(offset=-1))
        phase_advy = np.append(pa['Y'].iloc[0, -1] + tw.Q2, pa['Y'].values.diagonal(offset=-1))
        dphase_advy = np.append(dpa['Y'].iloc[0, -1], dpa['Y'].values.diagonal(offset=-1))
        phase_adv = tw[["S"]].copy()
        phase_adv['MUX'] = np.cumsum(phase_advx + dphase_advx) % 1 - .5
        phase_adv['MUY'] = np.cumsum(phase_advy + dphase_advy) % 1 - .5

        title = 'Phase'
        pstyle.set_style(self._plot_style, self._plot_manual)

        if combined:
            ax_dx = phase_adv.plot(x='S')
            ax_dx.set_title(title)
            pstyle.small_title(ax_dx)
            pstyle.set_name(title, ax_dx)
            pstyle.set_yaxis_label('phase', 'x,y', ax_dx, delta=False)
            ax_dy = ax_dx
        else:
            ax_dx = phase_adv.plot(x='S', y='MUX')
            ax_dx.set_title(title)
            pstyle.small_title(ax_dx)
            pstyle.set_name(title, ax_dx)
            pstyle.set_yaxis_label('phase', 'x', ax_dx, delta=False)

            ax_dy = phase_adv.plot(x='S', y='MUY')
            ax_dy.set_title(title)
            pstyle.small_title(ax_dy)
            pstyle.set_name(title, ax_dy)
            pstyle.set_yaxis_label('phase', 'y', ax_dy, delta=False)

        for ax in (ax_dx, ax_dy):
            self._nice_axes(ax)
            ax.legend()

    ################################
    #     Linear Chromaticity
    ################################

    def calc_linear_chromaticity(self):
        """ Eq. 31 in [2] """
        res = self._results_df

        if 'CHROMX' not in res:
            self._calc_chromatic_term()

        LOG.info("Calculating Linear Chromaticity")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            DQ1 = - 1/(4 * np.pi) * res['CHROMX'].dropna().sum(axis="index")
            DQ2 = 1/(4 * np.pi) * res['CHROMY'].dropna().sum(axis="index")

        self._results_df.DQ1 = DQ1
        self._results_df.DQ2 = DQ2
        LOG.debug("  Q'x: {:f}".format(DQ1))
        LOG.debug("  Q'y: {:f}".format(DQ2))

        self._log_added('DQ1', 'DQ2')

    def get_linear_chromaticity(self):
        """ available after calc_linear_chromaticity """
        return self._results_df.DQ1, self._results_df.DQ2

    ################################
    #     Chromatic Beating
    ################################

    def calc_chromatic_beating(self):
        """ Eq. 36 in [2] """
        tw = self.mad_twiss
        res = self._results_df
        phs_adv = self._phase_advance

        if 'CHROMX' not in res:
            self._calc_chromatic_term()

        LOG.info("Calculating Chromatic Beating")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            chromx = res['CHROMX'].dropna()
            chromy = res['CHROMY'].dropna()
            res['DBEATX'] = self._chromatic_beating(
                chromx,
                tau(phs_adv['X'].loc[chromx.index, :], tw.Q1),
                tw.Q1).transpose() - 1
            res['DBEATY'] = - self._chromatic_beating(
                chromy,
                tau(phs_adv['Y'].loc[chromy.index, :], tw.Q2),
                tw.Q2).transpose() - 1

        LOG.debug("  Pk2Pk chromatic beating DBEATX: {:g}".format(
            res['DBEATX'].max() - res['DBEATX'].min()))
        LOG.debug("  Pk2Pk chromatic beating DBEATY: {:g}".format(
            res['DBEATY'].max() - res['DBEATY'].min()))
        self._log_added('DBEATX', 'DBEATY')

    def get_chromatic_beating(self):
        """ available after calc_chromatic_beating """
        return self._results_df.loc[:, ["S", "DBEATX", "DBEATY"]]

    def plot_chromatic_beating(self, combined=True):
        """ available after calc_chromatic_beating """
        LOG.info("Plotting Chromatic Beating")
        chrom_beat = self.get_chromatic_beating().dropna()
        title = 'Chromatic Beating'
        pstyle.set_style(self._plot_style, self._plot_manual)

        if combined:
            ax_dx = chrom_beat.plot(x='S')
            ax_dx.set_title(title)
            pstyle.small_title(ax_dx)
            pstyle.set_name(title, ax_dx)
            pstyle.set_yaxis_label('dbetabeat', 'x,y', ax_dx)
            ax_dy = ax_dx
        else:
            ax_dx = chrom_beat.plot(x='S', y='DBEATX')
            ax_dx.set_title(title)
            pstyle.small_title(ax_dx)
            pstyle.set_name(title, ax_dx)
            pstyle.set_yaxis_label('dbetabeat', 'x', ax_dx)

            ax_dy = chrom_beat.plot(x='S', y='DBEATY')
            ax_dy.set_title(title)
            pstyle.small_title(ax_dy)
            pstyle.set_name(title, ax_dy)
            pstyle.set_yaxis_label('dbetabeat', 'y', ax_dy)

        for ax in (ax_dx, ax_dy):
            self._nice_axes(ax)
            ax.legend()

    @staticmethod
    def _chromatic_beating(chrom_term, tau, q):
        return 1 / (2 * np.sin(2 * np.pi * q)) * \
               np.cos(4 * np.pi * tau).mul(chrom_term, axis='index').sum(axis='index')

    ################################
    #  Chromatic Phase Adv Shifts (NA)
    ################################

    def calc_chromatic_phase_adv_shift(self):
        raise NotImplementedError('Chramatic Phase Advance Shift is not Implemented yet.')

    def get_chromatic_phase_adv_shift(self):
        raise NotImplementedError('Chramatic Phase Advance Shift is not Implemented yet.')

    def plot_chromatic_phase_adv_shift(self):
        raise NotImplementedError('Chramatic Phase Advance Shift is not Implemented yet.')

    ################################
    #  Second Order Dispersion (NA)
    ################################

    def calc_second_order_dispersion(self):
        """ Eq. 11, 41 - 43 in [2] """
        raise NotImplementedError('Second Order Dispersion is not Implemented yet.')

    def get_second_order_dispersion(self):
        raise NotImplementedError('Second Order Dispersion is not Implemented yet.')

    def plot_second_order_dispersion(self):
        raise NotImplementedError('Second Order Dispersion is not Implemented yet.')

    ################################
    #  Chromatic Coupling (NA)
    ################################

    def calc_chromatic_coupling(self):
        """ Eq. 47, 48, B77 in [2] """
        raise NotImplementedError('Chromatic Coupling is not Implemented yet.')

    def get_chromatic_coupling(self):
        raise NotImplementedError('Chromatic Coupling is not Implemented yet.')

    def plot_chromatic_coupling(self):
        raise NotImplementedError('Chromatic Coupling is not Implemented yet.')


    ################################
    #       Class Helpers
    ################################

    def _calc_chromatic_term(self):
        LOG.info("Calculating Chromatic Term.")
        res = self._results_df
        tw = self.mad_twiss

        #TODO: Decide weather DX, DY from Model or from calc

        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            mask = (tw['K1L'] != 0) | (tw['K2L'] != 0) | (tw['K2SL'] != 0)
            sum_term = tw.loc[mask, 'K1L'] - \
                       (tw.loc[mask, 'K2L'] * tw.loc[mask, 'DX']) + \
                       (tw.loc[mask, 'K2SL'] * tw.loc[mask, 'DY'])
            res['CHROMX'] = sum_term * tw.loc[mask, 'BETX']
            res['CHROMY'] = sum_term * tw.loc[mask, 'BETY']

        LOG.debug("Chromatic Term Calculated.")

    @staticmethod
    def _dphi(data, q):
        return data + np.where(data <= 0, q, 0)  # '<=' seems to be what MAD-X does

    @staticmethod
    def _tau(data, q):
        return data + np.where(data <= 0, q / 2, -q / 2)  # '<=' seems to be what MAD-X does

    def _nice_axes(self, ax):
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 3))
        pstyle.set_xaxis_label(ax)
        try:
            pstyle.set_xLimits(self.mad_twiss.SEQUENCE, ax)
        except pstyle.ArgumentError:
            pass
        if self._ip_pos is not None and len(self._ip_pos) > 0:
            pstyle.show_ir(self._ip_pos, ax)

    @staticmethod
    def _log_added(*args):
        if len(args) > 0:
            fields = "'" + "', '".join(args) + "'"
            LOG.debug("  Added fields to results: " + fields)




################################
#           MAIN
################################

def do_main():
    """ Main function launcher """
    model_2oct = os.path.abspath('tests/twiss_2oct_elements.dat')
    model_test = "/media/jdilly/Storage/Projects/2017_11_Analytical_Apporox_of_Beamparam/30_madx_tests/b1.twissKick"
    test_optics = TwissOptics(model_test)
    # test_optics.calc_rdts(['f1001', 'f1010'])
    # test_optics.calc_linear_dispersion()
    # test_optics.plot_linear_dispersion(combined=False)
    test_optics.calc_chromatic_beating()
    test_optics.plot_chromatic_beating()
    # test_optics.plot_phase_advance()

    plt.draw_all()
    plt.show()
    return

    nodk = "/media/jdilly/Storage/Projects/2017_11_Analytical_Apporox_of_Beamparam/30_madx_tests/b1.twissNoDK"
    dk = "/media/jdilly/Storage/Projects/2017_11_Analytical_Apporox_of_Beamparam/30_madx_tests/b1.twissDK"

    nodk_opt = TwissOptics(nodk)
    dk_opt = TwissOptics(dk)

    # nodk_opt.set_K(1, dk_opt.get_K(1))
    # nodk_opt.calc_phase_advance_shift()
    # nodk_opt.plot_phase_advance()
    # test_optics.calc_phase_advance_shift()


    dk_opt.calc_rdts(['f1001', 'f1010'])
    dk_opt.plot_rdts(['f1001', 'f1010'])


    plt.draw_all()
    plt.show()



if __name__ == '__main__':
    do_main()
