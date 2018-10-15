"""
Provides Classes to calculate optics from twiss parameters.

The calculation is based on formulas in [#AibaFirstbetabeatingmeasurement2009]_ and
[#FranchiAnalyticformulasrapid2017]_.

Only works properly for on-orbit twiss files.

 - Resonance Driving Terms: Eq. A8 in [#FranchiAnalyticformulasrapid2017]_
 - Linear Dispersion: Eq. 24 in [#FranchiAnalyticformulasrapid2017]_
 - Linear Chromaticity: Eq. 31 in [#FranchiAnalyticformulasrapid2017]_
 - Chromatic Beating: Eq. 36 in [#FranchiAnalyticformulasrapid2017]_

.. rubric:: References

.. [#AibaFirstbetabeatingmeasurement2009]
    M. Aiba et al.,
    First simultaneous measurement of sextupolar and octupolar resonance driving
    terms in a circular accelerator from turn-by-turn beam position monitor data,
    Phys. Rev. ST Accel. Beams 17, 074001 (2014).
    https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.17.074001

.. [#FranchiAnalyticformulasrapid2017]
    A. Franchi et al.,
    Analytic formulas for the rapid evaluation of the orbit response matrix
    and chromatic functions from lattice parameters in circular accelerators
    https://arxiv.org/abs/1711.06589

.. [#CalagaBetatroncouplingMerging2005]
    R. Calaga et al.,
    'Betatron coupling: Merging Hamiltonian and matrix approaches'
    Phys. Rev. ST Accel. Beams, vol. 8, no. 3, p. 034001, Mar. 2005.

"""

from math import factorial

import numpy as np
import pandas as pd

from twiss_optics.twiss_functions import assertion, get_all_rdts
from twiss_optics.twiss_functions import get_phase_advances, tau, dphi
from utils import logging_tools as logtool
from tfs_files import tfs_pandas as tfs
from utils.contexts import timeit
from utils.dict_tools import DotDict
from plotshop import plot_style as pstyle

LOG = logtool.get_logger(__name__)

PLOT_DEFAULTS = {
        "style": 'standard',
        "manual": {u'lines.linestyle': '-',
                   u'lines.marker': '',
                   }
}


################################
#        TwissOptics
################################


class TwissOptics(object):
    """ Class for calculating optics parameters from twiss-model.

    Args:
        model_path_or_df: Path to twissfile of model or DataFrame of model.
        quick_init: Initializes without calculating phase advances. Default: False
    """

    ################################
    #       init Functions
    ################################

    def __init__(self, model_path_or_df, quick_init=True):
        self.twiss_df = self._get_model_df(model_path_or_df)
        self._ip_pos = self._find_ip_positions()

        self._results_df = self._make_results_dataframe()

        self._phase_advance = None
        if not quick_init:
            self._phase_advance = self.get_phase_adv()

        self._plot_options = DotDict(PLOT_DEFAULTS)

    @staticmethod
    def _get_model_df(model_path_or_tfs):
        """ Check if DataFrame given, if not load model from file  """
        if isinstance(model_path_or_tfs, basestring):
            LOG.debug("Creating TwissOptics from '{:s}'".format(model_path_or_tfs))
            df = tfs.read_tfs(model_path_or_tfs, index="NAME")
        else:
            LOG.debug("Creating TwissOptics from input DataFrame")
            df = model_path_or_tfs
            if (len(df.index.values) == 0) or not isinstance(df.index.values[0], basestring):
                raise IndexError("Index of DataFrame needs to be the element names."
                                 "This does not seem to be the case.")
        return df

    def _find_ip_positions(self):
        """ Returns IP positions from Dataframe.
        Only needed for plotting, so if not present it will just skip it.
        Nice for LHC though.

        Load model into twiss_df first!
        """
        tw = self.twiss_df
        return tw.loc[tw.index.str.match(r"IP\d$", case=False), 'S']

    def _make_results_dataframe(self):
        """ Creating a dataframe used for storing results. """
        LOG.debug("Creating Results Dataframes.")
        results_df = tfs.TfsDataFrame(index=self.twiss_df.index)
        results_df["S"] = self.twiss_df["S"]
        return results_df

    def define_plot_style(self, **kwargs):
        """ Edit your desired plot style here """
        if not kwargs:
            options = DotDict(PLOT_DEFAULTS)
        else:
            options = DotDict(kwargs)
            if "style" not in options:
                options.style = PLOT_DEFAULTS["style"]
            if "manual" not in options:
                options.manual = PLOT_DEFAULTS["manual"]
        self._plot_options = options

    ################################
    #         Properties
    ################################

    def get_phase_adv(self):
        """ Wrapper for returning the matrix of phase-advances. """
        if self._phase_advance is None:
            self._phase_advance = get_phase_advances(self.twiss_df)
        return self._phase_advance

    def get_coupling(self, method='rdt'):
        """ Returns the coupling term.


        .. warning::
            This function changes sign of the real part of the RDTs compared to
            [#FranchiAnalyticformulasrapid2017]_ to be consistent with the RDT
            calculations from [#CalagaBetatroncouplingMerging2005]_ .

        Args:
            method: 'rdt' - Returns the values calculated by calc_rdts()
                    'cmatrix' - Returns the values calculated by calc_cmatrix()
        """
        if method == 'rdt':
            if "F1001" not in self._results_df or "F1010" not in self._results_df:
                self.calc_rdts(['F1001', 'F1010'])
                res_df = self._results_df.loc[:, ['S', 'F1001', 'F1010']]
                res_df.loc[:, "F1001"].real *= -1
                res_df.loc[:, "F1010"].real *= -1
            return res_df
        elif method == 'cmatrix':
            if "F1001_C" not in self._results_df:
                self.calc_cmatrix()
            res_df = self._results_df.loc[:, ['S', 'F1001_C', 'F1010_C']]
            return res_df.rename(columns=lambda x: x.replace("_C", ""))
        else:
            raise ValueError("method '{:s}' not recognized.".format(method))

    ################################
    #          C Matrix
    ################################

    def calc_cmatrix(self):
        """ Calculates C matrix and Coupling and Gamma from it.
        See [#CalagaBetatroncouplingMerging2005]_
        """
        tw = self.twiss_df
        res = self._results_df

        LOG.debug("Calculating CMatrix.")
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
        """ Return Gamma values. """
        if 'GAMMA_C' not in self._results_df:
            self.calc_cmatrix()
        res_df = self._results_df.loc[:, ['S', 'GAMMA_C']]
        return res_df.rename(columns=lambda x: x.replace("_C", ""))

    def get_cmatrix(self):
        """ Return the C-Matrix """
        if 'C11' not in self._results_df:
            self.calc_cmatrix()
        return self._results_df.loc[:, ['S', 'C11', 'C12', 'C21', 'C22']]

    ################################
    #   Resonance Driving Terms
    ################################

    def calc_rdts(self, order_or_rdts):
        """ Calculates the Resonance Driving Terms.
        
        Eq. A8 in [#FranchiAnalyticformulasrapid2017]_

        Args:
            order_or_rdts: int, string or list of strings
                If an int is given all Resonance Driving Terms up to this order
                will be calculated.
                The strings are assumed to be the desired driving term names, e.g. "F1001"
        """
        if isinstance(order_or_rdts, int):
            rdt_list = get_all_rdts(order_or_rdts)
        elif not isinstance(order_or_rdts, list):
            rdt_list = [order_or_rdts]
        else:
            rdt_list = order_or_rdts

        LOG.debug("Calculating RDTs: {:s}.".format(str(rdt_list)[1:-1]))
        with timeit(lambda t:
                    LOG.debug("  RDTs calculated in {:f}s".format(t))):

            i2pi = 2j * np.pi
            tw = self.twiss_df
            phs_adv = self.get_phase_adv()
            res = self._results_df

            for rdt in rdt_list:
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

                    try:
                        mask_in = tw[src] != 0
                        if sum(mask_in) == 0:
                            raise KeyError
                    except KeyError:
                        # either src is not in tw or all k's are zero.
                        LOG.warning("  All {:s} == 0. RDT '{:s}' will be zero.".format(src, rdt))
                        res.loc[:, rdt.upper()] = 0
                    else:
                        # the next three lines determine the main order of speed, hence
                        # - mask as much as possible
                        # - additions are faster than multiplications (-> applymap last)
                        phx = dphi(phs_adv['X'].loc[mask_in, :], tw.Q1)
                        phy = dphi(phs_adv['Y'].loc[mask_in, :], tw.Q2)
                        phase_term = ((j-k) * phx + (l-m) * phy).applymap(lambda p: np.exp(i2pi*p))

                        beta_term = tw.loc[mask_in, src] * \
                                    tw.loc[mask_in, 'BETX'] ** ((j+k) / 2.) * \
                                    tw.loc[mask_in, 'BETY'] ** ((l+m) / 2.)

                        res.loc[:, rdt.upper()] = sign * phase_term.multiply(
                            beta_term, axis="index").sum(axis=0).transpose() * denom

                        LOG.debug("  Average RDT amplitude |{:s}|: {:g}".format(rdt, np.mean(
                            np.abs(res.loc[:, rdt.upper()]))))

        self._log_added(*rdt_list)

    def get_rdts(self, rdt_names=None):
        """ Return Resonance Driving Terms. """
        if rdt_names:
            if not isinstance(rdt_names, list):
                rdt_names = [rdt_names]
            rdt_names = [rdt.upper() for rdt in rdt_names]
            new_rdts = [rdt for rdt in rdt_names if rdt not in self._results_df]
            if len(new_rdts) > 0:
                self.calc_rdts(new_rdts)
            return self._results_df.loc[:, ["S"] + rdt_names]
        else:
            return self._results_df.loc[:, self._results_df.columns.str.match(r'(S|F\d{4})$',
                                                                              case=False)]

    def plot_rdts(self, rdt_names=None, apply_fun=np.abs, combined=True):
        """ Plot Resonance Driving Terms """
        LOG.debug("Plotting Resonance Driving Terms")
        rdts = self.get_rdts(rdt_names)
        is_s = rdts.columns.str.match(r'S$', case=False)
        rdts = rdts.dropna()
        rdts.loc[:, ~is_s] = rdts.loc[:, ~is_s].applymap(apply_fun)
        pstyle.set_style(self._plot_options.style, self._plot_options.manual)

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

    def calc_linear_dispersion(self):
        """ Calculate the Linear Disperion.

        Eq. 24 in [#FranchiAnalyticformulasrapid2017]_
        """
        self._check_k_columns(["K0L", "K0SL", "K1SL"])
        tw = self.twiss_df
        phs_adv = self.get_phase_adv()
        res = self._results_df
        coeff_fun = self._linear_dispersion_coeff
        sum_fun = self._linear_dispersion_sum

        # Calculate
        LOG.debug("Calculate Linear Dispersion")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            # sources
            k0_mask = tw['K0L'] != 0
            k0s_mask = tw['K0SL'] != 0
            k1s_mask = tw['K1SL'] != 0

            mx_mask = k0_mask | k1s_mask  # magnets contributing to Dx,j (-> Dy,m)
            my_mask = k0s_mask | k1s_mask  # magnets contributing to Dy,j (-> Dx,m)

            if not any(mx_mask | my_mask):
                LOG.warning("  No linear dispersion contributions found. Values will be zero.")
                res['DX'] = 0.
                res['DY'] = 0.
                self._log_added('DX', 'DY')
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
            res.loc[:, 'DX'] = df.loc[:, 'COEFFX'] * \
                               sum_fun(tw.loc[mx_mask, 'K0L'],
                                       tw.loc[mx_mask, 'K1SL'],
                                       df.loc[mx_mask, 'DY'],
                                       tw.loc[mx_mask, 'BETX'],
                                       tau(phs_adv['X'].loc[mx_mask, :], tw.Q1)
                                       ).transpose()
            res.loc[:, 'DY'] = df.loc[:, 'COEFFY'] * \
                               sum_fun(-tw.loc[my_mask, 'K0SL'],  # MINUS!
                                       tw.loc[my_mask, 'K1SL'],
                                       df.loc[my_mask, 'DX'],
                                       tw.loc[my_mask, 'BETY'],
                                       tau(phs_adv['Y'].loc[my_mask, :], tw.Q2)
                                       ).transpose()

        LOG.debug("  Average linear dispersion Dx: {:g}".format(
                                                    np.mean(res['DX'])))
        LOG.debug("  Average linear dispersion Dy: {:g}".format(
                                                    np.mean(res['DY'])))
        self._log_added('DX', 'DY')

    def get_linear_dispersion(self):
        """ Return the Linear Dispersion.

        Available after calc_linear_dispersion!
        """
        if "DX" not in self._results_df or "DY" not in self._results_df:
            self.calc_linear_dispersion()
        return self._results_df.loc[:, ["S", "DX", "DY"]]

    def plot_linear_dispersion(self, combined=True):
        """ Plot the Linear Dispersion.

        Available after calc_linear_dispersion!

        Args:
            combined (bool): If 'True' plots x and y into the same axes.
        """
        LOG.debug("Plotting Linear Dispersion")
        lin_disp = self.get_linear_dispersion().dropna()
        title = 'Linear Dispersion'
        pstyle.set_style(self._plot_options.style, self._plot_options.manual)

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
        """ Helper to calculate the coefficient """
        return np.sqrt(beta) / (2 * np.sin(np.pi * q))

    @staticmethod
    def _linear_dispersion_sum(k, j, d, beta, tau):
        """ Helper to calculate the sum """
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
        """ Plots the phase advances between two consecutive elements

        Args:
            combined (bool): If 'True' plots x and y into the same axes.
        """
        raise NotImplementedError('Plotting Phase Advance Shift is not Implemented yet.')
        #TODO: reimplement the phase-advance shift calculations (if needed??)
        LOG.debug("Plotting Phase Advance")
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
        pstyle.set_style(self._plot_options.style, self._plot_options.manual)

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
        """ Calculate the Linear Chromaticity

        Eq. 31 in [#FranchiAnalyticformulasrapid2017]_
        """
        res = self._results_df

        if 'CHROMX' not in res:
            self._calc_chromatic_term()

        LOG.debug("Calculating Linear Chromaticity")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            DQ1 = - 1/(4 * np.pi) * res['CHROMX'].dropna().sum(axis="index")
            DQ2 = 1/(4 * np.pi) * res['CHROMY'].dropna().sum(axis="index")

        self._results_df.DQ1 = DQ1
        self._results_df.DQ2 = DQ2
        LOG.debug("  Q'x: {:f}".format(DQ1))
        LOG.debug("  Q'y: {:f}".format(DQ2))

        self._log_added('DQ1', 'DQ2')

    def get_linear_chromaticity(self):
        """ Return the Linear Chromaticity

        Available after calc_linear_chromaticity
        """
        try:
            return self._results_df.DQ1, self._results_df.DQ2
        except AttributeError:
            self.calc_linear_chromaticity()
            return self._results_df.DQ1, self._results_df.DQ2

    ################################
    #     Chromatic Beating
    ################################

    def calc_chromatic_beating(self):
        """ Calculate the Chromatic Beating

        Eq. 36 in [#FranchiAnalyticformulasrapid2017]_
        """
        tw = self.twiss_df
        res = self._results_df
        phs_adv = self.get_phase_adv()

        if 'CHROMX' not in res:
            self._calc_chromatic_term()

        LOG.debug("Calculating Chromatic Beating")
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
        """ Return the Chromatic Beating

         Available after calc_chromatic_beating
         """
        if "DBEATX" not in self._results_df or "DBEATY" not in self._results_df:
            self.calc_chromatic_beating()
        return self._results_df.loc[:, ["S", "DBEATX", "DBEATY"]]

    def plot_chromatic_beating(self, combined=True):
        """ Plot the Chromatic Beating

        Available after calc_chromatic_beating

        Args:
            combined (bool): If 'True' plots x and y into the same axes.
        """
        LOG.debug("Plotting Chromatic Beating")
        chrom_beat = self.get_chromatic_beating().dropna()
        title = 'Chromatic Beating'
        pstyle.set_style(self._plot_options.style, self._plot_options.manual)

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
        """ Eq. 11, 41 - 43 in [#FranchiAnalyticformulasrapid2017]_ """
        raise NotImplementedError('Second Order Dispersion is not Implemented yet.')

    def get_second_order_dispersion(self):
        raise NotImplementedError('Second Order Dispersion is not Implemented yet.')

    def plot_second_order_dispersion(self):
        raise NotImplementedError('Second Order Dispersion is not Implemented yet.')

    ################################
    #  Chromatic Coupling (NA)
    ################################

    def calc_chromatic_coupling(self):
        """ Eq. 47, 48, B77 in [#FranchiAnalyticformulasrapid2017]_ """
        raise NotImplementedError('Chromatic Coupling is not Implemented yet.')

    def get_chromatic_coupling(self):
        raise NotImplementedError('Chromatic Coupling is not Implemented yet.')

    def plot_chromatic_coupling(self):
        raise NotImplementedError('Chromatic Coupling is not Implemented yet.')

    ################################
    #       Class Helpers
    ################################

    def _calc_chromatic_term(self):
        """ Calculates the chromatic term which is common to all chromatic equations """
        LOG.debug("Calculating Chromatic Term.")
        self._check_k_columns(["K1L", "K2L", "K2SL"])
        res = self._results_df
        tw = self.twiss_df

        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            mask = (tw['K1L'] != 0) | (tw['K2L'] != 0) | (tw['K2SL'] != 0)
            if "DX" in tw and "DY" in tw:
                LOG.debug("Dispersion values found in model. Used for chromatic calculations")
                sum_term = tw.loc[mask, 'K1L'] - \
                           (tw.loc[mask, 'K2L'] * tw.loc[mask, 'DX']) + \
                           (tw.loc[mask, 'K2SL'] * tw.loc[mask, 'DY'])
            else:
                LOG.info("Dispersion values NOT found in model. Using analytic values.")
                if "DX" not in res or "DY" not in res:
                    self.calc_linear_dispersion()
                sum_term = tw.loc[mask, 'K1L'] - \
                           (tw.loc[mask, 'K2L'] * res.loc[mask, 'DX']) + \
                           (tw.loc[mask, 'K2SL'] * res.loc[mask, 'DY'])

            res['CHROMX'] = sum_term * tw.loc[mask, 'BETX']
            res['CHROMY'] = sum_term * tw.loc[mask, 'BETY']

        LOG.debug("Chromatic Term Calculated.")

    def _check_k_columns(self, k_columns):
        """ Check if K_columns are in twiss_df and add them all zero if not."""
        for k in k_columns:
            if k not in self.twiss_df:
                LOG.debug("Added {:s} with all zero to data-frame.".format(k))
                self.twiss_df[k] = 0.

    def _nice_axes(self, ax):
        """ Makes the axes look nicer """
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 3))
        pstyle.set_xaxis_label(ax)
        try:
            pstyle.set_xLimits(self.twiss_df.SEQUENCE, ax)
        except pstyle.ArgumentError:
            pass
        if self._ip_pos is not None and len(self._ip_pos) > 0:
            pstyle.show_ir(self._ip_pos, ax)

    @staticmethod
    def _log_added(*args):
        """ Logging Helper to log which fields were added """
        if len(args) > 0:
            fields = "'" + "', '".join(args) + "'"
            LOG.debug("  Added fields to results: " + fields)


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))

