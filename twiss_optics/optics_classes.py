from __future__ import print_function
import __init__
import os, sys
import numpy as np
import pandas as pd
from math import factorial
from Utilities import tfs_pandas
from Utilities.contexts import timeit
import matplotlib.pyplot as plt
from Utilities import logging_tools as logtool
import re


LOG = logtool.get_logger(__name__, os.devnull)


################################
#           CLASSES
################################

class TwissOptics(object):

    ################################
    #       init Functions
    ################################

    def __init__(self, model_path):
        LOG.debug("Creating TwissOptics from '{:s}'".format(model_path))
        self.mad_twiss = tfs_pandas.read_tfs(model_path).set_index('NAME')
        self._remove_nonnecessaries()
        self._add_result_dataframe()  #self._results_df
        self._add_element_types()  # self._elements
        self._add_phase_advances()
        # self._calc_rdt('f4000')
        # self.f4000 = self.mad_twiss['f4000']
        #
        # self.f4000.abs().plot()
        # plt.show()

    def _remove_nonnecessaries(self):
        LOG.debug("Removing non-necessaries:")
        LOG.debug("  Entries total: {:d}".format(self.mad_twiss.shape[0]))
        with timeit(lambda t:
                    LOG.debug("  Removed in {:f}s".format(t))):
            self.mad_twiss = self.mad_twiss.loc[regex_in("\A[MB]", self.mad_twiss.index)]
        LOG.debug("  Entries left: {:d}".format(self.mad_twiss.shape[0]))

    def _add_result_dataframe(self):
        LOG.debug("Adding results dataframe.")
        self._results_df = tfs_pandas.TfsDataFrame(
            index=self.mad_twiss.index
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
        """ Coupling Fucntion """
        return self._results_df.coupling

    @coupling.setter
    def coupling(self, new):
        self._results_df.coupling = new

    ################################
    #       Main Functions
    ################################

    def get_rdts(self, rdt_names=None):
        """ Return dataframe with rdts """
        if rdt_names:
            if not isinstance(rdt_names, list):
                rdt_names = [rdt_names]
            return self._results_df.loc[:, [rdt.upper() for rdt in rdt_names]]
        else:
            return self._results_df.loc[:, regex_in('\AF\d{4}', self._results_df.columns)]

    def calc_rdts(self, order):
        """ Calculate RDTs """
        if not isinstance(order, list):
            order = [order]

        LOG.debug("Calculating RDTs: {:s}.".format(str(order)[1:-1]))
        with timeit(lambda t:
                    LOG.debug("RDTs calculated in {:f}s".format(t))):

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
                        phx = self._dphi(phs_adv['X'].loc[k_mask, el_mask], 'X')
                        phy = self._dphi(phs_adv['Y'].loc[k_mask, el_mask], 'Y')
                        phase_term = ((j-k) * phx + (l-m) * phy).applymap(lambda p: np.exp(i2pi*p))

                        beta_term = tw.loc[k_mask, src] * \
                                    tw.loc[k_mask, 'BETX'] ** ((j+k) / 2.) * \
                                    tw.loc[k_mask, 'BETY'] ** ((l+m) / 2.)

                        res_df.loc[el_mask, rdt.upper()] = sign * np.sum(
                            phase_term.multiply(beta_term, axis="index")
                            , axis=0).transpose() * denom
                    else:
                        LOG.debug("All {:s} == 0. RDT '{:s}' will be zero.".format(src, rdt))
                        res_df.loc[el_mask, rdt.upper()] = 0

    ################################
    #       Class Helpers
    ################################

    def _dphi(self, data, plane):
        q = self.mad_twiss.Q1 if plane == 'X' else self.mad_twiss.Q2
        return data + np.where(data < 0, q, 0)

    def _tau(self, data, plane):
        q = self.mad_twiss.Q1 if plane == 'X' else self.mad_twiss.Q2
        return data + np.where(data < 0, q / 2, -q / 2)


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


################################
#           MAIN
################################

def do_main():
    ''' Main function launcher '''
    model_path = 'twiss_2octupoles.dat'
    two_octupoles = TwissOptics(model_path)
    two_octupoles.calc_rdt('f1000')


if __name__ == '__main__':
    do_main()
