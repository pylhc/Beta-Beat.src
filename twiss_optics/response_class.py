""" Provides Class to get response matrices from Twiss parameters.
    The calculation is based on formulas in [#fran], [#toma].

    Only works properly for on-orbit twiss files.

    Beta Beating Response:  Eq. A35 inserted into Eq. B45 in [#fran]
    Dispersion Response:    Eq. 25-27 in [#fran]
    Phase Advance Response: Eq. 28 in [#fran]
    Tune Response:          Eq. 7 in [#toma]

    For people reading the code:
    The response matrices are first calculated like:

    |  Elements of interest (j) --> ... |
    |Magnets (m)                        |
    |  |                                |
    |  v                                |
    |  .                                |
    |  .                                |
    |  .                                |
    |                                   |

    As this was avoided transposing all vectors in the beginning.
    At the end (of the calculation) the matrix is then transposed
    to fit the :math:'M \cdot \delta K' orientation.

    References:
        [#fran]  A. Franchi et al.,
                 'Analytic formulas for the rapid evaluation of the orbit response matrix
                 and chromatic functions from lattice parameters in circular accelerators'
                 NOT YET PUBLISHED

        [#toma]  R. Tomas, et al.,
                 'Review of linear optics measurement and correction for charged particle
                 accelerators.'
                 Physical Review Accelerators and Beams, 20(5), 54801. (2017)
                 https://doi.org/10.1103/PhysRevAccelBeams.20.054801
"""

import json

import numpy as np
import pandas as pd

from utils import logging_tools as logtool
from utils import tfs_pandas as tfs
from utils.contexts import timeit
from twiss_optics import sequence_parser
from twiss_optics.twiss_functions import get_phase_advances, tau, dphi
from twiss_optics.twiss_functions import regex_in, upper

LOG = logtool.get_logger(__name__)

DUMMY_ID = "DUMMY_PLACEHOLDER"


# Twiss Response Class ########################################################


class TwissResponse(object):
    """ Provides Response Matrices calculated from sequence, model and given variables.

    Args:
        varmap_or_seq_path: Path to sequence file. If there is a pre-parsed .varmap file in the same
            folder, it will use this one. (Hence, can also be the path to this file)
        model_or_path: Path to twiss-model file, or model
        variables: List of variable-names
        at_elements (str): Get response matrix for these elements. Can be:
            'bpms': All BPMS (Default)
            'bpms+': BPMS+ used magnets (== magnets defined by variables in varfile)
            'all': All BPMS and Magnets given in the model (Markers are removed)
    """

    ################################
    #            INIT
    ################################

    def __init__(self, varmap_or_seq_path, model_or_path, variables,
                 at_elements='bpms'):

        LOG.info("Calculating TwissResponse.")
        with timeit(lambda t: LOG.debug("  Total time TwissResponse: {:f}s".format(t))):
            # Get input
            self._twiss = self._get_model_twiss(model_or_path)
            self._var_to_el = self._get_variable_mapping(varmap_or_seq_path)
            self._elements_in = self._get_input_elements()
            self._elements_out = self._get_output_elements(at_elements)
            self._variables = variables

            # calculate all phase advances
            self._phase_advances = get_phase_advances(self._twiss)

            # All responses are calcluated as needed, see getters below!

            # slots for response matrices
            self._beta = None
            self._dispersion = None
            self._phase = None
            self._tune = None
            self._coupling = None

            # slots for normalized response matrices
            self._beta_norm = None  # i.e. beta-beating
            self._dispersion_norm = None

            # slots for mapped response matrices
            self._coupling_mapped = None
            self._beta_mapped = None
            self._dispersion_mapped = None
            self._phase_mapped = None
            self._tune_mapped = None

            # slots for normalized mapped response matrices
            self._beta_mapped_norm = None  # i.e. beta-beating
            self._dispersion_mapped_norm = None

    @staticmethod
    def _get_model_twiss(model_or_path):
        """ Load model, but keep only BPMs and Magnets """
        try:
            model = tfs.read_tfs(model_or_path, index="NAME")
        except TypeError:
            LOG.debug("Received model as DataFrame")
            model = model_or_path
        else:
            LOG.debug("Loaded Model from file '{:s}'".format(model_or_path))

        # Remove not needed Stuff
        LOG.debug("Removing non-necessary entries:")
        LOG.debug("  Entries total: {:d}".format(model.shape[0]))
        model = model.loc[regex_in(r"\A(M|BPM)", model.index), :]
        LOG.debug("  Entries left: {:d}".format(model.shape[0]))

        # Add Dummy for Phase Calculations
        model.loc[DUMMY_ID, ["S", "MUX", "MUY"]] = 0.0
        return model

    def _get_variable_mapping(self, varmap_or_seq_path):
        """ Get variable mapping as dictionary

        Define _variables first!
        """
        LOG.debug("Converting variables to magnet names.")
        variables = self._variables
        mapping = sequence_parser.load_or_parse_variable_mapping(varmap_or_seq_path, "dictionary")

        for order in ("K0L", "K0SL", "K1L", "K1SL"):
            if order not in mapping:
                mapping[order] = {}

        # check if all variables can be found
        check_var = [var for var in variables
                     if all(var not in mapping[order] for order in mapping)]
        if len(check_var) > 0:
            raise ValueError("Variables '{:s}' cannot be found in sequence!".format(
                ", ".join(check_var)
            ))

        # drop mapping for unused variables
        [mapping[order].pop(var) for order in mapping for var in mapping[order].keys()
         if var not in variables]

        return mapping

    def _get_input_elements(self):
        """ Return variable names of input elements.

        Define _var_to_el and _twiss first!
        """
        v2e = self._var_to_el
        tw = self._twiss

        el_in = dict.fromkeys(v2e.keys())
        for order in el_in:
            el_order = []
            for var in v2e[order]:
                el_order += upper(v2e[order][var].index)
            el_in[order] = tw.loc[list(set(el_order)), "S"].sort_values().index.tolist()
        return el_in

    def _get_output_elements(self, at_elements):
        """ Return name-array of elements to use for output.

        Define _elements_in first!
        """
        tw_idx = self._twiss.index

        if isinstance(at_elements, list):
            # elements specified
            if any(el not in tw_idx for el in at_elements):
                LOG.warning("One or more specified elements are not in the model.")
            return [idx for idx in tw_idx
                    if idx in at_elements]

        if at_elements == "bpms":
            # bpms only
            return [idx for idx in tw_idx
                    if idx.upper().startswith('B')]

        if at_elements == "bpms+":
            # bpms and the used magnets
            el_in = self._elements_in
            return [idx for idx in tw_idx
                    if (idx.upper().startswith('B')
                        or any(idx in el_in[order] for order in el_in))]

        if at_elements == "all":
            # all, obviously
            return [idx for idx in tw_idx if idx != DUMMY_ID]

    ################################
    #       Response Matrix
    ################################

    def _calc_coupling_response(self):
        """ Response Matrix for coupling.

        Eq. 10 in [#fran]
        """
        LOG.debug("Calculate Coupling Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}s".format(t))):
            tw = self._twiss
            adv = self._phase_advances
            el_out = self._elements_out
            k1s_el = self._elements_in["K1SL"]
            dcoupl = dict.fromkeys(["1001", "1010"])

            i2pi = 2j * np.pi
            phx = dphi(adv['X'].loc[k1s_el, el_out], tw.Q1).values
            phy = dphi(adv['Y'].loc[k1s_el, el_out], tw.Q2).values
            bet_term = np.sqrt(tw.loc[k1s_el, "BETX"].values * tw.loc[k1s_el, "BETY"].values)

            for plane in ["1001", "1010"]:
                phs_sign = -1 if plane == "1001" else 1
                dcoupl[plane] = tfs.TfsDataFrame(
                    bet_term[:, None] * np.exp(i2pi * (phx + phs_sign * phy)) /
                    (4 * (1 - np.exp(i2pi * (tw.Q1 + phs_sign * tw.Q2)))),
                    index=k1s_el, columns=el_out).transpose()
        return dcoupl

    def _calc_beta_response(self):
        """ Response Matrix for delta beta.

        Eq. A35 -> Eq. B45 in [#fran]
        """
        LOG.debug("Calculate Beta Response Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}s".format(t))):
            tw = self._twiss
            adv = self._phase_advances
            el_out = self._elements_out
            k1_el = self._elements_in["K1L"]
            dbeta = dict.fromkeys(["X", "Y"])

            for plane in ["X", "Y"]:
                col_beta = "BET" + plane
                q = tw.Q1 if plane == "X" else tw.Q2
                coeff_sign = -1 if plane == "X" else 1

                pi2tau = 2 * np.pi * tau(adv[plane].loc[k1_el, el_out], q)

                dbeta[plane] = tfs.TfsDataFrame(
                    tw.loc[el_out, col_beta].values[None, :] *
                    tw.loc[k1_el, col_beta].values[:, None] * np.cos(2 * pi2tau.values) *
                    (coeff_sign / (2 * np.sin(2 * np.pi * q))),
                    index=k1_el, columns=el_out).transpose()

        return dbeta

    def _calc_dispersion_response(self):
        """ Response Matrix for delta dispersion

            Eq. 25-27 in [#fran]
        """
        LOG.debug("Calculate Dispersion Response Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            tw = self._twiss
            adv = self._phase_advances
            el_out = self._elements_out
            els_in = self._elements_in

            disp_resp = dict.fromkeys(["X_K0L", "X_K1SL", "Y_K0SL", "Y_K1SL"])

            for plane in ["X", "Y"]:
                q = tw.Q1 if plane == "X" else tw.Q2
                type_plane = ("K0L" if plane == "X" else "K0SL", "K1SL")
                el_in_plane = (els_in[type_plane[0]], els_in[type_plane[1]])
                col_beta = "BET" + plane
                col_disp = "DY" if plane == "X" else "DX"

                if any((len(el_in_plane[0]), len(el_in_plane[1]))):
                    coeff = np.sqrt(tw.loc[el_out, col_beta].values) / (2 * np.sin(np.pi * q))

                for el_in, el_type in zip(el_in_plane, type_plane):
                    coeff_sign = -1 if el_type == "K0SL" else 1
                    out_str = "{p:s}_{t:s}".format(p=plane, t=el_type)

                    if len(el_in):
                        pi2tau = 2 * np.pi * tau(adv[plane].loc[el_in, el_out], q)
                        bet_term = np.sqrt(tw.loc[el_in, col_beta].values)
                        if el_type == "K1SL":
                            bet_term *= tw.loc[el_in, col_disp].values
                        m = coeff_sign * coeff[None, :] * bet_term[:, None] * np.cos(pi2tau)
                        disp_resp[out_str] = m.transpose()
                    else:
                        LOG.debug(
                            "  No '{:s}' variables found. ".format(el_type) +
                            "Dispersion Response '{:s}' will be empty.".format(out_str))
                        disp_resp[out_str] = tfs.TfsDataFrame(None, index=el_out)
        return disp_resp

    def _calc_phase_advance_response(self):
        """ Response Matrix for delta DPhi.

        Eq. 28 in [#fran]
        Reduced to only phase advances between consecutive elements,
        as the 3D-Matrix of all elements exceeds memory space
        (~11000^3 = 1331 Giga Elements)
        --> w = j-1:  DPhi(z,j) = DPhi(x, (j-1)->j)
        """
        LOG.debug("Calculate Phase Advance Response Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}s".format(t))):
            tw = self._twiss
            adv = self._phase_advances
            k1_el = self._elements_in["K1L"]

            el_out_all = [DUMMY_ID] + self._elements_out  # Add MU[XY] = 0.0 to the start
            el_out = el_out_all[1:]  # in these we are actually interested
            el_out_mm = el_out_all[0:-1]  # elements--

            if len(k1_el) > 0:
                dmu = dict.fromkeys(["X", "Y"])

                pi = tfs.TfsDataFrame(tw['S'][:, None] < tw['S'][None, :],  # pi(i,j) = s(i) < s(j)
                                      index=tw.index, columns=tw.index, dtype=int)

                pi_term = (pi.loc[k1_el, el_out].values -
                           pi.loc[k1_el, el_out_mm].values +
                           np.diag(pi.loc[el_out, el_out_mm].values)[None, :])

                for plane in ["X", "Y"]:
                    col_beta = "BET" + plane
                    q = tw.Q1 if plane == "X" else tw.Q2
                    coeff_sign = 1 if plane == "X" else -1

                    pi2tau = 2 * np.pi * tau(adv[plane].loc[k1_el, el_out_all], q)
                    brackets = (2 * pi_term +
                                ((np.sin(2 * pi2tau.loc[:, el_out].values) -
                                  np.sin(2 * pi2tau.loc[:, el_out_mm].values))
                                 / np.sin(2 * np.pi * q)
                                 ))
                    dmu[plane] = tfs.TfsDataFrame(
                        tw.loc[k1_el, col_beta].values[:, None] * brackets
                        * (coeff_sign / (8 * np.pi)),
                        index=k1_el, columns=el_out).transpose().apply(np.cumsum)
            else:
                LOG.debug("  No 'K1L' variables found. Phase Response will be empty.")
                dmu = {"X": tfs.TfsDataFrame(None, index=el_out),
                       "Y": tfs.TfsDataFrame(None, index=el_out)}

        return dmu

    def _calc_tune_response(self):
        """ Response vectors for Tune.

        Eq. 7 in [#toma]
        """
        LOG.debug("Calculate Tune Response Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}s".format(t))):
            tw = self._twiss
            k1_el = self._elements_in["K1L"]

            if len(k1_el) > 0:
                dtune = dict.fromkeys(["X", "Y"])

                dtune["X"] = 1/(4 * np.pi) * tw.loc[k1_el, ["BETX"]].transpose()
                dtune["X"].index = ["DQX"]

                dtune["Y"] = -1 / (4 * np.pi) * tw.loc[k1_el, ["BETY"]].transpose()
                dtune["Y"].index = ["DQY"]
            else:
                LOG.debug("  No 'K1L' variables found. Tune Response will be empty.")
                dtune = {"X": tfs.TfsDataFrame(None, index="DQX"),
                         "Y": tfs.TfsDataFrame(None, index="DQY")}

        return dtune

    ################################
    #       Normalizing
    ################################

    def _normalize_beta(self):
        """ Convert to Beta Beating """
        el_out = self._elements_out
        tw = self._twiss
        beta = self._beta

        self._beta_norm = dict.fromkeys(beta.keys())
        for plane in beta:
            col = "BET" + plane
            self._beta_norm[plane] = beta[plane].div(
                tw.loc[el_out, col], axis='index')

    def _normalize_dispersion(self):
        """ Convert to Normalized Dispersion """
        el_out = self._elements_out
        tw = self._twiss
        disp = self._dispersion

        self._dispersion_norm = dict.fromkeys(disp.keys())
        for plane in disp:
            col = "BET" + plane[0]
            self._dispersion_norm[plane] = disp[plane].div(
                np.sqrt(tw.loc[el_out, col]), axis='index')

    ################################
    #       Mapping
    ################################

    def _map_dispersion_response(self, normalized=False):
        """ Maps all dispersion matrices """
        var2k0 = self._var_to_el["K0L"]
        var2j0 = self._var_to_el["K0SL"]
        var2j1 = self._var_to_el["K1SL"]
        m2v = self._map_to_variables
        disp = self._dispersion_norm if normalized else self._dispersion

        return {
            "X_K0L": m2v(disp["X_K0L"], var2k0),
            "X_K1SL": m2v(disp["X_K1SL"], var2j1),
            "Y_K0SL": m2v(disp["Y_K0SL"], var2j0),
            "Y_K1SL": m2v(disp["Y_K1SL"], var2j1),
        }

    @staticmethod
    def _map_to_variables(df, mapping):
        """ Maps from magnets to variables using self._var_to_el.
            Could actually be done by matrix multiplication :math:'A \cdot var_to_el',
             yet, as var_to_el is very sparsely populated, looping is easier.

            Args:
                df: DataFrame or dictionary of DataFrames to map
                mapping: mapping to be applied (e.g. var_to_el[order])
            Returns:
                DataFrame or dictionary of mapped DataFrames
        """
        def map_fun(df, mapping):
            """ Actual mapping function """
            df_map = tfs.TfsDataFrame(index=df.index, columns=mapping.keys())
            for var, magnets in mapping.iteritems():
                df_map[var] = df.loc[:, upper(magnets.index)].mul(
                    magnets.values, axis="columns"
                ).sum(axis="columns")
            return df_map

        # convenience wrapper for dicts
        if isinstance(df, dict):
            mapped = dict.fromkeys(df.keys())
            for plane in mapped:
                mapped[plane] = map_fun(df[plane], mapping)
        else:
            mapped = map_fun(df, mapping)
        return mapped

    ################################
    #          Getters
    ################################

    def get_beta(self, mapped=True, normalized=True):
        """ Returns Response Matrix for Beta Beating """
        if not self._beta:
            self._beta = self._calc_beta_response()

        if mapped and not self._beta_mapped:
            self._beta_mapped = self._map_to_variables(self._beta,
                                                       self._var_to_el["K1L"])

        if normalized and not self._beta_norm:
            self._normalize_beta()
            if mapped and not self._beta_mapped_norm:
                self._beta_mapped_norm = self._map_to_variables(self._beta_norm,
                                                                self._var_to_el["K1L"])

        if mapped:
            return self._beta_mapped_norm if normalized else self._beta_mapped
        else:
            return self._beta_norm if normalized else self._beta

    def get_dispersion(self, mapped=True, normalized=True):
        """ Returns Response Matrix for Normalized Dispersion """
        # prepare what is needed:
        if not self._dispersion:
            self._dispersion = self._calc_dispersion_response()

        if mapped and not self._dispersion_mapped:
            self._dispersion_mapped = self._map_dispersion_response(normalized=False)

        if normalized and not self._dispersion_norm:
            self._dispersion_norm = self._normalize_dispersion_response()
            if mapped and not self._dispersion_mapped_norm:
                self._dispersion_mapped_norm = self._map_dispersion_response(normalized=True)

        # return what is needed
        if mapped:
            return self._dispersion_mapped_norm if normalized else self._dispersion_mapped
        else:
            return self._dispersion_norm if normalized else self._dispersion

    def get_phase(self, mapped=True):
        """ Returns Response Matrix for Total Phase """
        if not self._phase:
            self._phase = self._calc_phase_advance_response()
            self._phase_mapped = self._map_to_variables(self._phase, self._var_to_el["K1L"])

        if mapped:
            return self._phase_mapped
        else:
            return self._phase

    def get_tune(self, mapped=True):
        """ Returns Response Matrix for the Tunes """
        if not self._tune:
            self._tune = self._calc_tune_response()
            self._tune_mapped = self._map_to_variables(self._tune, self._var_to_el["K1L"])

        if mapped:
            return self._tune_mapped
        else:
            return self._tune

    def get_coupling(self, mapped=True):
        """ Returns Response Matrix for the Tunes """
        if not self._coupling:
            self._coupling = self._calc_coupling_response()
            self._coupling_mapped = self._map_to_variables(self._coupling, self._var_to_el["K1SL"])

        if mapped:
            return self._coupling_mapped
        else:
            return self._coupling

    def get_fullresponse(self, mapped=True):
        """ Returns all Response Matrices in a similar way as ``response_madx.py`` """
        # get all optical parametets
        tune = self.get_tune(mapped=mapped)
        beta = self.get_beta(mapped=mapped, normalized=False)
        disp = self.get_dispersion(mapped=mapped, normalized=False)
        phase = self.get_phase(mapped=mapped)
        couple = self.get_coupling(mapped=mapped)
        bbeat = self.get_beta(mapped=mapped, normalized=True)
        ndisp = self.get_dispersion(mapped=mapped, normalized=True)

        # merge tune to one
        q_df = tune["X"].append(tune["Y"])
        q_df.index = ["Q1", "Q2"]

        # return a big dictionary of things
        return {
            "BETX": beta["X"],
            "BETY": beta["Y"],
            "BBX": bbeat["X"],
            "BBY": bbeat["X"],
            "MUX": phase["Y"],
            "MUY": phase["Y"],
            "DX": response_add(disp["X_K0L"], disp["X_K1SL"]),
            "DY": response_add(disp["Y_K0SL"], disp["Y_K1SL"]),
            "NDX": response_add(ndisp["X_K0L"], ndisp["X_K1SL"]),
            "NDY": response_add(ndisp["Y_K0SL"], ndisp["Y_K1SL"]),
            "1001R": couple["1001"].apply(np.real),
            "1001I": couple["1001"].apply(np.imag),
            "1010R": couple["1010"].apply(np.real),
            "1010I": couple["1010"].apply(np.imag),
            "Q": q_df,
        }

    def get_variabel_names(self):
        return self._variables

    def get_variable_mapping(self, order=None):
        if order is None:
            return self._var_to_el
        else:
            return self._var_to_el[order]


# Associated Functions #########################################################


def get_delta(fullresp_or_tr, delta_k):
    """ Returns the deltas of :math:'response_matrix \cdot delta_k'.

    Args:
        fullresp_or_tr: Either the fullresponse dictionary or the TwissResponse Object
        delta_k: Pandas Series of variables and their delta-value

    Returns:
        TFS_DataFrame with elements as indices and the calculated deltas in the columns
    """
    if isinstance(fullresp_or_tr, TwissResponse):
        response = fullresp_or_tr.get_fullresponse()
    else:
        response = fullresp_or_tr

    columns = response.keys()
    index = response["BETX"].index

    delta_df = tfs.TfsDataFrame(None, index=index)
    for col in columns:
        # equivalent to .dot() but more efficient as delta_k is "sparse"
        if col == "Q":
            try:
                delta_q = (response[col].loc[:, delta_k.index] * delta_k
                           ).dropna(axis="columns").sum(axis="columns")
            except KeyError:
                # none of the delta_k are in DataFrame
                delta_q = pd.Series([0., 0.], index=["Q1", "Q2"])
            delta_df.headers["QX"] = delta_q["Q1"]
            delta_df.headers["QY"] = delta_q["Q2"]
        else:
            try:
                delta_df.loc[:, col] = (response[col].loc[:, delta_k.index] * delta_k
                                        ).dropna(axis="columns").sum(axis="columns")
            except KeyError:
                # none of the delta_k are in DataFrame
                delta_df.loc[:, col] = 0.
    return delta_df


def response_add(*args):
    """ Merges two or more Response Matrix DataFrames """
    base_df = args[0]
    for df in args[1:]:
        base_df = base_df.add(df, fill_value=0.)
    return base_df


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
