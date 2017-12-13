""" Provides Class to get response matrices from Twiss parameters.
    The calculation is based on formulas in [1,2].

    Only works properly for on-orbit twiss files.

    Beta Beating Response:  Eq. A35 inserted into Eq. B45 in [1]
    Dispersion Response*:    Eq. 26-27 in [1]
    Phase Advance Response: Eq. 28 in [1]
    Tune Response:          Eq. 7 in [2]

    References:
        [1]  A. Franchi et al.,
             Analytic formulas for the rapid evaluation of the orbit response matrix and chromatic
             functions from lattice parameters in circular accelerators
             NOT YET PUBLISHED

        [2]  R. Tomas, et al.,
             Review of linear optics measurement and correction for charged particle accelerators.
             Physical Review Accelerators and Beams, 20(5), 54801. (2017)
             https://doi.org/10.1103/PhysRevAccelBeams.20.054801

    * Because the dispersion is normalized here, the change comes not only from the change of the
    dispersion itself, but also from the change in the beta function.
    The response is linearised, see ./doc/normalized_dispersion_linearisation.pdf
"""

import os
import json
import numpy as np
from Utilities import logging_tools as logtool
from Utilities import tfs_pandas as tfs
from Utilities.contexts import timeit
from Utilities import iotools
from twiss_optics import sequence_parser

from twiss_optics.twiss_functions import get_phase_advances, tau, dphi
from twiss_optics.twiss_functions import assertion, regex_in, upper


LOG = logtool.get_logger(__name__)

EXCLUDE_CATEGORIES_DEFAULT = ["LQ", "MCBH", "MQX", "MQXT", "Q", "QIP15", "QIP2", "getListsByIR"]
DUMMY_ID = "DUMMY_PLACEHOLDER"

"""
=============================   Twiss Response Class   =============================
"""

class TwissResponse(object):
    """ Provides Response Matrices calculated from sequence, model and given variables.

    Args:
        seqfile_path: Path to sequence file. If there is a pre-parsed .varmap file in the same
            folder, it will use this one. (Hence, can also be the path to this file)
        modelfile_path: Path to twiss-model file
        varfile_path: Path to json file containing the variable-names
        exclude_categories: Names of the categories to exclude in varfile
        at_elements (str): Get response matrix for these elements. Can be:
            'bpms': All BPMS (Default)
            'bpms+': BPMS+ used magnets (== magnets defined by variables in varfile)
            'all': All BPMS and Magnets given in the model (Markers are removed)
    """

    ################################
    #            INIT
    ################################

    def __init__(self, seqfile_path, modelfile_path, varfile_path,
                 exclude_categories=EXCLUDE_CATEGORIES_DEFAULT,
                 at_elements='bpms'):
        LOG.info("Calculating TwissResponse.")
        with timeit(lambda t: LOG.debug("  Total time TwissResponse: {:f}s".format(t))):
            # Get input
            self._twiss = self._get_model_twiss(modelfile_path)
            self._variables = self._read_variables(varfile_path, exclude_categories)
            self._var_to_el = self._get_variable_mapping(seqfile_path)
            self._elements_in = self._get_input_elements()
            self._elements_out = self._get_output_elements(at_elements)

            # calculate all phase advances
            self._phase_advances = get_phase_advances(self._twiss)

            # get response matrices
            self._beta = self._calc_beta_response()
            self._dispersion = self._calc_dispersion_response()
            self._phase = self._calc_phase_advance_response()
            self._tune = self._calc_tune_response()

            # map response matrices to variabels
            self._beta_mapped = self._map_to_variables(self._beta, self._var_to_el["K1L"])
            self._dispersion_mapped = self._map_dispersion_response()
            self._phase_mapped = self._map_to_variables(self._phase, self._var_to_el["K1L"])
            self._tune_mapped = self._map_to_variables(self._tune, self._var_to_el["K1L"])


    @staticmethod
    def _get_model_twiss(modelfile_path):
        """ Load model, but keep only BPMs and Magnets """
        LOG.debug("Loading Model from file '{:s}'".format(modelfile_path))
        model = tfs.read_tfs(modelfile_path).set_index('NAME')

        # Remove not needed Stuff
        LOG.debug("Removing non-necessary entries:")
        LOG.debug("  Entries total: {:d}".format(model.shape[0]))
        model = model.loc[regex_in(r"\A(M|BPM)", model.index), :]
        LOG.debug("  Entries left: {:d}".format(model.shape[0]))

        # Add Dummy for Phase Calculations
        model.loc[DUMMY_ID, ["S", "MUX", "MUY"]] = 0.0
        return model

    @staticmethod
    def _read_variables(varfile_path, exclude_categories):
        """ Load variables list from json file """
        LOG.debug("Loading variables from file {:s}".format(varfile_path))

        with open(varfile_path, 'r') as varfile:
            var_dict = json.load(varfile)

        variables = []
        for category in var_dict.keys():
            if category not in exclude_categories:
                variables += var_dict[category]
        return list(set(variables))

    def _get_variable_mapping(self, seqfile_path):
        """ Get variable mapping as dictionary

        Define _variables first!
        """
        LOG.debug("Converting variables to magnet names.")
        variables = self._variables
        mapping = sequence_parser.load_or_parse_variable_mapping(seqfile_path, "dictionary")

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

    def _calc_beta_response(self):
        """ Response Matrix for delta betabeats.

        Eq. A35 -> Eq. B45 in [1]
        """
        LOG.debug("Calculate Beta Beating Response Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}s".format(t))):
            tw = self._twiss
            adv = self._phase_advances
            el_out = self._elements_out
            k1_el = self._elements_in["K1L"]
            dbetabeat = dict.fromkeys(["X", "Y"])

            dphix = dphi(adv["X"].loc[k1_el, el_out], tw.Q1)
            dphiy = dphi(adv["Y"].loc[k1_el, el_out], tw.Q2)

            dbetabeat["X"] = tfs.TfsDataFrame(
                -(1/(2*np.sin(2*np.pi*tw.Q1))) *
                  np.cos(2*np.pi*(2*dphix.values - tw.Q1)) *
                  tw.loc[k1_el, ["BETX"]].values, index=k1_el, columns=el_out).transpose()

            dbetabeat["Y"] = tfs.TfsDataFrame(
                (1 / (2 * np.sin(2 * np.pi * tw.Q2))) *
                np.cos(2 * np.pi * (2 * dphiy.values - tw.Q2)) *
                tw.loc[k1_el, ["BETY"]].values, index=k1_el, columns=el_out).transpose()

        return dbetabeat

    def _calc_dispersion_response(self):
        """ Response Matrix for delta normalized dispersion
            Eq. 26-27 in [1]
            Call after beta-beat calculation!
        """
        LOG.debug("Calculate Normalized Dispersion Response Matrix")
        with timeit(lambda t: LOG.debug("  Time needed: {:f}".format(t))):
            tw = self._twiss
            adv = self._phase_advances
            el_out = self._elements_out
            k0_el = self._elements_in["K0L"]
            k1_el = self._elements_in["K1L"]
            j0_el = self._elements_in["K0SL"]
            j1_el = self._elements_in["K1SL"]
            bbeat = self._beta

            disp_resp = dict.fromkeys(["X", "Y_J0", "Y_J1", "X_BB", "Y_BB"])

            # Matrix N in Eq. 26 (normalized)
            if len(k0_el) > 0:
                taux = tau(adv["X"].loc[k0_el, el_out], tw.Q1)
                n = 1/(2*np.sin(np.pi * tw.Q1)) * \
                    np.sqrt(tw.loc[k0_el, ["BETX"]].values) * np.cos(2 * np.pi * taux)
                disp_resp["X"] = n.transpose()
            else:
                LOG.debug("  No 'K0L' variables found. Dispersion Response 'X' will be empty.")
                disp_resp["X"] = tfs.TfsDataFrame(None, index=el_out)

            # Matrix S(J0) in Eq. 26 (normalized)
            if len(j0_el) > 0:
                tauy_j0 = tau(adv["Y"].loc[j0_el, el_out], tw.Q2)
                sj0 = -1 / (2 * np.sin(np.pi * tw.Q2)) * \
                  np.sqrt(tw.loc[j0_el, ["BETY"]].values) * np.cos(2 * np.pi * tauy_j0)
                disp_resp["Y_J0"] = sj0.transpose()
            else:
                LOG.debug("  No 'K0SL' variables found. Dispersion Response 'Y_J0' will be empty.")
                disp_resp["Y_J0"] = tfs.TfsDataFrame(None, index=el_out)

            # Matrix S(J1) in Eq. 26 (normalized)
            if len(j1_el) > 0:
                tauy_j1 = tau(adv["Y"].loc[j1_el, el_out], tw.Q2)
                sj1 = -1 / (2 * np.sin(np.pi * tw.Q2)) * tw.loc[j1_el, ["DX"]].values * \
                  np.sqrt(tw.loc[j1_el, ["BETY"]].values) * np.cos(2 * np.pi * tauy_j1)
                disp_resp["Y_J1"] = sj1.transpose()
            else:
                LOG.debug("  No 'K1SL' variables found. Dispersion Response 'Y_J1' will be empty.")
                disp_resp["Y_J1"] = tfs.TfsDataFrame(None, index=el_out)

            # Correction Terms rising from delta beta in normalization
            if len(k1_el) > 0:
                for plane in ["X", "Y"]:
                    disp_resp[plane+"_BB"] = tfs.TfsDataFrame(
                        -.5 * bbeat[plane].values * (tw.loc[k1_el, ["D"+plane]].values /
                                                     np.sqrt(tw.loc[k1_el, ["BET"+plane]].values)
                                                     ).transpose(),
                        index=el_out, columns=k1_el
                    )
            else:
                LOG.debug("  No 'K1L' variables found."
                          "Dispersion Responses 'X_BB' and 'Y_BB' will be empty.")
                disp_resp["X_BB"] = tfs.TfsDataFrame(None, index=el_out)
                disp_resp["Y_BB"] = tfs.TfsDataFrame(None, index=el_out)

        return disp_resp

    def _calc_phase_advance_response(self):
        """ Response Matrix for delta DPhi.

        Eq. 28 in [1]
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
                dadv = dict.fromkeys(["X", "Y"])

                pi = tfs.TfsDataFrame(tw['S'][:, None] < tw['S'][None, :],  # pi(i,j) = s(i) < s(j)
                                      index=tw.index, columns=tw.index, dtype=int)

                pi_term = (pi.loc[k1_el, el_out].values -
                           pi.loc[k1_el, el_out_mm].values +
                           np.diag(pi.loc[el_out, el_out_mm].values)[None, :])

                taux = 2 * np.pi * tau(adv['X'].loc[k1_el, el_out_all], tw.Q1)
                tauy = 2 * np.pi * tau(adv['Y'].loc[k1_el, el_out_all], tw.Q2)

                brackets_x = (2 * pi_term +
                              ((np.sin(2 * taux.loc[:, el_out].values) -
                                np.sin(2 * taux.loc[:, el_out_mm].values)) /
                               np.sin(2 * np.pi * tw.Q1)))

                # hint: np can only broadcast ndarrays, which is values of a DF but not of a Series
                dadv['X'] = tfs.TfsDataFrame(
                            brackets_x * tw.loc[k1_el, ['BETX']].values * (1 / (8 * np.pi)),
                            index=k1_el, columns=el_out)

                brackets_y = (2 * pi_term +
                              ((np.sin(2 * tauy.loc[:, el_out].values) -
                                np.sin(2 * tauy.loc[:, el_out_mm].values)) /
                               np.sin(2 * np.pi * tw.Q2)))

                # hint: np can only broadcast ndarrays, which is values of a DF but not of a Series
                dadv['Y'] = tfs.TfsDataFrame(
                            brackets_y * tw.loc[k1_el, ['BETY']].values * (-1 / (8 * np.pi)),
                            index=k1_el, columns=el_out)

                # switching to total phase and proper axis orientation
                dmu = {
                    "X": dadv['X'].transpose().apply(np.cumsum),
                    "Y": dadv['Y'].transpose().apply(np.cumsum),
                }
            else:
                LOG.debug("  No 'K1L' variables found. Phase Response will be empty.")
                dmu = {"X": tfs.TfsDataFrame(None, index=el_out),
                       "Y": tfs.TfsDataFrame(None, index=el_out)}

        return dmu

    def _calc_tune_response(self):
        """ Response vectors for Tune.

        Eq. 7 in [2]
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

    def _map_dispersion_response(self):
        """ Maps all dispersion matrices """
        var2k0 = self._var_to_el["K0L"]
        var2j0 = self._var_to_el["K0SL"]
        var2k1 = self._var_to_el["K1L"]
        var2j1 = self._var_to_el["K1SL"]
        m2v = self._map_to_variables
        disp = self._dispersion

        return {
            "X": m2v(disp["X"], var2k0),
            "X_BB": m2v(disp["X_BB"], var2k1),
            "Y_J0": m2v(disp["Y_J0"], var2j0),
            "Y_J1": m2v(disp["Y_J1"], var2j1),
            "Y_BB": m2v(disp["Y_BB"], var2k1),
        }

    @staticmethod
    def _map_to_variables(df, mapping):
        """ Maps from magnets to variables using self._var_to_el.
            Could actually be done by matrix multiplication
            A * var_to_el, yet, as var_to_el is very sparse,
            looping is easier.

            Args:
                df: DataFrame or dictionary of dataframes to map
                mapping: mapping to be applied (e.g. var_to_el[order])
            Returns:
                DataFrame or dictionary of mapped dataframes
        """
        def map_fun(df, mapping):
            """ Actual mapping """
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

    def get_beta_beat(self, mapped=True):
        """ Returns Response Matrix for Beta Beat """
        if mapped:
            return self._beta_mapped
        else:
            return self._beta

    def get_dispersion(self, mapped=True):
        """ Returns Response Matrix for Normalized Dispersion """
        if mapped:
            return self._dispersion_mapped
        else:
            return self._dispersion

    def get_phase(self, mapped=True):
        """ Returns Response Matrix for Total Phase """
        if mapped:
            return self._phase_mapped
        else:
            return self._phase

    def get_tune(self, mapped=True):
        """ Returns Response Matrix for the tunes """
        if mapped:
            return self._tune_mapped
        else:
            return self._tune

    def get_fullresponse(self):
        """ Returns all mapped Response Matrices """
        tune = self._tune_mapped
        beta = self._beta_mapped
        disp = self._dispersion_mapped
        phse = self._phase_mapped

        q_df = tune["X"].append(tune["Y"])
        q_df.index = ["Q1", "Q2"]

        return {
            "BBX": beta["X"],
            "BBY": beta["Y"],
            "MUX": phse["X"],
            "MUY": phse["Y"],
            "NDX": disp["X"].join(disp["X_BB"]),
            "NDY": disp["Y_J0"].join(disp["Y_J1"]).join(disp["Y_BB"]),
            "Q": q_df,
        }

    def get_variabel_names(self):
        return self._variables

    def get_variable_mapping(self, order=None):
        if order is None:
            return self._var_to_el
        else:
            return self._var_to_el[order]


"""
=============================   Main   =============================
"""


if __name__ == '__main__':
    root = iotools.get_absolute_path_to_betabeat_root()
    seq = os.path.join(root, 'twiss_optics', 'tests', 'lhcb1.seq')
    mod = os.path.join(root, 'twiss_optics', 'tests', 'twiss_dispersion.dat')
    var = os.path.join(root, 'MODEL', 'LHCB', 'fullresponse', 'LHCB1', 'AllLists.json')
    tr = TwissResponse(seq, mod, var)
    tr.get_fullresponse()
