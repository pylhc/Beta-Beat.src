import numpy as np
import sbs_math
import pandas as pd


PLANES = ("x", "y")


def get_all_propagables():
    return (Phase, BetaPhase, AlfaPhase)


def _buffered(function):
    def wrapper(self, plane):
        try:
            result = self._buffer[plane][function.__name__]
        except KeyError:
            result = function(self, plane)
            self._buffer[plane][function.__name__] = result
        return self._buffer[plane][function.__name__]
    return wrapper


class Propagable(object):

    def __init__(self, segment, meas):
        self._segment = segment
        self._meas = meas
        self._buffer = {"x": {}, "y": {}}
        self._segment_models = None

    @property
    def segment_models(self):
        if self._segment_models is None:
            raise ValueError("self.segment_models have not been set.")
        return self._segment_models

    @segment_models.setter
    def segment_models(self, segment_models):
        self._segment_models = segment_models

    def init_conds_dict(self):
        raise NotImplementedError("Abstract method")

    @_buffered
    def meas_front(self, plane):
        """Interpolation or measured deviations to front propagated model.
        """
        raise NotImplementedError("Abstract method")

    @_buffered
    def meas_back(self, plane):
        """Interpolation or measured deviations to back propagated model.
        """
        raise NotImplementedError("Abstract method")

    @_buffered
    def corr_front(self, plane):
        """Interpolation or corrected deviations to front propagated model.
        """
        raise NotImplementedError("Abstract method")

    @_buffered
    def corr_back(self, plane):
        """Interpolation or corrected deviations to back propagated model.
        """
        raise NotImplementedError("Abstract method")

    def write_to_file(self, output_dir):
        raise NotImplementedError("Abstract method")


class Phase(Propagable):

    _cols = ("NAME S MEASPHASE{0} STDERRPHASE{0} PROPPHASE{0} "
             "ERRPROPPHASE{0} CORPHASE{0} ERRCORPHASE{0} BACKPHASE{0} "
             "ERRBACKPHASE{0} BACKCORPHASE{0} ERRBACKCORPHASE{0}")

    def __init__(self, segment, meas):
        super(Phase, self).__init__(segment, meas)
        self.bet0, self.alf0, self.errbet0, self.erralf0 = {}, {}, {}, {}
        for plane in PLANES:
            self.bet0[plane], self.errbet0[plane] =\
                BetaPhase.get_at(self._segment.start, meas, plane)
            self.alf0[plane], self.erralf0[plane] =\
                AlfaPhase.get_at(self._segment.start, meas, plane)

    def init_conds_dict(self):
        # The phase is not necessary for the initial conditions.
        return {}

    @staticmethod
    def get_at(names, meas, plane):
        uplane = plane.upper()
        phase = meas.phasetot[plane].loc[names, "PHASE{}".format(uplane)]
        error = meas.phasetot[plane].loc[names, "STDPH{}".format(uplane)]
        return phase, error

    @_buffered
    def meas_front(self, plane):
        return self._comp_meas(plane, self._segment_models.front, 1)

    @_buffered
    def corr_front(self, plane):
        return self._comp_corr(plane,
                               self.segment_models.front,
                               self.segment_models.front_corrected)

    @_buffered
    def meas_back(self, plane):
        return self._comp_meas(plane, self._segment_models.back, -1)

    @_buffered
    def corr_back(self, plane):
        return self._comp_corr(plane,
                               self.segment_models.back,
                               self.segment_models.back_corrected)

    def write_to_file(self, seg_beats):
        for plane in ("x", "y"):
            uplane = plane.upper()
            names = _common_indices(self.segment_models.front.index,
                                    self._meas.phasetot[plane].index)
            dfm = pd.DataFrame(index=names,
                               columns=Phase._cols.format(uplane).split(" "))
            dfm.NAME = names
            dfm.S = self.segment_models.front.loc[names, "S"]
            meas_ph, err_meas_ph = Phase.get_at(names, self._meas, plane)
            dfm.loc[:, "MEASPHASE{}".format(uplane)] = meas_ph
            dfm.loc[:, "STDERRPHASE{}".format(uplane)] = err_meas_ph
            phs, err_phs = self.meas_front(plane)
            dfm.loc[:, "PROPPHASE{}".format(uplane)] = phs
            dfm.loc[:, "ERRPROPPHASE{}".format(uplane)] = err_phs
            phs, err_phs = self.corr_front(plane)
            dfm.loc[:, "CORPHASE{}".format(uplane)] = phs
            dfm.loc[:, "ERRCORPHASE{}".format(uplane)] = err_phs
            phs, err_phs = self.meas_back(plane)
            dfm.loc[:, "BACKPHASE{}".format(uplane)] = phs
            dfm.loc[:, "ERRBACKPHASE{}".format(uplane)] = err_phs
            phs, err_phs = self.corr_back(plane)
            dfm.loc[:, "BACKCORPHASE{}".format(uplane)] = phs
            dfm.loc[:, "ERRBACKCORPHASE{}".format(uplane)] = err_phs
            import matplotlib.pyplot as plt
            dfm.loc[:, "PROPPHASE{}".format(uplane)].plot()
            plt.show()
            seg_beats.phase[plane] = dfm

    def _comp_meas(self, plane, seg_model, sign):
        uplane = plane.upper()
        model_ph = seg_model.loc[:, "MU{}".format(uplane)]
        if not self._segment.element:
            meas_ph, meas_err = Phase.get_at(slice(None, None, None),
                                             self._meas, plane)
            names = _common_indices(seg_model.index, meas_ph.index)
            seg_meas_ph = sign * (meas_ph.loc[names] - meas_ph.loc[names[0]]) % 1.
            ph_beating = (seg_meas_ph[names] - model_ph[names]) % 1.
            ph_beating[ph_beating > 0.5] = ph_beating[ph_beating > 0.5] - 1
            bet0, errbet0 = self.bet0[plane], self.errbet0[plane]
            alf0, erralf0 = self.alf0[plane], self.erralf0[plane]
            prop_err = sbs_math.propagate_error_phase(errbet0, erralf0,
                                                      model_ph[names],
                                                      bet0, alf0)
            err_ph = _quadratic_add(meas_err, prop_err)
            return ph_beating, err_ph
        else:
            prop_ph = model_ph[model_ph.index[0]]
            bet0, errbet0 = self.bet0[plane], self.errbet0[plane]
            alf0, erralf0 = self.alf0[plane], self.erralf0[plane]
            prop_err = sbs_math.propagate_error_phase(errbet0, erralf0,
                                                      model_ph[model_ph.index[0]],
                                                      bet0, alf0)
            return prop_ph, prop_err

    def _comp_corr(self, plane, seg_model, seg_model_corr):
        uplane = plane.upper()
        model_ph = seg_model.loc[:, "MU{}".format(uplane)]
        corr_ph = seg_model_corr.loc[:, "MU{}".format(uplane)]
        if not self._segment.element:
            ph_beating = (corr_ph - model_ph) % 1.
            bet0, errbet0 = self.bet0[plane], self.errbet0[plane]
            alf0, erralf0 = self.alf0[plane], self.erralf0[plane]
            prop_err = sbs_math.propagate_error_phase(errbet0, erralf0,
                                                      model_ph,
                                                      bet0, alf0)
            return ph_beating, prop_err
        else:
            prop_ph = model_ph.iloc[0]
            bet0, errbet0 = self.bet0[plane], self.errbet0[plane]
            alf0, erralf0 = self.alf0[plane], self.erralf0[plane]
            prop_err = sbs_math.propagate_error_phase(errbet0, erralf0,
                                                      model_ph[model_ph.index[0]],
                                                      bet0, alf0)
            return prop_ph, prop_err


class BetaPhase(Propagable):

    _init_pattern = "bet{}_{}"

    def __init__(self, segment, meas):
        super(BetaPhase, self).__init__(segment, meas)

    def init_conds_dict(self):
        init_dict = {}
        for plane in PLANES:
            ini_cond, _ = BetaPhase.get_at(self._segment.start, self._meas, plane)
            ini_name = self._init_pattern.format(plane, "ini")
            init_dict[ini_name] = ini_cond
            end_cond, _ = BetaPhase.get_at(self._segment.end, self._meas, plane)
            end_name = self._init_pattern.format(plane, "end")
            init_dict[end_name] = end_cond
        return init_dict

    @staticmethod
    def get_at(names, meas, plane):
        uplane = plane.upper()
        beta = meas.beta[plane].loc[names, "BET{}".format(uplane)]
        error = meas.beta[plane].loc[names, "ERRBET{}".format(uplane)]
        return beta, error

    @_buffered
    def meas_front(self, plane):
        return self._comp_meas(plane, self.segment_models._front)

    @_buffered
    def corr_front(self, plane):
        pass

    @_buffered
    def meas_back(self, plane):
        pass

    @_buffered
    def corr_back(self, plane):
        pass

    def _comp_meas(self, plane, seg_model):
        uplane = plane.upper()
        model_beta = seg_model.loc[:, "BET{}".format(uplane)]
        model_ph = seg_model.loc[:, "MU{}".format(uplane)]
        if not self._segment.element:
            bphase, errbphase = BetaPhase.get_at(slice(None, None, None),
                                                 self._meas, plane)
            names = _common_indices(seg_model.index, bphase.index)
            model_beta = model_beta.loc[names]
            errbphase = errbphase / model_beta
            ph_beating = (bphase - model_beta) / model_beta
            bet0, errbet0 = self.bet0[plane], self.errbet0[plane]
            alf0, erralf0 = self.alf0[plane], self.erralf0[plane]
            prop_err = sbs_math.propagate_error_beta(errbet0, erralf0,
                                                     model_ph[names],
                                                     model_beta,
                                                     bet0, alf0)
            err_ph = _quadratic_add(errbphase, prop_err)
            return ph_beating, err_ph
        else:
            prop_beta = model_beta[model_beta.index[0]]
            bet0, errbet0 = self.bet0[plane], self.errbet0[plane]
            alf0, erralf0 = self.alf0[plane], self.erralf0[plane]
            prop_err = sbs_math.propagate_error_phase(errbet0, erralf0,
                                                      model_beta[model_beta.index[0]],
                                                      bet0, alf0)
            return prop_beta, prop_err


class AlfaPhase(Propagable):

    _init_pattern = "alf{}_{}"

    def __init__(self, segment, meas):
        super(AlfaPhase, self).__init__(segment, meas)

    def init_conds_dict(self):
        init_dict = {}
        for plane in PLANES:
            ini_cond, _ = AlfaPhase.get_at(self._segment.start, self._meas, plane)
            ini_name = self._init_pattern.format(plane, "ini")
            init_dict[ini_name] = ini_cond
            end_cond, _ = AlfaPhase.get_at(self._segment.end, self._meas, plane)
            end_name = self._init_pattern.format(plane, "end")
            init_dict[end_name] = end_cond
        return init_dict

    @staticmethod
    def get_at(names, meas, plane):
        uplane = plane.upper()
        beta = meas.beta[plane].loc[names, "ALF{}".format(uplane)]
        error = meas.beta[plane].loc[names, "ERRALF{}".format(uplane)]
        return beta, error

    @_buffered
    def meas_front(self, plane):
        pass

    @_buffered
    def corr_front(self, plane):
        pass

    @_buffered
    def meas_back(self, plane):
        pass

    @_buffered
    def corr_back(self, plane):
        pass


def _common_indices(*indices):
    """ Common indices with indicies[0] order
    """
    common = indices[0]
    for index in indices[1:]:
        common = common.intersection(index)
    return common


def _quadratic_add(*values):
    result = 0.
    for value in values:
        result += value ** 2
    return np.sqrt(result)
