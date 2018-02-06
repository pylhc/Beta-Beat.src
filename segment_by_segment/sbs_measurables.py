import os
import sys
from collections import namedtuple
from contextlib import contextmanager
import numpy as np
import sbs_math

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

from Utilities import tfs_pandas, logging_tools

LOGGER = logging_tools.get_logger(__name__)

PLANES = ("x", "y")


MeasValue = namedtuple("MeasValue", ("value", "error"))
MeasCoupling = namedtuple("MeasValue", ("abs", "error", "real", "imag"))
MeasRMatrix = namedtuple("MeasValue", ("r11", "r12", "r21", "r22"))


class Measurable(object):

    _init_pattern = None

    def __init__(self, segment, meas):
        self._segment = segment
        self._meas = meas
        self._is_available = False

    @staticmethod
    def get_at(name, meas):
        """
        Returns the value of this Measurable in the meas Measurement at name.
        """
        raise NotImplementedError("Abstract method")

    @property
    def is_available(self):
        return self._is_available

    def get_init_conds_dict(self):
        raise NotImplementedError("Abstract method")


# Obligatory measurables ###################################################
# This measurables must be available for segment_by_segment to work and will
# raise IOErrors otherwise.
############################################################################

class Phase(Measurable):

    cols = ("NAME S MEASPHASE{0} STDERRPHASE{0} PROPPHASE{0} "
            "ERRPROPPHASE{0} CORPHASE{0} ERRCORPHASE{0} BACKPHASE{0} "
            "ERRBACKPHASE{0} BACKCORPHASE{0} ERRBACKCORPHASE{0} "
            "MODEL_S")

    def __init__(self, segment, meas):
        super(Phase, self).__init__(segment, meas)
        self.alf_ini = Alfa.get_at(segment.start, meas)
        self.alf_end = Alfa.get_at(segment.end, meas)
        self.bet_ini = BetaPhase.get_at(segment.start, meas)
        self.bet_end = BetaPhase.get_at(segment.end, meas)
        self._is_available = True

    def get_init_conds_dict(self):
        # The phase is not used for initial conditions
        return {}

    def get_beatings_at(self, names, ref_model, backwards):
        results = {}
        for plane in PLANES:
            meas_ph = self._meas.phasetot[plane]
            names = _common_indices(ref_model.index, meas_ph.index)
            name0 = self._segment.start if not backwards else self._segment.end
            seg_meas_ph = (meas_ph.loc[names, "PHASE{}".format(plane)] -
                           meas_ph.loc[name0, "PHASE{}".format(plane)]) % 1.
            ph_beating = (seg_meas_ph -
                          ref_model.loc[names, "MU{}".format(plane)]) % 1.
            err_ph = self._get_error_at(names, ref_model, plane, backwards, meas_ph)
            results[plane] = (ph_beating, err_ph)
        return results

    def _get_error_at(self, names, ref_model, plane, backwards, meas_ph=None):
        if not backwards:
            bet0, errbet0 = self.bet_ini[plane]
            alf0, erralf0 = self.alf_ini[plane]
        else:
            bet0, errbet0 = self.bet_end[plane]
            alf0, erralf0 = self.alf_end[plane]
        dphi = ref_model.loc[names, "MU{}".format(plane)]
        prop_err_ph = sbs_math.propagate_error_phase(
            errbet0, erralf0, dphi, bet0, alf0
        )
        if meas_ph is None:
            return prop_err_ph
        stddphi = meas_ph.loc[names, "STDPH{}".format(plane)]
        return _quadratic_add(stddphi, prop_err_ph)


class BetaPhase(Measurable):

    _init_pattern = "bet{}_{}"

    def __init__(self, segment, meas):
        super(BetaPhase, self).__init__(segment, meas)
        self.alf_ini = Alfa.get_at(segment.start, meas)
        self.alf_end = Alfa.get_at(segment.end, meas)
        self.bet_ini = BetaPhase.get_at(segment.start, meas)
        self.bet_end = BetaPhase.get_at(segment.end, meas)
        self._is_available = True

    @staticmethod
    def get_at(name, meas, try_amp_kmod=False):
        result = {}
        for plane in PLANES:
            betas = []
            betas.append(BetaPhase._beta_and_error(name, meas.beta, plane))
            if not try_amp_kmod:
                result[plane] = betas[0]
                continue
            with _ignore(IOError):
                if name in meas.kmod_beta[plane].NAME:
                    betas.append(
                        BetaPhase._beta_and_error(name, meas.kmod_beta, plane)
                    )
            with _ignore(IOError):
                if name in meas.amp_beta[plane].NAME:
                    betas.append(
                        BetaPhase._beta_and_error(name, meas.amp_beta, plane)
                    )
            result[plane] = sorted(betas, key=lambda meas: meas.error)[0]
        return result

    def get_init_conds_dict(self):
        init_dict = {}
        for plane in PLANES:
            ini_cond = BetaPhase.get_at(self._segment.start, self._meas,
                                        try_amp_kmod=True)
            ini_name = self._init_pattern.format(plane, "ini")
            init_dict[ini_name] = ini_cond[plane].value
            end_cond = BetaPhase.get_at(self._segment.end, self._meas,
                                        try_amp_kmod=True)
            end_name = self._init_pattern.format(plane, "end")
            init_dict[end_name] = end_cond[plane].value
        return init_dict

    @staticmethod
    def _beta_and_error(name, values, plane):
        uplane = plane.upper()
        beta = values[plane].loc[name, "BET{}".format(uplane)]
        try:
            error = values[plane].loc[name, "ERRBET{}".format(uplane)]
        except KeyError:
            error = values[plane].loc[name, "STDBET{}".format(uplane)]
        return MeasValue(beta, error)

    def get_beatings_at(self, names, ref_model, backwards):
        results = {}
        for plane in PLANES:
            pass
        return results


class Alfa(Measurable):

    _init_pattern = "alf{}_{}"

    def __init__(self, segment, meas):
        super(Alfa, self).__init__(segment, meas)
        self.alf_ini = Alfa.get_at(segment.start, meas)
        self.alf_end = Alfa.get_at(segment.end, meas)
        self.bet_ini = BetaPhase.get_at(segment.start, meas)
        self.bet_end = BetaPhase.get_at(segment.end, meas)
        self._is_available = True

    @staticmethod
    def get_at(name, meas):
        result = {}
        for plane in PLANES:
            uplane = plane.upper()
            alfa = meas.beta[plane].loc[name, "ALF{}".format(uplane)]
            try:
                error = meas.beta[plane].loc[name, "ERRALF{}".format(uplane)]
            except KeyError:
                error = meas.beta[plane].loc[name, "STDALF{}".format(uplane)]
            result[plane] = MeasValue(alfa, error)
        return result

    def get_init_conds_dict(self):
        init_dict = {}
        for plane in PLANES:
            init_dict[self._init_pattern.format(plane, "ini")] =\
                self.alf_ini[plane].value
            # Negative sign for end alfa
            init_dict[self._init_pattern.format(plane, "end")] =\
                -self.alf_end[plane].value
        return init_dict


# Optional measurables #####################################################

class BetaKmod(Measurable):

    def __init__(self, segment, meas):
        super(BetaKmod, self).__init__(segment, meas)
        self.alf_ini = Alfa.get_at(segment.start, meas)
        self.alf_end = Alfa.get_at(segment.end, meas)
        try:
            self._meas.kmod_beta_x
            self._meas.kmod_beta_y
        except IOError:
            _report_file_not_found(BetaKmod)
            self._is_available = False
            return
        self._is_available = True

    @staticmethod
    def get_conditions_at(name, meas):
        result = {}
        for plane in PLANES:
            if name not in meas.kmod_beta[plane].NAME:
                continue
            result[plane] = BetaKmod._beta_and_error(name, meas.kmod_beta, plane)
        return result

    @staticmethod
    def _beta_and_error(name, values, plane):
        uplane = plane.upper()
        beta = values[plane].loc[name, "BET{}".format(uplane)]
        try:
            error = values[plane].loc[name, "ERRBET{}".format(uplane)]
        except KeyError:
            error = values[plane].loc[name, "STDBET{}".format(uplane)]
        return MeasValue(beta, error)

    def get_init_conds_dict(self):
        # Not used for initial conditions
        return {}


class BetaAmp(Measurable):

    def __init__(self, segment, meas):
        super(BetaAmp, self).__init__(segment, meas)
        self.alf_ini = Alfa.get_at(segment.start, meas)
        self.alf_end = Alfa.get_at(segment.end, meas)
        try:
            self._meas.amp_beta_x
            self._meas.amp_beta_y
        except IOError:
            _report_file_not_found(BetaAmp)
            self._is_available = False
            return
        self._is_available = True

    @staticmethod
    def get_at(name, meas):
        result = {}
        for plane in PLANES:
            if name not in meas.amp_beta[plane].NAME:
                continue
            result[plane] = BetaAmp._beta_and_error(name, meas.amp_beta, plane)
        return result

    @staticmethod
    def _beta_and_error(name, values, plane):
        uplane = plane.upper()
        beta = values[plane].loc[name, "BET{}".format(uplane)]
        try:
            error = values[plane].loc[name, "ERRBET{}".format(uplane)]
        except KeyError:
            error = values[plane].loc[name, "STDBET{}".format(uplane)]
        return MeasValue(beta, error)

    def get_init_conds_dict(self):
        # Not used for initial conditions
        return {}


class Dispersion(Measurable):

    _init_pattern = "d{}_{}"

    def __init__(self, segment, meas):
        super(Dispersion, self).__init__(segment, meas)
        try:
            self._meas.disp_x
            self._meas.disp_y
        except IOError:
            _report_file_not_found(Dispersion)
            self._is_available = False
            return
        self.alf_ini = Alfa.get_at(segment.start, meas)
        self.alf_end = Alfa.get_at(segment.end, meas)
        self.bet_ini = BetaPhase.get_at(segment.start, meas)
        self.bet_end = BetaPhase.get_at(segment.end, meas)
        self._is_available = True

    @staticmethod
    def get_at(name, meas):
        result = {}
        for plane in PLANES:
            uplane = plane.upper()
            try:
                disp = meas.disp[plane].loc[name, "D{}".format(uplane)]
                error = meas.disp[plane].loc[name, "STDD{}".format(uplane)]
            except KeyError:
                continue
            result[plane] = MeasValue(disp, error)
        if not result:
            return None
        return result

    def get_init_conds_dict(self):
        init_dict = {}
        for plane in PLANES:
            try:
                ini = Dispersion.get_at(self._segment.start, self._meas)
                end = Dispersion.get_at(self._segment.end, self._meas)
            except KeyError:
                return {}
            ini_name = self._init_pattern.format(plane, "ini")
            init_dict[ini_name] = ini[plane].value
            end_name = self._init_pattern.format(plane, "end")
            init_dict[end_name] = end[plane].value
        return init_dict


class DispersionP(Measurable):

    _init_pattern = "dp{}_{}"

    def __init__(self, segment, meas):
        super(DispersionP, self).__init__(segment, meas)
        try:
            self._meas.disp_x
            self._meas.disp_y
        except IOError:
            _report_file_not_found(DispersionP)
            self._is_available = False
            return
        self._is_available = True

    @staticmethod
    def get_at(name, meas):
        result = {}
        for plane in PLANES:
            uplane = plane.upper()
            try:
                disp = meas.disp[plane].loc[name, "DP{}".format(uplane)]
            except KeyError:
                continue
            result[plane] = MeasValue(disp, None)
        if not result:
            return None
        return result

    def get_init_conds_dict(self):
        init_dict = {}
        for plane in PLANES:
            try:
                ini = DispersionP.get_at(self._segment.start, self._meas)
                end = DispersionP.get_at(self._segment.end, self._meas)
            except KeyError:
                return {}
            ini_name = self._init_pattern.format(plane, "ini")
            init_dict[ini_name] = ini[plane].value
            # Negative sign for end dp
            end_name = self._init_pattern.format(plane, "end")
            init_dict[end_name] = -end[plane].value
        return init_dict


class CouplingTerms(Measurable):

    def __init__(self, segment, meas):
        super(CouplingTerms, self).__init__(segment, meas)
        try:
            self._meas.coupling
        except IOError:
            _report_file_not_found(CouplingTerms)
            self._is_available = False
            return
        self._is_available = True

    @staticmethod
    def get_at(name, meas):
        result = {}
        for term, suffix in (("1001", "1"), ("1010", "2")):
            abs_f = meas.coupling.loc[name, "F{}W".format(term)]
            eabs_f = meas.coupling.loc[name, "FWSTD{}".format(suffix)]
            f_real = meas.coupling.loc[name, "F{}R".format(term)]
            f_imag = meas.coupling.loc[name, "F{}I".format(term)]
            result[term] = MeasCoupling(abs_f, eabs_f, f_real, f_imag)
        return result

    def get_init_conds_dict(self):
        # The coupling f terms are not used for initial conditions
        return {}


class CouplingRMatrix(Measurable):

    def __init__(self, segment, meas):
        super(CouplingRMatrix, self).__init__(segment, meas)
        try:
            self._meas.coupling
        except IOError:
            _report_file_not_found(CouplingRMatrix)
            self._is_available = False
            return
        self._is_available = True

    @staticmethod
    def get_at(name, meas):
        r11, r12, r21, r22 = sbs_math.get_r_terms(name, meas)
        return MeasRMatrix(r11, r12, r21, r22)

    def get_init_conds_dict(self):
        init_dict = {}
        try:
            ini = CouplingRMatrix.get_at(self._segment.start, self._meas)
            end = CouplingRMatrix.get_at(self._segment.end, self._meas)
        except KeyError:
            return init_dict
        init_dict["ini_r11"] = ini.r11
        init_dict["ini_r12"] = ini.r12
        init_dict["ini_r21"] = ini.r21
        init_dict["ini_r22"] = ini.r22
        init_dict["end_r11"] = end.r11
        init_dict["end_r12"] = end.r12
        init_dict["end_r21"] = end.r21
        init_dict["end_r22"] = end.r22
        return init_dict


# Auxiliary functions #########################################################

def _get_beating(values, ref_values):
    beating = (values - ref_values) / ref_values
    beating.dropna()
    return beating


def _common_indices(*indices):
    common = indices[0]
    for index in indices[1:]:
        common = common.intersection(index)
    return common


def _quadratic_add(*values):
    result = 0.
    for value in values:
        result += value ** 2
    return np.sqrt(result)


@contextmanager
def _ignore(exception):
    try:
        yield
    except exception.__class__:
        pass


def _report_file_not_found(cls):
    LOGGER.debug(
        "{}: Needed files not found. Measurable disabled.".format(cls.__name__)
    )


MEASURABLES = (
    Phase,
    BetaPhase,
    BetaKmod,
    BetaAmp,
    Alfa,
    Dispersion,
    DispersionP,
    CouplingTerms,
    CouplingRMatrix,
)
