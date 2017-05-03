import os
from .matcher import Matcher
from .kmod_matcher import KmodMatcher
from Python_Classes4MAD import metaclass


class AmpMatcher(KmodMatcher):

    BETA_BEATING_CONSTR_WEIGHT = 1.

    @Matcher.override(KmodMatcher)
    def __init__(self, matcher_name, matcher_dict, match_path):
        Matcher.__init__(self, matcher_name, matcher_dict, match_path)
        self._sbs_amp_data = {}
        for beam in self.get_beams():
            self._sbs_amp_data[beam] = {}
            for plane in ["x", "y"]:
                sbs_amp_data_path = os.path.join(
                    self.get_match_data(beam).get_beam_match_sbs_path(),
                    'sbsampbetabeat' + plane + '_IP' + str(self._ip) + '.out'
                )
                self._sbs_amp_data[beam][plane] = metaclass.twiss(sbs_amp_data_path)

    @Matcher.override(KmodMatcher)
    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            this_amp_data = self._sbs_amp_data[beam][plane]
            for name in this_amp_data.NAME:
                index = this_amp_data.indx[name]
                beta_beating = getattr(this_amp_data, "BETABEATAMP" + plane.upper())[index]
                err_beta_beating = getattr(this_amp_data, "ERRBETABEATAMP" + plane.upper())[index]
                constr_string += self._get_constraint_instruction(
                    self._name + self._get_suffix() + plane + name,
                    beta_beating, err_beta_beating)

        return constr_string

    @Matcher.override(KmodMatcher)
    def _get_suffix(self):
        return ".ampbeating"
