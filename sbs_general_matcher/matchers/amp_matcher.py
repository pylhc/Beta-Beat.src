import os
from .matcher import Matcher
from .kmod_matcher import KmodMatcher
from Python_Classes4MAD import metaclass


class AmpMatcher(KmodMatcher):

    BETA_BEATING_CONSTR_WEIGHT = 1.

    @Matcher.override(KmodMatcher)
    def define_constraints(self):
        constr_string = ""
        for plane in ["x", "y"]:
            this_amp_data = self._get_amp_data(plane)
            for name in this_amp_data.NAME:
                index = this_amp_data.indx[name]
                beta_beating = getattr(this_amp_data, "BETABEATAMP" + plane.upper())[index]
                err_beta_beating = getattr(this_amp_data, "ERRBETABEATAMP" + plane.upper())[index]
                constr_string += self._get_constraint_instruction(
                    self.name + self._get_suffix() + plane + name,
                    beta_beating, err_beta_beating)

        return constr_string

    @Matcher.override(KmodMatcher)
    def _get_suffix(self):
        return ".ampbeating"

    def _get_amp_data(self, plane):
        sbs_amp_data_path = os.path.join(
            os.path.join(self.matcher_path, "sbs"),
            'sbsampbetabeat' + plane + '_' + self.segment.label + '.out'
        )
        return metaclass.twiss(sbs_amp_data_path)
