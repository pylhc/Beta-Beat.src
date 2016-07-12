import os
from .matcher import Matcher
from .kmod_matcher import KmodMatcher
from Python_Classes4MAD import metaclass


class AmpMatcher(KmodMatcher):

    @Matcher.override(KmodMatcher)
    def __init__(self, matcher_name, matcher_dict, match_path):
        super(AmpMatcher, self).__init__(matcher_name, matcher_dict, match_path)
        self._sbs_kmod_data = {}
        for beam in self.get_beams():
            self._sbs_kmod_data[beam] = {}
            for plane in ["x", "y"]:
                sbs_kmod_data_path = os.path.join(
                    self.get_match_data(beam).get_beam_match_sbs_path(),
                    'sbsampbetabeat' + plane + '_IP' + str(self._ip) + '.out'
                )
                self._sbs_kmod_data[beam][plane] = metaclass.twiss(sbs_kmod_data_path)

    @Matcher.override(KmodMatcher)
    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            this_kmod_data = self._sbs_kmod_data[beam][plane]
            for name in this_kmod_data.NAME:
                index = this_kmod_data.indx[name]
                beta_beating = getattr(this_kmod_data, "BETABEATAMP" + plane.upper())[index]
                err_beta_beating = getattr(this_kmod_data, "ERRBETABEATAMP" + plane.upper())[index]
                s = this_kmod_data.S[index]
                weight = self.get_constraint_weight(beta_beating, err_beta_beating,
                                                    lambda value: abs(value) <= 0.25)
                if weight is None:
                    weight = 1.0

                constr_string += '    constraint, weight = ' + str(weight) + ' , '
                constr_string += 'expr =  ' + self._name + self._get_suffix() + plane + name + ' = ' + str(beta_beating) + '; '
                constr_string += '!   S = ' + str(s)
                constr_string += ';\n'
        return constr_string

    @Matcher.override(KmodMatcher)
    def _get_suffix(self):
        return ".ampbeating"
