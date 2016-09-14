import os
from .matcher import Matcher
from .phase_matcher import PhaseMatcher
from Python_Classes4MAD import metaclass


class KmodMatcher(PhaseMatcher):

    BETA_BEATING_CONSTR_WEIGHT = 1.

    @Matcher.override(PhaseMatcher)
    def __init__(self, matcher_name, matcher_dict, match_path):
        super(KmodMatcher, self).__init__(matcher_name, matcher_dict, match_path)
        self._sbs_kmod_data = {}
        for beam in self.get_beams():
            self._sbs_kmod_data[beam] = {}
            for plane in ["x", "y"]:
                sbs_kmod_data_path = os.path.join(
                    self.get_match_data(beam).get_beam_match_sbs_path(),
                    'sbskmodbetabeat' + plane + '_IP' + str(self._ip) + '.out'
                )
                self._sbs_kmod_data[beam][plane] = metaclass.twiss(sbs_kmod_data_path)

    @Matcher.override(PhaseMatcher)
    def define_aux_vars(self):
        beatings_str = ""
        for beam in self.get_beams():
            for plane in ["x", "y"]:
                for name in self._sbs_kmod_data[beam][plane].NAME:
                    beatings_str += self._name + self._get_suffix() + plane + name + ' := '
                    beatings_str += "(table(twiss, " + name + ", bet" + plane + ")"
                    beatings_str += " - table(" + self._get_nominal_table_name(beam) + ", " + name + ", bet" + plane + ")) /\n"
                    beatings_str += "table(" + self._get_nominal_table_name(beam) + ", " + name + ", bet" + plane + ");\n"

        variables_s_str = ""
        for variable in self.get_all_variables():
            variables_s_str += self.get_name() + '.' + variable + '_0' + ' = ' + variable + ';\n'

        return PhaseMatcher.DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ_B1": "lhcb1_" + self.get_front_or_back() + "_" + self.get_name(),
            "SEQ_B2": "lhcb2_" + self.get_front_or_back() + "_" + self.get_name(),
            "INIT_VALS_B1": "b1_" + self.get_ini_end() + "_" + self.get_name(),
            "INIT_VALS_B2": "b2_" + self.get_ini_end() + "_" + self.get_name(),
            "B1_TABLE_NAME": self._get_nominal_table_name(1),
            "B2_TABLE_NAME": self._get_nominal_table_name(2),
            "PHASES": beatings_str,
            "S_VARIABLES": variables_s_str,
        }

    @Matcher.override(PhaseMatcher)
    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            this_kmod_data = self._sbs_kmod_data[beam][plane]
            for name in this_kmod_data.NAME:
                index = this_kmod_data.indx[name]
                beta_beating = getattr(this_kmod_data, "BETABEAT" + plane.upper())[index]
                err_beta_beating = getattr(this_kmod_data, "ERRBETABEAT" + plane.upper())[index]
                constr_string += self._get_constraint_instruction(
                    self._name + self._get_suffix() + plane + name,
                    beta_beating, err_beta_beating)

        return constr_string

    def _get_suffix(self):
        return ".kmodbeating"
