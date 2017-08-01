import os
from Python_Classes4MAD import metaclass
from .matcher import Matcher
from .phase_matcher import PhaseMatcher


class KmodMatcher(PhaseMatcher):

    BETA_BEATING_CONSTR_WEIGHT = 1.

    @Matcher.override(PhaseMatcher)
    def define_aux_vars(self):
        beatings_str = ""
        beam = self.get_segment().get_beam()
        for plane in ["x", "y"]:
            for name in self._get_kmod_data(plane).NAME:
                beatings_str += self._name + self._get_suffix() + plane + name + ' := '
                beatings_str += "(table(twiss, " + name + ", bet" + plane + ")"
                beatings_str += " - table(" + self._get_nominal_table_name(beam) + ", " + name + ", bet" + plane + ")) /\n"
                beatings_str += "table(" + self._get_nominal_table_name(beam) + ", " + name + ", bet" + plane + ");\n"

        variables_s_str = ""
        for variable in self.get_variables():
            variables_s_str += self.get_name() + '.' + variable + '_0' + ' = ' + variable + ';\n'

        return PhaseMatcher.DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ": "lhcb" + str(beam) + "_" + self._front_or_back + "_" + self._name,
            "INIT_VALS": "b" + str(beam) + "_" + self._ini_end + "_" + self._name,
            "TABLE_NAME": self._get_nominal_table_name(),
            "PHASES": beatings_str,
            "S_VARIABLES": variables_s_str,
        }

    @Matcher.override(PhaseMatcher)
    def define_constraints(self):
        constr_string = ""
        for plane in ["x", "y"]:
            this_kmod_data = self._get_kmod_data(plane)
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

    def _get_kmod_data(self, plane):
        sbs_kmod_data_path = os.path.join(
            os.path.join(self._matcher_path, "sbs"),
            'sbskmodbetabeat' + plane + '_' + self._segment.label + '.out'
        )
        return metaclass.twiss(sbs_kmod_data_path)
