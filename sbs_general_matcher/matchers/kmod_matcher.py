import os
from Python_Classes4MAD import metaclass
from .matcher import Matcher
from .phase_matcher import PhaseMatcher


class KmodMatcher(PhaseMatcher):

    BETA_BEATING_CONSTR_WEIGHT = 1.

    BETA_BEATING_TMPL = (
        "{varname} := ((table(twiss, {bpm_name}, bet{plane}) - table({nominal_table_name}, {bpm_name}, bet{plane})) / "
        "table({nominal_table_name}, {bpm_name}, bet{plane})) / {error};"
    )

    @Matcher.override(PhaseMatcher)
    def define_aux_vars(self):
        beatings_str = ""
        for plane in ["x", "y"]:
            this_kmod_data = self._get_kmod_data(plane)
            for name in self._get_kmod_data(plane).NAME:
                index = this_kmod_data.indx[name]
                err_beta_beating = 1.
                if self.use_errors:
                    err_beta_beating = getattr(this_kmod_data, "ERRBETABEAT" + plane.upper())[index]
                beatings_str += KmodMatcher.BETA_BEATING_TMPL.format(
                    varname=self.name + self._get_suffix() + plane + name,
                    bpm_name = name,
                    nominal_table_name=self._get_nominal_table_name(),
                    plane=plane,
                    error=err_beta_beating,
                ) + "\n"
        variables_s_str = ""
        for variable in self.get_variables():
            variables_s_str += self.name + '.' + variable + '_0' + ' = ' + variable + ';\n'

        aux_vars_str = PhaseMatcher.SEGMENT_TWISS_TMPL.format(
            seq=self.get_sequence_name(),
            init_vals=self.get_initvals_name(),
            table_name=self._get_nominal_table_name(),
        )
        aux_vars_str += beatings_str
        aux_vars_str += variables_s_str
        return aux_vars_str

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
                    self.name + self._get_suffix() + plane + name,
                    beta_beating, err_beta_beating)

        return constr_string

    def _get_suffix(self):
        return ".kmodbeating"

    def _get_kmod_data(self, plane):
        sbs_kmod_data_path = os.path.join(
            os.path.join(self.matcher_path, "sbs"),
            'sbskmodbetabeat' + plane + '_' + self.segment.label + '.out'
        )
        return metaclass.twiss(sbs_kmod_data_path)
