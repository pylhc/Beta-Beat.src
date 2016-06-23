import os
from .matcher import Matcher
from Python_Classes4MAD import metaclass

DEF_CONSTR_AUX_VALUES_TEMPLATE = """
use, period=%(SEQ_B1)s;
twiss, beta0=%(INIT_VALS_B1)s, chrom, table=%(B1_TABLE_NAME)s;

use, period=%(SEQ_B2)s;
twiss, beta0=%(INIT_VALS_B2)s, chrom, table=%(B2_TABLE_NAME)s;

%(S_VARIABLES)s
%(D_VARIABLES)s
"""

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
ALL_LISTS = os.path.join(CURRENT_PATH, '..', '..', 'MODEL', 'LHCB', 'fullresponse')


class CouplingMatcher(Matcher):

    def __init__(self, name, ip,
                 measurement_data_b1_path, measurement_data_b2_path,
                 match_path,
                 use_errors, front_or_back,
                 exclude_constr_string, exclude_vars_string,
                 all_lists=ALL_LISTS):
        super(CouplingMatcher, self).__init__(
            name, ip,
            measurement_data_b1_path, measurement_data_b2_path,
            match_path,
            use_errors, front_or_back,
            exclude_constr_string, exclude_vars_string, all_lists
        )
        self.f_terms_strings_b1 = self._get_f_terms_strings(1)
        self.f_terms_strings_b2 = self._get_f_terms_strings(2)

    def define_aux_values(self):
        variables_s_str = ""
        variables_d_str = ""

        for variable in self.variables_common[str(self.ip)]:
            variables_s_str += self.name + '.' + variable + '_0' + ' = ' + variable + ';\n'
            variables_d_str += variable + ' := ' + self.name + "." + variable + '_0' + ' + d' + variable + ';\n'

        return DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ_B1": "lhcb1_" + self.front_or_back + "_" + self.name,
            "SEQ_B2": "lhcb2_" + self.front_or_back + "_" + self.name,
            "INIT_VALS_B1": "b1_" + self.ini_end + "_" + self.name,
            "INIT_VALS_B2": "b2_" + self.ini_end + "_" + self.name,
            "B1_TABLE_NAME": self._get_nominal_table_name(1),
            "B2_TABLE_NAME": self._get_nominal_table_name(2),
            "S_VARIABLES": variables_s_str,
            "D_VARIABLES": variables_d_str,
        }

    def _get_nominal_table_name(self, beam):
        return self.name + ".twiss.b" + str(beam)

    def define_variables(self):
        def_variables_string = ""
        for variable in self.variables_common[str(self.ip)]:
            if variable not in self.excluded_variables_list:
                def_variables_string += '    vary, name=d' + variable
                def_variables_string += ', step:=1e-4;\n'
        return def_variables_string

    def define_constraints(self, beam):
        constr_string = ""
        sbs_data = metaclass.twiss(
            os.path.join(self.get_match_data_for_beam(beam).beam_match_sbs_path,
                         'sbscouple_IP' + str(self.ip) + '.out')
        )

        for index in range(0, len(sbs_data.NAME)):
            name = sbs_data.NAME[index]
            if name not in self.excluded_constraints_list:
                name = sbs_data.NAME[index]
                s = sbs_data.S[index]

                f1001r = sbs_data.F1001REMEAS[index]
                f1001i = sbs_data.F1001IMMEAS[index]

                constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                constr_string += 'expr = ' + self.name + "." + name + '_f1001r = ' + str(f1001r) + '; \n'
                constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                constr_string += 'expr = ' + self.name + "." + name + '_f1001i = ' + str(f1001i) + '; \n'
                constr_string += '!   S = ' + str(s)
                constr_string += ';\n'
        return constr_string

    def update_constraints_values(self, beam):
        update_constraints_str = ""
        update_constraints_str += getattr(self, "f_terms_strings_b" + str(beam))
        return update_constraints_str

    def update_variables_definition(self):
        update_vars_str = ""
        for variable in self.variables_common[str(self.ip)]:
            update_vars_str += "        " + variable + ' := ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return update_vars_str

    def generate_changeparameters(self):
        changeparameters_str = ""
        for variable in self.variables_common[str(self.ip)]:
            if variable not in self.excluded_variables_list:
                changeparameters_str += 'select,flag=save,pattern=\"d' + variable + '\";\n'
        return changeparameters_str

    def apply_correction(self):
        apply_correction_str = ""
        for variable in self.variables_common[str(self.ip)]:
            if variable not in self.excluded_variables_list:
                apply_correction_str += variable + ' = ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return apply_correction_str

    def _get_f_terms_strings(self, beam):
        f_terms_string = ""
        sbs_data = metaclass.twiss(
            os.path.join(self.get_match_data_for_beam(beam).beam_match_sbs_path,
                         'sbscouple_IP' + str(self.ip) + '.out')
        )
        for bpm_name in sbs_data.NAME:
            f_terms_string += "exec, get_f_terms_for(twiss, " + bpm_name + ");\n"
            f_terms_string += self.name + "." + bpm_name + "_f1001r = " + bpm_name + "_f1001r;\n"
            f_terms_string += self.name + "." + bpm_name + "_f1001i = " + bpm_name + "_f1001i;\n"
            f_terms_string += self.name + "." + bpm_name + "_f1010r = " + bpm_name + "_f1010r;\n"
            f_terms_string += self.name + "." + bpm_name + "_f1010i = " + bpm_name + "_f1010i;\n"
        return f_terms_string

    @classmethod
    def from_matcher_dict(cls, matcher_name, matcher_dict, match_path):
        instance = Matcher.from_matcher_dict(matcher_name, matcher_dict, match_path)
        instance.__class__ = cls
        instance.f_terms_strings_b1 = instance._get_f_terms_strings(1)
        instance.f_terms_strings_b2 = instance._get_f_terms_strings(2)
        return instance
