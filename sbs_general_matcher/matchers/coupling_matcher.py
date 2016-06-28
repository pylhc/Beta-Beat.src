import os
import json
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


class CouplingMatcher(Matcher):

    ALL_LISTS = os.path.join(CURRENT_PATH, '..', '..', 'MODEL', 'LHCB', 'fullresponse')

    @Matcher.override(Matcher)
    def __init__(self, matcher_name, matcher_dict, match_path):
        super(CouplingMatcher, self).__init__(matcher_name, matcher_dict, match_path)
        self.f_terms_strings_b1 = self._get_f_terms_strings(1)
        self.f_terms_strings_b2 = self._get_f_terms_strings(2)
        if "all_lists" in matcher_dict:
            all_lists = matcher_dict["all_lists"]
        else:
            all_lists = CouplingMatcher.ALL_LISTS
        all_list_file = json.load(
            file(os.path.join(all_lists, "LHCB1", "AllLists_couple.json"), 'r')
        )
        self._variables_common = all_list_file['getListsByIR'][str(self.get_ip())]

    @Matcher.override(Matcher)
    def define_aux_vars(self):
        variables_s_str = ""
        variables_d_str = ""

        for variable in self.get_variables():
            variables_s_str += self.get_name() + '.' + variable + '_0' + ' = ' + variable + ';\n'
            variables_d_str += variable + ' := ' + self.get_name() + "." + variable + '_0' + ' + d' + variable + ';\n'

        return DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ_B1": "lhcb1_" + self.get_front_or_back() + "_" + self.get_name(),
            "SEQ_B2": "lhcb2_" + self.get_front_or_back() + "_" + self.get_name(),
            "INIT_VALS_B1": "b1_" + self.get_ini_end() + "_" + self.get_name(),
            "INIT_VALS_B2": "b2_" + self.get_ini_end() + "_" + self.get_name(),
            "B1_TABLE_NAME": self._get_nominal_table_name(1),
            "B2_TABLE_NAME": self._get_nominal_table_name(2),
            "S_VARIABLES": variables_s_str,
            "D_VARIABLES": variables_d_str,
        }

    def _get_nominal_table_name(self, beam):
        return self.get_name() + ".twiss.b" + str(beam)

    @Matcher.override(Matcher)
    def get_variables(self):
        variable_list = []
        for variable in self._get_common_variables():
            if variable not in self._excluded_variables_list:
                variable_list.append(variable)
        return variable_list

    def _get_common_variables(self):
        return self._variables_common

    @Matcher.override(Matcher)
    def define_constraints(self, beam):
        constr_string = ""
        sbs_data = metaclass.twiss(
            os.path.join(self.get_match_data(beam).get_beam_match_sbs_path(),
                         'sbscouple_IP' + str(self.get_ip()) + '.out')
        )

        for index in range(0, len(sbs_data.NAME)):
            name = sbs_data.NAME[index]
            if name not in self._excluded_constraints_list:
                name = sbs_data.NAME[index]
                s = sbs_data.S[index]

                f1001r = sbs_data.F1001REMEAS[index]
                f1001i = sbs_data.F1001IMMEAS[index]

                constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                constr_string += 'expr = ' + self.get_name() + "." + name + '_f1001r = ' + str(f1001r) + '; \n'
                constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                constr_string += 'expr = ' + self.get_name() + "." + name + '_f1001i = ' + str(f1001i) + '; \n'
                constr_string += '!   S = ' + str(s)
                constr_string += ';\n'
        return constr_string

    @Matcher.override(Matcher)
    def update_constraints_values(self, beam):
        update_constraints_str = ""
        update_constraints_str += getattr(self, "f_terms_strings_b" + str(beam))
        return update_constraints_str

    @Matcher.override(Matcher)
    def update_variables_definition(self):
        update_vars_str = ""
        for variable in self.get_variables():
            update_vars_str += "        " + variable + ' := ' + self.get_name() + "." + variable + '_0 + d' + variable + ';\n'
        return update_vars_str

    @Matcher.override(Matcher)
    def generate_changeparameters(self):
        changeparameters_str = ""
        for variable in self.variables_common[str(self.ip)]:
            if variable not in self.excluded_variables_list:
                changeparameters_str += 'select,flag=save,pattern=\"d' + variable + '\";\n'
        return changeparameters_str

    @Matcher.override(Matcher)
    def apply_correction(self):
        apply_correction_str = ""
        for variable in self.variables_common[str(self.ip)]:
            if variable not in self.excluded_variables_list:
                apply_correction_str += variable + ' = ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return apply_correction_str

    def _get_f_terms_strings(self, beam):
        f_terms_string = ""
        sbs_data = metaclass.twiss(
            os.path.join(self.get_match_data(beam).get_beam_match_sbs_path(),
                         'sbscouple_IP' + str(self.get_ip()) + '.out')
        )
        for bpm_name in sbs_data.NAME:
            f_terms_string += "exec, get_f_terms_for(twiss, " + bpm_name + ");\n"
            f_terms_string += self.get_name() + "." + bpm_name + "_f1001r = " + bpm_name + "_f1001r;\n"
            f_terms_string += self.get_name() + "." + bpm_name + "_f1001i = " + bpm_name + "_f1001i;\n"
            f_terms_string += self.get_name() + "." + bpm_name + "_f1010r = " + bpm_name + "_f1010r;\n"
            f_terms_string += self.get_name() + "." + bpm_name + "_f1010i = " + bpm_name + "_f1010i;\n"
        return f_terms_string
