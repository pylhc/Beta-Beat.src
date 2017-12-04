import os
from .matcher import Matcher
from Python_Classes4MAD import metaclass
import math

DEF_CONSTR_AUX_VALUES_TEMPLATE = """
use, period=%(SEQ)s;
twiss, beta0=%(INIT_VALS)s, chrom, table=%(TABLE_NAME)s;

%(S_VARIABLES)s
%(D_VARIABLES)s
"""

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
USE_ABS = True
USE_F1010 = False


class CouplingMatcher(Matcher):

    COUP_CORR_CLASSES = ["MQSX"]

    @Matcher.override(Matcher)
    def get_variables(self, exclude=True):
        variables = self.segment.get_segment_vars(
            classes=CouplingMatcher.COUP_CORR_CLASSES,
        )
        if exclude:
            variables = [
                var for var in variables
                if var not in self.excluded_variables
            ]
        return variables

    @Matcher.override(Matcher)
    def define_aux_vars(self):
        variables_s_str = ""
        variables_d_str = ""

        for variable in self.get_variables():
            variables_s_str += self.name + '.' + variable + '_0' + ' = ' + variable + ';\n'
            variables_d_str += variable + ' := ' + self.name + "." + variable + '_0' + ' + d' + variable + ';\n'

        return DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ": self.get_sequence_name(),
            "INIT_VALS": self.get_initvals_name(),
            "TABLE_NAME": self._get_nominal_table_name(),
            "S_VARIABLES": variables_s_str,
            "D_VARIABLES": variables_d_str,
        }

    @Matcher.override(Matcher)
    def define_constraints(self):
        constr_string = ""
        sbs_data = metaclass.twiss(
            os.path.join(os.path.join(self.matcher_path, "sbs"),
                         'sbscouple_' + self.segment.label + '.out')
        )

        for index in range(0, len(sbs_data.NAME)):
            name = sbs_data.NAME[index]
            if name not in self.excluded_constraints:
                name = sbs_data.NAME[index]
                s = sbs_data.S[index]

                f1001r = sbs_data.F1001REMEAS[index]
                f1001i = sbs_data.F1001IMMEAS[index]
                f1001abs = math.sqrt(f1001r ** 2 + f1001i ** 2)
                f1010r = sbs_data.F1010REMEAS[index]
                f1010i = sbs_data.F1010IMMEAS[index]
                f1010abs = math.sqrt(f1010r ** 2 + f1010i ** 2)

                if (USE_ABS):
                    constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                    constr_string += 'expr = ' + self.name + "." + name + '_f1001abs = ' + str(f1001abs) + '; \n'
                    if USE_F1010:
                        constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                        constr_string += 'expr = ' + self.name + "." + name + '_f1010abs = ' + str(f1010abs) + '; \n'
                else:
                    constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                    constr_string += 'expr = ' + self.name + "." + name + '_f1001r = ' + str(f1001r) + '; \n'
                    constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                    constr_string += 'expr = ' + self.name + "." + name + '_f1001i = ' + str(f1001i) + '; \n'
                    if USE_F1010:
                        constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                        constr_string += 'expr = ' + self.name + "." + name + '_f1010r = ' + str(f1010r) + '; \n'
                        constr_string += '   constraint, weight = ' + str(1.0) + ' , '
                        constr_string += 'expr = ' + self.name + "." + name + '_f1010i = ' + str(f1010i) + '; \n'

                constr_string += '!   S = ' + str(s)
                constr_string += ';\n'
        return constr_string

    @Matcher.override(Matcher)
    def update_constraints_values(self):
        update_constraints_str = ""
        update_constraints_str += self._get_f_terms_strings()
        return update_constraints_str

    @Matcher.override(Matcher)
    def update_variables_definition(self):
        update_vars_str = ""
        for variable in self.get_variables():
            update_vars_str += "        " + variable + ' := ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return update_vars_str

    @Matcher.override(Matcher)
    def generate_changeparameters(self):
        changeparameters_str = ""
        for variable in self.get_variables():
            changeparameters_str += 'select,flag=save,pattern=\"d' + variable + '\";\n'
        return changeparameters_str

    @Matcher.override(Matcher)
    def apply_correction(self):
        apply_correction_str = ""
        for variable in self.get_variables():
            apply_correction_str += variable + ' = ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return apply_correction_str

    def _get_f_terms_strings(self):
        f_terms_string = ""
        sbs_data = metaclass.twiss(
            os.path.join(os.path.join(self.matcher_path, "sbs"),
                         'sbscouple_' + self.segment.label + '.out')
        )
        for bpm_name in sbs_data.NAME:
            f_terms_string += "exec, get_f_terms_for(twiss, " + bpm_name + ");\n"
            if (USE_ABS):
                f_terms_string += self.name + "." + bpm_name + "_f1001abs = " + bpm_name + "_f1001abs;\n"
                f_terms_string += self.name + "." + bpm_name + "_f1010abs = " + bpm_name + "_f1010abs;\n"
            else:
                f_terms_string += self.name + "." + bpm_name + "_f1001r = " + bpm_name + "_f1001r;\n"
                f_terms_string += self.name + "." + bpm_name + "_f1001i = " + bpm_name + "_f1001i;\n"
                f_terms_string += self.name + "." + bpm_name + "_f1010r = " + bpm_name + "_f1010r;\n"
                f_terms_string += self.name + "." + bpm_name + "_f1010i = " + bpm_name + "_f1010i;\n"
        return f_terms_string
