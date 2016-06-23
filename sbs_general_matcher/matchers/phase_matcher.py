import os
from .matcher import Matcher
from Python_Classes4MAD import metaclass

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class PhaseMatcher(Matcher):

    DEF_CONSTR_AUX_VALUES_TEMPLATE = """
    use, period=%(SEQ_B1)s;
    twiss, beta0=%(INIT_VALS_B1)s, chrom, table=%(B1_TABLE_NAME)s;

    use, period=%(SEQ_B2)s;
    twiss, beta0=%(INIT_VALS_B2)s, chrom, table=%(B2_TABLE_NAME)s;

    %(PHASES)s

    %(S_VARIABLES)s
    %(D_VARIABLES)s
    """

    def define_aux_values(self):
        phases_str = ""
        variables_s_str = ""
        variables_d_str = ""

        sign = ""
        if "b" in self.front_or_back:
            sign = "-"
        for beam in [1, 2]:
            for plane in ["x", "y"]:
                sbs_data_path = os.path.join(
                    self.get_match_data_for_beam(beam).beam_match_sbs_path,
                    'sbsphase' + plane + 't_IP' + str(self.ip) + '.out'
                )
                sbs_data = metaclass.twiss(sbs_data_path)
                for name in sbs_data.NAME:
                    phases_str += self.name + '.dmu' + plane + name + ' := '
                    phases_str += sign + "(table(twiss, " + name + ", mu" + plane + ") - "
                    phases_str += "table(" + self._get_nominal_table_name(beam) + ", " + name + ", mu" + plane + "));\n"

        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self.ip)]:
                variables_s_str += self.name + '.' + variable + '_0' + ' = ' + variable + ';\n'
                variables_d_str += variable + ' := ' + self.name + "." + variable + '_0' + ' + d' + variable + ';\n'

        return PhaseMatcher.DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ_B1": "lhcb1_" + self.front_or_back + "_" + self.name,
            "SEQ_B2": "lhcb2_" + self.front_or_back + "_" + self.name,
            "INIT_VALS_B1": "b1_" + self.ini_end + "_" + self.name,
            "INIT_VALS_B2": "b2_" + self.ini_end + "_" + self.name,
            "B1_TABLE_NAME": self._get_nominal_table_name(1),
            "B2_TABLE_NAME": self._get_nominal_table_name(2),
            "PHASES": phases_str,
            "S_VARIABLES": variables_s_str,
            "D_VARIABLES": variables_d_str,
        }

    def _get_nominal_table_name(self, beam):
        return self.name + ".twiss.b" + str(beam)

    def define_variables(self):
        def_variables_string = ""
        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self.ip)]:
                if variable not in self.excluded_variables_list:
                    def_variables_string += '    vary, name=d' + variable
                    def_variables_string += ', step:=1e-4;\n'
        return def_variables_string

    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            sbs_data = metaclass.twiss(
                os.path.join(self.get_match_data_for_beam(beam).beam_match_sbs_path,
                             'sbsphase' + plane + 't_IP' + str(self.ip) + '.out')
            )

            is_back = "b" in self.front_or_back
            for index in range(0, len(sbs_data.NAME)):
                name = sbs_data.NAME[index]
                if name not in self.excluded_constraints_list:
                    if is_back is not True:
                        phase = sbs_data.PROPPHASEX[index] if plane == "x" else sbs_data.PROPPHASEY[index]
                        error = sbs_data.ERRPROPPHASEX[index] if plane == "x" else sbs_data.ERRPROPPHASEY[index]
                    else:
                        phase = sbs_data.BACKPHASEX[index] if plane == "x" else sbs_data.BACKPHASEY[index]
                        error = sbs_data.ERRBACKPHASEX[index] if plane == "x" else sbs_data.ERRBACKPHASEY[index]
                    s = sbs_data.S[index]

                    weight = self.get_constraint_weight(phase, error, lambda value: abs(value) <= 0.25)

                    weight = 1.0
                    constr_string += '    constraint, weight = ' + str(weight) + ' , '
                    constr_string += 'expr =  ' + self.name + '.dmu' + plane + name + ' = ' + str(phase) + '; '

                    constr_string += '!   S = ' + str(s)
                    constr_string += ';\n'
        return constr_string

    def update_constraints_values(self, beam):
        return ""

    def update_variables_definition(self):
        update_vars_str = ""
        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self.ip)]:
                update_vars_str += "        " + variable + ' := ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return update_vars_str

    def generate_changeparameters(self):
        changeparameters_str = ""
        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self.ip)]:
                if variable not in self.excluded_variables_list:
                    changeparameters_str += 'select,flag=save,pattern=\"d' + variable + '\";\n'
        return changeparameters_str

    def apply_correction(self):
        apply_correction_str = ""
        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self.ip)]:
                if variable not in self.excluded_variables_list:
                    apply_correction_str += variable + ' = ' + self.name + "." + variable + '_0 + d' + variable + ';\n'
        return apply_correction_str
