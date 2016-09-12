import os
import json
from .matcher import Matcher
from Python_Classes4MAD import metaclass

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class PhaseMatcher(Matcher):

    ALL_LISTS = os.path.join(CURRENT_PATH, '..', '..', 'MODEL', 'LHCB', 'fullresponse')

    @Matcher.override(Matcher)
    def __init__(self, matcher_name, matcher_dict, match_path):
        super(PhaseMatcher, self).__init__(matcher_name, matcher_dict, match_path)
        self._variables = {}
        all_lists = None
        if "all_lists" in matcher_dict:
            all_lists = matcher_dict["all_lists"]
        else:
            all_lists = PhaseMatcher.ALL_LISTS
        all_list_file_beam1 = json.load(
            file(os.path.join(all_lists, "LHCB1", "AllLists.json"), 'r')
        )
        all_list_file_beam2 = json.load(
            file(os.path.join(all_lists, "LHCB2", "AllLists.json"), 'r')
        )
        self._variables[1] = all_list_file_beam1['getListsByIR'][1][str(self.get_ip())]
        self._variables[2] = all_list_file_beam2['getListsByIR'][1][str(self.get_ip())]
        self._variables_common = all_list_file_beam1['getListsByIR'][0][str(self.get_ip())]

    DEF_CONSTR_AUX_VALUES_TEMPLATE = """
    use, period=%(SEQ_B1)s;
    twiss, beta0=%(INIT_VALS_B1)s, chrom, table=%(B1_TABLE_NAME)s;

    use, period=%(SEQ_B2)s;
    twiss, beta0=%(INIT_VALS_B2)s, chrom, table=%(B2_TABLE_NAME)s;

    %(PHASES)s

    %(S_VARIABLES)s
    """

    @Matcher.override(Matcher)
    def define_aux_vars(self):
        phases_str = ""
        variables_s_str = ""

        sign = ""
        if "b" in self._front_or_back:
            sign = "-"
        for beam in self.get_beams():
            for plane in ["x", "y"]:
                sbs_data_path = os.path.join(
                    self.get_match_data(beam).get_beam_match_sbs_path(),
                    'sbsphase' + plane + 't_IP' + str(self._ip) + '.out'
                )
                sbs_data = metaclass.twiss(sbs_data_path)
                for name in sbs_data.NAME:
                    phases_str += self._name + '.dmu' + plane + name + ' := '
                    phases_str += sign + "(table(twiss, " + name + ", mu" + plane + ") - "
                    phases_str += "table(" + self._get_nominal_table_name(beam) + ", " + name + ", mu" + plane + "));\n"

        for variable in self.get_all_variables():
            variables_s_str += self._name + '.' + variable + '_0' + ' = ' + variable + ';\n'

        return PhaseMatcher.DEF_CONSTR_AUX_VALUES_TEMPLATE % {
            "SEQ_B1": "lhcb1_" + self._front_or_back + "_" + self._name,
            "SEQ_B2": "lhcb2_" + self._front_or_back + "_" + self._name,
            "INIT_VALS_B1": "b1_" + self._ini_end + "_" + self._name,
            "INIT_VALS_B2": "b2_" + self._ini_end + "_" + self._name,
            "B1_TABLE_NAME": self._get_nominal_table_name(1),
            "B2_TABLE_NAME": self._get_nominal_table_name(2),
            "PHASES": phases_str,
            "S_VARIABLES": variables_s_str,
        }

    def _get_nominal_table_name(self, beam):
        return self._name + ".twiss.b" + str(beam)

    @Matcher.override(Matcher)
    def get_all_variables(self):
        variable_list = []
        for variable in self.get_common_variables():
            if variable not in self._excluded_variables_list:
                variable_list.append(variable)
        for beam in self.get_beams():
            for variable in self.get_variables_for_beam(beam):
                if variable not in self._excluded_variables_list:
                    variable_list.append(variable)
        return variable_list

    @Matcher.override(Matcher)
    def get_variables_for_beam(self, beam):
        return self._variables[beam]

    @Matcher.override(Matcher)
    def get_common_variables(self):
        return self._variables_common

    @Matcher.override(Matcher)
    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            sbs_data = metaclass.twiss(
                os.path.join(self.get_match_data(beam).get_beam_match_sbs_path(),
                             'sbsphase' + plane + 't_IP' + str(self._ip) + '.out')
            )

            is_back = "b" in self._front_or_back
            for index in range(0, len(sbs_data.NAME)):
                name = sbs_data.NAME[index]
                if name not in self._excluded_constraints_list:
                    if is_back is not True:
                        phase = sbs_data.PROPPHASEX[index] if plane == "x" else sbs_data.PROPPHASEY[index]
                        error = sbs_data.ERRPROPPHASEX[index] if plane == "x" else sbs_data.ERRPROPPHASEY[index]
                    else:
                        phase = sbs_data.BACKPHASEX[index] if plane == "x" else sbs_data.BACKPHASEY[index]
                        error = sbs_data.ERRBACKPHASEX[index] if plane == "x" else sbs_data.ERRBACKPHASEY[index]
                    s = sbs_data.S[index]

                    if self._use_errors:
                        constr_string += '    constraint, weight = 1.0, '
                        constr_string += 'expr =  ' + self._name + '.dmu' + plane + name + ' > ' + str(phase - error) + ';'
                        constr_string += '!   S = ' + str(s) + ';\n'

                        constr_string += '    constraint, weight = 1.0, '
                        constr_string += 'expr =  ' + self._name + '.dmu' + plane + name + ' < ' + str(phase + error) + ';'
                        constr_string += '!   S = ' + str(s) + ';\n'
                    else:
                        constr_string += '    constraint, weight = 1.0, '
                        constr_string += 'expr =  ' + self._name + '.dmu' + plane + name + ' = ' + str(phase) + ';'
                        constr_string += '!   S = ' + str(s) + ';\n'
        return constr_string

    @Matcher.override(Matcher)
    def update_constraints_values(self, beam):
        return ""

    @Matcher.override(Matcher)
    def update_variables_definition(self):
        update_vars_str = ""
        for variable in self.get_all_variables():
            update_vars_str += "        " + variable + ' := ' + self._name + "." + variable + '_0 + d' + variable + ';\n'
        return update_vars_str

    @Matcher.override(Matcher)
    def generate_changeparameters(self):
        changeparameters_str = ""
        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self._ip)]:
                if variable not in self.excluded_variables_list:
                    changeparameters_str += 'select,flag=save,pattern=\"d' + variable + '\";\n'
        return changeparameters_str

    @Matcher.override(Matcher)
    def apply_correction(self):
        apply_correction_str = ""
        for variable_list in [self.variables_beam1, self.variables_beam2,
                              self.variables_common]:
            for variable in variable_list[str(self._ip)]:
                if variable not in self.excluded_variables_list:
                    apply_correction_str += variable + ' = ' + self._name + "." + variable + '_0 + d' + variable + ';\n'
        return apply_correction_str
