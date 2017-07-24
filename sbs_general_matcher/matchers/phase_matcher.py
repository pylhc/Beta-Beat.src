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
    """

    BETA_CORR_CLASSES = ["MQT", "MQM", "MQY", "MQX", "MQXT"]

    @Matcher.override(Matcher)
    def get_variables(self, exclude=True):
        variables = self._segment.get_segment_vars(
            classes=PhaseMatcher.BETA_CORR_CLASSES,
        )
        if exclude:
            variables = [
                var for var in variables
                if var not in self._excluded_variables_list
            ]
        return variables

    @Matcher.override(Matcher)
    def define_aux_vars(self):
        phases_str = ""
        variables_s_str = ""

        sign = ""
        if "b" in self._front_or_back:
            sign = "-"
        for plane in ["x", "y"]:
            sbs_data_path = os.path.join(
                os.path.join(self.get_matcher_path(), "sbs"),
                'sbsphase' + plane + 't_' + str(self._segment.label) + '.out'
            )
            sbs_data = metaclass.twiss(sbs_data_path)
            for name in sbs_data.NAME:
                phases_str += self._name + '.dmu' + plane + name + ' := '
                phases_str += sign + "(table(twiss, " + name + ", mu" + plane + ") - "
                phases_str += "table(" + self._get_nominal_table_name() + ", " + name + ", mu" + plane + "));\n"

        for variable in self.get_variables():
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

    def _get_nominal_table_name(self, beam=None):
        if beam is None:
            beam = self._segment.get_beam()
        return self._name + ".twiss.b" + str(beam)

    @Matcher.override(Matcher)
    def define_constraints(self):
        constr_string = ""
        for plane in ["x", "y"]:
            sbs_data = metaclass.twiss(
                os.path.join(self.get_matcher_path(), "sbs",
                             'sbsphase' + plane + 't_' + str(self._segment.label) + '.out')
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
                    constr_string += self._get_constraint_instruction(
                        self._name + '.dmu' + plane + name,
                        phase, error)
        return constr_string

    @Matcher.override(Matcher)
    def update_constraints_values(self):
        return ""

    @Matcher.override(Matcher)
    def update_variables_definition(self):
        update_vars_str = ""
        for variable in self.get_variables():
            update_vars_str += "        " + variable + ' := ' + self._name + "." + variable + '_0 + d' + variable + ';\n'
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
            apply_correction_str += variable + ' = ' + self._name + "." + variable + '_0 + d' + variable + ';\n'
        return apply_correction_str
