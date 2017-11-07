import os
from .matcher import Matcher
from Python_Classes4MAD import metaclass

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class PhaseMatcher(Matcher):

    @Matcher.override(Matcher)
    def get_variables(self, exclude=True):
        variables = self.segment.get_segment_vars(
            classes=self.var_classes,
        )
        if exclude:
            variables = [
                var for var in variables
                if var not in self.excluded_variables
            ]
        return variables

    PH_ERR_TMPL = (
        "{matcher_name}.dmu{plane}{name} := "
        "{sign}((table(twiss, {name}, mu{plane}) - "
        "table({nominal_table_name}, {name}, mu{plane})));\n"
    )

    SEGMENT_TWISS_TMPL = (
        "use, period={seq};\n"
        "twiss, beta0={init_vals}, chrom, table={table_name};\n"
    )

    @Matcher.override(Matcher)
    def define_aux_vars(self):
        phases_str = ""
        sign = "" if "f" in self.propagation else "-"
        for plane in ["x", "y"]:
            sbs_data_path = os.path.join(
                os.path.join(self.matcher_path, "sbs"),
                'sbsphase' + plane + 't_' + str(self.segment.label) + '.out'
            )
            sbs_data = metaclass.twiss(sbs_data_path)
            for name in sbs_data.NAME:
                phases_str += PhaseMatcher.PH_ERR_TMPL.format(
                    matcher_name=self.name, sign=sign, name=name, plane=plane,
                    nominal_table_name=self._get_nominal_table_name(),
                )

        variables_s_str = ""
        for variable in self.get_variables():
            variables_s_str += self.name + '.' + variable + '_0' + ' = ' + variable + ';\n'

        aux_vars_str = PhaseMatcher.SEGMENT_TWISS_TMPL.format(
            seq=self.get_sequence_name(),
            init_vals=self.get_initvals_name(),
            table_name=self._get_nominal_table_name(),
        )
        aux_vars_str += phases_str
        aux_vars_str += variables_s_str
        return aux_vars_str

    @Matcher.override(Matcher)
    def define_constraints(self):
        constr_string = ""
        for plane in ["x", "y"]:
            sbs_data = metaclass.twiss(
                os.path.join(self.matcher_path, "sbs",
                             'sbsphase' + plane + 't_' + str(self.segment.label) + '.out')
            )

            is_back = "b" in self.propagation
            for index in range(0, len(sbs_data.NAME)):
                name = sbs_data.NAME[index]
                if name not in self.excluded_constraints:
                    if is_back is not True:
                        phase = sbs_data.PROPPHASEX[index] if plane == "x" else sbs_data.PROPPHASEY[index]
                        error = sbs_data.ERRPROPPHASEX[index] if plane == "x" else sbs_data.ERRPROPPHASEY[index]
                    else:
                        phase = sbs_data.BACKPHASEX[index] if plane == "x" else sbs_data.BACKPHASEY[index]
                        error = sbs_data.ERRBACKPHASEX[index] if plane == "x" else sbs_data.ERRBACKPHASEY[index]
                    constr_string += self._get_constraint_instruction(
                        self.name + '.dmu' + plane + name,
                        phase, error)
        return constr_string

    @Matcher.override(Matcher)
    def update_constraints_values(self):
        return ""

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
