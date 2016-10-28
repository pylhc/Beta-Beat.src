import os


class TemplateProcessor(object):

    DEFAULT_TEMPLATE_PATH = os.path.join(
        "..", "..", "madx", "templates", "lhc_super_matcher.madx"
    )

    START_MATCH = """
match, use_macro;
"""
    END_MATCH = """
lmdif, tolerance:=1e-24, calls:=1200;
endmatch;
"""
    SAVE_CHANGEPARAMETERS = """
save, file="%(MATCH_PATH)s/changeparameters.madx";
"""
    CALL_CHANGEPARAMETERS = """
call, file="%(MATCH_PATH)s/changeparameters.madx";
"""

    def __init__(self, matchers_list, match_path, lhc_mode, minimize,
                 madx_templates_runner):
        self._matchers_list = matchers_list
        self._match_path = match_path
        self._lhc_mode = lhc_mode
        self._minimize = minimize
        self._madx_templates_runner = madx_templates_runner
        self._set_up_collections()

    def _set_up_collections(self):
        self._variables = set()
        self._extract_sequences_list = []
        self._set_initial_values_list = []
        self._aux_var_definition_list = []
        self._define_variables_list = []
        self._set_matching_macros_list = []
        self._gen_changeparameters_list = []
        self._run_corrected_twiss_list = []

    def run(self):
        self._process_matchers()
        self._madx_templates_runner.lhc_super_matcher_madx(
            self._lhc_mode,
            "\n".join(self._extract_sequences_list),
            "\n".join(self._set_initial_values_list),
            "\n".join(self._aux_var_definition_list),
            TemplateProcessor.START_MATCH,
            "\n".join(self._define_variables_list),
            "\n".join(self._set_matching_macros_list),
            TemplateProcessor.END_MATCH,
            "\n".join(self._gen_changeparameters_list),
            TemplateProcessor.SAVE_CHANGEPARAMETERS % {"MATCH_PATH": self._match_path},
            "\n".join(self._run_corrected_twiss_list),
        )

    def run_just_twiss(self):
        for matcher in self._matchers_list:
            self._extract_sequences(matcher)
            self._set_initial_values(matcher)
            self._define_aux_vars(matcher)
            self._run_corrected_twiss(matcher)
        self._madx_templates_runner.lhc_super_matcher_madx(
            self._lhc_mode,
            "\n".join(self._extract_sequences_list),
            "\n".join(self._set_initial_values_list),
            "\n".join(self._aux_var_definition_list),
            "",
            "\n".join(self._define_variables_list),
            "\n".join(self._set_matching_macros_list),
            "",
            "\n".join(self._gen_changeparameters_list),
            TemplateProcessor.CALL_CHANGEPARAMETERS % {"MATCH_PATH": self._match_path},
            "\n".join(self._run_corrected_twiss_list),
        )

    def _process_matchers(self):
        for matcher in self._matchers_list:
            self._extract_sequences(matcher)
            self._set_initial_values(matcher)
            self._define_aux_vars(matcher)
            self._collect_variables(matcher)
            self._set_matching_macros(matcher)
            self._generate_changeparameters(matcher)
            self._run_corrected_twiss(matcher)
        self._define_variables()

    EXTRACT_SEQUENCES_TEMPLATE = """
exec, extract_segment_sequence(LHCB%(BEAM_NUM)s, %(FRONT_SEQ)s, %(BACK_SEQ)s, %(START_FROM)s, %(END_AT)s);
exec, beam_LHCB%(BEAM_NUM)s(%(FRONT_SEQ)s);
exec, beam_LHCB%(BEAM_NUM)s(%(BACK_SEQ)s);
"""

    def _extract_sequences(self, matcher):
        extract_sequences_str = ""
        for beam in matcher.get_beams():
            extract_sequences_str += TemplateProcessor.EXTRACT_SEQUENCES_TEMPLATE % {
                "BEAM_NUM": str(beam),
                "FRONT_SEQ": "lhcb" + str(beam) + "_f_" + matcher.get_name(),
                "BACK_SEQ": "lhcb" + str(beam) + "_b_" + matcher.get_name(),
                "START_FROM": matcher.get_match_data(beam).get_range_start_name(),
                "END_AT": matcher.get_match_data(beam).get_range_end_name(),
            }
        self._extract_sequences_list.append(extract_sequences_str)

    SET_INITIAL_VALUES_TEMPLATE = """
option, -echo, -info;
call, file = "%(MODIFIERS_OPTICS)s";
option, echo, info;
exec, save_initial_and_final_values(
    LHCB%(BEAM_NUM)s,
    %(START_FROM)s,
    %(END_AT)s,
    "%(PATH)s/measurement_%(LABEL)s.madx",
    %(B_INI)s,
    %(B_END)s
);
exec, twiss_segment(%(FRONT_SEQ)s, "%(PATH)s/twiss_%(LABEL)s.dat", %(B_INI)s);
exec, twiss_segment(%(BACK_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_back.dat", %(B_END)s);
"""

    def _set_initial_values(self, matcher):
        set_initial_values_str = ""
        for beam in matcher.get_beams():
            set_initial_values_str += TemplateProcessor.SET_INITIAL_VALUES_TEMPLATE % {
                "MODIFIERS_OPTICS": matcher.get_match_data(beam).get_modifiers(),
                "BEAM_NUM": str(beam),
                "PATH": matcher.get_match_data(beam).get_beam_match_sbs_path(),
                "LABEL": "IP" + str(matcher.get_ip()),
                "START_FROM": matcher.get_match_data(beam).get_range_start_name(),
                "END_AT": matcher.get_match_data(beam).get_range_end_name(),
                "FRONT_SEQ": "lhcb" + str(beam) + "_f_" + matcher.get_name(),
                "BACK_SEQ": "lhcb" + str(beam) + "_b_" + matcher.get_name(),
                "B_INI": "b" + str(beam) + "_ini_" + matcher.get_name(),
                "B_END": "b" + str(beam) + "_end_" + matcher.get_name()
            }
        self._set_initial_values_list.append(set_initial_values_str)

    def _define_aux_vars(self, matcher):
        self._aux_var_definition_list.append(matcher.define_aux_vars())

    def _collect_variables(self, matcher):
        if type(matcher.get_all_variables()) is list:
            for variable in matcher.get_all_variables():
                self._variables.add(variable)
        elif type(matcher.get_all_variables()) is str:
            self._define_variables_list.append(matcher.get_all_variables())

    def _define_variables(self):
        def_variables_string = ""
        for variable in self._variables:
            def_variables_string += '    vary, name=d' + variable
            def_variables_string += ', step := 1e-5;\n'
        self._define_variables_list.append(def_variables_string)

    MATCHING_MACRO_TEMPLATE = """
%(MACRO_NAME)s: macro = {
use, period=%(SEQ)s;

option, -echo, -info;
call, file = "%(MODIFIERS_OPTICS)s";
%(UPDATE_VARIABLES)s

twiss, beta0=%(B_INI_END)s, chrom;
%(UPDATE_CONSTRAINTS)s
option, echo, info;
print, text = "=========== step for %(MACRO_NAME)s ===========";
};
%(DEFINE_CONSTRAINTS)s
"""

    def _set_matching_macros(self, matcher):
        matching_macros = ""
        for beam in matcher.get_beams():
            define_constr_str = matcher.define_constraints(beam)
            if self._minimize:
                define_constr_str += self._minimize_variables(matcher)
            matching_macros += TemplateProcessor.MATCHING_MACRO_TEMPLATE % {
                "MODIFIERS_OPTICS": matcher.get_match_data(beam).get_modifiers(),
                "MACRO_NAME": "macro_" + matcher.get_name() + "_b" + str(beam),
                "BEAM_NUM": str(beam),
                "SEQ": "lhcb" + str(beam) + "_" + matcher.get_front_or_back() + "_" + matcher.get_name(),
                "UPDATE_CONSTRAINTS": matcher.update_constraints_values(beam),
                "UPDATE_VARIABLES": matcher.update_variables_definition(),
                "B_INI_END": "b" + str(beam) + "_" + matcher.get_ini_end() + "_" + matcher.get_name(),
                "DEFINE_CONSTRAINTS": define_constr_str
            }
            matching_macros += "\n"
        self._set_matching_macros_list.append(matching_macros)

    def _minimize_variables(self, matcher):
        variables = matcher.get_all_variables()
        minimize_vars_str = ""
        minimize_vars_str += '    constraint, weight = 1.0, expr = sqrt('
        for variable in variables:
            minimize_vars_str += variable + "^2 +"
        minimize_vars_str = minimize_vars_str[:-1]
        minimize_vars_str += ') = 0.0; \n'
        return minimize_vars_str

    def _generate_changeparameters(self, matcher):
        changeparameters_str = ""
        for variable in matcher.get_all_variables():
            changeparameters_str += 'select,flag=save,pattern=\"d' + variable + '\";\n'
        self._gen_changeparameters_list.append(changeparameters_str)

    RUN_CORRECTED_TWISS_TEMPLATE = """
option, -echo, -info;
call, file = "%(MODIFIERS_OPTICS)s";
%(APPLY_CORRECTIONS)s
option, echo, info;
exec, twiss_segment(%(FRONT_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_cor.dat", %(B_INI)s);
exec, twiss_segment(%(BACK_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_cor_back.dat", %(B_END)s);
"""

    def _run_corrected_twiss(self, matcher):
        run_corrected_twiss_str = ""

        apply_correction_str = ""
        for variable in matcher.get_all_variables():
            apply_correction_str += variable + ' = ' + matcher.get_name() + "." + variable + '_0 + d' + variable + ';\n'

        for beam in matcher.get_beams():
            run_corrected_twiss_str += TemplateProcessor.RUN_CORRECTED_TWISS_TEMPLATE % {
                "MODIFIERS_OPTICS": matcher.get_match_data(beam).get_modifiers(),
                "APPLY_CORRECTIONS": apply_correction_str,
                "BEAM_NUM": str(beam),
                "PATH": matcher.get_match_data(beam).get_beam_match_sbs_path(),
                "LABEL": "IP" + str(matcher.get_ip()),
                "FRONT_SEQ": "lhcb" + str(beam) + "_f_" + matcher.get_name(),
                "BACK_SEQ": "lhcb" + str(beam) + "_b_" + matcher.get_name(),
                "B_INI": "b" + str(beam) + "_ini_" + matcher.get_name(),
                "B_END": "b" + str(beam) + "_end_" + matcher.get_name()
            }
        self._run_corrected_twiss_list.append(run_corrected_twiss_str)
