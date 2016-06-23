import os
from .phase_matcher import PhaseMatcher
from Python_Classes4MAD import metaclass

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
ALL_LISTS = os.path.join(CURRENT_PATH, '..', '..', 'MODEL', 'LHCB', 'fullresponse')

SET_INITIAL_VALUES_TEMPLATE = """
option, -echo, -info;
call, file = "%(MODIFIERS_OPTICS)s";
option, echo, info;
exec, save_initial_and_final_values(LHCB%(BEAM_NUM)s, %(START_FROM)s, %(END_AT)s, "%(PATH)s/measurement_from_amplitude_%(LABEL)s.madx", %(B_INI)s, %(B_END)s);
exec, twiss_segment(%(FRONT_SEQ)s, "%(PATH)s/twiss_%(LABEL)s.dat", %(B_INI)s);
exec, twiss_segment(%(BACK_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_back.dat", %(B_END)s);
"""


class BetaMatcher(PhaseMatcher):

    def define_aux_values(self):
        variables_s_str = ""
        variables_d_str = ""

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
            "PHASES": "",
            "S_VARIABLES": variables_s_str,
            "D_VARIABLES": variables_d_str,
        }

    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            calibrated_betas = metaclass.twiss(self.get_calibrated_betas_path(beam, plane))

            for index in range(0, len(calibrated_betas.NAME)):
                name = calibrated_betas.NAME[index]
                if name not in self.excluded_constraints_list:
                    calibrated_beta = calibrated_betas.BETX[index] if plane == "x" else calibrated_betas.BETY[index]
                    calibrated_beta_error = calibrated_betas.BETXSTD[index] if plane == "x" else calibrated_betas.BETYSTD[index]
                    s = calibrated_betas.S[index]

                    weight = self.get_constraint_weight(calibrated_beta, calibrated_beta_error, lambda value: abs(value) <= 0.25)

                    weight = 1.0
                    constr_string += '    constraint, weight = ' + str(weight) + ' , '
                    constr_string += 'expr =  table(twiss, ' + name + ', bet' + plane + ') = '
                    constr_string += str(calibrated_beta) + ';\n'

                    constr_string += '!   S = ' + str(s)
                    constr_string += ';\n'
        return constr_string

    def get_calibrated_betas_path(self, beam, plane):
        measurement_path = self.get_match_data_for_beam(beam).beam_match_path
        return os.path.join(measurement_path, "getampbeta" + plane + ".out")
