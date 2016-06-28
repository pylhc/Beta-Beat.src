from .matcher import Matcher
from .phase_matcher import PhaseMatcher
from Python_Classes4MAD import metaclass


class IpMatcher(PhaseMatcher):

    @Matcher.override(PhaseMatcher)
    def __init__(self, matcher_name, matcher_dict, match_path):
        super(IpMatcher, self).__init__(matcher_name, matcher_dict, match_path)
        IpMatcher._check_attribute(matcher_name, matcher_dict, "ip_data_b1")
        IpMatcher._check_attribute(matcher_name, matcher_dict, "ip_data_b2")
        self.ip_data_b1 = metaclass.twiss(matcher_dict["ip_data_b1"])
        self.ip_data_b2 = metaclass.twiss(matcher_dict["ip_data_b2"])

    @Matcher.override(PhaseMatcher)
    def define_aux_vars(self):
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

    @Matcher.override(PhaseMatcher)
    def define_constraints(self, beam):
        constr_string = ""
        for plane in ["x", "y"]:
            ip_data = getattr(self, "ip_data_b" + str(beam))
            ip_index = ip_data.indx["IP" + str(self.ip)]

            measured_beta = getattr(ip_data, "BET" + plane.upper())[ip_index]
            measured_beta_err = getattr(ip_data, "ERRBET" + plane.upper())[ip_index]
            beta_weight = self.get_constraint_weight(measured_beta, measured_beta_err, lambda x: True)
            beta_weight = 1.0
            constr_string += '    constraint, weight = ' + str(beta_weight) + ', '
            constr_string += 'expr =  table(twiss, IP' + str(self.ip) + ', bet' + plane + ') = '
            constr_string += str(measured_beta) + ';\n'

            measured_alpha = getattr(ip_data, "ALF" + plane.upper())[ip_index]
            measured_alpha_err = getattr(ip_data, "ERRALF" + plane.upper())[ip_index]
            alpha_weight = self.get_constraint_weight(measured_alpha, measured_alpha_err, lambda x: True)
            alpha_weight = 1.0
            constr_string += '    constraint, weight = ' + str(alpha_weight) + ', '
            constr_string += 'expr =  table(twiss, IP' + str(self.ip) + ', alf' + plane + ') = '
            constr_string += str(measured_alpha) + ';\n'
        return constr_string
