import os
import sys
import shutil
import json
from Python_Classes4MAD import metaclass
from Utilities import iotools

EXTRACT_SEQUENCES_TEMPLATE = """
exec, extract_segment_sequence(LHCB%(BEAM_NUM)s, %(FRONT_SEQ)s, %(BACK_SEQ)s, %(START_FROM)s, %(END_AT)s);
exec, beam_LHCB%(BEAM_NUM)s(%(FRONT_SEQ)s);
exec, beam_LHCB%(BEAM_NUM)s(%(BACK_SEQ)s);
"""

SET_INITIAL_VALUES_TEMPLATE = """
option, -echo, -info;
call, file = "%(MODIFIERS_OPTICS)s";
option, echo, info;
exec, save_initial_and_final_values(LHCB%(BEAM_NUM)s, %(START_FROM)s, %(END_AT)s, "%(PATH)s/measurement_%(LABEL)s.madx", %(B_INI)s, %(B_END)s);
exec, twiss_segment(%(FRONT_SEQ)s, "%(PATH)s/twiss_%(LABEL)s.dat", %(B_INI)s);
exec, twiss_segment(%(BACK_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_back.dat", %(B_END)s);
"""

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

RUN_CORRECTED_TWISS_TEMPLATE = """
option, -echo, -info;
call, file = "%(MODIFIERS_OPTICS)s";
%(APPLY_CORRECTIONS)s
option, echo, info;
exec, twiss_segment(%(FRONT_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_cor.dat", %(B_INI)s);
exec, twiss_segment(%(BACK_SEQ)s, "%(PATH)s/twiss_%(LABEL)s_cor_back.dat", %(B_END)s);
"""

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class Matcher(object):

    ERROR_CONSTRAINT_FACTOR = 1000
    MAX_WEIGHT = 4
    MAX_REL_DELTA_LIMIT = 0.1
    ALL_LISTS = os.path.join(CURRENT_PATH, '..', '..', 'MODEL', 'LHCB', 'fullresponse')

    def __init__(self, name, ip,
                 measurement_data_b1_path, measurement_data_b2_path,
                 match_path,
                 use_errors, front_or_back,
                 exclude_constr_string, exclude_vars_string, all_lists=None):
        self.name = name
        self.match_data_b1 = MatchData(name, ip, measurement_data_b1_path, match_path, 1)
        self.match_data_b2 = MatchData(name, ip, measurement_data_b2_path, match_path, 2)

        self.ip = ip
        self.use_errors = use_errors

        self.excluded_constraints_list = self._parse_exclude_string(
            exclude_constr_string
        )
        self.excluded_variables_list = self._parse_exclude_string(
            exclude_vars_string
        )

        self.front_or_back = front_or_back.lower()
        assert self.front_or_back in ["front", "back", "f", "b"]
        self.front_or_back = self.front_or_back[0]
        self.ini_end = "ini" if self.front_or_back == "f" else "end"

        if all_lists is None:
            all_lists = Matcher.ALL_LISTS
        self.variables_beam1 = json.load(
            file(os.path.join(all_lists, "LHCB1", "AllLists.json"), 'r')
        )['getListsByIR'][1]
        self.variables_common, self.variables_beam2 = json.load(
            file(os.path.join(all_lists, "LHCB2", "AllLists.json"), 'r')
        )['getListsByIR']

    def define_aux_values(self):
        """Returns the MAD-X string to define the auxiliary values to use
        during the matching"""
        raise NotImplementedError

    def define_variables(self):
        """Returns the MAD-X string to define the matching variables for this
        matcher"""
        raise NotImplementedError

    def define_constraints(self, beam):
        """Returns two MAD-X strings to define the matching constraints for this
        matcher for beam 1 and 2"""
        raise NotImplementedError

    def update_constraints_values(self, beam):
        """Returns the MAD-X string that updates the value of the constraints to
        let MAD-X reevaluate in every iteration."""
        raise NotImplementedError

    def update_variables_definition(self, beam):
        """Returns the MAD-X string that updates the definition of the variables
        that may have been override by the modifiers file."""
        raise NotImplementedError

    def generate_changeparameters(self):
        """Returns the MAD-X string that selects the variables to dump to the
        changeparameters file."""
        raise NotImplementedError

    def apply_correction(self):
        """Returns the MAD-X string that applies the final correction to the
        variables, in order to get the corrected twiss files"""
        raise NotImplementedError

    def extract_sequences(self):
        extract_sequences_str = ""
        for beam in [1, 2]:
            extract_sequences_str += EXTRACT_SEQUENCES_TEMPLATE % {
                "BEAM_NUM": str(beam),
                "FRONT_SEQ": "lhcb" + str(beam) + "_f_" + self.name,
                "BACK_SEQ": "lhcb" + str(beam) + "_b_" + self.name,
                "START_FROM": getattr(self, "match_data_b" + str(beam)).range_start_name,
                "END_AT": getattr(self, "match_data_b" + str(beam)).range_end_name,
            }
        return extract_sequences_str

    def set_initial_values(self):
        set_initial_values_str = ""
        for beam in [1, 2]:
            set_initial_values_str += SET_INITIAL_VALUES_TEMPLATE % {
                "MODIFIERS_OPTICS": getattr(self, "match_data_b" + str(beam)).modifiers,
                "BEAM_NUM": str(beam),
                "PATH": getattr(self, "match_data_b" + str(beam)).beam_match_sbs_path,
                "LABEL": "IP" + str(self.ip),
                "START_FROM": getattr(self, "match_data_b" + str(beam)).range_start_name,
                "END_AT": getattr(self, "match_data_b" + str(beam)).range_end_name,
                "FRONT_SEQ": "lhcb" + str(beam) + "_f_" + self.name,
                "BACK_SEQ": "lhcb" + str(beam) + "_b_" + self.name,
                "B_INI": "b" + str(beam) + "_ini_" + self.name,
                "B_END": "b" + str(beam) + "_end_" + self.name
            }
        return set_initial_values_str

    def set_matching_macros(self):
        matching_macros = ""
        for beam in [1, 2]:
            matching_macros += MATCHING_MACRO_TEMPLATE % {
                "MODIFIERS_OPTICS": getattr(self, "match_data_b" + str(beam)).modifiers,
                "MACRO_NAME": "macro_" + self.name + "_b" + str(beam),
                "BEAM_NUM": str(beam),
                "SEQ": "lhcb" + str(beam) + "_" + self.front_or_back + "_" + self.name,
                "UPDATE_CONSTRAINTS": self.update_constraints_values(beam),
                "UPDATE_VARIABLES": self.update_variables_definition(),
                "B_INI_END": "b" + str(beam) + "_" + self.ini_end + "_" + self.name,
                "DEFINE_CONSTRAINTS": self.define_constraints(beam)
            }
            matching_macros += "\n"
        return matching_macros

    def run_corrected_twiss(self):
        run_corrected_twiss_str = ""
        for beam in [1, 2]:
            run_corrected_twiss_str += RUN_CORRECTED_TWISS_TEMPLATE % {
                "MODIFIERS_OPTICS": getattr(self, "match_data_b" + str(beam)).modifiers,
                "APPLY_CORRECTIONS": self.apply_correction(),
                "BEAM_NUM": str(beam),
                "PATH": getattr(self, "match_data_b" + str(beam)).beam_match_sbs_path,
                "LABEL": "IP" + str(self.ip),
                "FRONT_SEQ": "lhcb" + str(beam) + "_f_" + self.name,
                "BACK_SEQ": "lhcb" + str(beam) + "_b_" + self.name,
                "B_INI": "b" + str(beam) + "_ini_" + self.name,
                "B_END": "b" + str(beam) + "_end_" + self.name
            }
        return run_corrected_twiss_str

    @classmethod
    def from_matcher_dict(cls, matcher_name, matcher_dict, match_path):
        for attribute_name in ["type", "ip", "beam1_path", "beam2_path",
                               "use_errors", "propagation"]:
            Matcher._check_attribute(matcher_name, matcher_dict, attribute_name)

        exclude_constr_string = ""
        if "exclude_constraints" in matcher_dict:
            exclude_constr_string = matcher_dict["exclude_constraints"]
        exclude_vars_string = ""
        if "exclude_variables" in matcher_dict:
            exclude_vars_string = matcher_dict["exclude_variables"]
        all_lists = None
        if "all_lists" in matcher_dict:
            all_lists = matcher_dict["all_lists"]
        print "Successfully read matcher " + matcher_name
        instance = Matcher(
            matcher_name, matcher_dict["ip"],
            str(matcher_dict["beam1_path"]), str(matcher_dict["beam2_path"]),
            str(match_path),
            matcher_dict["use_errors"], matcher_dict["propagation"],
            exclude_constr_string, exclude_vars_string,
            all_lists=all_lists
        )
        instance.__class__ = cls
        return instance

    @staticmethod
    def _check_attribute(base_dict_name, base_dict, attribute_name):
        if attribute_name not in base_dict:
            print >> sys.stderr, 'Cannot find ' + attribute_name + ' attribute in ' + base_dict_name + '. Aborting.'
            sys.exit(-1)

    def _parse_exclude_string(self, exclude_string):
        if not exclude_string == "":
            exclude_list = [var_name.strip()
                            for var_name in exclude_string.strip('"').split(",")]
        else:
            exclude_list = []
        return exclude_list

    def get_match_data_for_beam(self, beam):
        if beam == 1:
            return self.match_data_b1
        elif beam == 2:
            return self.match_data_b2
        else:
            return None

    def get_constraint_weight(self, value, error, filter_function):
        if not filter_function(value):
            weight = 1e-6
        elif not self.use_errors:
            weight = 1.
        else:
            weight = 1. / ((Matcher.ERROR_CONSTRAINT_FACTOR * error) ** 2)
            if weight > Matcher.MAX_WEIGHT:
                weight = Matcher.MAX_WEIGHT


class MatchData():

    def __init__(self, name, ip, measurement_data_path, match_path, beam):
        self.beam_match_path = os.path.join(match_path, "Beam" + str(beam) + "_" + name)
        self.beam_match_sbs_path = os.path.join(self.beam_match_path, "sbs")

        print "Copying measurement files for beam", str(beam),\
              "into match folder..."
        iotools.create_dirs(self.beam_match_sbs_path)
        self._copy_measurement_files(ip,
                                     measurement_data_path)

        print "Getting matching range for beam" + str(beam) + "..."
        ((self.range_start_s,
          self.range_start_name),
        (self.range_end_s,
         self.range_end_name)) = MatchData._get_match_bpm_range(
            os.path.join(self.beam_match_sbs_path, "sbsphasext_IP" + str(ip) + ".out")
        )
        print "Matching range for Beam " + str(beam),\
              self.range_start_name, self.range_end_name

        self.modifiers = os.path.join(self.beam_match_sbs_path, "modifiers.madx")
        assert os.path.isfile(self.modifiers)

    def _copy_measurement_files(self, ip, measurement_path):
        # GetLLM output files:
        MatchData._copy_files_with_extension(measurement_path,
                                             self.beam_match_path, ".out")
        # SbS output files for the given IP:
        MatchData._copy_files_which_contains(os.path.join(measurement_path, "sbs"),
                                             self.beam_match_sbs_path,
                                             "IP" + str(ip))
        # SbS MAD-X files (not necessary but useful):
        MatchData._copy_files_with_extension(os.path.join(measurement_path, "sbs"),
                                             self.beam_match_sbs_path,
                                             ".madx")

    @staticmethod
    def _copy_files_with_extension(src, dest, ext):
        MatchData._copy_files_with_filter(src, dest,
                                          lambda file_name: file_name.endswith(ext))

    @staticmethod
    def _copy_files_which_contains(src, dest, substring):
        MatchData._copy_files_with_filter(src, dest,
                                          lambda file_name: substring in file_name)

    @staticmethod
    def _copy_files_with_filter(src, dest, filter_function):
        src_files = MatchData._get_filtered_file_list(src, filter_function)
        for file_name in src_files:
            full_file_name = os.path.join(src, file_name)
            shutil.copy(full_file_name, dest)

    @staticmethod
    def _get_filtered_file_list(src, filter_function):
        filtered_file_list = []
        original_file_list = os.listdir(src)
        for file_name in original_file_list:
            if os.path.isfile(os.path.join(src, file_name)) and filter_function(file_name):
                filtered_file_list.append(file_name)
        return filtered_file_list

    @staticmethod
    def _get_match_bpm_range(file_path):
        twiss_data = metaclass.twiss(file_path)
        bpms_with_distances_list = zip(twiss_data.S, twiss_data.NAME)
        bpms_with_distances_list.sort()
        return bpms_with_distances_list[0], bpms_with_distances_list[-1]
