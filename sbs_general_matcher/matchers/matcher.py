import os
import logging
import shutil
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
LOGGER = logging.getLogger(__name__)


class Matcher(object):

    def __init__(self, matcher_name, matcher_dict, match_path):
        for attribute_name in ["ip", "use_errors", "propagation"]:
            Matcher._check_attribute(matcher_name, matcher_dict, attribute_name)

        self._name = matcher_name
        self._match_path = match_path
        self._ip = matcher_dict["ip"]
        self._use_errors = matcher_dict["use_errors"]
        self._front_or_back = matcher_dict["propagation"].lower()

        beams_paths = {}
        for beam in [1, 2]:
            try:
                Matcher._check_attribute(matcher_name, matcher_dict,
                                         "beam" + str(beam) + "_path")
                beams_paths[beam] = str(matcher_dict["beam" + str(beam) + "_path"])
            except InputError:
                pass
        if len(beams_paths) == 0:
            raise InputError("No beam1_path or beam2_path defined for matcher " + matcher_name)

        exclude_constr_string = ""
        if "exclude_constraints" in matcher_dict:
            exclude_constr_string = matcher_dict["exclude_constraints"]
        exclude_vars_string = ""
        if "exclude_variables" in matcher_dict:
            exclude_vars_string = matcher_dict["exclude_variables"]
        LOGGER.info("Successfully read matcher " + matcher_name)

        self._match_data = {}
        for beam in beams_paths:
            self._match_data[beam] = MatchData(matcher_name, self._ip,
                                               beams_paths[beam],
                                               match_path, beam)

        self._excluded_constraints_list = self._parse_exclude_string(
            exclude_constr_string
        )
        self._excluded_variables_list = self._parse_exclude_string(
            exclude_vars_string
        )

        assert self._front_or_back in ["front", "back", "f", "b"]
        self._front_or_back = self._front_or_back[0]
        self._ini_end = "ini" if self._front_or_back == "f" else "end"

    def define_aux_vars(self):
        """Returns the MAD-X string to define the auxiliary values to use
        during the matching"""
        raise NotImplementedError

    def get_all_variables(self):
        """Returns either a list of variable names or a MAD-X string defining
        the variables. If a list is given, this function can be used by other
        matchers to read the variable names and the variables won't be duplicated
        in the resulting MAD-X script."""
        return []

    def get_common_variables(self):
        """Returns either a list of variable names or a MAD-X string defining
        the variables. If a list is given, this function can be used by other
        matchers to read the variable names and the variables won't be duplicated
        in the resulting MAD-X script."""
        return []

    def get_variables_for_beam(self, beam):
        """Returns either a list of variable names or a MAD-X string defining
        the variables. If a list is given, this function can be used by other
        matchers to read the variable names and the variables won't be duplicated
        in the resulting MAD-X script."""
        return []

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

    def get_name(self):
        return self._name

    def get_match_path(self):
        return self._match_path

    def get_ip(self):
        return self._ip

    def get_beams(self):
        return self._match_data.keys()

    def get_match_data(self, beam):
        return self._match_data[beam]

    def get_front_or_back(self):
        return self._front_or_back

    def get_ini_end(self):
        return self._ini_end

    def set_exclude_variables(self, excluded_variables_list):
        self._excluded_variables_list = excluded_variables_list

    def set_disabled_constraints(self, disbled_constraints):
        self._excluded_constraints_list = disbled_constraints

    @staticmethod
    def _check_attribute(base_dict_name, base_dict, attribute_name):
        if attribute_name not in base_dict:
            raise InputError('Cannot find ' + attribute_name + ' attribute in ' + base_dict_name + '. Aborting.')

    def _get_constraint_instruction(self, constr_name, value, error, sigmas=1.):
        if self._use_errors:
            upper_bound = str(value + error * sigmas)
            lower_bound = str(value - error * sigmas)
            constr_string = '    constraint, weight = 1.0, '
            constr_string += 'expr =  ' + constr_name + ' < ' + upper_bound + ';\n'
            constr_string += '    constraint, weight = 1.0, '
            constr_string += 'expr =  ' + constr_name + ' > ' + lower_bound + ';\n'
        else:
            constr_string = '    constraint, weight = 1.0, '
            constr_string += 'expr =  ' + constr_name + ' = ' + str(value) + ';\n'
        return constr_string

    def _parse_exclude_string(self, exclude_string):
        if not exclude_string == "":
            exclude_list = [var_name.strip()
                            for var_name in exclude_string.strip('"').split(",")]
        else:
            exclude_list = []
        return exclude_list

    @staticmethod
    def override(parent_cls):
        """
        Decorator that checks if the decorated method overrides some method in the given parent class.
        Similar to java @override annotation.
        To use it:
        @Matcher.override(Superclass)
        def some_method():
            pass
        """
        def override_decorator(method):
            if not (method.__name__ in dir(parent_cls)):
                raise TypeError(method.__name__ + " must override a method from class " + parent_cls.__name__ + ".")
            else:
                method._doc__ = getattr(parent_cls, method.__name__).__doc__
            return method
        return override_decorator


class InputError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MatchData():

    def __init__(self, name, ip, measurement_data_path, match_path, beam):
        self._beam_match_path = os.path.join(match_path, "Beam" + str(beam) + "_" + name)
        self._beam_match_sbs_path = os.path.join(self._beam_match_path, "sbs")

        LOGGER.info("Copying measurement files for beam " + str(beam) +
                    " into match folder...")
        iotools.create_dirs(self._beam_match_sbs_path)
        self._copy_measurement_files(ip,
                                     measurement_data_path)

        LOGGER.info("Getting matching range for beam " + str(beam) + "...")
        ((self._range_start_s,
          self._range_start_name),
         (self._range_end_s,
          self._range_end_name)) = MatchData._get_match_bpm_range(
            os.path.join(self._beam_match_sbs_path, "sbsphasext_IP" + str(ip) + ".out")
        )
        LOGGER.info("Matching range for Beam " + str(beam) + ": " +
                    self._range_start_name + " " + self._range_end_name)

        self._modifiers = os.path.join(self._beam_match_sbs_path, "modifiers.madx")
        assert os.path.isfile(self._modifiers)

    def get_beam_match_path(self):
        return self._beam_match_path

    def get_beam_match_sbs_path(self):
        return self._beam_match_sbs_path

    def get_range_start_name(self):
        return self._range_start_name

    def get_range_start_s(self):
        return self._range_start_s

    def get_range_end_name(self):
        return self._range_end_name

    def get_range_end_s(self):
        return self._range_end_s

    def get_modifiers(self):
        return self._modifiers

    def _copy_measurement_files(self, ip, measurement_path):
        # GetLLM output files:
        MatchData._copy_files_with_extension(measurement_path,
                                             self._beam_match_path, ".out")
        # SbS output files for the given IP:
        MatchData._copy_files_which_contains(os.path.join(measurement_path, "sbs"),
                                             self._beam_match_sbs_path,
                                             "IP" + str(ip))
        # SbS MAD-X files (not necessary but useful):
        MatchData._copy_files_with_extension(os.path.join(measurement_path, "sbs"),
                                             self._beam_match_sbs_path,
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
