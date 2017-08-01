import os
import logging
import shutil
from Python_Classes4MAD import metaclass
from Utilities import iotools
from model import manager

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

    def __init__(self, lhc_mode, beam, matcher_name, matcher_dict,
                 matcher_variables, match_path):
        self._name = matcher_name
        self._main_match_path = match_path
        self._label = matcher_dict["label"]
        self._use_errors = matcher_dict["use_errors"]
        self._front_or_back = matcher_dict["propagation"].lower()
        self._variables_cls = matcher_variables
        self._matcher_path = os.path.join(
            match_path,
            self._name
        )
        measurement_path = matcher_dict["path"]
        match_path = os.path.join(
            match_path,
            self._name
        )
        Matcher._copy_measurement_files(
            self._label, measurement_path, match_path
        )
        self._segment = Matcher._get_segment(
            lhc_mode, beam, self._matcher_path, self._label
        )

        self._excluded_constraints_list = []
        self._excluded_variables_list = []
        if "exclude_constraints" in matcher_dict:
            self._excluded_constraints_list = (
                matcher_dict["exclude_constraints"].strip('"').split(",")
            )
        if "exclude_variables" in matcher_dict:
            self._excluded_variables_list = (
                matcher_dict["exclude_variables"].strip('"').split(",")
            )

        assert self._front_or_back in ["front", "back", "f", "b"]
        self._front_or_back = self._front_or_back[0]
        self._ini_end = "ini" if self._front_or_back == "f" else "end"

        LOGGER.info("Successfully read matcher " + matcher_name)

    def get_name(self):
        return self._name

    def get_main_match_path(self):
        return self._main_match_path

    def get_matcher_path(self):
        return self._matcher_path

    def get_segment(self):
        return self._segment

    def get_front_or_back(self):
        return self._front_or_back

    def get_variables(self, exclude=True):
        """
        Returns the variables names to use to match. If exclude is true,
        it will not return the variables in the excluded variables list.
        """
        raise NotImplementedError

    def set_exclude_variables(self, excluded_variables_list):
        self._excluded_variables_list = excluded_variables_list

    def set_disabled_constraints(self, disbled_constraints):
        self._excluded_constraints_list = disbled_constraints

    def define_aux_vars(self):
        """Returns the MAD-X string to define the auxiliary values to use
        during the matching"""
        raise NotImplementedError

    def define_constraints(self):
        """Returns two MAD-X strings to define the matching constraints for this
        matcher."""
        raise NotImplementedError

    def update_constraints_values(self):
        """Returns the MAD-X string that updates the value of the constraints to
        let MAD-X reevaluate in every iteration."""
        raise NotImplementedError

    def update_variables_definition(self):
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

    def _get_constraint_instruction(self, constr_name,
                                    value, error, sigmas=1.):
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

    @staticmethod
    def override(parent_cls):
        """
        Decorator that checks if the decorated method overrides some method
        in the given parent class.
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
                method.__doc__ = getattr(parent_cls, method.__name__).__doc__
            return method
        return override_decorator

    @staticmethod
    def _get_segment(lhc_mode, beam, match_path, label):
        LOGGER.info("Getting matching range for beam " + str(beam) + "...")
        (_, range_start), (_, range_end) = _get_match_bpm_range(
            os.path.join(os.path.join(match_path, "sbs"),
                         "sbsphasext_" + label + ".out")
        )
        LOGGER.info("Matching range for Beam " + str(beam) + ": " +
                    range_start + " " + range_end)
        accel_cls = manager.get_accel_class(
            "lhc", lhc_mode=lhc_mode, beam=beam
        )()
        optics_file = os.path.join(
            os.path.join(match_path, "sbs"), "modifiers.madx"
        )
        segment = accel_cls.get_segment(label,
                                        range_start,
                                        range_end,
                                        optics_file)
        return segment

    @staticmethod
    def _copy_measurement_files(label, measurement_path, match_math):
        iotools.create_dirs(match_math)
        iotools.create_dirs(os.path.join(match_math, "sbs"))
        # GetLLM output files:
        _copy_files_with_extension(measurement_path,
                                   match_math, ".out")
        # SbS output files for the given label:
        _copy_files_which_contains(
            os.path.join(measurement_path, "sbs"),
            os.path.join(match_math, "sbs"),
            label
        )
        # SbS MAD-X files (not necessary but useful):
        _copy_files_with_extension(
            os.path.join(measurement_path, "sbs"),
            os.path.join(match_math, "sbs"),
            ".madx"
        )


class InputError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def _copy_files_with_extension(src, dest, ext):
    _copy_files_with_filter(
        src, dest,
        lambda file_name: file_name.endswith(ext)
    )


def _copy_files_which_contains(src, dest, substring):
    _copy_files_with_filter(
        src, dest,
        lambda file_name: substring in file_name
    )


def _copy_files_with_filter(src, dest, filter_function):
    src_files = _get_filtered_file_list(src, filter_function)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        shutil.copy(full_file_name, dest)


def _get_filtered_file_list(src, filter_function):
    filtered_file_list = []
    original_file_list = os.listdir(src)
    for file_name in original_file_list:
        if (os.path.isfile(os.path.join(src, file_name)) and
                filter_function(file_name)):
            filtered_file_list.append(file_name)
    return filtered_file_list


def _get_match_bpm_range(file_path):
    twiss_data = metaclass.twiss(file_path)
    bpms_with_distances_list = zip(twiss_data.S, twiss_data.NAME)
    bpms_with_distances_list.sort()
    return bpms_with_distances_list[0], bpms_with_distances_list[-1]
