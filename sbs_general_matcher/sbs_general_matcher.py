from __future__ import print_function
import __init__  # @UnusedImport
import os
import sys
import json
import argparse
import logging
import log_handler
from matchers import (matcher,
                      phase_matcher,
                      coupling_matcher,
                      kmod_matcher,
                      amp_matcher)
from template_manager.template_processor import TemplateProcessor
from SegmentBySegment import SegmentBySegment
from madx import madx_templates_runner


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

MATCHER_TYPES = {
    "phase": phase_matcher.PhaseMatcher,
    "coupling": coupling_matcher.CouplingMatcher,
    "kmod": kmod_matcher.KmodMatcher,
    "amp": amp_matcher.AmpMatcher,
}


LOGGER = logging.getLogger(__name__)


def main(input_file_path):
    LOGGER.info("+++ Segment-by-segment general matcher +++")

    input_data = InputData.init_from_json_file(input_file_path)
    log_handler.add_file_handler(LOGGER, input_data.match_path)
    run_full_madx_matching(input_data)


def _run_madx_matching(input_data, just_twiss=False):
    _try_to_delete_twiss_cors(input_data)
    madx_templates = madx_templates_runner.MadxTemplates(
        log_file=os.path.join(input_data.match_path,
                              "match_madx_log.out"),
        output_file=os.path.join(input_data.match_path,
                                 "resolved_madx_match.madx")
    )
    template_processor = TemplateProcessor(input_data.matchers,
                                           input_data.match_path,
                                           input_data.lhc_mode,
                                           input_data.minimize,
                                           madx_templates)
    if not just_twiss:
        template_processor.run()
    else:
        template_processor.run_just_twiss()
    _check_twiss_cors(input_data)


def run_full_madx_matching(input_data):
    _run_madx_matching(input_data)
    _write_sbs_data_for_matchers(input_data)
    _build_changeparameters_file(input_data)


def run_twiss_and_sbs(input_data):
    _run_madx_matching(input_data, just_twiss=True)
    _write_sbs_data_for_matchers(input_data)
    _build_changeparameters_file(input_data)


def _manipulate_twiss_cors(input_data, function):
    for this_matcher in input_data.matchers:
        for beam in this_matcher.get_beams():
            matcher_path = this_matcher.get_match_data(
                beam).get_beam_match_path()
            matcher_path_sbs = os.path.join(matcher_path, "sbs")
            if os.path.isdir(matcher_path_sbs):
                twiss_cor_path = os.path.join(
                    matcher_path_sbs,
                    "twiss_IP" + str(this_matcher.get_ip()) + "_cor.dat"
                )
                function(twiss_cor_path)
                twiss_cor_back_path = os.path.join(
                    matcher_path_sbs,
                    "twiss_IP" + str(this_matcher.get_ip()) + "_cor_back.dat"
                )
                function(twiss_cor_back_path)


def _try_to_delete_twiss_cors(input_data):
    def function(twiss_cor_path):
        if os.path.isfile(twiss_cor_path):
            os.remove(twiss_cor_path)
    _manipulate_twiss_cors(input_data, function)


def _check_twiss_cors(input_data):
    def function(twiss_cor_path):
        if not os.path.isfile(twiss_cor_path):
            raise TwissFailedError(twiss_cor_path)
    _manipulate_twiss_cors(input_data, function)


class TwissFailedError(Exception):
    def __init__(self, twiss_path):
        super(TwissFailedError, self).__init__(
            "Twiss failed for twiss file: " + twiss_path
        )


def _write_sbs_data_for_matchers(input_data):
    for this_matcher in input_data.matchers:
        for beam in this_matcher.get_beams():
            _write_sbs_data(
                beam,
                str(this_matcher.get_ip()),
                this_matcher.get_match_data(beam).get_beam_match_path(),
                this_matcher.get_match_data(beam).get_range_start_name(),
            )


def _write_sbs_data(beam, ip, temporary_path, range_start_name):
    save_path = os.path.join(temporary_path, "sbs")
    input_data = SegmentBySegment._InputData(temporary_path)
    prop_models = SegmentBySegment._PropagatedModels(save_path, "IP" + str(ip))
    SegmentBySegment.getAndWriteData("IP" + ip, input_data, None, prop_models, save_path, False, False, True, False, "LHCB" + str(beam), None, None, None)


def _build_changeparameters_file(input_data):
    original_changeparameters_file = os.path.join(input_data.match_path, "changeparameters.madx")
    changeparameters_match_file = os.path.join(input_data.match_path, "changeparameters_match.madx")
    with open(original_changeparameters_file, "r") as original_changeparameters_data:
        with open(changeparameters_match_file, "w") as changeparameters_match_data:
            for original_line in original_changeparameters_data:
                parts = original_line.split("=")
                variable_name = parts[0].replace("d", "", 1).strip()
                variable_value = -float(parts[1].replace(";", "").strip())
                sign = " - " if variable_value < 0.0 else " + "
                changeparameters_match_data.write(
                    variable_name + " = " +
                    variable_name + sign + str(abs(variable_value)) + ";\n"
                )
            changeparameters_match_data.write("return;\n")


class InputData():

    @staticmethod
    def init_from_json_file(input_file_path):
        instance = InputData()
        with open(input_file_path, "r") as input_file:
            input_data = InputData._byteify(json.load(input_file))
            instance._check_and_assign_attribute(input_data, "lhc_mode")
            instance._check_and_assign_attribute(input_data, "match_path")
            instance.matchers = []
            if "minimize" in input_data:
                instance.minimize = input_data["minimize"]
            if "matchers" in input_data:
                instance._get_matchers_list(input_data)
        return instance

    @staticmethod
    def init_from_matchers_list(lhc_mode, match_path, minimize, matchers_list):
        instance = InputData()
        instance.lhc_mode = lhc_mode
        instance.match_path = match_path
        instance.matchers = matchers_list
        instance.minimize = minimize
        return instance

    def _get_matchers_list(self, input_data):
        raw_matchers_list = input_data["matchers"]
        for matcher_name, matcher_data in raw_matchers_list.iteritems():
            matcher.Matcher._check_attribute(matcher_name,
                                             matcher_data,
                                             "type")
            matcher_type = matcher_data["type"]
            MatcherClass = MATCHER_TYPES.get(matcher_type, None)
            if MatcherClass is None:
                raise ValueError('Unknown matcher type: ' + matcher_type +
                                 ' must be in: ' + str(MATCHER_TYPES.keys()))
            self.matchers.append(
                MatcherClass(matcher_name, matcher_data, self.match_path)
            )

    def _check_and_assign_attribute(self, input_data, attribute_name):
        matcher.Matcher._check_attribute("input data",
                                         input_data,
                                         attribute_name)
        setattr(self, attribute_name, input_data[attribute_name])

    # This transforms annoying unicode string into common byte string
    @staticmethod
    def _byteify(input_data):
        if isinstance(input_data, dict):
            return dict([(InputData._byteify(key), InputData._byteify(value))
                         for key, value in input_data.iteritems()])
        elif isinstance(input_data, list):
            return [InputData._byteify(element) for element in input_data]
        elif isinstance(input_data, unicode):
            return input_data.encode('utf-8')
        else:
            return input_data


def _run_gui(lhc_mode=None, match_path=None, input_dir=None):
    try:
        from gui import gui
    except ImportError as e:
        LOGGER.debug("ImportError importing GUI", exc_info=1)
        LOGGER.info("Cannot start GUI using the current Python installation:")
        LOGGER.info(str(e))
        LOGGER.info("Launching OMC Anaconda Python...")
        _run_gui_anaconda()
        return
    gui.main(lhc_mode, match_path, input_dir)


def _run_gui_anaconda():
    from subprocess import call
    if not sys.platform == "darwin":  # This is Mac
        if "win" in sys.platform:
            LOGGER.error("There is not Windows version of Anaconda in OMC.\
                         Aborting.")
            return
    interpreter = os.path.join("/afs", "cern.ch", "work", "o", "omc",
                               "anaconda", "bin", "python")
    command = sys.argv
    command.insert(0, interpreter)
    call(command)


def _parse_args():
    if len(sys.argv) >= 2:
        first_arg = sys.argv[1]
        if first_arg == "gui":
            LOGGER.info("Running GUI...")
            parser = argparse.ArgumentParser()
            parser.add_argument(
                "gui", help="Run GUI mode.",
                type=str,
            )
            parser.add_argument(
                "--lhc_mode", help="LHC mode.",
                type=str,
            )
            parser.add_argument(
                "--input", help="Beta-beating GUI output directory.",
                type=str,
            )
            parser.add_argument(
                "--match", help="Match ouput directory.",
                type=str,
            )
            args = parser.parse_args()
            _run_gui(args.lhc_mode, args.match, args.input)
        elif os.path.isfile(first_arg):
            LOGGER.info("Given input is a file, matching from JSON file...")
            main(os.path.abspath(first_arg))
    elif len(sys.argv) == 1:
        LOGGER.info("No given input, running GUI...")
        _run_gui()


if __name__ == "__main__":
    log_handler.set_up_console_logger(logging.getLogger(""))
    _parse_args()
