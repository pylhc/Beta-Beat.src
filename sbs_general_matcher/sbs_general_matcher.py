import __init__  # @UnusedImport
import os
import sys
import json
from matchers import matcher, phase_matcher, coupling_matcher, kmod_matcher, beta_matcher
from template_manager.template_processor import TemplateProcessor
from SegmentBySegment import SegmentBySegment
from madx import madx_templates_runner


CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

MATCHER_TYPES = {
    "phase": phase_matcher.PhaseMatcher,
    "coupling": coupling_matcher.CouplingMatcher,
    "kmod": kmod_matcher.KmodMatcher,
    "beta": beta_matcher.BetaMatcher
}


def main(input_file_path):

    print "+++ Segment-by-segment general matcher +++"

    print "Preparing MAD-X script..."
    input_data = InputData(input_file_path)
    _run_madx_matching(input_data)
    for matcher in input_data.matchers:
        for beam in matcher.get_beams():
            _write_sbs_data(beam, str(matcher.get_ip()),
                            matcher.get_match_data(beam).get_beam_match_path(),
                            matcher.get_match_data(beam).get_range_start_name(),
                            )
    _build_changeparameters_file(input_data)


def _run_madx_matching(input_data):
    madx_templates = madx_templates_runner.MadxTemplates(
        log_file=os.path.join(input_data.match_path, "match_madx_log.out"),
        output_file=os.path.join(input_data.match_path, "resolved_madx_match.madx")
    )
    template_processor = TemplateProcessor(input_data.matchers,
                                           input_data.match_path,
                                           input_data.lhc_mode,
                                           madx_templates)
    template_processor.run()


def _write_sbs_data(beam, ip, temporary_path, range_start_name):
    save_path = os.path.join(temporary_path, "sbs")
    input_data = SegmentBySegment._InputData(temporary_path)
    prop_models = SegmentBySegment._PropagatedModels(save_path, "IP" + str(ip))
    SegmentBySegment.getAndWriteData("IP" + ip, input_data, None, prop_models, save_path, False, False, True, False, "LHCB" + str(beam), None, None, None)


def _build_changeparameters_file(input_data):
    original_changeparameters_file = open(os.path.join(input_data.match_path, "changeparameters.madx"), "r")
    changeparameters_match_file = open(os.path.join(input_data.match_path, "changeparameters_match.madx"), "w")
    for original_line in original_changeparameters_file.readlines():
        parts = original_line.split("=")
        variable_name = parts[0].replace("d", "", 1).strip()
        variable_value = -float(parts[1].replace(";", "").strip())
        if variable_value < 0.0:
            print >> changeparameters_match_file, variable_name, " = ", variable_name, " - ", abs(variable_value), ";"
        else:
            print >> changeparameters_match_file, variable_name, " = ", variable_name, " + ", abs(variable_value), ";"
    print >> changeparameters_match_file, "return;"


class InputData():
    def __init__(self, input_file_path):
        with open(input_file_path, "r") as input_file:
            input_data = InputData._byteify(json.load(input_file))
            self._check_and_assign_attribute(input_data, "lhc_mode")
            self._check_and_assign_attribute(input_data, "match_path")
            self.matchers = []
            if "matchers" in input_data:
                self._get_matchers_list(input_data)

    def _get_matchers_list(self, input_data):
        raw_matchers_list = input_data["matchers"]
        for matcher_name, matcher_data in raw_matchers_list.iteritems():
            matcher.Matcher._check_attribute(matcher_name, matcher_data, "type")
            matcher_type = matcher_data["type"]
            MatcherClass = MATCHER_TYPES.get(matcher_type, None)
            if MatcherClass is None:
                print >> sys.stderr, 'Unknown matcher type: ' + matcher_type +\
                                     ' must be in: ' + str(MATCHER_TYPES.keys())
                sys.exit(-1)
            self.matchers.append(MatcherClass(matcher_name, matcher_data, self.match_path))

    def _check_and_assign_attribute(self, input_data, attribute_name):
            matcher.Matcher._check_attribute("input data", input_data, attribute_name)
            setattr(self, attribute_name, input_data[attribute_name])

    # This transforms annoying unicode string into common byte string
    @staticmethod
    def _byteify(input_data):
        if isinstance(input_data, dict):
            return dict([(InputData._byteify(key), InputData._byteify(value)) for key, value in input_data.iteritems()])
        elif isinstance(input_data, list):
            return [InputData._byteify(element) for element in input_data]
        elif isinstance(input_data, unicode):
            return input_data.encode('utf-8')
        else:
            return input_data


if __name__ == "__main__":
    if len(sys.argv) == 2:
        _, _input_path = sys.argv
        main(os.path.abspath(_input_path))
    elif len(sys.argv) == 1:
        from gui import gui
        gui.main()
