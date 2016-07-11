import sys
from matcher_view_default import MatcherViewDefault, MatcherControllerDefault
from ..matchers_models.matcher_model_amp import MatcherModelAmp


class MatcherViewAmp(MatcherViewDefault):

    def __init__(self, controller, parent=None):
        super(MatcherViewAmp, self).__init__(controller, parent)
        self.setWindowTitle("New amp matcher")


class MatcherControllerAmp(MatcherControllerDefault):

    def _get_matcher_model(self, main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation):
        return MatcherModelAmp(main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation)

    def _get_matcher_prefix(self):
        return "amp"


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
