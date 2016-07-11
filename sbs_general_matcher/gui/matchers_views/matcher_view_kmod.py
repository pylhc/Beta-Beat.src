import sys
from matcher_view_default import MatcherViewDefault, MatcherControllerDefault
from ..matchers_models.matcher_model_kmod import MatcherModelKmod


class MatcherViewKmod(MatcherViewDefault):

    def __init__(self, controller, parent=None):
        super(MatcherViewKmod, self).__init__(controller, parent)
        self.setWindowTitle("New kmod matcher")


class MatcherControllerKmod(MatcherControllerDefault):

    def _get_matcher_model(self, main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation):
        return MatcherModelKmod(main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation)

    def _get_matcher_prefix(self):
        return "kmod"


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
