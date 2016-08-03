import sys
from matcher_view_default import MatcherViewDefault, MatcherControllerDefault
from ..matchers_models.matcher_model_coupling import MatcherModelCoupling


class MatcherViewCoupling(MatcherViewDefault):

    def __init__(self, controller, parent=None):
        super(MatcherViewCoupling, self).__init__(controller, parent)
        self.setWindowTitle("New coupling matcher")


class MatcherControllerCoupling(MatcherControllerDefault):

    def _get_matcher_model(self, main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation):
        return MatcherModelCoupling(main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation)

    def _get_matcher_prefix(self):
        return "coupling"

if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
