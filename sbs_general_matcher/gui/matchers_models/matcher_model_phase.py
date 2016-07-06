import sys
from matcher_model_default import MatcherModelDefault
from matchers.phase_matcher import PhaseMatcher


class MatcherModelPhase(MatcherModelDefault):

    def get_matcher_dict(self):
        matcher_dict = {}
        matcher_dict["type"] = "phase"
        matcher_dict["ip"] = self._ip
        matcher_dict["beam1_path"] = self._beam1_path
        matcher_dict["beam2_path"] = self._beam2_path
        matcher_dict["use_errors"] = self._use_errors
        matcher_dict["propagation"] = self._propagation
        return matcher_dict

    def create_matcher(self, match_path):
        self._matcher = PhaseMatcher(self._name, self.get_matcher_dict(), match_path)


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
