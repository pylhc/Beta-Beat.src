import sys
import os


class MatcherModelDefault(object):

    def __init__(self, controller, name, beam1_path, beam2_path, ip, use_errors, propagation):
        self._controller = controller
        self._name = name
        self._beam1_path = beam1_path
        self._beam2_path = beam2_path
        self._ip = ip
        self._use_errors = use_errors
        self._propagation = propagation
        self._matcher = None

    def get_name(self):
        return self._name

    def get_beam1_path(self):
        return self._beam1_path

    def get_beam2_path(self):
        return self._beam2_path

    def get_beam1_output_path(self):
        return os.path.join(self._controller.get_match_path(), "Beam1_" + self._name)

    def get_beam2_output_path(self):
        return os.path.join(self._controller.get_match_path(), "Beam2_" + self._name)

    def get_ip(self):
        return self._ip

    def get_use_errors(self):
        return self._use_errors

    def get_propagation(self):
        return self._propagation

    def get_matcher_dict(self):
        matcher_dict = {}
        matcher_dict["ip"] = self._ip
        matcher_dict["beam1_path"] = self._beam1_path
        matcher_dict["beam2_path"] = self._beam2_path
        matcher_dict["use_errors"] = self._use_errors
        matcher_dict["propagation"] = self._propagation
        return matcher_dict

    def create_matcher(self, match_path):
        raise NotImplementedError

    def get_plotter(self, figures):
        raise NotImplementedError

    def get_matcher(self):
        return self._matcher


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
