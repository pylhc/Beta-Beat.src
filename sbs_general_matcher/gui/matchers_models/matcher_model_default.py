import sys


class MatcherModelDefault(object):

    def __init__(self, name, beam1_path, beam2_path, ip, use_errors, propagation):
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

    def get_ip(self):
        return self._ip

    def get_use_errors(self):
        return self._use_errors

    def get_propagation(self):
        return self._propagation

    def get_matcher_dict(self):
        raise NotImplementedError

    def create_matcher(self, match_path):
        raise NotImplementedError

    def get_matcher(self):
        return self._matcher


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
