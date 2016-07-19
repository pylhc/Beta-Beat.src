import sys
import os
import shutil
from Python_Classes4MAD import metaclass


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
        self._elements_positions = None

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

    def get_ignore_vars_list(self):
        return self._ignore_vars_list

    def set_ignore_vars_list(self, ignore_vars_list):
        if self._matcher is not None:
            self._matcher.set_exclude_variables(ignore_vars_list)
        else:
            raise ValueError("Cannot set excluded variables of not created matcher")

    def set_disabled_constraints(self, disabled_constraints):
        if self._matcher is not None:
            self._matcher.set_disabled_constraints(disabled_constraints)
        else:
            raise ValueError("Cannot set diabled constraints of not created matcher")

    def get_variables_for_beam(self):
        return {1: self._matcher.get_variables_for_beam(1),
                2: self._matcher.get_variables_for_beam(2)}

    def get_common_variables(self):
        return self._matcher.get_common_variables()

    def get_matcher_dict(self):
        matcher_dict = {}
        matcher_dict["ip"] = self._ip
        if self._beam1_path is not None:
            matcher_dict["beam1_path"] = self._beam1_path
        if self._beam2_path is not None:
            matcher_dict["beam2_path"] = self._beam2_path
        matcher_dict["use_errors"] = self._use_errors
        matcher_dict["propagation"] = self._propagation
        return matcher_dict

    def create_matcher(self, match_path):
        raise NotImplementedError

    def delete_matcher(self):
        shutil.rmtree(self.get_beam1_output_path())
        shutil.rmtree(self.get_beam2_output_path())

    def get_plotter(self, figures):
        raise NotImplementedError

    def get_matcher(self):
        return self._matcher

    def get_elements_positions(self, beam):
        if self._elements_positions is None:
            self._read_elements_positions()
        return self._elements_positions[beam]

    def _read_elements_positions(self):
        self._elements_positions = {}
        for beam in self._matcher.get_beams():
            self._elements_positions[beam] = {}
            segment_model = metaclass.twiss(
                os.path.join(getattr(self, "get_beam" + str(beam) + "_output_path")(),
                             "sbs", "twiss_IP" + str(self.get_ip()) + ".dat")
            )
            for index in range(len(segment_model.NAME)):
                self._elements_positions[beam][segment_model.S[index]] = segment_model.NAME[index]


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
