import sys
import os
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from ..widgets import FileSelectionComboWidget


class MatcherViewDefault(QtWidgets.QDialog):

    IPS = ["1", "2", "3", "4", "5", "6", "7", "8"]

    def __init__(self, controller, parent=None):
        super(MatcherViewDefault, self).__init__(parent)
        self._controller = controller
        self._selected_matcher = None
        self._build_gui()

    def _build_gui(self):
        self.resize(655, 184)
        main_layout = QtWidgets.QVBoxLayout()
        self.setLayout(main_layout)

        selectors_layout = QtWidgets.QVBoxLayout()
        self._beam1_file_selector = FileSelectionComboWidget(label_text="Beam 1 path:")
        selectors_layout.addWidget(self._beam1_file_selector)
        self._beam2_file_selector = FileSelectionComboWidget(label_text="Beam 2 path:")
        selectors_layout.addWidget(self._beam2_file_selector)
        main_layout.addLayout(selectors_layout)

        other_controls_layout = QtWidgets.QHBoxLayout()
        self._ip_combo_box = QtWidgets.QComboBox()
        self._ip_combo_box.addItem("Select IP")
        self._ip_combo_box.addItems(MatcherViewDefault.IPS)
        other_controls_layout.addWidget(self._ip_combo_box)
        self._use_errors_checkbox = QtWidgets.QCheckBox("Use errorbars")
        other_controls_layout.addWidget(self._use_errors_checkbox, Qt.AlignRight)
        self._back_prop_checkbox = QtWidgets.QCheckBox("Back propagation")
        other_controls_layout.addWidget(self._back_prop_checkbox, Qt.AlignRight)
        main_layout.addLayout(other_controls_layout)

        buttons_layout = QtWidgets.QHBoxLayout()
        add_button = QtWidgets.QPushButton("Add", self)
        add_button.clicked.connect(self._controller._accept_action)
        buttons_layout.addWidget(add_button)
        cancel_button = QtWidgets.QPushButton("Cancel", self)
        cancel_button.clicked.connect(self._controller._cancel_action)
        buttons_layout.addWidget(cancel_button)
        main_layout.addLayout(buttons_layout)

    def get_ip(self):
        try:
            return int(self._ip_combo_box.currentText())
        except ValueError:
            return None

    def get_beam1_path(self):
        selected_file = str(self._beam1_file_selector.get_selected_file())
        if selected_file == "":
            return None
        return selected_file

    def get_beam2_path(self):
        selected_file = str(self._beam2_file_selector.get_selected_file())
        if selected_file == "":
            return None
        return selected_file

    def reload_list(self, beam, dir_list):
        if beam == 1:
            self._beam1_file_selector.reload_items(dir_list)
        elif beam == 2:
            self._beam2_file_selector.reload_items(dir_list)

    def get_use_errors(self):
        return self._use_errors_checkbox.isChecked()

    def get_back_prop(self):
        return self._back_prop_checkbox.isChecked()

    def show_error(self, text):
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setText(text)
        msg.setWindowTitle("Input error")
        msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msg.exec_()

    def plot(self, dto):
        raise NotImplementedError


class MatcherControllerDefault(object):

    def __init__(self, main_controller):
        self._main_controller = main_controller
        self._view = MatcherViewDefault(self)
        self._matcher_model = None
        self.load_possible_measurements()

    def show_view(self):
        return self._view.exec_()

    def clone_matcher(self, matcher_to_clone):
        ip = matcher_to_clone.get_ip()
        name = self._get_new_name(ip)
        beam1_path = matcher_to_clone.get_beam1_path()
        beam2_path = matcher_to_clone.get_beam2_path()
        use_errors = matcher_to_clone.get_use_errors()
        propagation = matcher_to_clone.get_propagation()

        self._matcher_model = self._get_matcher_model(
            self._main_controller, name,
            beam1_path, beam2_path, ip, use_errors, propagation
        )

    def load_possible_measurements(self):
        for beam in [1, 2]:
            possible_measurements = self._main_controller.get_posible_measurements(beam)
            self._view.reload_list(beam, possible_measurements)

    def _accept_action(self):
        if not self._check_input():
            return
        name = self._get_new_name(self._view.get_ip())
        ip = self._view.get_ip()
        beam1_path = self._view.get_beam1_path()
        beam2_path = self._view.get_beam2_path()
        use_errors = self._view.get_use_errors()
        if not self._view.get_back_prop():
            propagation = "front"
        else:
            propagation = "back"

        self._matcher_model = self._get_matcher_model(
            self._main_controller, name,
            beam1_path, beam2_path, ip, use_errors, propagation
        )
        self._view.accept()

    def get_selected_matcher_model(self):
        return self._matcher_model

    def _get_matcher_prefix(self):
        raise NotImplementedError

    def _get_new_name(self, ip):
        counter = 1
        suffix = ""
        if not self._view.get_back_prop():
            propagation = "f"
        else:
            propagation = "b"

        def new_name(suffix):
            return (self._get_matcher_prefix() + "_ip" + str(ip) + "_" +
                    propagation + suffix)

        while not self._main_controller.is_this_matcher_name_ok(new_name(suffix)):
            counter += 1
            suffix = "_" + str(counter)
        return new_name(suffix)

    def _cancel_action(self):
        return self._view.reject()

    def _check_input(self):
        ip = self._view.get_ip()
        if ip not in [1, 2, 3, 4, 5, 6, 7, 8]:
            self._view.show_error("Please select a valid IP.")
            return False
        beam1_path = self._view.get_beam1_path()
        beam2_path = self._view.get_beam2_path()
        if beam1_path is None and beam2_path is None:
            self._view.show_error("Al least one beam results directory has to be given.")
            return False
        if not beam1_path is None and not os.path.isdir(beam1_path):
            self._view.show_error("The selected beam 1 path must be an existing directory, containing GetLLM results.")
            return False
        if not beam2_path is None and not os.path.isdir(beam2_path):
            self._view.show_error("The selected beam 2 path must be an existing directory, containing GetLLM results.")
            return False
        use_errors = self._view.get_use_errors()
        if not type(use_errors) is bool:
            self._view.show_error("Use errors must be a boolean.")
            return False
        back_prop = self._view.get_back_prop()
        if not type(back_prop) is bool:
            self._view.show_error("Back propagation must be a boolean.")
            return False
        return True

    def _get_matcher_model(self, main_controller, name, beam1_path, beam2_path, ip, use_errors, propagation):
        raise NotImplementedError


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
