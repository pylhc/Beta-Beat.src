import sys
import os
from PyQt4 import QtGui
from PyQt4.QtCore import Qt
from ..widgets import FileSelectionDialogWidget


class MatcherViewDefault(QtGui.QDialog):

    IPS = ["1", "2", "5", "8"]

    def __init__(self, controller, parent=None):
        super(MatcherViewDefault, self).__init__(parent)
        self._controller = controller
        self._selected_matcher = None
        self._build_gui()

    def _build_gui(self):
        self.resize(655, 184)
        main_layout = QtGui.QVBoxLayout()
        self.setLayout(main_layout)

        selectors_layout = QtGui.QVBoxLayout()
        self._beam1_file_selector = FileSelectionDialogWidget(label_text="Beam 1 path:")
        selectors_layout.addWidget(self._beam1_file_selector)
        self._beam2_file_selector = FileSelectionDialogWidget(label_text="Beam 2 path:")
        selectors_layout.addWidget(self._beam2_file_selector)
        main_layout.addLayout(selectors_layout)

        other_controls_layout = QtGui.QHBoxLayout()
        self._ip_combo_box = QtGui.QComboBox()
        self._ip_combo_box.addItem("Select IP")
        self._ip_combo_box.addItems(MatcherViewDefault.IPS)
        other_controls_layout.addWidget(self._ip_combo_box)
        self._use_errors_checkbox = QtGui.QCheckBox("Use errorbars")
        other_controls_layout.addWidget(self._use_errors_checkbox, Qt.AlignRight)
        self._back_prop_checkbox = QtGui.QCheckBox("Back propagation")
        other_controls_layout.addWidget(self._back_prop_checkbox, Qt.AlignRight)
        main_layout.addLayout(other_controls_layout)

        buttons_layout = QtGui.QHBoxLayout()
        add_button = QtGui.QPushButton("Add", self)
        add_button.clicked.connect(self._controller._accept_action)
        buttons_layout.addWidget(add_button)
        cancel_button = QtGui.QPushButton("Cancel", self)
        cancel_button.clicked.connect(self._controller._cancel_action)
        buttons_layout.addWidget(cancel_button)
        main_layout.addLayout(buttons_layout)

    def get_ip(self):
        try:
            return int(self._ip_combo_box.currentText())
        except ValueError:
            return None

    def get_beam1_path(self):
        selected_file = self._beam1_file_selector.get_selected_file()
        if selected_file == "":
            return None
        return str(selected_file)

    def get_beam2_path(self):
        selected_file = self._beam2_file_selector.get_selected_file()
        if selected_file == "":
            return None
        return str(selected_file)

    def get_use_errors(self):
        return self._use_errors_checkbox.isChecked()

    def get_back_prop(self):
        return self._back_prop_checkbox.isChecked()

    def show_error(self, text):
        msg = QtGui.QMessageBox()
        msg.setIcon(QtGui.QMessageBox.Critical)
        msg.setText(text)
        msg.setWindowTitle("Input error")
        msg.setStandardButtons(QtGui.QMessageBox.Ok)
        msg.exec_()

    def plot(self, dto):
        raise NotImplementedError


class MatcherControllerDefault(object):

    def __init__(self, main_controller):
        self._main_controller = main_controller
        self._view = MatcherViewDefault(self)
        self._matcher_model = None

    def show_view(self):
        return self._view.exec_()

    def _accept_action(self):
        if not self._check_input():
            return
        name = self._get_new_name()
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

    def _get_new_name(self):
        counter = 1
        suffix = ""
        if not self._view.get_back_prop():
            propagation = "f"
        else:
            propagation = "b"

        def new_name(suffix):
            return (self._get_matcher_prefix() + "_ip" + str(self._view.get_ip()) + "_" +
                    propagation + suffix)

        while not self._main_controller.is_this_matcher_name_ok(new_name(suffix)):
            counter += 1
            suffix = "_" + str(counter)
        return new_name(suffix)

    def _cancel_action(self):
        return self._view.reject()

    def _check_input(self):
        ip = self._view.get_ip()
        if ip not in [1, 2, 5, 8]:
            self._view.show_error("Please select a valid IP.")
            return False
        beam1_path = self._view.get_beam1_path()
        if beam1_path is None and beam1_path is None:
            self._view.show_error("Al least one beam results directory has to be given.")
            return False
        if not beam1_path is None and not os.path.isdir(beam1_path):
            self._view.show_error("The selected beam 1 path must be an existing directory, containing GetLLM results.")
            return False
        beam2_path = self._view.get_beam2_path()
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
