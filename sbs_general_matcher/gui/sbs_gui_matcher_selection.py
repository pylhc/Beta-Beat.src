import os
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from matchers_models.matcher_model_coupling import MatcherModelCoupling
from matchers_models.matcher_model_phase import MatcherModelPhase
from matchers_models.matcher_model_amp import MatcherModelAmp
from matchers_models.matcher_model_kmod import MatcherModelKmod
from widgets import FileSelectionComboWidget


class SbSGuiMatcherTypeSelection(QtWidgets.QDialog):
    """
    View that shows the different matcher types available and opens the
    matcher definition view.
    """

    WINDOW_TITLE = "Select matcher to add"
    MATCHER_TYPES = {
        "phase": MatcherModelPhase,
        "kmod": MatcherModelKmod,
        "beta from amplitude": MatcherModelAmp,
        "coupling": MatcherModelCoupling,
    }

    def __init__(self, main_controller, clone_matcher=None, parent=None):
        super(SbSGuiMatcherTypeSelection, self).__init__(parent)
        self._main_controller = main_controller
        self._clone_matcher = clone_matcher
        self._selected_matcher = None
        self._build_gui()

    def _build_gui(self):
        self.setWindowTitle(SbSGuiMatcherTypeSelection.WINDOW_TITLE)
        self.resize(320, 240)

        dialog_layout = QtWidgets.QVBoxLayout(self)
        buttons_layout = QtWidgets.QHBoxLayout()

        self._matchers_list_widget = SbSGuiMatcherTypeSelection._get_matchers_list_widget()
        dialog_layout.addWidget(self._matchers_list_widget)

        if self._clone_matcher is None:
            add_button = QtWidgets.QPushButton("Add", self)

            add_button.clicked.connect(self._accept_action)
            buttons_layout.addWidget(add_button)
        else:
            clone_button = QtWidgets.QPushButton("Clone", self)

            clone_button.clicked.connect(self._clone_action)
            buttons_layout.addWidget(clone_button)

        cancel_button = QtWidgets.QPushButton("Cancel", self)
        cancel_button.clicked.connect(self.reject)
        buttons_layout.addWidget(cancel_button)

        dialog_layout.addLayout(buttons_layout, stretch=1)

    @staticmethod
    def _get_matchers_list_widget():
        matchers_list_widget = QtWidgets.QListWidget()
        for key in SbSGuiMatcherTypeSelection.MATCHER_TYPES:
            matchers_list_widget.addItem(key)
        return matchers_list_widget

    def _accept_action(self):
        self._selected_matcher = str(self._matchers_list_widget.currentItem().text())
        selected_model_cls = SbSGuiMatcherTypeSelection.MATCHER_TYPES[self._selected_matcher]
        matcher_controller_instance = MatcherControllerDefault(
            self._selected_matcher,
            selected_model_cls,
            self._main_controller
        )
        result_code = matcher_controller_instance.show_view()
        if result_code == QtWidgets.QDialog.Accepted:
            self._selected_matcher = matcher_controller_instance.get_selected_matcher_model()
            self.accept()

    def _clone_action(self):
        self._selected_matcher = str(self._matchers_list_widget.currentItem().text())
        selected_model_cls = SbSGuiMatcherTypeSelection.MATCHER_TYPES[self._selected_matcher]
        matcher_controller_instance = MatcherControllerDefault(
            self._selected_matcher,
            selected_model_cls,
            self._main_controller
        )
        matcher_controller_instance.clone_matcher(self._clone_matcher)
        self._selected_matcher = matcher_controller_instance.get_selected_matcher_model()
        self.accept()

    def get_selected_matcher(self):
        return self._selected_matcher

# End SbSGuiMatcherTypeSelection ##############################################


class MatcherViewDefault(QtWidgets.QDialog):

    def __init__(self, name, controller, parent=None):
        super(MatcherViewDefault, self).__init__(parent)
        self._controller = controller
        self._name = name
        self._build_gui()

    def _build_gui(self):
        self.resize(705, 184)
        main_layout = QtWidgets.QVBoxLayout()
        self.setLayout(main_layout)
        self.setWindowTitle("New " + self._name + " matcher")

        selectors_layout = QtWidgets.QHBoxLayout()
        self._beam_combo_box = QtWidgets.QComboBox()
        for beam in ("Beam 1", "Beam 2"):
            self._beam_combo_box.addItem(beam)
        self._beam_combo_box.setCurrentIndex(0)
        self._beam_combo_box.currentIndexChanged.connect(
            lambda index: self._controller.load_possible_measurements(
                index + 1
            )
        )
        selectors_layout.addWidget(self._beam_combo_box)
        selectors_layout.addWidget(QtWidgets.QLabel("->"))
        self._file_selector = FileSelectionComboWidget(label_text="Path:")
        self._file_selector.on_item_change = self._controller.reload_segments
        selectors_layout.addWidget(self._file_selector)

        main_layout.addLayout(selectors_layout)

        other_controls_layout = QtWidgets.QHBoxLayout()
        self._segment_combo_box = QtWidgets.QComboBox()
        other_controls_layout.addWidget(self._segment_combo_box)
        self._use_errors_checkbox = QtWidgets.QCheckBox("Use errorbars")
        other_controls_layout.addWidget(self._use_errors_checkbox,
                                        Qt.AlignRight)
        self._back_prop_checkbox = QtWidgets.QCheckBox("Back propagation")
        other_controls_layout.addWidget(self._back_prop_checkbox,
                                        Qt.AlignRight)
        main_layout.addLayout(other_controls_layout)

        buttons_layout = QtWidgets.QHBoxLayout()
        add_button = QtWidgets.QPushButton("Add", self)
        add_button.clicked.connect(self._controller._accept_action)
        buttons_layout.addWidget(add_button)
        cancel_button = QtWidgets.QPushButton("Cancel", self)
        cancel_button.clicked.connect(self._controller._cancel_action)
        buttons_layout.addWidget(cancel_button)
        main_layout.addLayout(buttons_layout)

    def get_selected_beam(self):
        return self._beam_combo_box.currentIndex() + 1

    def get_segment_label(self):
        if self._segment_combo_box.count() == 0:
            return None
        return str(self._segment_combo_box.currentText())

    def get_selected_path(self):
        selected_file = str(self._file_selector.get_selected_file())
        if selected_file == "":
            return None
        return selected_file

    def reload_meas_list(self, dir_list):
        self._file_selector.reload_items(dir_list)

    def reload_segment_list(self, segments):
        self._segment_combo_box.clear()
        for segment in segments:
            self._segment_combo_box.addItem(segment)
        self._segment_combo_box.setCurrentIndex(0)

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

    def __init__(self, name, model_class, main_controller):
        self._main_controller = main_controller
        self._view = MatcherViewDefault(name, self)
        self.load_possible_measurements(1)
        self._name = name
        self._model_cls = model_class
        self._matcher_model = None

    def show_view(self):
        return self._view.exec_()

    def clone_matcher(self, matcher_to_clone):
        label = matcher_to_clone.get_segment_label()
        name = self._get_new_name(label)
        meas_path = matcher_to_clone.get_meas_path()
        use_errors = matcher_to_clone.get_use_errors()
        propagation = matcher_to_clone.get_propagation()

        self._matcher_model = self._get_matcher_model(
            self._main_controller.get_match_path(), name, beam,
            meas_path, label, use_errors, propagation
        )

    def reload_segments(self, selected_meas):
        selected_meas_sbs = os.path.join(selected_meas, "sbs")
        segments = []
        if not os.path.isdir(selected_meas_sbs):
            self._view.reload_segment_list(segments)
            return
        for filename in os.listdir(selected_meas_sbs):
            if (filename.startswith("twiss_") and
                    filename.endswith(".dat") and
                    "cor" not in filename and
                    "back" not in filename):
                segments.append(
                    filename.replace("twiss_", "").replace(".dat", "")
                )
        self._view.reload_segment_list(segments)

    def load_possible_measurements(self, beam):
        possible_measurements = self._main_controller.get_posible_measurements(
            beam
        )
        self._view.reload_meas_list(possible_measurements)
        self._view.reload_segment_list([])

    def _accept_action(self):
        if not self._check_input():
            return
        beam = self._view.get_selected_beam()
        name = self._get_new_name(self._view.get_segment_label(), beam)
        label = self._view.get_segment_label()
        meas_path = self._view.get_selected_path()
        use_errors = self._view.get_use_errors()
        if not self._view.get_back_prop():
            propagation = "front"
        else:
            propagation = "back"

        self._matcher_model = self._get_matcher_model(
            self._main_controller.get_match_path(), name, beam,
            meas_path, label, use_errors, propagation
        )
        self._view.accept()

    def get_selected_matcher_model(self):
        return self._matcher_model

    def _get_new_name(self, label, beam):
        counter = 1
        suffix = ""
        if not self._view.get_back_prop():
            propagation = "f"
        else:
            propagation = "b"

        def new_name(suffix):
            return (self._name + "_" + label + "_" +
                    propagation + "_beam" + str(beam) + suffix)

        while not self._main_controller.is_this_matcher_name_ok(new_name(suffix)):
            counter += 1
            suffix = "_" + str(counter)
        return new_name(suffix)

    def _cancel_action(self):
        return self._view.reject()

    def _check_input(self):
        segment_label = self._view.get_segment_label()
        if segment_label is None:
            self._view.show_error("Please select a valid segment.")
            return False
        meas_path = self._view.get_selected_path()
        if meas_path is None:
            self._view.show_error("Please select a measurement directory.")
            return False
        if meas_path is not None and not os.path.isdir(meas_path):
            self._view.show_error(
                "The selected path must be an existing directory, "
                "containing GetLLM results."
            )
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

    def _get_matcher_model(self, match_path, name, beam, meas_path,
                           label, use_errors, propagation):
        return self._model_cls(match_path, name, beam,
                               meas_path, label,
                               use_errors, propagation)
