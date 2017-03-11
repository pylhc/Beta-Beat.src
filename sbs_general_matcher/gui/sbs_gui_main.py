import sys
import os
import subprocess
import logging
from PyQt5 import QtWidgets
from PyQt5.QtCore import QThread, Qt, QFileSystemWatcher, pyqtSignal
from contextlib import contextmanager
from sbs_gui_matcher_selection import SbSGuiMatcherSelection
from widgets import InitialConfigPopup
from sbs_gui_match_result_view import SbSGuiMatchResultController
import sbs_general_matcher


LOGGER = logging.getLogger(__name__)


class SbSGuiMain(QtWidgets.QMainWindow):

    WINDOW_TITLE = "Segment-by-segment general matcher GUI"

    def __init__(self, controller, parent=None):
        super(SbSGuiMain, self).__init__(parent)

        self._controller = controller
        self._active_background_dialog = None
        self._build_gui()

    def _build_gui(self):

        self.setWindowTitle(SbSGuiMain.WINDOW_TITLE)
        screen_shape = QtWidgets.QDesktopWidget().screenGeometry()
        self.resize(2 * screen_shape.width() / 3,
                    2 * screen_shape.height() / 3)

        self._main_widget = SbSGuiMain.SbSMainWidget(self._controller)
        self.setCentralWidget(self._main_widget)

        main_menu = self.menuBar()
        main_menu.setNativeMenuBar(False)
        matchers_menu = main_menu.addMenu("Matchers")
        matchers_menu.addAction(self._get_new_matcher_action())
        matchers_menu.addAction(self._get_clone_matcher_action())
        matchers_menu.addAction(self._get_remove_matcher_action())

    def add_tab(self, name, widget):
        self._main_widget._matchers_tabs_widget.addTab(widget, name)

    def remove_tab(self, index):
        self._main_widget._matchers_tabs_widget.removeTab(index)

    def get_selected_matcher_index(self):
        return self._main_widget._matchers_tabs_widget.currentIndex()

    def show_background_task_dialog(self, message):
        self._active_background_dialog = SbSGuiMain.BackgroundTaskDialog(
            message, parent=self
        )
        self._active_background_dialog.setModal(True)
        self._active_background_dialog.setVisible(True)

    def hide_background_task_dialog(self):
        if self._active_background_dialog is not None:
            self._active_background_dialog.setVisible(False)
            self._active_background_dialog = None

    def show_error_dialog(self, title, message):
        message_box = QtWidgets.QMessageBox(
            QtWidgets.QMessageBox.Critical,
            title,
            message,
            QtWidgets.QMessageBox.Ok,
            self
        )
        message_box.exec_()

    def is_minimize_selected(self):
        return self._main_widget.is_minimize_selected()

    def _get_new_matcher_action(self):
        new_matcher_action = QtWidgets.QAction("New matcher...", self)
        new_matcher_action.triggered.connect(self._controller.new_matcher)
        return new_matcher_action

    def _get_remove_matcher_action(self):
        remove_matcher_action = QtWidgets.QAction("Remove matcher", self)
        remove_matcher_action.triggered.connect(self._controller.remove_matcher)
        return remove_matcher_action

    def _get_clone_matcher_action(self):
        clone_matcher_action = QtWidgets.QAction("Clone matcher", self)
        clone_matcher_action.triggered.connect(self._controller.clone_matcher)
        return clone_matcher_action

    class SbSMainWidget(QtWidgets.QWidget):
        def __init__(self, controller, parent=None):
            super(SbSGuiMain.SbSMainWidget, self).__init__(parent)
            self._controller = controller
            self._build_gui()

        def is_minimize_selected(self):
            return self._minimize_checkbox.isChecked()

        def _build_gui(self):
            main_layout = QtWidgets.QVBoxLayout()
            self.setLayout(main_layout)

            self._matchers_tabs_widget = QtWidgets.QTabWidget()
            main_layout.addWidget(self._matchers_tabs_widget)

            lower_panel_layout = QtWidgets.QHBoxLayout()
            buttons_layout = QtWidgets.QVBoxLayout()
            global_options_layout = QtWidgets.QVBoxLayout()
            lower_panel_layout.addLayout(buttons_layout, stretch=3)
            lower_panel_layout.addLayout(global_options_layout)
            main_layout.addLayout(lower_panel_layout)

            self._minimize_checkbox = QtWidgets.QCheckBox("Minimize variables")
            global_options_layout.addWidget(self._minimize_checkbox)

            run_button = QtWidgets.QPushButton("Run matching")
            run_button.clicked.connect(self._controller.run_matching)
            buttons_layout.addWidget(run_button)

            edit_corr_button = QtWidgets.QPushButton("Edit corrections file")
            edit_corr_button.clicked.connect(self._controller.edit_corrections_file)
            buttons_layout.addWidget(edit_corr_button)

    class BackgroundTaskDialog(QtWidgets.QMessageBox):
        def __init__(self, message, parent=None):
            super(SbSGuiMain.BackgroundTaskDialog, self).__init__(
                QtWidgets.QMessageBox.NoIcon,
                "Please wait...",
                message
            )
            self.setWindowFlags(Qt.CustomizeWindowHint)
            self.setStandardButtons(QtWidgets.QMessageBox.NoButton)
            self.resize(420, 240)


class SbSGuiMainController(object):

    def __init__(self):
        self._view = SbSGuiMain(self)
        self._match_path = None
        self._possible_measurements = {1: [], 2: []}
        self._input_dir = None
        self._matchers_tabs = []
        self._current_thread = None
        self._active_watcher = None

    @staticmethod
    def ask_for_initial_config(lhc_mode, match_path):
        initial_config_popup = InitialConfigPopup(lhc_mode, match_path)
        initial_config_popup.setWindowTitle("Please choose an output path")
        result_code = initial_config_popup.exec_()
        if result_code == QtWidgets.QDialog.Accepted:
            return initial_config_popup.get_selected_lhc_mode(), initial_config_popup.get_selected_file()
        else:
            return None, None

    def set_match_path(self, match_path):
        self._match_path = match_path
        self._corrections_file = os.path.join(match_path,
                                              "changeparameters.madx")

    def set_lhc_mode(self, lhc_mode):
        self._lhc_mode = lhc_mode

    def set_input_dir(self, input_dir):
        if input_dir is None:
            return
        self._input_dir = os.path.abspath(input_dir)
        self._find_measurements()

    def _find_measurements(self):
        self._possible_measurements = {1: [], 2: []}
        if self._input_dir is None:
            return
        for lhcb1or2 in os.listdir(self._input_dir):
            if lhcb1or2 == "LHCB1":
                beam = 1
            elif lhcb1or2 == "LHCB2":
                beam = 2
            else:
                continue
            results_path = os.path.join(self._input_dir, lhcb1or2,
                                        "Results")
            if not os.path.isdir(results_path):
                continue
            for individual_result in os.listdir(results_path):
                individual_result_path = os.path.join(results_path,
                                                      individual_result)
                if os.path.isdir(individual_result_path):
                    self._possible_measurements[beam].append(
                        individual_result_path
                    )

    def get_posible_measurements(self, beam):
        return self._possible_measurements[beam]

    def get_match_path(self):
        return self._match_path

    def show_view(self):
        self._view.show()

    def new_matcher(self):
        self._new_matcher_from_chooser()

    def clone_matcher(self):
        if len(self._matchers_tabs) == 0:
            return
        index = self._view.get_selected_matcher_index()
        self._new_matcher_from_chooser(self._matchers_tabs[index].model)

    def _new_matcher_from_chooser(self, matcher_to_clone=None):
        sbs_gui_matcher_selection_dialog = SbSGuiMatcherSelection(
            self, clone_matcher=matcher_to_clone
        )
        result_code = sbs_gui_matcher_selection_dialog.exec_()
        if result_code == QtWidgets.QDialog.Accepted:
            selected_matcher_model = sbs_gui_matcher_selection_dialog.get_selected_matcher()
            with self._heavy_task("Copying files..."):
                selected_matcher_model.create_matcher(self._match_path)
            variables_for_beam = selected_matcher_model.get_variables_for_beam()
            variables_common = selected_matcher_model.get_common_variables()
            tab_controller = SbSGuiMatchResultController(
                variables_for_beam,
                variables_common,
                lambda beam: selected_matcher_model.get_elements_positions(beam)
            )
            self._matchers_tabs.append(SbSGuiMainController.Tab(selected_matcher_model, tab_controller))
            self._view.add_tab(selected_matcher_model.get_name(), tab_controller.get_view())

    @contextmanager
    def _heavy_task(self, message):
        self._view.show_background_task_dialog(message)
        try:
            yield
        except Exception as e:
            LOGGER.exception(str(e))
            self._view.show_error_dialog("Error", str(e))
        finally:
            self._view.hide_background_task_dialog()

    def is_this_matcher_name_ok(self, matcher_name):
        for matcher_tab in self._matchers_tabs:
            model_name = matcher_tab.model.get_name()
            if matcher_name == model_name:
                return False
        return True

    def remove_matcher(self):
        if len(self._matchers_tabs) == 0:
            return
        index = self._view.get_selected_matcher_index()
        with self._heavy_task("Deleting matcher..."):
            self._matchers_tabs[index].model.delete_matcher()
            del(self._matchers_tabs[index])
            self._view.remove_tab(index)

    def run_matching(self, just_twiss=False):
        matchers_list = []
        for matcher_tab in self._matchers_tabs:
            matcher_tab.model.set_ignore_vars_list(
                matcher_tab.results_controller.get_unselected_variables()
            )
            matcher_tab.model.set_disabled_constraints(
                matcher_tab.results_controller.get_disabled_constraints()
            )
            matchers_list.append(matcher_tab.model.get_matcher())
        minimize = self._view.is_minimize_selected()
        input_data = sbs_general_matcher.InputData.init_from_matchers_list(
            self._lhc_mode, self._match_path, minimize, matchers_list
        )

        if not just_twiss:
            def background_task():
                had_active_watcher = False
                if self._active_watcher is not None:
                    had_active_watcher = True
                    self._active_watcher.removePath(self._match_path)
                    self._active_watcher = None
                sbs_general_matcher.run_full_madx_matching(input_data)
                if had_active_watcher:
                    self._watch_dir(self._match_path)
        else:
            def background_task():
                sbs_general_matcher.run_twiss_and_sbs(input_data)

        self._current_thread = BackgroundThread(
            self._view,
            background_task,
            message="Running matching...",
            on_end_function=self._on_match_end,
            on_exception_function=self._on_match_exception
        )
        self._current_thread.start()

    def _on_match_end(self):
        for index in range(len(self._matchers_tabs)):
            matcher_model = self._matchers_tabs[index].model
            results_controller = self._matchers_tabs[index].results_controller
            figures = results_controller.get_figures()
            matcher_model.get_plotter(figures).plot()
            results_controller.update_variables(
                matcher_model.get_match_results()
            )
        self._current_thread = None

    def _on_match_exception(self, message):
        self._current_thread = None

    def edit_corrections_file(self):
        if not os.path.isfile(self._corrections_file):
            open(self._corrections_file, "a").close()  # Create empty file
        self._launch_text_editor(self._corrections_file)
        if (self._active_watcher is not None and
                self._corrections_file in self._active_watcher.files()):
            return
        self._watch_dir(self._match_path)

    def _watch_dir(self, directory):
        self._active_watcher = QFileSystemWatcher([directory])
        self._active_watcher.directoryChanged.connect(self._match_dir_changed)

    def _launch_text_editor(self, file_path):
        if sys.platform.startswith('darwin'):
            subprocess.call(('open', file_path))
        elif os.name == 'nt':
            os.startfile(file_path)
        elif os.name == 'posix':
            subprocess.call(('xdg-open', file_path))

    def _match_dir_changed(self, path):
        if os.path.samefile(path, self._match_path):
            if self._current_thread is None:
                self.run_matching(just_twiss=True)

    class Tab(object):
        def __init__(self, matcher_model, matcher_results_controller):
            self.model = matcher_model
            self.results_controller = matcher_results_controller


class BackgroundThread(QThread):

    on_exception = pyqtSignal([str])

    def __init__(self, view, function, message=None,
                 on_end_function=None, on_exception_function=None):
        QThread.__init__(self)
        self._view = view
        self._function = function
        self._message = message
        self._on_end_function = on_end_function
        self._on_exception_function = on_exception_function

    def run(self):
        try:
            self._function()
        except Exception as e:
            LOGGER.exception(str(e))
            self.on_exception.emit(str(e))

    def start(self):
        self.finished.connect(self._on_end)
        self.on_exception.connect(self._on_exception)
        super(BackgroundThread, self).start()
        self._view.show_background_task_dialog(self._message)

    def _on_end(self):
        self._view.hide_background_task_dialog()
        self._on_end_function()

    def _on_exception(self, exception_message):
        self._view.hide_background_task_dialog()
        self._view.show_error_dialog("Error", exception_message)
        self._on_exception_function(exception_message)


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
