import sys
import os
import subprocess
from PyQt4 import QtGui
from PyQt4.QtCore import QThread, Qt, QFileSystemWatcher, pyqtSignal
from sbs_gui_matcher_selection import SbSGuiMatcherSelection
from widgets import InitialConfigPopup
from sbs_gui_match_result_view import SbSGuiMatchResultController
import sbs_general_matcher


class SbSGuiMain(QtGui.QMainWindow):

    WINDOW_TITLE = "Segment-by-segment general matcher GUI"

    def __init__(self, controller, parent=None):
        super(SbSGuiMain, self).__init__(parent)

        self._controller = controller
        self._active_background_dialog = None
        self._build_gui()

    def _build_gui(self):

        self.setWindowTitle(SbSGuiMain.WINDOW_TITLE)
        screen_shape = QtGui.QDesktopWidget().screenGeometry()
        self.resize(2 * screen_shape.width() / 3, 2 * screen_shape.height() / 3)

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
        self._active_background_dialog = SbSGuiMain.BackgroundTaskDialog(message)
        self._active_background_dialog.setModal(True)
        self._active_background_dialog.setVisible(True)

    def hide_background_task_dialog(self):
        self._active_background_dialog.setVisible(False)
        self._active_background_dialog = None

    def show_error_dialog(self, title, message):
        message_box = QtGui.QMessageBox(
            QtGui.QMessageBox.Critical,
            title,
            message,
            QtGui.QMessageBox.Ok,
            self
        )
        message_box.exec_()

    def _get_new_matcher_action(self):
        new_matcher_action = QtGui.QAction("New matcher...", self)
        new_matcher_action.triggered.connect(self._controller.new_matcher)
        return new_matcher_action

    def _get_remove_matcher_action(self):
        remove_matcher_action = QtGui.QAction("Remove matcher", self)
        remove_matcher_action.triggered.connect(self._controller.remove_matcher)
        return remove_matcher_action

    def _get_clone_matcher_action(self):
        clone_matcher_action = QtGui.QAction("Clone matcher", self)
        clone_matcher_action.triggered.connect(self._controller.clone_matcher)
        return clone_matcher_action

    class SbSMainWidget(QtGui.QWidget):
        def __init__(self, controller, parent=None):
            super(SbSGuiMain.SbSMainWidget, self).__init__(parent)
            self._controller = controller
            self._build_gui()

        def _build_gui(self):
            main_layout = QtGui.QVBoxLayout()
            self.setLayout(main_layout)

            self._matchers_tabs_widget = QtGui.QTabWidget()
            main_layout.addWidget(self._matchers_tabs_widget)

            run_button = QtGui.QPushButton("Run matching")
            run_button.clicked.connect(self._controller.run_matching)
            main_layout.addWidget(run_button)

            edit_corr_button = QtGui.QPushButton("Edit corrections file")
            edit_corr_button.clicked.connect(self._controller.edit_corrections_file)
            main_layout.addWidget(edit_corr_button)

    class BackgroundTaskDialog(QtGui.QMessageBox):
        def __init__(self, message, parent=None):
            super(SbSGuiMain.BackgroundTaskDialog, self).__init__(
                QtGui.QMessageBox.NoIcon,
                "Please wait...",
                message
            )
            self.setWindowFlags(Qt.CustomizeWindowHint)
            self.setStandardButtons(QtGui.QMessageBox.NoButton)
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
        if result_code == QtGui.QDialog.Accepted:
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
            results_path = os.path.join(self._input_dir, lhcb1or2, "Results")
            if not os.path.isdir(results_path):
                continue
            for individual_result in os.listdir(results_path):
                individual_result_path = os.path.join(results_path,
                                                      individual_result)
                self._possible_measurements[beam].append(individual_result_path)

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
        if result_code == QtGui.QDialog.Accepted:
            selected_matcher_model = sbs_gui_matcher_selection_dialog.get_selected_matcher()
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
        input_data = sbs_general_matcher.InputData.init_from_matchers_list(
            self._lhc_mode, self._match_path, matchers_list
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

        self._current_thread = SbSGuiMainController.BackgroundThread(background_task)
        self._current_thread.finished.connect(self._on_match_end)
        self._current_thread.on_exception.connect(self._on_match_exception)
        self._current_thread.start()
        self._view.show_background_task_dialog("Running matching...")

    def _on_match_end(self):
        for index in range(len(self._matchers_tabs)):
            matcher_model = self._matchers_tabs[index].model
            results_controller = self._matchers_tabs[index].results_controller
            figures = results_controller.get_figures()
            matcher_model.get_plotter(figures).plot()
            results_controller.update_variables(matcher_model.get_match_results())
        self._view.hide_background_task_dialog()
        self._current_thread = None

    def _on_match_exception(self, message):
        self._view.show_error_dialog("Error", message)
        self._view.hide_background_task_dialog()
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

        def __init__(self, function):
            QThread.__init__(self)
            self._function = function

        def run(self):
            try:
                self._function()
            except Exception as e:
                self.on_exception.emit(str(e))


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
