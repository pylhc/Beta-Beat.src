import sys
from PyQt4 import QtGui
from sbs_gui_matcher_selection import SbSGuiMatcherSelection
from widgets import InitialConfigPopup
from sbs_gui_match_result_view import SbSGuiMatchResultController
import sbs_general_matcher


class SbSGuiMain(QtGui.QMainWindow):

    WINDOW_TITLE = "Segment-by-segment general matcher GUI"

    def __init__(self, controller, parent=None):
        super(SbSGuiMain, self).__init__(parent)

        self._controller = controller
        self._build_gui()

    def _build_gui(self):

        self.setWindowTitle(SbSGuiMain.WINDOW_TITLE)
        screen_shape = QtGui.QDesktopWidget().screenGeometry()
        self.resize(2 * screen_shape.width() / 3, 2 * screen_shape.height() / 3)

        self._main_widget = _SbSMainWidget(self._controller)
        self.setCentralWidget(self._main_widget)

        main_menu = self.menuBar()
        main_menu.setNativeMenuBar(False)
        matchers_menu = main_menu.addMenu("Matchers")
        matchers_menu.addAction(self._get_new_matcher_action())
        matchers_menu.addAction(self._get_remove_matcher_action())

    def add_tab(self, name, widget):
        self._main_widget._matchers_tabs_widget.addTab(widget, name)

    def remove_tab(self, index):
        self._main_widget._matchers_tabs_widget.removeTab(index)

    def get_selected_matcher_index(self):
        return self._main_widget._matchers_tabs_widget.currentIndex()

    def _get_new_matcher_action(self):
        new_matcher_action = QtGui.QAction("New matcher...", self)
        new_matcher_action.triggered.connect(self._controller.new_matcher)
        return new_matcher_action

    def _get_remove_matcher_action(self):
        remove_matcher_action = QtGui.QAction("Remove matcher", self)
        remove_matcher_action.triggered.connect(self._controller.remove_matcher)
        return remove_matcher_action


class _SbSMainWidget(QtGui.QWidget):
    def __init__(self, controller, parent=None):
        super(_SbSMainWidget, self).__init__(parent)
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


class SbSGuiMainController(object):

    def __init__(self):
        self._view = SbSGuiMain(self)
        self._match_path = None
        self._matchers_tabs = []

    @staticmethod
    def ask_for_initial_config():
        initial_config_popup = InitialConfigPopup()
        initial_config_popup.setWindowTitle("Please choose an output path")
        result_code = initial_config_popup.exec_()
        if result_code == QtGui.QDialog.Accepted:
            return initial_config_popup.get_selected_lhc_mode(), initial_config_popup.get_selected_file()
        else:
            return None, None

    def set_match_path(self, match_path):
        self._match_path = match_path

    def set_lhc_mode(self, lhc_mode):
        self._lhc_mode = lhc_mode

    def get_match_path(self):
        return self._match_path

    def show_view(self):
        self._view.show()

    def new_matcher(self):
        sbs_gui_matcher_selection_dialog = SbSGuiMatcherSelection(self)
        result_code = sbs_gui_matcher_selection_dialog.exec_()
        if result_code == QtGui.QDialog.Accepted:
            selected_matcher_model = sbs_gui_matcher_selection_dialog.get_selected_matcher()
            selected_matcher_model.create_matcher(self._match_path)
            variables_for_beam = selected_matcher_model.get_variables_for_beam()
            variables_common = selected_matcher_model.get_common_variables()
            tab_controller = SbSGuiMatchResultController(variables_for_beam, variables_common)
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

    def run_matching(self):
        matchers_list = []
        for index in range(len(self._matchers_tabs)):
            matcher_model = self._matchers_tabs[index].model

        for matcher_tab in self._matchers_tabs:
            matcher_tab.model.set_ignore_vars_list(
                matcher_tab.results_controller.get_unselected_variables()
            )
            matchers_list.append(matcher_tab.model.get_matcher())
        input_data = sbs_general_matcher.InputData.init_from_matchers_list(
            self._lhc_mode, self._match_path, matchers_list
        )
        sbs_general_matcher.run_full_madx_matching(input_data)
        for index in range(len(self._matchers_tabs)):
            matcher_model = self._matchers_tabs[index].model
            figures = self._matchers_tabs[index].results_controller.get_figures()
            matcher_model.get_plotter(figures).plot()

    class Tab(object):
        def __init__(self, matcher_model, matcher_results_controller):
            self.model = matcher_model
            self.results_controller = matcher_results_controller


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
