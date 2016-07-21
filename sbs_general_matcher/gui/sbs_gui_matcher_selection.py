import sys
from PyQt4 import QtGui
from matchers_views.matcher_view_phase import MatcherControllerPhase
from matchers_views.matcher_view_kmod import MatcherControllerKmod
from matchers_views.matcher_view_amp import MatcherControllerAmp
from matchers_views.matcher_view_coupling import MatcherControllerCoupling


class SbSGuiMatcherSelection(QtGui.QDialog):

    WINDOW_TITLE = "Select matcher to add"
    MATCHER_TYPES = {
        "phase": MatcherControllerPhase,
        "kmod": MatcherControllerKmod,
        "beta from amplitude": MatcherControllerAmp,
        "coupling": MatcherControllerCoupling,
    }

    def __init__(self, main_controller, parent=None):
        super(SbSGuiMatcherSelection, self).__init__(parent)
        self._main_controller = main_controller
        self._build_gui()
        self._selected_matcher = None

    def _build_gui(self):
        self.setWindowTitle(SbSGuiMatcherSelection.WINDOW_TITLE)
        self.resize(320, 240)

        dialog_layout = QtGui.QVBoxLayout(self)
        buttons_layout = QtGui.QHBoxLayout()

        self._matchers_list_widget = SbSGuiMatcherSelection._get_matchers_list_widget()
        dialog_layout.addWidget(self._matchers_list_widget)

        add_button = QtGui.QPushButton("Add", self)

        add_button.clicked.connect(self._accept_action)
        buttons_layout.addWidget(add_button)

        cancel_button = QtGui.QPushButton("Cancel", self)
        cancel_button.clicked.connect(self.reject)
        buttons_layout.addWidget(cancel_button)

        dialog_layout.addLayout(buttons_layout, stretch=1)

    @staticmethod
    def _get_matchers_list_widget():
        matchers_list_widget = QtGui.QListWidget()
        for key in SbSGuiMatcherSelection.MATCHER_TYPES:
            matchers_list_widget.addItem(key)
        return matchers_list_widget

    def _accept_action(self):
        self._selected_matcher = str(self._matchers_list_widget.currentItem().text())
        matcher_controller_instance = SbSGuiMatcherSelection.MATCHER_TYPES[self._selected_matcher](self._main_controller)
        result_code = matcher_controller_instance.show_view()
        if result_code == QtGui.QDialog.Accepted:
            self._selected_matcher = matcher_controller_instance.get_selected_matcher_model()
            self.accept()

    def get_selected_matcher(self):
        return self._selected_matcher


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
