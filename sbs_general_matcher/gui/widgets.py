import sys
from PyQt4 import QtGui


class FileSelectionDialogWidget(QtGui.QWidget):
    def __init__(self, label_text="", parent=None):
        super(FileSelectionDialogWidget, self).__init__(parent)
        layout = QtGui.QHBoxLayout(self)
        label = QtGui.QLabel(label_text)
        layout.addWidget(label)
        self._text_area = QtGui.QLineEdit()
        layout.addWidget(self._text_area)
        self._select_button = QtGui.QPushButton("Choose...")
        self._select_button.clicked.connect(self._select_file_action)
        layout.addWidget(self._select_button)

    def get_selected_file(self):
        return self._text_area.text()

    def _select_file_action(self):
        file_selection_dialog = QtGui.QFileDialog()
        self._text_area.setText(file_selection_dialog.getExistingDirectory())


class InitialConfigPopup(QtGui.QDialog):

    LHC_MODES = ["lhc_runII_2016", "lhc_runII_2016_ats", "lhc_runII", "lhc_runI", "hllhc"]

    def __init__(self, parent=None):
        super(InitialConfigPopup, self).__init__(parent)
        self.resize(655, 90)
        main_layout = QtGui.QVBoxLayout(self)

        self._lhc_mode_combo = QtGui.QComboBox()
        self._lhc_mode_combo.addItems(InitialConfigPopup.LHC_MODES)
        main_layout.addWidget(self._lhc_mode_combo)

        self._file_selector = FileSelectionDialogWidget()
        main_layout.addWidget(self._file_selector)

        buttons_layout = QtGui.QHBoxLayout()
        accept_button = QtGui.QPushButton("Accept")
        accept_button.clicked.connect(self.accept)
        buttons_layout.addWidget(accept_button)
        cancel_button = QtGui.QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        buttons_layout.addWidget(cancel_button)
        main_layout.addLayout(buttons_layout)

    def get_selected_file(self):
        return str(self._file_selector.get_selected_file())

    def get_selected_lhc_mode(self):
        return str(self._lhc_mode_combo.currentText())


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
