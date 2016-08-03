import sys
import constants
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

    def set_selected_file(self, text):
        self._text_area.setText(text)

    def _select_file_action(self):
        file_selection_dialog = QtGui.QFileDialog()
        self._text_area.setText(file_selection_dialog.getExistingDirectory())


class FileSelectionComboWidget(QtGui.QWidget):

    added_items = []

    def __init__(self, label_text="", initial_list=[], parent=None):
        super(FileSelectionComboWidget, self).__init__(parent)
        layout = QtGui.QHBoxLayout(self)
        label = QtGui.QLabel(label_text)
        layout.addWidget(label)
        self._items_combo = QtGui.QComboBox()
        layout.addWidget(self._items_combo)
        self._items_combo.addItems(initial_list + self.added_items)
        self._select_button = QtGui.QPushButton("Add...")
        self._select_button.clicked.connect(self._select_file_action)
        layout.addWidget(self._select_button)

    def get_selected_file(self):
        return self._items_combo.currentText()

    def add_item(self, text):
        self.added_items.append(text)
        self._items_combo.insertItem(0, text)
        self._items_combo.setCurrentIndex(0)

    def reload_items(self, text_list):
        self._items_combo.clear()
        self._items_combo.addItem("")
        for text in text_list:
            self._items_combo.addItem(text)
        self._items_combo.setCurrentIndex(0)

    def _select_file_action(self):
        file_selection_dialog = QtGui.QFileDialog()
        self.add_item(file_selection_dialog.getExistingDirectory())


class InitialConfigPopup(QtGui.QDialog):

    def __init__(self, lhc_mode, match_path, parent=None):
        super(InitialConfigPopup, self).__init__(parent)
        self.resize(655, 90)
        main_layout = QtGui.QVBoxLayout(self)

        self._lhc_mode_combo = QtGui.QComboBox()
        self._lhc_mode_combo.addItems(constants.LHC_MODES)
        if lhc_mode is not None:
            if lhc_mode not in constants.LHC_MODES:
                raise ValueError("Invalid lhc mode, must be one of " +
                                 str(constants.LHC_MODES))
            else:
                self._lhc_mode_combo.setCurrentIndex(constants.LHC_MODES.index(lhc_mode))
        main_layout.addWidget(self._lhc_mode_combo)

        self._file_selector = FileSelectionDialogWidget()
        main_layout.addWidget(self._file_selector)
        if match_path is not None:
            self._file_selector.set_selected_file(match_path)

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
