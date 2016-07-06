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


class FileSelectionPopup(QtGui.QDialog):
    def __init__(self, parent=None):
        super(FileSelectionPopup, self).__init__(parent)
        self.resize(655, 90)
        main_layout = QtGui.QVBoxLayout(self)
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
        return self._file_selector.get_selected_file()

if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
