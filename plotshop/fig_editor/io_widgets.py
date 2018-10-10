from PyQt5 import QtWidgets, QtGui, QtCore

# Column Selector ############################################################

class HLine(QtWidgets.QFrame):
    """ A horizontal line for use in a widget """
    def __init__(self):
        super(HLine, self).__init__()
        self.setFrameShape(self.HLine)
        self.setFrameShadow(self.Sunken)


class ColumnSelectorDialog(QtWidgets.QDialog):
    """ Dialogbox for selecting columns.

    The selected columns can be requested by get_selected_columns() after closing.
    """

    def __init__(self, columns, single_line=False, parent=None):
        super(ColumnSelectorDialog, self).__init__(parent)
        self._accepted = False  # for checking which way the window was closed
        self._available_columns = columns
        self._single_line = single_line

        if single_line:
            self._pre_cols = ["Label"]
        else:
            self._pre_cols = ["Legend"]
        self._pre_cols += ["Columns: ", "--NONE--"]

        layout = QtWidgets.QVBoxLayout()

        self._lines = []
        self._lines_layout = QtWidgets.QGridLayout()
        layout.addLayout(self._lines_layout)

        buttons = QtWidgets.QHBoxLayout()
        self._common_x = QtWidgets.QCheckBox('Common X-Column', self)
        buttons.addWidget(self._common_x)
        self._label_by_y = QtWidgets.QCheckBox('{} by Y-Column'.format(self._pre_cols[0]), self)
        self._label_by_y.stateChanged.connect(self._on_label_by_y_change)
        buttons.addWidget(self._label_by_y)

        self._accept = QtWidgets.QPushButton('OK', self)
        self._accept.clicked.connect(self._on_accept)
        buttons.addWidget(self._accept)

        self._abort = QtWidgets.QPushButton('Cancel', self)
        self._abort.clicked.connect(self._on_abort)
        buttons.addWidget(self._abort)

        layout.addLayout(buttons)
        self.setLayout(layout)

        self._add_column_header()
        self._add_new_line()
        self.exec_()

    # Layout #################################################################

    def _add_column_header(self):
        grid = self._lines_layout
        all_cols = self._pre_cols + self._available_columns
        for idx, col in enumerate(all_cols):
            txt = QtWidgets.QLabel(col)
            grid.addWidget(txt, 0, idx)

        self._add_separator()

    def _add_separator(self):
        grid = self._lines_layout
        all_cols = self._pre_cols + self._available_columns
        row = grid.rowCount()
        grid.addWidget(HLine(), row, 0, 1, len(all_cols))

    def _add_new_line(self):
        grid = self._lines_layout
        idx_line = len(self._lines)
        line_dict = {}
        row = grid.rowCount()
        col_offset = len(self._pre_cols)

        label = QtWidgets.QLineEdit()
        label.setText("Line {:d}".format(idx_line + 1))
        label.textEdited[str].connect(lambda str, il=idx_line: self._on_label_change(il, str))
        line_dict["label"] = label

        grid.addWidget(label, row + 1, 0)

        for idx_type, type in enumerate(["X-Values", "Y-Values", "Error-Values"]):
            # add data type
            grid.addWidget(QtWidgets.QLabel(type), row + idx_type, 1)

            # add buttons
            button_group = QtWidgets.QButtonGroup(self)
            button_group.buttonClicked[int].connect(
                lambda idx, il=idx_line, it=idx_type: self._on_box_checked(il, it, idx)
            )

            button_group = self._new_button(button_group, id=-10,
                                            row=row+idx_type, col=2,
                                            tooltip="Uncheck", checked=True)

            for idx_col, col in enumerate(self._available_columns):
                button_group = self._new_button(button_group, id=idx_col,
                                                row=row+idx_type, col=col_offset+idx_col,
                                                tooltip="Set {} to {}".format(type, col))
            line_dict[type[0]] = button_group
        self._lines.append(line_dict)
        self._add_separator()

    # Listeners ##############################################################

    def _on_label_by_y_change(self, checked):
        if checked:
            for idx_line, line in enumerate(self._lines):
                self._set_y_label_to_column(idx_line, line["Y"].checkedId())

    def _on_label_change(self, idx_line, label):
        if not self._single_line and idx_line == len(self._lines) - 1:
            self._add_new_line()

    def _on_box_checked(self, idx_line, idx_type, idx_button):
        if not self._single_line and idx_line == len(self._lines) - 1:
            self._add_new_line()

        if self._common_x.isChecked():
            if idx_type == 0:
                self._set_all_x_to(idx_button)
            elif idx_type == 1:
                line = self._lines[idx_line]
                if line["X"].checkedId() < 0:
                    common_x_id = self._lines[0]["X"].checkedId()
                    line["X"].button(common_x_id).setChecked(True)

        if self._label_by_y.isChecked() and idx_type == 1:
            self._set_y_label_to_column(idx_line, idx_button)

    def _on_accept(self):
        self._accepted = True
        self.close()

    def _on_abort(self):
        self.close()

    def closeEvent(self, event):
        """ When the window is closed.

        Done this way to count closing by x-button as abort
        """
        if not self._accepted:
            self._lines = []

    # Other ##################################################################

    def _new_button(self, group, id, row, col, tooltip, checked=False):
        grid = self._lines_layout
        box = QtWidgets.QCheckBox()
        box.setToolTip(tooltip)
        box.setChecked(checked)
        group.addButton(box)
        group.setId(box, id)
        grid.addWidget(box, row, col)
        grid.setAlignment(box, QtCore.Qt.AlignHCenter)
        return group

    def _set_all_x_to(self, idx):
        for line in self._lines[:-1]:
            line["X"].button(idx).setChecked(True)

    def _set_y_label_to_column(self, idx_line, idx_button):
        if idx_button >= 0:
            self._lines[idx_line]["label"].setText(self._available_columns[idx_button])

    # Public #################################################################

    def get_selected_columns(self):
        lines = []
        cols = self._available_columns
        for line in self._lines:
            idx_x = line["X"].checkedId()
            idx_y = line["Y"].checkedId()
            idx_e = line["E"].checkedId()
            if idx_x >= 0 and idx_y >= 0:
                label = line["label"].text()
                error = ""
                if idx_e >= 0:
                    error = cols[idx_e]

                lines.append({
                    "l": label,
                    "x": cols[idx_x],
                    "y": cols[idx_y],
                    "e": error,
                })
        return lines
