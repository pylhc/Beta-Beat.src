import sys
from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


class SbSGuiMatchResultView(QtGui.QWidget):

    def __init__(self, controller, variables_for_beam, variables_common, parent=None):
        super(SbSGuiMatchResultView, self).__init__(parent)

        self._controller = controller
        self._variables_for_beam = variables_for_beam
        self._variables_common = variables_common
        self._build_gui()

    def _build_gui(self):
        main_layout = QtGui.QHBoxLayout()

        beam1_layout = QtGui.QVBoxLayout()
        beam1_frame = _BorderedGroupBox("Beam 1")
        beam1_frame.setLayout(beam1_layout)

        self._beam1_upper_figure = plt.figure()
        beam1_layout.addLayout(self._get_new_canvas_layout(self._beam1_upper_figure))
        self._beam1_lower_figure = plt.figure()
        beam1_layout.addLayout(self._get_new_canvas_layout(self._beam1_lower_figure))
        main_layout.addWidget(beam1_frame)

        beam2_layout = QtGui.QVBoxLayout()
        beam2_frame = _BorderedGroupBox("Beam 2")
        beam2_frame.setLayout(beam2_layout)

        self._beam2_upper_figure = plt.figure()
        beam2_layout.addLayout(self._get_new_canvas_layout(self._beam2_upper_figure))
        self._beam2_lower_figure = plt.figure()
        beam2_layout.addLayout(self._get_new_canvas_layout(self._beam2_lower_figure))
        main_layout.addWidget(beam2_frame)

        variables_layout = QtGui.QVBoxLayout()
        variables_frame = _BorderedGroupBox("Variables")
        variables_frame.setLayout(variables_layout)

        self._beam1_vars_layout = QtGui.QVBoxLayout()
        if not len(self._variables_for_beam[1]) == 0:
            beam1_widget = QtGui.QWidget()
            beam1_widget.setLayout(self._beam1_vars_layout)
            beam1_scroll = QtGui.QScrollArea()
            beam1_scroll.setWidgetResizable(True)
            beam1_scroll.setWidget(beam1_widget)
            beam1_vars_frame = _BorderedGroupBox("Beam 1")
            beam1_vars_frame_layout = QtGui.QHBoxLayout()
            beam1_vars_frame_layout.addWidget(beam1_scroll)
            beam1_vars_frame.setLayout(beam1_vars_frame_layout)
            variables_layout.addWidget(beam1_vars_frame)
            for variable in self._variables_for_beam[1]:
                self._beam1_vars_layout.addWidget(QtGui.QCheckBox(variable))

        self._beam2_vars_layout = QtGui.QVBoxLayout()
        if not len(self._variables_for_beam[2]) == 0:
            beam2_widget = QtGui.QWidget()
            beam2_widget.setLayout(self._beam2_vars_layout)
            beam2_scroll = QtGui.QScrollArea()
            beam2_scroll.setWidgetResizable(True)
            beam2_scroll.setWidget(beam2_widget)
            beam2_vars_frame = _BorderedGroupBox("Beam 2")
            beam2_vars_frame_layout = QtGui.QHBoxLayout()
            beam2_vars_frame_layout.addWidget(beam2_scroll)
            beam2_vars_frame.setLayout(beam2_vars_frame_layout)
            variables_layout.addWidget(beam2_vars_frame)
            for variable in self._variables_for_beam[2]:
                self._beam2_vars_layout.addWidget(QtGui.QCheckBox(variable))

        self._common_vars_layout = QtGui.QVBoxLayout()
        if not len(self._variables_common) == 0:
            common_widget = QtGui.QWidget()
            common_widget.setLayout(self._common_vars_layout)
            common_scroll = QtGui.QScrollArea()
            common_scroll.setWidgetResizable(True)
            common_scroll.setWidget(common_widget)
            common_vars_frame = _BorderedGroupBox("Common")
            common_vars_frame_layout = QtGui.QHBoxLayout()
            common_vars_frame_layout.addWidget(common_scroll)
            common_vars_frame.setLayout(common_vars_frame_layout)
            variables_layout.addWidget(common_vars_frame)
            for variable in self._variables_common:
                self._common_vars_layout.addWidget(QtGui.QCheckBox(variable))

        select_all_checkbox = QtGui.QCheckBox("Toggle select all")
        select_all_checkbox.stateChanged.connect(self._toogle_select_all)
        variables_layout.addWidget(select_all_checkbox)

        main_layout.addWidget(variables_frame)

        self.setLayout(main_layout)

    def _get_new_canvas_layout(self, figure):
        layout = QtGui.QVBoxLayout()
        canvas = FigureCanvas(figure)
        toolbar = NavigationToolbar(canvas, self)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        return layout

    def get_figures(self):
        return ((self._beam1_upper_figure, self._beam1_lower_figure),
                (self._beam2_upper_figure, self._beam2_lower_figure))

    def get_unselected_variables(self):
        unselected_vars = []

        def add_text_to_list(checkbox):
            if not checkbox.isChecked():
                unselected_vars.append(str(checkbox.text()))

        self._loop_through_checkboxes(add_text_to_list)
        return unselected_vars

    def _toogle_select_all(self, state):
        checked = bool(state)

        def toggle_checkbox(checkbox):
            checkbox.setChecked(checked)

        self._loop_through_checkboxes(toggle_checkbox)

    def _loop_through_checkboxes(self, function):
        for layout in [self._beam1_vars_layout,
                       self._beam2_vars_layout,
                       self._common_vars_layout]:
            for index in range(layout.count()):
                checkbox = layout.itemAt(index).widget()
                if type(checkbox) is QtGui.QCheckBox:
                    function(checkbox)


class _BorderedGroupBox(QtGui.QGroupBox):

    GROUP_BOX_STYLE = """
        QGroupBox {
            border: 1px solid gray;
            border-radius: 3px;
        }
        QGroupBox::title {
            background-color: transparent;
            subcontrol-position: top left;
            padding:2 13px;
        }
    """

    def __init__(self, label, parent=None):
        super(_BorderedGroupBox, self).__init__(label, parent)
        # self.setStyleSheet(_BorderedGroupBox.GROUP_BOX_STYLE)


class SbSGuiMatchResultController(object):

    def __init__(self, variables_for_beam, variables_common):
        self._view = SbSGuiMatchResultView(self, variables_for_beam, variables_common)

    def get_view(self):
        return self._view

    def get_unselected_variables(self):
        return self._view.get_unselected_variables()

    def get_figures(self):
        return self._view.get_figures()


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
