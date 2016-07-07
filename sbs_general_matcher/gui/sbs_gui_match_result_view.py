import sys
from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


class SbSGuiMatchResultView(QtGui.QWidget):

    def __init__(self, controller, parent=None):
        super(SbSGuiMatchResultView, self).__init__(parent)

        self._controller = controller
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
        beam1_vars_frame = _BorderedGroupBox("Beam 1")
        beam1_vars_frame.setLayout(self._beam1_vars_layout)
        variables_layout.addWidget(beam1_vars_frame)

        self._beam2_vars_layout = QtGui.QVBoxLayout()
        beam2_vars_frame = _BorderedGroupBox("Beam 2")
        beam2_vars_frame.setLayout(self._beam2_vars_layout)
        variables_layout.addWidget(beam2_vars_frame)

        self._common_vars_layout = QtGui.QVBoxLayout()
        common_vars_frame = _BorderedGroupBox("Common")
        common_vars_frame.setLayout(self._common_vars_layout)
        variables_layout.addWidget(common_vars_frame)

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
        self.setStyleSheet(_BorderedGroupBox.GROUP_BOX_STYLE)


class SbSGuiMatchResultController(object):

    def __init__(self):
        self._view = SbSGuiMatchResultView(self)

    def get_view(self):
        return self._view


if __name__ == "__main__":
    print >> sys.stderr, "This module is meant to be imported."
    sys.exit(-1)
