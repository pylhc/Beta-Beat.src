import logging
import os
import six
from PyQt5 import QtWidgets, QtCore, QtGui

import matplotlib
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT


from gui_utils import get_icon
import options_figure


# Log Dialog Window ##########################################################

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
DATE_FORMAT = '%d/%m/%Y %H:%M:%S'


class LogDialog(QtWidgets.QDialog, logging.StreamHandler):
    """ Logging dialog window, created by Jaime """

    _update_signal = QtCore.pyqtSignal(str)

    def __init__(self, fmt=LOG_FORMAT, datefmt=DATE_FORMAT, parent=None):
        QtWidgets.QDialog.__init__(self, parent=parent)
        logging.Handler.__init__(self)

        self.resize(855, 655)
        layout = QtWidgets.QVBoxLayout()
        self._log_text = QtWidgets.QPlainTextEdit(parent)
        self._log_text.setReadOnly(True)
        layout.addWidget(self._log_text)
        self.setLayout(layout)

        self.setFormatter(logging.Formatter(fmt, datefmt))
        self.setLevel(logging.DEBUG)

        self._update_signal.connect(self._update_log, QtCore.Qt.QueuedConnection)

    def _update_log(self, msg):
        self._log_text.appendPlainText(msg)

    def emit(self, record):
        msg = self.format(record)
        self._update_signal.emit(msg)


class LogStatusBar(QtWidgets.QStatusBar, logging.StreamHandler):
    """ Logging dialog window, created by Jaime """

    _update_signal = QtCore.pyqtSignal(str)

    def __init__(self, fmt=LOG_FORMAT, datefmt=DATE_FORMAT, parent=None):
        QtWidgets.QStatusBar.__init__(self, parent=parent)
        logging.Handler.__init__(self)

        self._status_text = QtWidgets.QLabel()
        self.addWidget(self._status_text, 1)

        self.setFormatter(logging.Formatter(fmt, datefmt))
        self.setLevel(logging.INFO)

        self._update_signal.connect(self.showMessage, QtCore.Qt.QueuedConnection)

    def emit(self, record):
        msg = self.format(record)
        self._update_signal.emit(msg)


# Navigation Toolbar #########################################################

ICON_SIZE_NAVTOOLBAR = 24


class NavigationToolbar(NavigationToolbar2QT):
    """ Customized Navigation Toolbar """

    LEGEND_LOCATIONS = [1, 2, 3, 4, None]

    def __init__(self, canvas, save_fun=None, load_fun=None, export_fun=None, import_fun=None, parent=None):
        self.load_data = load_fun
        self.export_plot = export_fun
        self.import_plot = import_fun

        self.toolitems = list(self.toolitems)

        if save_fun is not None:
            self.save_figure = save_fun  # otherwise defined by super-class

        if self.load_data is not None:
            self.toolitems.append((
                    "Load Data", "Plot data from tfs files.",
                    "folder_open", "load_data"
                ))

        if self.export_plot is not None:
            self.toolitems.append((
                    "Export Plot", "Export plot to .dly file format.",
                    "export", "export_plot"
                ))

        if self.import_plot is not None:
            self.toolitems.append((
                    "Import Plot", "Import plot from .dly file format.",
                    "import", "import_plot"
                ))

        super(NavigationToolbar, self).__init__(canvas, parent)
        self._legend_location_index = 0
        self.figure = canvas.figure
        self.axes = canvas.figure.gca()

    def _init_toolbar(self):
        """ Called from the super function """
        self.basedir = os.path.join(matplotlib.rcParams['datapath'], 'images')

        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.addSeparator()
            else:
                try:
                    icon = get_icon(image_file)
                except IOError:
                    icon = self._icon(image_file + '.png')

                a = self.addAction(icon,
                                   text, getattr(self, callback))
                self._actions[callback] = a
                if callback in ['zoom', 'pan']:
                    a.setCheckable(True)
                if tooltip_text is not None:
                    a.setToolTip(tooltip_text)
                if text == 'Subplots':
                    a = self.addAction(self._icon("qt4_editor_options.png"),
                                       'Customize', self.edit_parameters)
                    a.setToolTip('Edit axis, curve and image parameters')

        self.buttons = {}

        # Add the x,y location widget at the right side of the toolbar
        # The stretch factor is 1 which means any resizing of the toolbar
        # will resize this label instead of the buttons.
        if self.coordinates:
            self.locLabel = QtWidgets.QLabel("", self)
            self.locLabel.setAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignTop)
            self.locLabel.setSizePolicy(
                QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                      QtWidgets.QSizePolicy.Ignored))
            labelAction = self.addWidget(self.locLabel)
            labelAction.setVisible(True)

        # reference holder for subplots_adjust window
        self.adj_window = None

        # Aesthetic adjustments - we need to set these explicitly in PyQt5
        # otherwise the layout looks different - but we don't want to set it if
        # not using HiDPI icons otherwise they look worse than before.

        self.setIconSize(QtCore.QSize(ICON_SIZE_NAVTOOLBAR, ICON_SIZE_NAVTOOLBAR))
        self.layout().setSpacing(12)

    def move_legend(self, *args):
        self._legend_location_index = (
            self._legend_location_index + 1
        ) % len(NavigationToolbar.LEGEND_LOCATIONS)
        for axes in self.axes:
            legend = axes.get_legend()
            if legend is not None:
                loc = NavigationToolbar.LEGEND_LOCATIONS[
                    self._legend_location_index
                ]
                if loc is None:
                    legend.set_visible(False)
                else:
                    legend.set_visible(True)
                    legend._set_loc(loc)
        self.draw()

    def edit_parameters(self):
        allaxes = self.canvas.figure.get_axes()
        if not allaxes:
            QtWidgets.QMessageBox.warning(
                self.parent, "Error", "There are no axes to edit.")
            return
        elif len(allaxes) == 1:
            axes, = allaxes
        else:
            titles = []
            for axes in allaxes:
                name = (axes.get_title() or
                        " - ".join(filter(None, [axes.get_xlabel(),
                                                 axes.get_ylabel()])) or
                        "<anonymous {} (id: {:#x})>".format(
                            type(axes).__name__, id(axes)))
                titles.append(name)
            item, ok = QtWidgets.QInputDialog.getItem(
                self.parent, 'Customize', 'Select axes:', titles, 0, False)
            if ok:
                axes = allaxes[titles.index(six.text_type(item))]
            else:
                return

        options_figure.figure_edit(axes, self)
