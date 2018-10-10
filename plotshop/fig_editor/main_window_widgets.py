import logging
import os
import time
from PyQt5 import QtWidgets, QtCore

import matplotlib
import six
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

import options_figure
import options_utils
from gui_utils import get_icon

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


class NavigationToolbar(NavigationToolbar2QT):
    """ Customized Navigation Toolbar """

    ICON_SIZE = 24

    def __init__(self, canvas,
                 save_fun=None, load_fun=None, export_fun=None, import_fun=None, parent=None):
        self.load_data = load_fun
        self.export_plot = export_fun
        self.import_plot = import_fun

        self.toolitems = list(self.toolitems)

        if hasattr(canvas, "move_legend") and hasattr(canvas, "update_legend"):
            self.update_legend = canvas.update_legend
            self.move_legend = canvas.move_legend

            self.toolitems.insert(7, (None, None, None, None))

            self.toolitems.insert(8, (
                "Move Legend", "Move the legend location to predefined settings.",
                "arrows", "move_legend"
            ))

            self.toolitems.insert(9, (
                "Update Legend", "Update the legend.",
                "refresh", "update_legend"
            ))

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

        self.setIconSize(QtCore.QSize(NavigationToolbar.ICON_SIZE, NavigationToolbar.ICON_SIZE))
        self.layout().setSpacing(12)

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


# Figure Canvas ################################################################


class FigureCanvasExt(FigureCanvas):
    """ Extended FigureCanvas.

    TODO: Multiaxes
    - Legend positions can be easily saved and restored

    """

    LEGEND_LOCATIONS = [0, 1, 2, 3, 4, None]

    def __init__(self, figure):
        super(FigureCanvasExt, self).__init__(figure)
        self._legend_location_index = 0
        self._legend_locations = list(FigureCanvasExt.LEGEND_LOCATIONS)

    def move_legend(self):
        self._legend_location_index = (
                                              self._legend_location_index + 1
                                      ) % len(self._legend_locations)
        self._set_legend_loc()

    def update_legend(self):
        axes = self.figure.gca()
        options_utils.regenerate_legend(axes, force_new=True)
        self._save_legend_loc()
        self.draw()

    def _set_legend_loc(self):
        legend = self.figure.gca().get_legend()
        if legend is not None:
            loc = self._legend_locations[
                self._legend_location_index
            ]
            if loc is None:
                legend.set_visible(False)
            else:
                legend.set_visible(True)
                legend._set_loc(loc)

            self.draw()

    def _save_legend_loc(self):
        legend = self.figure.gca().get_legend()
        if legend:
            self._legend_location_index = 0
            self._legend_locations = ([legend._get_loc()] +
                                      list(FigureCanvasExt.LEGEND_LOCATIONS))

    def update_figure(self, figure):
        """ Change Figure for this canvas """
        figure.canvas = self
        self.figure = figure
        self._save_legend_loc()
        self._set_legend_loc()
        self.draw()

    def set_pickers(self, tol=0.5):
        """ Set all artists to send pick events when they are clicked """
        axes = self.figure.get_axes()
        for ax in axes:
            for l in ax.lines:
                width = max(1.1*l.get_markersize(), 1.5*l.get_linewidth())
                l.set_picker(width)

            for a in [ax.xaxis, ax.yaxis]:
                a.set_picker(tol)
                for c in a.get_children():
                    c.set_picker(True)

            for t in ax.texts:
                t.set_picker(True)

            legend = ax.get_legend()
            if legend is not None:
                legend.set_picker(True)

            ax.title.set_picker(True)
            # ax.set_picker(True)

    def unset_pickers(self):
        for c in self.figure.get_children():
            try:
                c.set_picker(False)
                c.set_draggable(False)
            except AttributeError:
                pass


# Figure Canvas ################################################################


class DragHandler(object):
    """ Handles the dragging of objects around the screen """
    time_tol = 0.1  # time between drag'n'drop to be registered (in sec)

    def __init__(self, artist, mouseevent):
        self.artist = artist
        self.time = time.time()
        self.pick_pos_data = [mouseevent.xdata, mouseevent.ydata]
        self.pick_pos_disp = [mouseevent.x, mouseevent.y]
        self.artist_pos_data = artist.get_position()
        self.artist_pos_disp = artist._get_xy_display()

    def move(self, mouseevent):
        # if (time.time() - self.time) < self.time_tol:
        #     return

        mouse_pos_disp = [mouseevent.x, mouseevent.y]
        mouse_pos_data = [mouseevent.xdata, mouseevent.ydata]

        if any(pos is None for pos in self.pick_pos_data + mouse_pos_data):
            # something was outside of the plotting area -> use figure coodinates
            return  # TODO: DOES NOT WORK AS EXPECTED
            to_pos = self._get_new_position(self.artist_pos_disp,
                                            mouse_pos_disp,
                                            self.pick_pos_disp)
            self.artist.transform = self.artist.figure.transFigure
            to_pos = self.artist.transform.inverted().transform(to_pos)
        else:
            # keep everything inside the plotting area
            to_pos = self._get_new_position(self.artist_pos_data,
                                            mouse_pos_data,
                                            self.pick_pos_data)

        self.artist.set_position(to_pos)

    def _get_new_position(self, from_pos, mouse_pos, pick_pos):
        return [fp + mp - pp for fp, mp, pp in zip(from_pos, mouse_pos, pick_pos)]