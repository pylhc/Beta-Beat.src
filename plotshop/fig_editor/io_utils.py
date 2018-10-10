import os
import six
import matplotlib
import logging

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog

from io_widgets import ColumnSelectorDialog
from tfs_files import tfs_pandas as tfs
from plotshop import plot_tfs

LOG = logging.getLogger(__name__)


# Public Methods #############################################################


def load_tfs():
    """ Plots a tfs-file.

    TODO:
        * check all tfs files for common columns -> make user choose which column
        * changemarkers tickbox
    """
    LOG.debug("Load Tfs clicked.")
    paths = QFileDialog().getOpenFileNames(None, 'Load file(s)', '')[0]
    fig = None
    if paths:
        LOG.info("Files chosen: {:s}".format(", ".join(paths)))
        if len(paths) > 1:
            # load all files and check for common columns
            df_list, common_cols = _get_all_tfs_and_common_columns(paths)
            column_selector = ColumnSelectorDialog(common_cols, single_line=True)
        elif len(paths) == 1:
            # load only one tfs
            LOG.debug("Loading only one file")
            try:
                df = tfs.read_tfs(paths[0])
            except tfs.TfsFormatError:
                LOG.error("File '{}' is not of TFS format!".format(paths[0]))
            else:
                column_selector = ColumnSelectorDialog(df.columns.tolist())
                selected = column_selector.get_selected_columns()
                if selected:
                    fig = plot_tfs.plot_single_file(
                        files=paths,
                        x_cols=[s["x"] for s in selected],
                        y_cols=[s["y"] for s in selected],
                        e_cols=[s["e"] for s in selected],
                        labels=[s["l"] for s in selected],
                        no_show=True,
                    )
        return fig
    LOG.debug("No files chosen.")
    return None


def load_dly():
    """ Load data and layout from .dly file """
    path = QFileDialog().getOpenFileName(None, 'Load Figure', '', '*.dly')[0]
    if path:
        LOG.info("Importing file '{}'".format(path))
        with open(path, "rb") as f:
            pass # TODO
    return None


def save_dly(fig):
    """ Save data and layout into .dly file """
    path = QFileDialog().getSaveFileName(None, "Save Figure", "", "*.dly")[0]
    if path:
        with open(path, "wb") as f:
            pass  #TODO

def save_layout(figure):
    """ Saves the current state of the layout.

    Into ini-file
    * Save labels
    * Save Title
    * Save x_lim
    * Save y_lim
    * Gridstatus
    * axis-scaling
    * the subplot parameters


    Returns:

    """
    pass


def load_layout(figure):
    """ Load ini-file with layout properties

    """


def save_figure(figure):
    """ Save figure

     From backend_qt5.NavigationToolbar2QT
     """
    filetypes = figure.canvas.get_supported_filetypes_grouped()
    sorted_filetypes = sorted(six.iteritems(filetypes))
    default_filetype = figure.canvas.get_default_filetype()

    startpath = os.path.expanduser(
        matplotlib.rcParams['savefig.directory'])
    start = os.path.join(startpath, figure.canvas.get_default_filename())
    filters = []
    selected_filter = None

    for name, exts in sorted_filetypes:
        exts_list = " ".join(['*.%s' % ext for ext in exts])
        filter = '%s (%s)' % (name, exts_list)
        if default_filetype in exts:
            selected_filter = filter
        filters.append(filter)
    filters = ';;'.join(filters)

    fname, filter = QtWidgets.QFileDialog.getSaveFileName(
        caption="Choose a filename to save to",
        directory=start,
        filter=filters,
        initialFilter=selected_filter,
    )

    if fname:
        # Save dir for next time, unless empty str (i.e., use cwd).
        if startpath != "":
            matplotlib.rcParams['savefig.directory'] = (
                os.path.dirname(six.text_type(fname)))
        try:
            figure.savefig(six.text_type(fname))
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                None, "Error saving file", six.text_type(e),
                QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)




# Private Methods #############################################################


def _get_all_tfs_and_common_columns(paths):
    tfs_list = [tfs.read_tfs(p) for p in paths]
    cols = tfs_list[0].columns
    for t in tfs_list[1:]:
        cols = cols.intersection(t.columns)
    return tfs_list, cols
