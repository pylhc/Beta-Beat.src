import logging

from PyQt5.QtWidgets import QFileDialog

from io_widgets import ColumnSelectorDialog
from tfs_files import tfs_pandas as tfs
from utils.plotting import plot_tfs

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

def save_layout(self):
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


def load_layout(self):
    """ Load ini-file with layout properties

    """

# Private Methods #############################################################


def _get_all_tfs_and_common_columns(paths):
    tfs_list = [tfs.read_tfs(p) for p in paths]
    cols = tfs_list[0].columns
    for t in tfs_list[1:]:
        cols = cols.intersection(t.columns)
    return tfs_list, cols
