"""
Module plotshop.plot_tfs
-------------------------

Wrapper to easily plot tfs-files. With entrypoint functionality.
"""


import sys
from os.path import abspath, join, pardir
sys.path.append(abspath(join(__file__, pardir, pardir)))

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt, gridspec, rcParams
import numpy as np

from utils.entrypoint import EntryPointParameters, entrypoint
from plotshop import plot_style as ps
from tfs_files import tfs_pandas as tfs
from utils import logging_tools


LOG = logging_tools.get_logger(__name__)


# Constants, Style and Arguments #############################################


def get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--files",
        help="Twiss files to plot",
        name="files",
        required=True,
        nargs="+",
        type=basestring,
    )
    params.add_parameter(
        flags=["-y", "--y_cols"],
        help="List of column names to plot (e.g. BETX, BETY)",
        name="y_cols",
        required=True,
        type=basestring,
        nargs="+",
    )
    params.add_parameter(
        flags=["-x", "--x_cols"],
        help="List of column names to use as x-values.",
        name="x_cols",
        type=basestring,
        nargs="+",
    )
    params.add_parameter(
        flags=["-e", "--e_cols"],
        help="List of parameters to get error values from.",
        name="e_cols",
        type=basestring,
        nargs="+",
    )
    params.add_parameter(
        flags="--labels",
        help="Y-Lables for the plots, default: y_col.",
        name="labels",
        type=basestring,
        nargs="+",
    )
    params.add_parameter(
        flags="--source_names",
        help="Names for the sources for the plots, default: filenames.",
        name="source_names",
        type=basestring,
        nargs="+",
    )
    params.add_parameter(
        flags="--output",
        help="Base-Name of the output files. _'y_col'.pdf will be attached.",
        name="output",
        type=basestring,
    )
    params.add_parameter(
        flags="--changemarker",
        help="Changes marker for each line in the plot.",
        action="store_true",
        name="change_marker",
    )
    params.add_parameter(
        flags="--nolegend",
        help="Deactivates the legend.",
        action="store_true",
        name="no_legend",
    )
    params.add_parameter(
        flags="--noshow",
        help="Suppresses opening plotting windows.",
        action="store_true",
        name="no_show",
    )
    params.add_parameter(
        flags="--xy",
        help="Plots X and Y for the give parameters into one figure (two axes).",
        action="store_true",
        name="xy",
    )
    params.add_parameter(
        flags="--autoscale",
        help="Scales the plot, so that this percentage of points is inside the picture.",
        type=float,
        name="auto_scale",
    )
    params.add_parameter(
        flags="--singlefig",
        help="Plots into one figure. (Works only with one file so far).",
        action="store_true",
        name="single_fig",
    )
    return params


IR_POS_DEFAULT = {
    "LHCB1": {
        'IP1': 23519.36962,
        'IP2': 192.923,
        'IP3': 3525.207216,
        'IP4': 6857.491433,
        'IP5': 10189.77565,
        'IP6': 13522.21223,
        'IP7': 16854.64882,
        'IP8': 20175.8654,
    },
    "LHCB2": {
        'IP1': 3195.252584,
        'IP2': 6527.5368,
        'IP3': 9859.973384,
        'IP4': 13192.40997,
        'IP5': 16524.84655,
        'IP6': 19857.13077,
        'IP7': 23189.41498,
        'IP8': 26510.4792,
    }
}

MANUAL_STYLE = {
    # differences to the standard style
    u'lines.markersize': 5.0,
    u'lines.linestyle': u'',
}

ERROR_ALPHA = 1.  # Set errorbar transparency
MAX_LEGENDLENGTH = 78  # maximum length of legend letters before linebreak
COMPLEX_NAMES = [p+ext for p in ["1001", "1010"] for ext in "RI"]  # Endings of columns that contain complex data

# Main invocation ############################################################


@entrypoint(get_params(), strict=True)
def plot(opt):
    """ Plots data from different twiss-input files into one plot.

    Keyword Args:
        Required
        files (basestring): Twiss files to plot
                            **Flags**: --files
        y_cols (basestring): List of column names to plot (e.g. BETX, BETY)
                             **Flags**: ['-y', '--y_cols']
        Optional
        auto_scale (float): Scales the plot, so that this percentage of
                            points is inside the picture.
                            **Flags**: --autoscale
        change_marker: Changes marker for each line in the plot.
                       **Flags**: --changemarker
                       **Action**: ``store_true``
        e_cols (basestring): List of parameters to get error values from.
                             **Flags**: ['-e', '--e_cols']
        labels (basestring): Y-Lables for the plots, default: y_col.
                             **Flags**: --labels
        no_legend: Deactivates the legend.
                   **Flags**: --nolegend
                   **Action**: ``store_true``
        no_show: Suppresses opening plotting windows.
                 **Flags**: --noshow
                 **Action**: ``store_true``
        output (basestring): Base-Name of the output files. _'y_col'.pdf will be attached.
                             **Flags**: --output
        single_fig (bool): Plots into one figure. (Works only with one file so far).
                           **Flags**: --singlefig
                           **Action**: ``store_true``

        source_names (basestring): Names for the sources for the plots, default: filenames.
                                   **Flags**: --source_names
        x_cols (basestring): List of column names to use as x-values.
                             **Flags**: ['-x', '--x_cols']
        xy (bool): Plots X and Y for the give parameters into one figure (two axes).
                   **Flags**: --xy
                   **Action**: ``store_true``
    """
    LOG.debug("Starting plotting of tfs files: {:s}".format(", ".join(opt.files)))

    # preparations
    opt = _check_opt(opt)

    # extract data
    twiss_data = _get_data(opt.files)

    # plotting
    figs = _create_plots(opt.x_cols, opt.y_cols, opt.e_cols, twiss_data,
                         opt.source_names, opt.labels,
                         opt.xy, opt.change_marker, opt.no_legend, opt.auto_scale, opt.single_fig)

    # exports
    if opt.output:
        _export_plots(figs, opt.output)

    if not opt.no_show:
        plt.show()

    if opt.single_fig:
        return figs[0]
    return figs


# Private Functions ##########################################################


def _get_data(files):
    """ Load all data from files """
    return [tfs.read_tfs(f) for f in files]


class _LoopGenerator:
    """ Takes care of the loop order and collects the figures. """

    def __call__(self):
        # just here so the IDE does not complain, reassigned below
        pass

    def __init__(self, x_cols, y_cols, e_cols, datas, source_names, labels, xy, single_fig):
        # self.gs = _get_gridspec(xy)
        self.figs = None
        self.x_cols = x_cols
        self.y_cols = y_cols
        self.e_cols = e_cols
        self.datas = datas
        self.source_names = source_names
        self.labels = labels
        self.xy = xy
        self.single_fig = single_fig

        if single_fig:
            self.__call__ = self._do_sigle_fig
        else:
            self.__call__ = self._do_multi_fig

    def get_figs(self):
        return self.figs

    @staticmethod
    def _get_current_axes(axs, idx_plot):
        try:
            return axs[idx_plot]
        except TypeError:
            return axs

    def _get_fig(self):
        try:
            return plt.subplots(1+self.xy, 1, constrained_layout=False)  # matplotlib>=2.2
        except TypeError:
            return plt.subplots(1+self.xy, 1)

    def _do_sigle_fig(self):
        """ single figure: loop over the two axes as outer loop """
        fig, axs = self._get_fig()
        fig.canvas.set_window_title("File '{:s}'".format(self.source_names[0]))
        self.figs = [fig]

        # define y_label
        y_label = ""
        if len(self.labels) == 1:
            y_label = self.labels[0]
        for idx_plot in range(1 + self.xy):
            ax = self._get_current_axes(axs, idx_plot)
            for idx_col in range(len(self.x_cols)):
                yield (ax, idx_plot, idx_col, self.datas[0],
                       self.x_cols[idx_col], self.y_cols[idx_col], self.e_cols[idx_col],
                       self.labels[idx_col], y_label,
                       idx_col == (len(self.x_cols)-1),
                       )

    def _do_multi_fig(self):
        """ a figure for each parameter: loop over the columns as outer loop """
        self.figs = {}
        for idx_col in range(len(self.x_cols)):
            fig, axs = self._get_fig()
            p_title = self.labels[idx_col]
            if self.xy:
                p_title += "_dualPlane"
            fig.canvas.set_window_title("Parameter '{:s}'".format(p_title))
            self.figs[p_title] = fig
            for idx_plot in range(1 + self.xy):
                ax = self._get_current_axes(axs, idx_plot)
                for idx_data in range(len(self.datas)):
                    yield (ax, idx_plot, idx_data, self.datas[idx_data],
                           self.x_cols[idx_col], self.y_cols[idx_col], self.e_cols[idx_col],
                           self.source_names[idx_data], self.labels[idx_col],
                           idx_data == (len(self.datas)-1)
                           )

    def get_ncols(self):
        """ Returns the number of columns for the legend.

        Done here, as this class divides single-plot from multiplot anyway
        """
        names = self.labels if self.single_fig else self.source_names
        return max([int(np.floor(MAX_LEGENDLENGTH/max([len(l) for l in names]))), 1])


def _create_plots(x_cols, y_cols, e_cols, datas, source_names, labels,
                  xy, change_marker, no_legend, auto_scale, single_fig=False):
    # create layout
    ps.set_style("standard", MANUAL_STYLE)
    ir_positions, x_is_position = _get_ir_positions(datas, x_cols)

    y_lims = None
    the_loop = _LoopGenerator(x_cols, y_cols, e_cols, datas, source_names, labels, xy, single_fig)
    for ax, idx_plot, idx, data, x_col, y_col, e_col, legend, y_label, last_line in the_loop():
        # plot data
        y_label_from_col, y_plane, y_col, e_col, chromatic = _get_names_and_columns(idx_plot, xy,
                                                                                     y_col, e_col)

        x_val, y_val, e_val = _get_column_data(data, x_col, y_col, e_col)

        ebar = ax.errorbar(x_val, y_val, yerr=e_val,
                           ls=rcParams[u"lines.linestyle"], fmt=get_marker(idx, change_marker),
                           label=legend)

        _change_ebar_alpha(ebar)

        if auto_scale:
            current_y_lims = _get_auto_scale(y_val, auto_scale)
            if y_lims is None:
                y_lims = current_y_lims
            else:
                y_lims = [min(y_lims[0], current_y_lims[0]),
                          max(y_lims[1], current_y_lims[1])]
            if last_line:
                ax.set_ylim(*y_lims)

        # things to do only once
        if last_line:
            # setting the y_label
            if y_label is None:
                _set_ylabel(ax, y_col, y_label_from_col, y_plane, chromatic)
            else:
                y_label_from_label = ""
                if y_label:
                    y_label_from_label, y_plane, _, _, chromatic = _get_names_and_columns(
                        idx_plot, xy, y_label, "")
                _set_ylabel(ax, y_label, y_label_from_label, y_plane, chromatic)

            # setting x limits
            if x_is_position:
                try:
                    ps.set_xLimits(data.SEQUENCE, ax)
                except (AttributeError, ps.ArgumentError):
                    pass

            # setting visibility, ir-markers and label
            if xy and idx_plot == 0:
                ax.axes.get_xaxis().set_visible(False)
                if x_is_position and ir_positions:
                    ps.show_ir(ir_positions, ax, mode='lines')
            else:
                if x_is_position:
                    ps.set_xaxis_label(ax)
                    if ir_positions:
                        ps.show_ir(ir_positions, ax, mode='outside')

            if not no_legend and idx_plot == 0:
                ps.make_top_legend(ax, the_loop.get_ncols())

    return the_loop.get_figs()


def _export_plots(figs, output):
    """ Export all created figures to PDF """
    for param in figs:
        pdf_path = "{:s}_{:s}.pdf".format(output, param)
        mpdf = PdfPages(pdf_path)
        fig = figs[param]

        try:
            mpdf.savefig(fig, bbox_inches='tight')
            LOG.debug("Exported GetLLM results to PDF '{:s}'".format(pdf_path))
        finally:
            mpdf.close()


# Helper #####################################################################


def _check_opt(opt):
    """ Sanity checks for the opt structure """
    if opt.source_names is None:
        opt.source_names = opt.files
    elif len(opt.source_names) != len(opt.files):
        raise AttributeError("The number of legends (source_names) and number of files differ!")

    if opt.single_fig and len(opt.files) > 1:
        raise ValueError("Single figure plotting mode works only with one file!")

    if opt.labels is None:
        opt.labels = [None] * len(opt.y_cols)
    elif len(opt.labels) != len(opt.y_cols):
        raise AttributeError("The number of labels and number of y columns differ!")

    if opt.e_cols is None:
        opt.e_cols = [""] * len(opt.y_cols)
    elif len(opt.e_cols) != len(opt.y_cols):
        raise AttributeError("The number of error columns and number of y columns differ!")

    if opt.x_cols is None:
        opt.x_cols = ["S"] * len(opt.y_cols)
    elif len(opt.x_cols) != len(opt.y_cols):
        raise AttributeError("The number of x columns and number of y columns differ!")

    return opt


def get_marker(idx, change):
    """ Return the marker used """
    if change:
        return ps.MarkerList.get_marker(idx)
    else:
        return rcParams['lines.marker']


def _get_auto_scale(y_val, scaling):
    """ Find the y-limits so that scaling% of the points are visible """
    y_sorted = sorted(y_val)
    n_points = len(y_val)
    y_min = y_sorted[int(((1 - scaling/100.) / 2.) * n_points)]
    y_max = y_sorted[int(((1 + scaling/100.) / 2.) * n_points)]
    return y_min, y_max


def _find_ir_pos(all_data):
    """ Return the middle positions of the interaction regions """
    ip_names = ["IP" + str(i) for i in range(1, 9)]
    for data in all_data:
        try:
            ip_pos = data.loc[ip_names, 'S'].values
        except KeyError:
            try:
                # loading failed, use defaults
                return IR_POS_DEFAULT[data.SEQUENCE]
            except AttributeError:
                # continue looking
                pass
        else:
            return dict(zip(ip_names, ip_pos))

    # did not find ips or defaults
    return {}


def _get_ir_positions(all_data, x_cols):
    """ Check if x is position around the ring and return ir positions if possible """
    ir_pos = None
    x_is_pos = all([xc == "S" for xc in x_cols])
    if x_is_pos:
        ir_pos = _find_ir_pos(all_data)
    return ir_pos, x_is_pos


def _get_gridspec(xy):
    """ Single or dual plot """
    if xy:
        return gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    return gridspec.GridSpec(1, 1, height_ratios=[1])


def _get_column_data(data, x_col, y_col, e_col):
    """ Extract column data """
    x_val = data[x_col]
    y_val = data[y_col]
    try:
        e_val = data[e_col]
    except KeyError:
        e_val = None
    return x_val, y_val, e_val


def _get_names_and_columns(idx_plot, xy, y_col, e_col):
    """ Names and columns """
    chromatic = False
    if xy:
        if y_col[-5:] in COMPLEX_NAMES:
            plane_map = {0: "R", 1: "I"}
            y_name = plane_map[idx_plot]
            if "C" == y_col[0]:
                y_plane_name = y_col[1:]
                chromatic = True
            else:
                y_plane_name = y_col
        else:
            plane_map = {0: "X", 1: "Y"}
            y_name = y_col
            y_plane_name = plane_map[idx_plot]
        y_col_full = y_col + plane_map[idx_plot]
        e_col_full = e_col + plane_map[idx_plot]
    else:
        if y_col[-5:] in COMPLEX_NAMES:
            y_name = y_col[-1]
            if "C" == y_col[0]:
                y_plane_name = y_col[1:-1]
                chromatic = True
            else:
                y_plane_name = y_col[:-1]
        else:
            y_name = y_col[:-1]
            y_plane_name = y_col[-1]

        y_col_full = y_col
        e_col_full = e_col
    return y_name, y_plane_name, y_col_full, e_col_full, chromatic


def _change_ebar_alpha(ebar):
    """ loop through bars (ebar[1]) and caps (ebar[2]) and set the alpha value """
    for bars_or_caps in ebar[1:]:
        for bar_or_cap in bars_or_caps:
            bar_or_cap.set_alpha(ERROR_ALPHA)


def _set_ylabel(ax, default, y_label, y_plane, chromatic):
    """ Tries to set a mapped y label, otherwise the default """
    try:
        ps.set_yaxis_label(_map_proper_name(y_label),
                           y_plane, ax, chromcoup=chromatic)
    except (KeyError, ps.ArgumentError):
        ax.set_ylabel(default)


def _map_proper_name(name):
    """ Maps to a name understood by plotstyle. """
    return {
        "BET": "beta",
        "BB": "betabeat",
        "D": "dispersion",
        "ND": "norm_dispersion",
        "MU": "phase",
        "X": "co",
        "Y": "co",
        "PHASE": "phase",
        "I": "imag",
        "R": "real",
    }[name.upper()]

# Script Mode ################################################################


if __name__ == "__main__":
    plot()
