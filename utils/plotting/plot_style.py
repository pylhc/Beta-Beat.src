"""
Helper functions to make the most awesome* plots out there.

* please feel free to add more stuff
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from utils import tfs_pandas as tfs


class ArgumentError(Exception):
    pass


# commented values are part of matplotlib 2.1 but not 1.5
_PRESENTATION_PARAMS = {
    # u'axes.autolimit_mode': u'data',
    u'backend': u'pdf',
    u'axes.edgecolor': u'k',
    u'axes.facecolor': u'w',
    u'axes.grid': True,
    u'axes.grid.axis': u'both',
    u'axes.grid.which': u'major',
    u'axes.labelcolor': u'k',
    # u'axes.labelpad': 4.0,
    u'axes.labelsize': u'x-large',
    u'axes.labelweight': u'bold',
    u'axes.linewidth': 1.8,
    # u'axes.titlepad': 16.0,
    u'axes.titlesize': u'xx-large',
    u'axes.titleweight': u'bold',
    u'figure.edgecolor': u'w',
    u'figure.facecolor': u'w',
    u'figure.figsize': [10.24, 7.68],
    u'figure.frameon': True,
    u'figure.titlesize': u'xx-large',
    u'figure.titleweight': u'bold',
    u'font.size': 15.0,
    u'font.stretch': u'normal',
    u'font.weight': u'bold',
    u'font.family': 'sans-serif',
    u'font.serif': ['Computer Modern'],
    u'grid.alpha': 1.0,
    u'grid.color': u'#b0b0b0',
    u'grid.linestyle': u'-',
    u'grid.linewidth': 1,
    u'legend.edgecolor': u'0.8',
    u'legend.facecolor': u'inherit',
    u'legend.fancybox': True,
    u'legend.fontsize': u'medium',
    u'legend.framealpha': 0.9,
    u'legend.frameon': True,
    u'legend.handleheight': 0.7,
    u'legend.handlelength': 2.0,
    u'legend.handletextpad': 0.8,
    u'legend.labelspacing': 0.5,
    u'legend.loc': u'best',
    u'legend.markerscale': 1.2,
    u'legend.numpoints': 1,
    u'legend.scatterpoints': 1,
    u'legend.shadow': False,
    u'lines.antialiased': True,
    # u'lines.color': u'C0',
    u'lines.linestyle': u'-',
    u'lines.linewidth': 1,
    u'lines.marker': u'.',
    u'lines.markeredgewidth': 1.0,
    u'lines.markersize': 14.0,
    u'lines.solid_capstyle': u'projecting',
    u'lines.solid_joinstyle': u'round',
    u'markers.fillstyle': u'full',
    u'text.antialiased': True,
    u'text.color': u'k',
    # u'xtick.alignment': u'center',
    # u'xtick.bottom': True,
    # u'xtick.color': u'k',
    # u'xtick.direction': u'out',
    # u'xtick.labelsize': u'medium',
    # u'xtick.major.bottom': True,
    # u'xtick.major.pad': 3.5,
    # u'xtick.major.size': 3.5,
    # u'xtick.major.top': True,
    # u'xtick.major.width': 1.2,
    # u'xtick.minor.bottom': True,
    # u'xtick.minor.pad': 3.4,
    # u'xtick.minor.size': 2.0,
    # u'xtick.minor.top': True,
    # u'xtick.minor.visible': False,
    # u'xtick.minor.width': 1,
    # u'xtick.top': False,
    # u'ytick.alignment': u'center_baseline',
    # u'ytick.color': u'k',
    # u'ytick.direction': u'out',
    # u'ytick.labelsize': u'medium',
    # u'ytick.left': True,
    # u'ytick.major.left': True,
    # u'ytick.major.pad': 3.5,
    # u'ytick.major.right': True,
    # u'ytick.major.size': 3.5,
    # u'ytick.major.width': 1.2,
    # u'ytick.minor.left': True,
    # u'ytick.minor.pad': 3.4,
    # u'ytick.minor.right': True,
    # u'ytick.minor.size': 2.0,
    # u'ytick.minor.visible': False,
    # u'ytick.minor.width': 1,
    # # u'ytick.right': False
}


_STANDARD_PARAMS = {
    # u'axes.autolimit_mode': u'data',
    u'axes.edgecolor': u'k',
    u'axes.facecolor': u'w',
    u'axes.grid': True,
    u'axes.grid.axis': u'both',
    u'axes.grid.which': u'major',
    u'axes.labelcolor': u'k',
    # u'axes.labelpad': 4.0,
    u'axes.labelsize': u'medium',
    u'axes.labelweight': u'bold',
    u'axes.linewidth': 1.5,
    # u'axes.titlepad': 6.0,
    u'axes.titlesize': u'x-large',
    u'axes.titleweight': u'bold',
    u'figure.edgecolor': u'w',
    u'figure.facecolor': u'w',
    u'figure.figsize': [10.24, 7.68],
    u'figure.frameon': True,
    u'figure.titlesize': u'large',
    u'figure.titleweight': u'normal',
    u'font.size': 15.0,
    u'font.stretch': u'normal',
    u'font.weight': u'bold',
    u'font.family': 'sans-serif',
    u'font.serif': ['Computer Modern'],
    u'grid.alpha': 1.0,
    u'grid.color': u'#b0b0b0',
    u'grid.linestyle': u'-',
    u'grid.linewidth': 0.6,
    u'legend.edgecolor': u'0.8',
    u'legend.facecolor': u'inherit',
    u'legend.fancybox': True,
    u'legend.fontsize': u'medium',
    u'legend.framealpha': 0.8,
    u'legend.frameon': False,
    u'legend.handleheight': 0.7,
    u'legend.handlelength': 2.0,
    u'legend.handletextpad': 0.8,
    u'legend.labelspacing': 0.5,
    u'legend.loc': u'best',
    u'legend.markerscale': .8,
    u'legend.numpoints': 1,
    u'legend.scatterpoints': 1,
    u'legend.shadow': False,
    u'lines.antialiased': True,
    # u'lines.color': u'C0',
    u'lines.linestyle': u'-',
    u'lines.linewidth': 1.5,
    u'lines.marker': u'o',
    u'lines.markeredgewidth': 1.0,
    u'lines.markersize': 8.0,
    u'lines.solid_capstyle': u'projecting',
    u'lines.solid_joinstyle': u'round',
    u'markers.fillstyle': u'none',
    u'text.antialiased': True,
    u'text.color': u'k',
    # u'xtick.alignment': u'center',
    # u'xtick.bottom': True,
    # u'xtick.color': u'k',
    # u'xtick.direction': u'out',
    # u'xtick.labelsize': u'medium',
    # u'xtick.major.bottom': True,
    # u'xtick.major.pad': 3.5,
    # u'xtick.major.size': 3.5,
    # u'xtick.major.top': True,
    # u'xtick.major.width': 0.8,
    # u'xtick.minor.bottom': True,
    # u'xtick.minor.pad': 3.4,
    # u'xtick.minor.size': 2.0,
    # u'xtick.minor.top': True,
    # u'xtick.minor.visible': False,
    # u'xtick.minor.width': 0.6,
    # u'xtick.top': False,
    # u'ytick.alignment': u'center_baseline',
    # u'ytick.color': u'k',
    # u'ytick.direction': u'out',
    # u'ytick.labelsize': u'medium',
    # u'ytick.left': True,
    # u'ytick.major.left': True,
    # u'ytick.major.pad': 3.5,
    # u'ytick.major.right': True,
    # u'ytick.major.size': 3.5,
    # u'ytick.major.width': 0.8,
    # u'ytick.minor.left': True,
    # u'ytick.minor.pad': 3.4,
    # u'ytick.minor.right': True,
    # u'ytick.minor.size': 2.0,
    # u'ytick.minor.visible': False,
    # u'ytick.minor.width': 0.6,
    # u'ytick.right': False
}


# Style ######################################################################


def set_style(style='standard', manual=None):
    """Sets the style for all following plots.

    Args:
        style: Choose Style, either 'standard' or 'presentation'
        manual: Dict of manual parameters to update. Convention: "REMOVE_ENTRY" removes entry
    """
    if style == 'standard':
        params = _STANDARD_PARAMS.copy()
    elif style == 'presentation':
        params = _PRESENTATION_PARAMS.copy()
    else:
        raise ArgumentError("Style '" + style + "' not found.")

    if manual:
        for key in manual.keys():
            if manual[key] == "REMOVE_ENTRY":
                params.pop(key)
            else:
                params[key] = manual[key]

    matplotlib.rcParams.update(params)


# Tools ######################################################################


def sync2d(axs, ax_str='xy', ax_lim=()):
    """
    Synchronizes the limits for the given axes

    Args:
        axs: list of axes or figures, or figure with multiple axes
        ax_lim: predefined limits (list or list of lists)
        ax_str: string 'x','y' or 'xy' defining the axes to sync
    """
    if isinstance(axs, list):
        if isinstance(axs[0], matplotlib.figure.Figure):
            # axs is list of figures: get all axes and call sync2D
            sync2d([ax for fig in axs for ax in fig.axes], ax_str)

        elif isinstance(axs[0], matplotlib.axes.Axes) and len(axs) > 1:
            # axs is list of axes
            if 'x' in ax_str:
                if len(ax_lim) == 0:
                    # find x limits
                    x_min = min([ax.get_xlim()[0] for ax in axs])
                    x_max = max([ax.get_xlim()[1] for ax in axs])
                    x_lim = [x_min, x_max]
                else:
                    # use defined limits
                    if len(ax_lim[0]) == 1:
                        x_lim = ax_lim
                    else:
                        x_lim = ax_lim[0]

                for ax in axs:
                    ax.set_xlim(x_lim)

            if 'y' in ax_str:
                if len(ax_lim) == 0:
                    # find x limits
                    y_min = min([ax.get_ylim()[0] for ax in axs])
                    y_max = max([ax.get_ylim()[1] for ax in axs])
                    y_lim = [y_min, y_max]
                else:
                    # use defined limits
                    if len(ax_lim[0]) == 1:
                        y_lim = ax_lim
                    else:
                        y_lim = ax_lim[1]

                for ax in axs:
                    ax.set_ylim(y_lim)

    elif isinstance(axs, matplotlib.figure.Figure):
        # axs is one figure: call sync2D with axes from figure
        sync2d(axs.axes, ax_str)

    else:
        raise TypeError(__file__[:-3] + '.sync2d input is of unknown type (' + type(axs) + ')')


def set_xLimits(accel, ax=None):
    """
    Sets the x-limits to the regularly used ones

    Args:
        accel: Name of the Accelerator
        ax:  Axes to put the label on (default: gca())
    """
    if not ax:
        ax = plt.gca()

    if accel.startswith("LHCB"):
        ax.set_xlim(-200, 27000)
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(base=1000.0))
        ax.xaxis.set_major_locator(mtick.MultipleLocator(base=5000.0))
    elif accel.startswith("ESRF"):
        ax.set_xlim(-5, 850)
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(base=50.0))
        ax.xaxis.set_major_locator(mtick.MultipleLocator(base=100.0))
    else:
        raise ArgumentError("Accelerator '" + accel + "' unknown.")


class MarkerList(object):
    """ Create a list of predefined markers """
    # markers = ["s", "o", ">", "D", "v", "*", "h", "^", "p", "X", "<", "P"]  # matplotlib 2.++
    markers = ["s", "o", ">", "D", "v", "*", "h", "^", "p", "<"]

    def __init__(self):
        self.idx = 0

    @classmethod
    def get_marker(cls, marker_num):
        """ Return marker of index marker_num

         Args:
             marker_num (int): return maker at this position in list (mod len(list))
        """
        return cls.markers[marker_num % len(cls.markers)]

    def get_next_marker(self):
        """ Return the next marker in the list (circularly wrapped) """
        marker = self.get_marker(self.idx)
        self.idx += 1
        return marker


# Labels #####################################################################


# List of common y-labels. Sorry for the ugly.
_ylabels = {
    "beta":               r'$\beta_{{{0}}} [m]$',
    "betabeat":           r'$\Delta \beta_{{{0}}} / \beta_{{{0}}}$',
    "betabeat_permile":   r'$\Delta \beta_{{{0}}} / \beta_{{{0}}} [$'u'\u2030'r'$]$',
    "dbeta":              r"$\beta'_{{{0}}} [m]$",
    "dbetabeat":          r'$1/\beta_{{{0}}} \cdot \partial\beta_{{{0}}} / \partial\delta_{{{0}}}$',
    "norm_dispersion":    r'$\frac{{\Delta D_{{{0}}}}}{{\beta_{{{0}}}}} [m]$',
    "norm_dispersion_mu": r'$\frac{{D_{{{0}}}}}{{\sqrt{{\beta_{{{0}}}}}}} [\mu m^{{1/2}}]$',
    "phase":              r'$\phi_{{{0}}} [2\pi]$',
    "phasetot":           r'$\phi_{{{0}}} [2\pi]$',
    "phase_milli":        r'$\phi_{{{0}}} [2\pi\cdot10^{{-3}}]$',
    "dispersion":         r'$D_{{{0}}} [m]$',
    "dispersion_mm":      r'$D_{{{0}}} [mm]$',
    "co":                 r'${0} [mm]$',
    "tune":               r'$Q_{{{0}}} [Hz]$',
    "nattune":            r'$Nat Q_{{{0}}} [Hz]$',
    "chromamp":           r'$W_{{{0}}}$',
    "real":               r'$re({0})$',
    "imag":               r'$im({0})$',
    "absolute":           r'$|{0}|$',
}


def set_yaxis_label(param, plane, ax=None, delta=False, chromcoup=False):  # plot x and plot y
    """ Set y-axis labels.

    Args:
        param: One of the ylabels above
        plane: Usually x or y, but can be any string actually to be placed into the label ({0})
        ax: Axes to put the label on (default: gca())
        delta: If True adds a Delta before the label (default: False)
    """
    if not ax:
        ax = plt.gca()
    try:
        label = _ylabels[param].format(plane)
    except KeyError:
        raise ArgumentError("Label '" + param + "' not found.")

    if delta:
        label = r'$\Delta ' + label[1:]

    if chromcoup:
        label = label[:-1] + r'/\Delta\delta$'

    ax.set_ylabel(label)


def set_xaxis_label(ax=None):
    """ Sets the standard x-axis label

    Args:
        ax: Axes to put the label on (default: gca())
    """
    if not ax:
        ax = plt.gca()
    ax.set_xlabel(r'Longitudinal location [m]')


def show_ir(ip_dict, ax=None, mode='inside'):
    """ Plots the interaction regions into the background of the plot.

    Args:
        ip_dict: dict, dataframe or series containing "IPLABEL" : IP_POSITION
        ax:  Axes to put the irs on (default: gca())
        mode: 'inside', 'outside' + 'nolines' or just 'lines'
    """
    if not ax:
        ax = plt.gca()

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    lines = 'nolines' not in mode
    inside = 'inside' in mode
    lines_only = 'inside' not in mode and 'outside' not in mode and 'lines' in mode

    if isinstance(ip_dict, (pd.DataFrame, pd.Series)):
        if isinstance(ip_dict, pd.DataFrame):
            ip_dict = ip_dict.iloc[:, 0]
        d = {}
        for ip in ip_dict.index:
            d[ip] = ip_dict.loc[ip]
        ip_dict = d

    for ip in ip_dict.keys():
        if xlim[0] <= ip_dict[ip] <= xlim[1]:
            xpos = ip_dict[ip]

            if lines:
                ax.axvline(xpos, linestyle=':', color='grey', marker='', zorder=0)

            if not lines_only:
                if inside:
                    ypos = ylim[0] + (ylim[1] - ylim[0]) * 0.01
                    ax.text(xpos, ypos, ip, color='grey', ha='center')
                else:
                    ax.text(xpos, ylim[1] * 1.03, ip, ha='center')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def get_ip_positions(path):
    """ Returns a dict of IP positions from tfs-file of path.

    Args:
        path (str): Path to the tfs-file containing IP-positions
    """
    df = tfs.read_tfs(path).set_index('NAME')
    ip_names = ["IP" + str(i) for i in range(1, 9)]
    ip_pos = df.loc[ip_names, 'S'].values
    return dict(zip(ip_names, ip_pos))


def set_name(name, fig_or_ax=None):
    """ Sets the name of the figure or axes

    Args:
        name (str): Sting to set as name.
        fig_or_ax: Figure or Axes to to use.
            If 'None' takes current figure. (Default: None)
    """
    if not fig_or_ax:
        fig_or_ax = plt.gcf()
    try:
        fig_or_ax.figure.canvas.set_window_title(name)
    except AttributeError:
        fig_or_ax.canvas.set_window_title(name)


def get_name(fig_or_ax=None):
    """ Returns the name of the figure or axes

    Args:
        fig_or_ax: Figure or Axes to to use.
            If 'None' takes current figure. (Default: None)
    """
    if not fig_or_ax:
        fig_or_ax = plt.gcf()
    try:
        return fig_or_ax.figure.canvas.get_window_title()
    except AttributeError:
        return fig_or_ax.canvas.get_window_title()


def set_annotation(text, ax=None):
    """ Writes an annotation on the top right of the axes

    Args:
        text: The annotation
        ax: Axes to set annotation on. If 'None' takes current Axes. (Default: None)
    """
    if not ax:
        ax = plt.gca()

    annotation = get_annotation(ax, by_reference=True)

    if annotation is None:
        ax.text(1.0, 1.0, text,
                verticalalignment='bottom',
                horizontalalignment='right',
                transform=ax.transAxes,
                label='plot_style_annotation')
    else:
        annotation.set_text(text)


def get_annotation(ax=None, by_reference=False):
    """ Returns the annotation set by set_annotation()

    Args:
        ax: Axes to get annotation from. If 'None' takes current Axes. (Default: None)
        by_reference (bool): If true returns the reference to the annotation,
            otherwise the text as string. (Default: False)
    """
    if not ax:
        ax = plt.gca()

    for c in ax.get_children():
        if c.get_label() == 'plot_style_annotation':
            if by_reference:
                return c
            else:
                return c.get_text()
    return None


def small_title(ax=None):
    """ Alternative to annotation, which lets you use the title-functions

    Args:
        ax: Axes to use. If 'None' takes current Axes. (Default: None)
    """
    if not ax:
        ax = plt.gca()

    # could not get set_title() to work properly, so one parameter at a time
    ax.title.set_position([1.0, 1.02])
    ax.title.set_transform(ax.transAxes)
    ax.title.set_fontsize(matplotlib.rcParams['font.size'])
    ax.title.set_fontweight(matplotlib.rcParams['font.weight'])
    ax.title.set_verticalalignment('bottom')
    ax.title.set_horizontalalignment('right')

