"""
Helper functions to make the most awesome* plots out there.

* please feel free to add more stuff
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


class ArgumentError(Exception):
    pass


_PRESENTATION_PARAMS = {
    u'axes.autolimit_mode': u'data',
    u'axes.axisbelow': u'line',
    u'backend': u'pdf',
    u'axes.edgecolor': u'k',
    u'axes.facecolor': u'w',
    u'axes.grid': True,
    u'axes.grid.axis': u'both',
    u'axes.grid.which': u'major',
    u'axes.labelcolor': u'k',
    u'axes.labelpad': 4.0,
    u'axes.labelsize': u'x-large',
    u'axes.labelweight': u'bold',
    u'axes.linewidth': 1.8,
    u'axes.titlepad': 6.0,
    u'axes.titlesize': u'xx-large',
    u'axes.titleweight': u'bold',
    u'figure.edgecolor': u'w',
    u'figure.facecolor': u'w',
    u'figure.figsize': [10.24, 7.68],
    u'figure.frameon': True,
    u'figure.titlesize': u'xx-large',
    u'figure.titleweight': u'bold',
    u'font.size': 14.0,
    u'font.stretch': u'normal',
    u'font.weight': u'bold',
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
    u'lines.color': u'C0',
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
    u'xtick.alignment': u'center',
    u'xtick.bottom': True,
    u'xtick.color': u'k',
    u'xtick.direction': u'out',
    u'xtick.labelsize': u'medium',
    u'xtick.major.bottom': True,
    u'xtick.major.pad': 3.5,
    u'xtick.major.size': 3.5,
    u'xtick.major.top': True,
    u'xtick.major.width': 1.2,
    u'xtick.minor.bottom': True,
    u'xtick.minor.pad': 3.4,
    u'xtick.minor.size': 2.0,
    u'xtick.minor.top': True,
    u'xtick.minor.visible': False,
    u'xtick.minor.width': 1,
    u'xtick.top': False,
    u'ytick.alignment': u'center_baseline',
    u'ytick.color': u'k',
    u'ytick.direction': u'out',
    u'ytick.labelsize': u'medium',
    u'ytick.left': True,
    u'ytick.major.left': True,
    u'ytick.major.pad': 3.5,
    u'ytick.major.right': True,
    u'ytick.major.size': 3.5,
    u'ytick.major.width': 1.2,
    u'ytick.minor.left': True,
    u'ytick.minor.pad': 3.4,
    u'ytick.minor.right': True,
    u'ytick.minor.size': 2.0,
    u'ytick.minor.visible': False,
    u'ytick.minor.width': 1,
    u'ytick.right': False
}


_STANDARD_PARAMS = {
    u'axes.autolimit_mode': u'data',
    u'axes.axisbelow': u'line',
    u'axes.edgecolor': u'k',
    u'axes.facecolor': u'w',
    u'axes.grid': True,
    u'axes.grid.axis': u'both',
    u'axes.grid.which': u'major',
    u'axes.labelcolor': u'k',
    u'axes.labelpad': 4.0,
    u'axes.labelsize': u'medium',
    u'axes.labelweight': u'bold',
    u'axes.linewidth': 1.5,
    u'axes.titlepad': 6.0,
    u'axes.titlesize': u'x-large',
    u'axes.titleweight': u'bold',
    u'figure.edgecolor': u'w',
    u'figure.facecolor': u'w',
    u'figure.figsize': [10.24, 7.68],
    u'figure.frameon': True,
    u'figure.titlesize': u'large',
    u'figure.titleweight': u'normal',
    u'font.size': 12.0,
    u'font.stretch': u'normal',
    u'font.weight': u'bold',
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
    u'lines.color': u'C0',
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
    u'xtick.alignment': u'center',
    u'xtick.bottom': True,
    u'xtick.color': u'k',
    u'xtick.direction': u'out',
    u'xtick.labelsize': u'medium',
    u'xtick.major.bottom': True,
    u'xtick.major.pad': 3.5,
    u'xtick.major.size': 3.5,
    u'xtick.major.top': True,
    u'xtick.major.width': 0.8,
    u'xtick.minor.bottom': True,
    u'xtick.minor.pad': 3.4,
    u'xtick.minor.size': 2.0,
    u'xtick.minor.top': True,
    u'xtick.minor.visible': False,
    u'xtick.minor.width': 0.6,
    u'xtick.top': False,
    u'ytick.alignment': u'center_baseline',
    u'ytick.color': u'k',
    u'ytick.direction': u'out',
    u'ytick.labelsize': u'medium',
    u'ytick.left': True,
    u'ytick.major.left': True,
    u'ytick.major.pad': 3.5,
    u'ytick.major.right': True,
    u'ytick.major.size': 3.5,
    u'ytick.major.width': 0.8,
    u'ytick.minor.left': True,
    u'ytick.minor.pad': 3.4,
    u'ytick.minor.right': True,
    u'ytick.minor.size': 2.0,
    u'ytick.minor.visible': False,
    u'ytick.minor.width': 0.6,
    u'ytick.right': False
}


"""
=============================   Style   =============================
"""


def set_style(style='standard', manual=None):
    """Sets the style for all following plots.

    :param style: Choose Style, either 'standard' or 'presentation'
    :param manual: Dict of manual parameters to update. Convention: "REMOVE_ENTRY" removes entry
    :return:
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


"""
=============================   Tools   =============================
"""


def sync2d(axs, ax_str='xy', ax_lim=()):
    """
    Syncronizes the limits for the given axes

    :param axs: list of axes or figures, or figure with multiple axes
    :param ax_lim: predefined limits (list or list of lists)
    :param ax_str: string 'x','y' or 'xy' defining the axes to sync
    :return:
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


"""
=============================   Labels   =============================
"""

# List of common y-labels. Sorry for the ugly.
_ylabels = {
    "beta":               r'$\beta_{{{0}}} [m]$',
    "betabeat":           r'$\Delta \beta_{{{0}}} / \beta_{{{0}}}$',
    "betabeat_permile":   r'$\Delta \beta_{{{0}}} / \beta_{{{0}}} [$'u'\u2030'r'$]$',
    "norm_dispersion":    r'$ND_{{{0}}} [m^{{1/2}}]$',
    "norm_dispersion_mu": r'$\frac{{D_{{{0}}}}}{{\sqrt{{\beta_{{{0}}}}}}} [\mu m^{{1/2}}]$',
    "phase":              r'$\phi_{{{0}}} [2\pi]$',
    "phasetot":           r'$\phi_{{{0}}} [2\pi]$',
    "phase_milli":        r'$\phi_{{{0}}} [2\pi\cdot10^{{-3}}]$',
    "dispersion":         r'$D_{{{0}}} [m]$',
    "dispersion_mm":      r'$D_{{{0}}} [mm]$',
    "co":                 r'{{{0}}} [mm]',
    "tune":               r'$Q_{{{0}}} [Hz]$',
    "nattune":            r'$Nat Q_{{{0}}} [Hz]$',
    "real":               r'$re({0})$',
    "imag":               r'$im({0})$',
    "absolute":           r'$|{0}|$',
}


def set_yaxis_label(param, plane, ax=None, delta=False):  # plot x and plot y
    """
    Set y-axis labels.

    :param param: One of the ylabels above
    :param plane: Usually x or y, but can be any string actually to be placed into the label ({0})
    :param ax: Axes to put the label on (default: gca())
    :param delta: If True adds a Delta before the label (default: False)
    :return:
    """
    if not ax:
        ax = plt.gca()
    try:
        label = _ylabels[param].format(plane)
    except KeyError:
        raise ArgumentError("Label '" + param + "' not found.")

    if delta:
        label = r'$\Delta ' + label[1:]
    ax.set_ylabel(label)


def set_xaxis_label(ax=None):
    """
    Sets the standard x-axis label

    :param ax: Axes to put the label on (default: gca())
    :return:
    """
    if not ax:
        ax = plt.gca()
    ax.set_xlabel(r'Longitudinal location [m]')


def show_ir(ip_dict, ax=None):
    """
    Plots the interaction regions into the background of the plot.

    :param ip_dict: dict, dataframe or series containing "IPLABEL" : IP_POSITION
    :param ax:  Axes to put the irs on (default: gca())
    :return:
    """
    if not ax:
        ax = plt.gca()

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

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
            ypos = ylim[0] + (ylim[1] - ylim[0]) * 0.01
            ax.text(xpos, ypos, ip, color='grey', ha='center')
            ax.axvline(xpos, linestyle=':', color='grey', marker='', zorder=0)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def set_xLimits(accel, ax=None):
    """
    Sets the x-limits to the regularly used ones

    :param accel: Name of the Accelerator
    :param ax:  Axes to put the label on (default: gca())
    :return:
    """
    if not ax:
        ax = plt.gca()

    if accel[0:4] == "LHCB":
        ax.set_xlim(-200, 27000)
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(base=1000.0))
        ax.xaxis.set_major_locator(mtick.MultipleLocator(base=5000.0))
    elif accel[0:4] == "ESRF":
        ax.set_xlim(-5, 850)
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(base=50.0))
        ax.xaxis.set_major_locator(mtick.MultipleLocator(base=100.0))
    else:
        raise ArgumentError("Accelerator '" + accel + "' unknown.")
