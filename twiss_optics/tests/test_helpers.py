import os
import numpy as np
from Utilities import logging_tools
from Utilities import iotools
from Utilities import tfs_pandas as tfs
from Utilities.dict_tools import DotDict
from Utilities.plotting import plot_style as pstyle
from collections import OrderedDict
import matplotlib.pyplot as plt

LOG = logging_tools.get_logger(__name__)

"""
====================== Error Functions ======================
"""

# Relative Error measurements to get rid of spikes where value is close to zero
# see https://en.wikipedia.org/wiki/Relative_change_and_difference and references


def error_rel(x, y):
    return (x - y) / y


def error_relmean(x, y):
    return 2 * (x - y) / (x + y)


def error_relabs(x, y):
    return np.abs(x - y) / np.abs(y)


def error_relmax(x, y):
    return (x - y) / np.max([np.abs(x), np.abs(y)], axis=0)


def error_relabsmax(x, y):
    return np.abs(x - y) / np.max([np.abs(x), np.abs(y)], axis=0)


def error_relabsmean(x, y):
    return 2 * np.abs(x - y) / (np.abs(x) + np.abs(y))


def error_log(x, y):
    return 100 * np.log(x/y)


def error_plusone(x, y):
    return (x - y) / (y + 1)


def error_abs(x, y):
    return np.abs(x-y)


def error_meanabs(x, y):
    mean = np.mean(np.abs(y))
    if mean != 0:
        return (x - y) / mean
    else:
        return np.ones(x.shape)


"""
====================== Config ======================
"""

opt = DotDict({
    'plot': {
        'show': False,
        'save': True,
        'style': 'presentation',
        'manual': {
                   u'lines.linestyle': '-',
                   u'lines.marker': '',
                   u'lines.markersize': 5,
                   u'lines.linewidth': 4,
                   u'markers.fillstyle': 'none',
                  },
        'second_linestyle': '--',
        'second_marker': '',
    },
    'error_fun': error_meanabs
})


name_map = {
        "BET": "beta",
        "MU": "phase",
        "D": "dispersion",
        "Q": "tune",
    }

"""
====================== Plot Helpers ======================
"""


def plot_df_comparison(df_one, df_two, title, planes, ylabel,
                       accel=None, plane_labels=None, data_labels=None):
    pstyle.set_style(opt.plot.style, opt.plot.manual)
    if data_labels is None:
        data_labels = ["Analytical", "MADX"]

    if "S" in df_one:
        x_axis = df_one.loc[:, "S"]
    elif "S" in df_two:
        x_axis = df_two.loc[:, "S"]
    else:
        raise IOError("'S' is not available in either input data-frame.")

    for iplane, plane in enumerate(planes):
        if plane_labels is None:
            plane_label = plane[-1]
        else:
            plane_label = plane_labels[iplane]

        # normal plot
        fig, ax = plt.subplots()
        title_str = "{:s} {:s}".format(title, plane)
        ax.set_title(title_str)
        pstyle.set_name(title_str, ax)
        pstyle.set_yaxis_label(ylabel, plane_label, ax)
        ax.plot(x_axis, df_one[plane], label=data_labels[0])
        ax.plot(x_axis, df_two[plane], label=data_labels[1],
                linestyle=opt.plot.second_linestyle,
                marker=opt.plot.second_marker)
        nice_plots(fig, ax, accel)

        # error values
        error_values = opt.error_fun(df_one[plane], df_two[plane])

        # error plot
        fig_err, ax_err = plt.subplots()
        title_err = "{:s} {:s} - Error".format(title, plane)
        ax_err.set_title(title_err)
        pstyle.set_name(title_err, ax_err)
        pstyle.set_yaxis_label(ylabel, plane_label, ax_err, delta=True)
        ax_err.plot(x_axis, error_values, label='error')

        # rms
        err_mean = np.mean(error_values)
        err_rms = rms(error_values)
        LOG.info('  {p:s} Error - Mean: {m:g},  RMS: {r:g}'.format(p=plane, m=err_mean, r=err_rms))
        pstyle.set_annotation('Mean: {m:g},  RMS: {r:g}'.format(m=err_mean, r=err_rms), ax_err)
        nice_plots(fig_err, ax_err, accel)


def plot_single_magnet(results, title):
    pstyle.set_style(opt.plot.style, opt.plot.manual)
    _stat_plot(results["madx_mean"], title + " - TwR vs MADX MEAN")
    _stat_plot(results["madx_rms"], title + " - TwR vs MADX RMS", log=True)
    _stat_plot(results["madx_resp_mean"], title + " - TwR vs MaR MEAN")
    _stat_plot(results["madx_resp_rms"], title + " - TwR vs MaR RMS", log=True)


def _stat_plot(df, title, log=False):
    """ Plot over statistics with exponents """
    fig, ax = plt.subplots()
    if log:
        ax.set_yscale("log", nonposx='clip')

    ax.set_title(title)
    pstyle.set_name(title, ax)
    ax.set_ylabel(r"Relative Error $\frac{\delta Twiss - \delta MADX}{mean(|\delta MADX|)}$")
    ax.set_xlabel(r"Approximate $log_{10}(\delta k)$ ")

    x_axis = df.index.tolist()
    parameters = sorted(set([col.replace("_min", "").replace("_max", "") for col in df.columns]))
    for param in parameters:
        errb = [df[param] - df[param + "_min"],
                df[param + "_max"] - df[param]]
        ax.errorbar(x_axis, df[param], yerr=errb, label=r"$\delta$" + param,
                    markersize=8, capsize=10,
                    elinewidth=opt.plot.manual["lines.linewidth"]/2)[-1][0].set_linestyle('--')

    ax.legend()
    if not log:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 3))
    plt.draw()
    save_fig(fig)


def nice_plots(fig, ax, accel=None):
    ax.legend()
    try:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 3))
    except AttributeError:
        pass
    pstyle.set_xaxis_label(ax)
    if accel:
        try:
            pstyle.set_xLimits(accel, ax)
        except pstyle.ArgumentError:
            pass
    plt.draw()
    save_fig(fig)


def save_fig(fig, formats=('svg', 'pdf')):
    if not opt.plot.save:
        return

    iotools.create_dirs(os.path.abspath(os.path.join('.', 'plots')))
    name = pstyle.get_name(fig).replace("'", ""
                              ).replace('"', ''
                              ).replace(' ', '_'
                              ).replace("$", ""
                              ).replace("\\", "")
    for f in formats:
        fig.savefig(os.path.abspath(os.path.join('.', 'plots', name + '.' + f)),
                    format=f, dpi=800)



"""
====================== Other Helpers ======================
"""


def rms(x):
    return np.sqrt(np.mean(np.square(x)))
