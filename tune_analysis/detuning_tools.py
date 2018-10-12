"""
Some tools for amplitude detuning, mainly plotting.

Important Convention:
    The beta-parameter in the ODR models go upwards with order, i.e.
    |  beta[0] = y-Axis offset
    |  beta[1] = slope
    |  beta[2] = quadratic term
    |  etc.

"""
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.odr import RealData, Model, ODR

import constants as const
from utils import logging_tools
from plotshop import plot_style as ps

LOG = logging_tools.get_logger(__name__)


# Linear ODR ###################################################################


def linear_model(beta, x):
    """ Return a linear model ``beta[0] + beta[1] * x``.

    Args:
        beta: beta[0] = y-offset
              beta[1] = slope
        x: x-value
    """
    return beta[0] + beta[1] * x


def do_linear_odr(x, y, x_err, y_err):
    """ Returns linear odr fit.

    Important: In contrast to "normal" ODR, this function uses x_err and y_err as weights
    instead of their squares.

    Args:
        x: Series of x data
        y: Series of y data
        x_err: Series of x data errors
        y_err: Series of y data errors

    Returns: Linear odr fit. Betas see ``linear_model()``.
    """
    lin_model = Model(linear_model)
    data = RealData(x, y, sx=x_err, sy=y_err)
    odr_fit = ODR(data, lin_model, beta0=[0., 1.]).run()
    print_odr_result(LOG.debug, odr_fit)
    return odr_fit


def print_odr_result(printer, odr_out):
        """ Logs the odr output results.

        Adapted from odr_output pretty print.
        """
        printer('Beta: {}'.format(odr_out.beta).replace("\n", ""))
        printer('Beta Std Error: {}'.format(odr_out.sd_beta).replace("\n", ""))
        printer('Beta Covariance: {}'.format(odr_out.cov_beta).replace("\n", ""))
        if hasattr(odr_out, 'info'):
            printer('Residual Variance: {}'.format(odr_out.res_var).replace("\n", ""))
            printer('Inverse Condition #: {}'.format(odr_out.inv_condnum).replace("\n", ""))
            printer('Reason(s) for Halting:')
            for r in odr_out.stopreason:
                printer('  {}'.format(r).replace("\n", ""))


def plot_linear_odr(ax, odr_fit, lim):
    """ Adds a linear odr fit to axes.
    """
    x_fit = np.linspace(lim[0], lim[1], 2)
    line_fit = odr_fit.beta[1] * x_fit
    ax.plot(x_fit, line_fit, marker="", linestyle='--', color='k',
            label='${:.4f}\, \pm\, {:.4f}$'.format(odr_fit.beta[1], odr_fit.sd_beta[1]))


# Data Extraction ##############################################################


def get_ampdet_data_from_kickac(kickac_df, action_plane, tune_plane):
    """ Extract the data needed for odr and plotting from the kickac dataframe.

    Args:
        kickac_df: Dataframe containing the data
        action_plane: Plane of the action
        tune_plane: Plane of the tune

    Returns:
        Dictionary containing x,y, x_err and y_err

    """
    columns = {"x": const.get_action_col(action_plane),
               "x_err": const.get_action_err_col(action_plane),
               "y": const.get_natq_corr_col(tune_plane),
               "y_err": const.get_total_natq_std_col(tune_plane),
               }
    data = {key: kickac_df.loc[:, columns[key]] for key in columns.keys()}
    return data


# General Plotting #############################################################


def plot_detuning(x, y, x_err, y_err, labels, x_min=None, x_max=None, y_min=None, y_max=None,
                  odr_fit=None, odr_plot=plot_linear_odr, output=None, show=True):
    """ Plot amplitude detuning.

    Args:
        x: Action data.
        y: Tune data.
        x_err: Action error.
        y_err: Tune error.
        x_min: Lower action range to plot.
        x_max: Upper action range to plot.
        y_min: Lower tune range to plot.
        y_max: Upper tune range to plot.
        odr_fit: results of the odr-fit (e.g. see do_linear_odr)
        odr_plot: function to plot odr_fit (e.g. see plot_linear_odr)
        labels: Dict of labels to use for the data ("line"), the x-axis ("x") and the y-axis ("y")
        output: Output file of the plot.
        show: Show the plot in window.

    Returns:
        Plotted Figure
    """
    ps.set_style("standard",
                 {u"lines.marker": u"o",
                  u"lines.linestyle": u"",
                  u'figure.figsize': [9.5, 4],
                  }
                 )

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x_min = 0 if x_min is None else x_min
    x_max = max(x + x_err)*1.01 if x_max is None else x_max
    offset = 0
    if odr_fit:
        odr_plot(ax, odr_fit, lim=[x_min, x_max])
        offset = odr_fit.beta[0]

    ax.errorbar(x, y-offset, xerr=x_err, yerr=y_err, label=labels.get("line", None))

    default_labels = const.get_paired_lables("{}", "{}")
    ax.set_xlabel(labels.get("x", default_labels[0]))
    ax.set_ylabel(labels.get("y", default_labels[1]))

    ax.set_xlim(left=x_min, right=x_max)
    ax.set_ylim(bottom=y_min, top=y_max)

    plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=2,)
    fig.tight_layout()

    if output:
        fig.savefig(output)
        ps.set_name(os.path.basename(output))

    if show:
        plt.draw()

    return fig


# Script Mode #################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
