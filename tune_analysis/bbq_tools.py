import datetime
import os
import pytz

import matplotlib.dates as mdates
import numpy as np
from matplotlib import pyplot as plt, gridspec
from matplotlib.ticker import FormatStrFormatter

from parameter_config import *
from utils import logging_tools
from utils.contexts import suppress_warnings
from utils.plotting import plot_style as ps

TIME_COL = get_time_col()
PLANES = get_planes()

COL_MAV = get_mav_col
COL_IN_MAV = get_used_in_mav_col
COL_BBQ = get_bbq_col

LOG = logging_tools.get_logger(__name__)


def get_moving_average(data_series, length=20,
                       min_val=None, max_val=None, fine_length=None, fine_cut=None):
    """ Get a moving average of the ``data_series`` over ``length`` entries.
    The data can be filtered beforehand.
    The values are shifted, so that the averaged value takes ceil((length-1)/2) values previous
    and floor((length-1)/2) following values into account.

    Args:
        data_series: Series of data
        length: length of the averaging window
        min_val: minimum value (for filtering)
        max_val: maximum value (for filtering)
        fine_length: length of the averaging window for fine cleaning
        fine_cut: allowed deviation for fine cleaning

    Returns: filtered and averaged Series and the mask used for filtering data.
    """
    LOG.debug("Calculating BBQ moving average of length {:d}.".format(length))

    if bool(fine_length) != bool(fine_cut):
        raise NotImplementedError("To activate fine cleaning, both "
                                  "'fine_window' and 'fine_cut' are needed.")

    if min_val is not None:
        min_mask = data_series <= min_val
    else:
        min_mask = np.zeros(data_series.size, dtype=bool)

    if max_val is not None:
        max_mask = data_series >= max_val
    else:
        max_mask = np.zeros(data_series.size, dtype=bool)

    cut_mask = min_mask | max_mask
    data_mav = _get_interpolated_moving_average(data_series, cut_mask, length)

    if fine_length is not None:
        min_mask = data_series <= (data_mav - fine_cut)
        max_mask = data_series >= (data_mav + fine_cut)
        cut_mask = min_mask | max_mask
        data_mav = _get_interpolated_moving_average(data_series, cut_mask, fine_length)

    return data_mav, cut_mask


def add_to_kickac_df(kickac_df, bbq_series, column):
    """ Add bbq values from series to kickac dataframe into column.

    Args:
        kickac_df: kickac dataframe
                  (needs to contain column "TIME_COL" or has time as index)
        bbq_series: series of bbq data with time as index
        column: column name to add the data into

    Returns: modified kickac dataframe

    """
    time_indx = kickac_df.index
    if TIME_COL in kickac_df:
        time_indx = kickac_df[TIME_COL]

    values = []
    for time in time_indx:
        values.append(bbq_series.iloc[bbq_series.index.get_loc(time, method="nearest")])
    kickac_df[column] = values
    return kickac_df


def plot_bbq_data(bbq_df,
                  interval=None, xmin=None, xmax=None, ymin=None, ymax=None,
                  output=None, show=True, two_plots=False):
    """ Plot BBQ data.

    Args:
        bbq_df: BBQ Dataframe with moving average columns
        interval: start and end time of used interval, will be marked with red bars
        xmin: Lower x limit (time)
        xmax: Upper x limit (time)
        ymin: Lower y limit (tune)
        ymax: Upper y limit (tune)
        output: Path to the output file
        show: Shows plot if `True`
        two_plots: Plots each tune in it's own axes if `True`

    Returns:
        Plotted figure

    """
    LOG.debug("Plotting BBQ data.")

    ps.set_style("standard", {u"lines.marker": u""})

    fig = plt.figure()

    if two_plots:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
        ax = [fig.add_subplot(gs[1]), fig.add_subplot(gs[0])]
    else:
        gs = gridspec.GridSpec(1, 1, height_ratios=[1])
        ax = fig.add_subplot(gs[0])
        ax = [ax, ax]

    tz = get_experiment_timezone()
    bbq_df.index = [datetime.datetime.fromtimestamp(time, tz=tz) for time in bbq_df.index]

    for idx, plane in enumerate(PLANES):
        color = ps.get_mpl_color(idx)
        mask = bbq_df[COL_IN_MAV(plane)]

        with suppress_warnings(UserWarning):  # caused by _nolegend_
            bbq_df.plot(
                y=COL_BBQ(plane), ax=ax[idx], color=color, alpha=.2,
                label="_nolegend_"
            )
        bbq_df.loc[mask, :].plot(
            y=COL_BBQ(plane), ax=ax[idx], color=color, alpha=.4,
            label="$Q_{:s}$ filtered".format(plane.lower())
        )
        bbq_df.plot(
            y=COL_MAV(plane), ax=ax[idx], color=color,
            label="$Q_{:s}$ moving av.".format(plane.lower())
        )

        if ymin is None and two_plots:
            ax[idx].set_ylim(bottom=min(bbq_df.loc[mask, COL_BBQ(plane)]))

        if ymax is None and two_plots:
            ax[idx].set_ylim(top=max(bbq_df.loc[mask, COL_BBQ(plane)]))

    # things to add/do only once if there is only one plot
    for idx in range(1+two_plots):
        if interval:
            ax[idx].axvline(x=interval[0], color="red")
            ax[idx].axvline(x=interval[1], color="red")

        ax[idx].set_ylabel('Tune')
        ax[idx].set_ylim(bottom=ymin, top=ymax)
        ax[idx].yaxis.set_major_formatter(FormatStrFormatter('%.5f'))

        ax[idx].set_xlim(left=xmin, right=xmax)
        ax[idx].set_xlabel('Time')
        ax[idx].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

        # don't show labels on upper plot (if two plots)
        if idx:
            # use the visibility to allow cursor x-position to be shown
            ax[idx].tick_params(labelbottom=False)
            ax[idx].xaxis.get_label().set_visible(False)

    plt.tight_layout()

    if output:
        fig.savefig(output)
        ps.set_name(os.path.basename(output))

    if show:
        plt.draw()

    return fig


# Private methods ############################################################

def _get_interpolated_moving_average(data_series, clean_mask, length):
    """ Returns the moving average of data series with a window of length and interpolated NaNs"""
    data_mav = data_series.copy()
    data_mav[clean_mask] = np.NaN

    # 'interpolate' fills nan based on index/values of neighbours
    data_mav = data_mav.interpolate("index").fillna(method="bfill").fillna(method="ffill")

    shift = -int((length-1)/2)  # Shift average to middle value
    return data_mav.rolling(length).mean().shift(shift).fillna(
        method="bfill").fillna(method="ffill")


# Script Mode #################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
