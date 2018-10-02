"""
Entrypoint for amplitude detuning analysis.

This module provides functionality to run amplitude detuning analysis with
additionally getting BBQ data from timber, averaging and filtering this data and
subtracting it from the measurement data.

Furthermore, the orthogonal distance regression is utilized to get a
linear fit from the measurements.

Also, plotting functionality is integrated, for the amplitude detuning as well as for the bbq data.


:author: Joschua Dilly
"""

import datetime
import os

import matplotlib.pyplot as plt

from tune_analysis import bbq_tools, timber_extract, detuning_tools
import tune_analysis.constants as ta_const
from utils import logging_tools
from tfs_files import tfs_pandas as tfs
from utils.dict_tools import ParameterError
from utils.entrypoint import entrypoint, EntryPointParameters

# Globals ####################################################################

# Column Names
COL_TIME = ta_const.get_time_col
COL_BBQ = ta_const.get_bbq_col
COL_MAV = ta_const.get_mav_col
COL_IN_MAV = ta_const.get_used_in_mav_col
COL_NATQ = ta_const.get_natq_col
COL_CORRECTED = ta_const.get_natq_corr_col

PLANES = ta_const.get_planes()
TIMBER_KEY = ta_const.get_timber_bbq_key


LOG = logging_tools.get_logger(__name__)

# Get Parameters #############################################################


def _get_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--label",
        help="Label to identify this run.",
        name="label",
        type=str,
    )
    params.add_parameter(
        flags="--beam",
        help="Which beam to use.",
        name="beam",
        required=True,
        type=int,
    )
    params.add_parameter(
        flags="--orientation",
        help="Orientation we are in. 'H' or 'V'.",
        name="orientation",
        required=True,
        type=str,
    )
    params.add_parameter(
        flags="--timberin",
        help="Fill number of desired data or path to presaved tfs-file",
        name="timber_in",
        required=True,
    )
    params.add_parameter(
        flags="--timberout",
        help="Output location to save fill as tfs-file",
        name="timber_out",
        type=str,
    )
    params.add_parameter(
        flags="--bbqout",
        help="Output location to save bbq data as tfs-file",
        name="bbq_out",
        type=str,
    )
    params.add_parameter(
        flags="--kickac",
        help="Location of the kickac file",
        name="kickac_path",
        type=str,
        required=True,
    )
    params.add_parameter(
        flags="--kickacout",
        help="If given, writes out the modified kickac file",
        name="kickac_out",
        type=str,
    )
    params.add_parameter(
        flags="--window",
        help="Length of the moving average window. (# data points)",
        name="window_length",
        type=int,
        default=20,
    )

    # cleaning method one:
    params.add_parameter(
        flags="--tunex",
        help="Horizontal Tune. For BBQ cleaning.",
        name="tune_x",
        type=float,
    )
    params.add_parameter(
        flags="--tuney",
        help="Vertical Tune. For BBQ cleaning.",
        name="tune_y",
        type=float,
    )
    params.add_parameter(
        flags="--tunecut",
        help="Cuts for the tune. For BBQ cleaning.",
        name="tune_cut",
        type=float,
    )
    # cleaning method two:
    params.add_parameter(
        flags="--tunexmin",
        help="Horizontal Tune minimum. For BBQ cleaning.",
        name="tune_x_min",
        type=float,
    )
    params.add_parameter(
        flags="--tunexmax",
        help="Horizontal Tune minimum. For BBQ cleaning.",
        name="tune_x_max",
        type=float,
    )
    params.add_parameter(
        flags="--tuneymin",
        help="Vertical  Tune minimum. For BBQ cleaning.",
        name="tune_y_min",
        type=float,
    )
    params.add_parameter(
        flags="--tuneymax",
        help="Vertical Tune minimum. For BBQ cleaning.",
        name="tune_y_max",
        type=float,
    )

    # fine cleaning
    params.add_parameter(
        flags="--finewindow",
        help="Length of the moving average window. (# data points)",
        name="fine_window",
        type=int,
    )
    params.add_parameter(
        flags="--finecut",
        help="Cut (i.e. tolerance) of the tune for the fine cleaning.",
        name="fine_cut",
        type=float,
    )

    # Plotting
    params.add_parameter(
        flags="--bbqplot",
        help="Save the bbq plot here.",
        name="bbq_plot_out",
        type=str,
    )
    params.add_parameter(
        flags="--bbqplotshow",
        help="Show the bbq plot.",
        name="bbq_plot_show",
        action="store_true",
    )
    params.add_parameter(
        flags="--bbqplottwo",
        help="Two plots for the bbq plot.",
        name="bbq_plot_two",
        action="store_true",
    )
    params.add_parameter(
        flags="--bbqplotfull",
        help="Plot the full bqq data with interval as lines.",
        name="bbq_plot_full",
        action="store_true",
    )
    params.add_parameter(
        flags="--ampdetplot",
        help="Save the amplitude detuning plot here.",
        name="ampdet_plot_out",
        type=str,
    )
    params.add_parameter(
        flags="--ampdetplotshow",
        help="Show the amplitude detuning plot.",
        name="ampdet_plot_show",
        action="store_true",
    )
    params.add_parameter(
        flags="--ampdetplotymin",
        help="Minimum tune (y-axis) in amplitude detuning plot.",
        name="ampdet_plot_ymin",
        type=float,
    )
    params.add_parameter(
        flags="--ampdetplotymax",
        help="Maximum tune (y-axis) in amplitude detuning plot.",
        name="ampdet_plot_ymax",
        type=float,
    )
    params.add_parameter(
        flags="--ampdetplotxmin",
        help="Minimum action (x-axis) in amplitude detuning plot.",
        name="ampdet_plot_xmin",
        type=float,
    )
    params.add_parameter(
        flags="--ampdetplotxmax",
        help="Maximum action (x-axis) in amplitude detuning plot.",
        name="ampdet_plot_xmax",
        type=float,
    )

    # Debug
    params.add_parameter(
        flags="--debug",
        help="Activates Debug mode",
        name="debug",
        action="store_true",
    )
    params.add_parameter(
        flags="--logfile",
        help="Logfile if debug mode is active.",
        name="logfile",
        type=str,
    )

    return params


def _get_plot_params():
    params = EntryPointParameters()
    params.add_parameter(
        flags="--in",
        help="BBQ data as data frame or tfs file.",
        name="input",
        required=True,
    )
    params.add_parameter(
        flags="--out",
        help="Save figure to this location.",
        name="output",
        type=str,
    )
    params.add_parameter(
        flags="--show",
        help="Show plot.",
        name="show",
        action="store_true"
    )
    params.add_parameter(
        flags="--xmin",
        help="Lower x-axis limit. (yyyy-mm-dd HH:mm:ss.mmm)",
        name="xmin",
        type=str,
    )
    params.add_parameter(
        flags="--ymin",
        help="Lower y-axis limit.",
        name="ymin",
        type=float,
    )
    params.add_parameter(
        flags="--xmax",
        help="Upper x-axis limit. (yyyy-mm-dd HH:mm:ss.mmm)",
        name="xmax",
        type=str,
    )
    params.add_parameter(
        flags="--ymax",
        help="Upper y-axis limit.",
        name="ymax",
        type=float,
    )
    params.add_parameter(
        flags="--interval",
        help="x_axis interval that was used in calculations.",
        name="interval",
        type=str,
        nargs=2,
    )
    params.add_parameter(
        flags="--two",
        help="Plot two axis into the figure.",
        name="two_plots",
        action="store_true",
    )
    return params


# Main #########################################################################


@entrypoint(_get_params(), strict=True)
def analyse_with_bbq_corrections(opt):
    """ Create amplitude detuning analysis with BBQ correction from timber data.

    Keyword Args:
        ampdet_plot_out (str): Save the amplitude detuning plot here.
                          **Flags**: --ampdetplot
        ampdet_plot_show: Show the amplitude detuning plot.
                          **Flags**: --ampdetplotshow
                          **Action**: ``store_true``
        ampdet_plot_xmax (float): Maximum action (x-axis) in amplitude detuning plot.
                          **Flags**: --ampdetplotxmax
        ampdet_plot_xmin (float): Minimum action (x-axis) in amplitude detuning plot.
                                  **Flags**: --ampdetplotxmin
        ampdet_plot_ymax (float): Maximum tune (y-axis) in amplitude detuning plot.
                                  **Flags**: --ampdetplotymax
        ampdet_plot_ymin (float): Minimum tune (y-axis) in amplitude detuning plot.
                                  **Flags**: --ampdetplotymin
        bbq_out (str): Output location to save bbq data as tfs-file
                       **Flags**: --bbqout
        bbq_plot_full: Plot the full bqq data with interval as lines.
                       **Flags**: --bbqplotfull
                       **Action**: ``store_true``
        bbq_plot_out (str): Save the bbq plot here.
                            **Flags**: --bbqplot
        bbq_plot_show: Show the bbq plot.
                       **Flags**: --bbqplotshow
                       **Action**: ``store_true``
        bbq_plot_two: Two plots for the bbq plot.
                      **Flags**: --bbqplottwo
                      **Action**: ``store_true``
        debug: Activates Debug mode
               **Flags**: --debug
               **Action**: ``store_true``
        fine_cut (float): Cut (i.e. tolerance) of the tune for the fine cleaning.
                          **Flags**: --finecut
        fine_window (int): Length of the moving average window. (# data points)
                           **Flags**: --finewindow
        kickac_out (str): If given, writes out the modified kickac file
                          **Flags**: --kickacout
        label (str): Label to identify this run.
                     **Flags**: --label
        logfile (str): Logfile if debug mode is active.
                       **Flags**: --logfile
        timber_out (str): Output location to save fill as tfs-file
                          **Flags**: --timberout
        tune_cut (float): Cuts for the tune. For BBQ cleaning.
                          **Flags**: --tunecut
        tune_x (float): Horizontal Tune. For BBQ cleaning.
                        **Flags**: --tunex
        tune_x_max (float): Horizontal Tune minimum. For BBQ cleaning.
                            **Flags**: --tunexmax
        tune_x_min (float): Horizontal Tune minimum. For BBQ cleaning.
                            **Flags**: --tunexmin
        tune_y (float): Vertical Tune. For BBQ cleaning.
                        **Flags**: --tuney
        tune_y_max (float): Vertical Tune minimum. For BBQ cleaning.
                            **Flags**: --tuneymax
        tune_y_min (float): Vertical  Tune minimum. For BBQ cleaning.
                            **Flags**: --tuneymin
        window_length (int): Length of the moving average window. (# data points)
                             **Flags**: --window
                             **Default**: ``20``

     """
    LOG.info("Starting Amplitude Detuning Analysis")
    with logging_tools.DebugMode(active=opt.debug, log_file=opt.logfile):
        opt = _check_analyse_opt(opt)
        figs = {}

        # get data
        bbq_df = _get_timber_data(opt.beam, opt.timber_in, opt.timber_out)
        kickac_df = tfs.read_tfs(opt.kickac_path, index=COL_TIME())
        x_interval = _get_approx_bbq_interval(bbq_df, kickac_df.index, opt.window_length)

        # add moving average to kickac
        kickac_df, bbq_df = _add_moving_average(kickac_df, bbq_df,
                                                **opt.get_subdict([
                                                    "window_length",
                                                    "tune_x_min", "tune_x_max",
                                                    "tune_y_min", "tune_y_max",
                                                    "fine_cut", "fine_window"]
                                                )
                                                )

        # add corrected values to kickac
        kickac_df = _add_corrected_natural_tunes(kickac_df)

        # output kickac and bbq data
        if opt.kickac_out:
            tfs.write_tfs(opt.kickac_out, kickac_df, save_index=COL_TIME())

        if opt.bbq_out:
            tfs.write_tfs(opt.bbq_out, bbq_df.loc[x_interval[0]:x_interval[1]],
                          save_index=COL_TIME())

        if opt.bbq_plot_out or opt.bbq_plot_show:
            if opt.bbq_plot_full:
                figs["bbq"] = bbq_tools.plot_bbq_data(
                    bbq_df,
                    output=opt.bbq_plot_out,
                    show=opt.bbq_plot_show,
                    two_plots=opt.bbq_plot_two,
                    interval=[str(datetime.datetime.fromtimestamp(xint)) for xint in x_interval],
                )
            else:
                figs["bbq"] = bbq_tools.plot_bbq_data(
                    bbq_df.loc[x_interval[0]:x_interval[1]],
                    output=opt.bbq_plot_out,
                    show=opt.bbq_plot_show,
                    two_plots=opt.bbq_plot_two,
                )

        # amplitude detuning analysis
        plane = ta_const.get_plane_from_orientation(opt.orientation)
        for other_plane in PLANES:
            labels = ta_const.get_paired_lables(plane, other_plane)
            id_str = "J{:s}_Q{:s}".format(plane.upper(), other_plane.upper())

            # get proper data
            columns = ta_const.get_paired_columns(plane, other_plane)
            data = {key: kickac_df.loc[:, columns[key]] for key in columns.keys()}

            # plotting
            try:
                output = os.path.splitext(opt.ampdet_plot_out)
                output = "{:s}_{:s}{:s}".format(output[0], id_str, output[1])
            except AttributeError:
                output = None

            figs[id_str] = detuning_tools.plot_detuning(
                odr_plot=detuning_tools.linear_odr_plot,
                labels={"x": labels[0], "y": labels[1], "line": opt.label},
                output=output,
                show=opt.ampdet_plot_show,
                x_min=opt.ampdet_plot_xmin,
                x_max=opt.ampdet_plot_xmax,
                y_min=opt.ampdet_plot_ymin,
                y_max=opt.ampdet_plot_ymax,
                **data
            )

    if opt.bbq_plot_show or opt.ampdet_plot_show:
        plt.show()

    return figs


@entrypoint(_get_plot_params(), strict=True)
def plot_bbq_data(opt):
    """ Plot BBQ wrapper.

    Keyword Args:
        interval (str): x_axis interval that was used in calculations.
                        **Flags**: --interval
        output (str): Save figure to this location.
                      **Flags**: --out
        show: Show plot.
              **Flags**: --show
              **Action**: ``store_true``
        two_plots: Plot two axis into the figure.
                   **Flags**: --two
                   **Action**: ``store_true``
        xmax (str): Upper x-axis limit. (yyyy-mm-dd HH:mm:ss.mmm)
                    **Flags**: --xmax
        xmin (str): Lower x-axis limit. (yyyy-mm-dd HH:mm:ss.mmm)
                    **Flags**: --xmin
        ymax (float): Upper y-axis limit.
                      **Flags**: --ymax
        ymin (float): Lower y-axis limit.
                      **Flags**: --ymin
    """
    LOG.info("Plotting BBQ.")
    if isinstance(opt.input, basestring):
        bbq_df = tfs.read_tfs(opt.input, index=COL_TIME())
    else:
        bbq_df = opt.input
    opt.pop("input")

    bbq_tools.plot_bbq_data(bbq_df, **opt)

    if opt.show:
        plt.show()


# Private Functions ############################################################


def _check_analyse_opt(opt):
    """ Perform manual checks on opt-sturcture """
    LOG.debug("Checking Options.")

    # for label
    if opt.label is None:
        opt.label = "Amplitude Detuning for Beam {:d}".format(opt.beam)

    # check if cleaning is properly specified
    if (any([opt.tune_x, opt.tune_y, opt.tune_cut])
            and any([opt.tune_x_min, opt.tune_x_max, opt.tune_y_min, opt.tune_y_max])
    ):
        raise ParameterError("Choose either the method of cleaning BBQ"
                             "with tunes and cut or with min and max values")

    for plane in PLANES:
        tune = "tune_{:s}".format(plane.lower())
        if opt[tune]:
            if opt.tune_cut is None:
                raise ParameterError("Tune cut is needed for cleaning tune.")
            opt["{:s}_min".format(tune)] = opt[tune] - opt.tune_cut
            opt["{:s}_max".format(tune)] = opt[tune] + opt.tune_cut

    if bool(opt.fine_cut) != bool(opt.fine_window):
        raise ParameterError("To activate fine cleaning, both fine cut and fine window need"
                             "to be specified")
    return opt


def _get_approx_bbq_interval(bbq_df, time_array, window_length):
    """ Get data in approximate time interval,
    for averaging based on window length and kickac interval """
    bbq_tmp = bbq_df.dropna()

    i_start = max(bbq_tmp.index.get_loc(time_array[0], method='nearest') - int(window_length/2.),
                  0
                  )
    i_end = min(bbq_tmp.index.get_loc(time_array[-1], method='nearest') + int(window_length/2.),
                len(bbq_tmp.index)
                )

    return bbq_tmp.index[i_start], bbq_tmp.index[i_end]


def _get_timber_data(beam, input, output):
    """ Return Timber data from input """
    LOG.debug("Getting timber data from '{}'".format(input))
    try:
        fill_number = int(input)
    except ValueError:
        fill = tfs.read_tfs(input, index=COL_TIME())
        fill.drop([COL_MAV(p) for p in PLANES if COL_MAV(p) in fill.columns],
                  axis='columns')
    else:
        timber_keys = [TIMBER_KEY(plane, beam) for plane in PLANES]
        bbq_cols = [COL_BBQ(plane) for plane in PLANES]

        fill = timber_extract.lhc_fill_to_tfs(fill_number,
                                              keys=timber_keys,
                                              names=dict(zip(timber_keys, bbq_cols)))

        if output:
            tfs.write_tfs(output, fill, save_index=COL_TIME())

    return fill


def _add_moving_average(kickac_df, bbq_df, **kwargs):
    """ Adds the moving average of the bbq data to kickac_df and bbq_df. """
    LOG.debug("Calculating moving average.")
    for plane in PLANES:
        tune = "tune_{:s}".format(plane.lower())
        bbq_mav, mask = bbq_tools.get_moving_average(bbq_df[COL_BBQ(plane)],
                                                     length=kwargs["window_length"],
                                                     min_val=kwargs["{}_min".format(tune)],
                                                     max_val=kwargs["{}_max".format(tune)],
                                                     fine_length=kwargs["fine_window"],
                                                     fine_cut=kwargs["fine_cut"],
                                                     )
        bbq_df[COL_MAV(plane)] = bbq_mav
        bbq_df[COL_IN_MAV(plane)] = ~mask
        kickac_df = bbq_tools.add_to_kickac_df(kickac_df, bbq_mav, COL_MAV(plane))
    return kickac_df, bbq_df


def _add_corrected_natural_tunes(kickac_df):
    """ Adds the corrected natural tunes to kickac """
    for plane in PLANES:
        kickac_df[COL_CORRECTED(plane)] = \
            kickac_df[COL_NATQ(plane)] - kickac_df[COL_MAV(plane)]
    return kickac_df


# Script Mode ##################################################################


if __name__ == '__main__':
    analyse_with_bbq_corrections()
