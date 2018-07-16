import numpy as np
from parameter_config import get_time_col


TIME_COL = get_time_col()


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

    if fine_length:
        min_mask = data_series <= (data_mav - fine_cut)
        max_mask = data_series >= (data_mav + fine_cut)
        cut_mask = min_mask | max_mask
        data_mav = _get_interpolated_moving_average(data_mav, cut_mask, fine_length)

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


# Private methods ############################################################

def _get_interpolated_moving_average(data_series, clean_mask, length):
    """ Returns the moving average of data series with a window of length and interpolated NaNs"""
    data_mav = data_series.copy()
    data_mav[clean_mask] = np.NaN

    # 'interpolate' fills nan based on index/values of neighbours
    data_mav = data_mav.interpolate("index")

    shift = -int((length-1)/2)  # Shift average to middle value
    return data_mav.rolling(length).mean().shift(shift).interpolate("index")


# Script Mode #################################################################

if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
