"""
Module tfs_files.diff_files
-------------------------------

Functions to get the difference between two tfs files or dataframes.
This is very general, i.e. not as results oriented as ``getdiff.py``.
"""
from tfs_files import tfs_pandas
import pandas as pd
import numpy as np


def get_diff_two_dataframes(df1, df2, diff_columns, is_error=None, prefix="", index=None, keep_colums=(), out_file=None):
    """ Get the difference of common elements of specific columns between two dataframes.

    Merges on index.

        Args:
            df1 (DataFrame or Path): First dataframe, Minuend
            df2 (DataFrame or Path): Second dataframe, Subtrahend
            diff_columns (list of stings): List of columns to get the difference of
            is_error (list of booleans): defines if the column in question is an error column
            prefix (str): Prefix for difference columns (default: "")
            index (str): index column - most likely needed when reading/writing files
            keep_colums (list of strings): additional columns to keep in the returned dataframe
            out_file (Path): if given, writes the result into this file

        Returns:
            DataFrame containing difference columns and kept columns.
    """
    # convert from files to dataframes
    df1 = _get_dataframe(df1, index)
    df2 = _get_dataframe(df2, index)

    # check input
    _check_for_missing_columns(df1, df2, diff_columns)
    _check_for_missing_columns(df1, df2, keep_colums)

    # merge dataframes
    merged = pd.merge(df1, df2, how='inner',
                      left_index=True, right_index=True,
                      suffixes=('_df1', '_df2'))

    # calculate difference
    if is_error is None:
        is_error = [False] * len(diff_columns)
    else:
        if len(is_error) != len(diff_columns):
            raise ValueError(
                "The length of the is_error switch needs to correspond to the columns."
            )

    for idx, col in enumerate(diff_columns):
        if is_error[idx]:
            merged[prefix + col] = .5 * np.sqrt(np.square(merged[col + '_df1'])
                                                + np.square(merged[col + '_df2']))
        else:
            merged[prefix + col] = merged[col + '_df1'] - merged[col + '_df2']

    # copy columns to be kept
    for col in keep_colums:
        for suffix in ["", "_df1", "_df2"]:
            try:
                merged[col] = merged[col + suffix]
            except KeyError:
                pass
            else:
                break

    # prepare output
    merged = merged.loc[:, keep_colums + [prefix + c for c in diff_columns]]
    if out_file:
        tfs_pandas.write_tfs(out_file, merged, save_index=index)
    return merged


# Helpers #####################################################################


def _get_dataframe(tfs_df, index):
    try:
        return tfs_pandas.read_tfs(tfs_df, index=index)
    except TypeError:
        return tfs_df


def _check_for_missing_columns(df1, df2, columns):
    missing_columns = [col for col in columns for df in [df1, df2] if col not in df.columns]
    if any(missing_columns):
        raise KeyError(
            "The following columns can not be found in either dataframe: {:}".format(
                list(set(missing_columns)))
        )
