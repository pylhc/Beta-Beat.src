from __future__ import print_function
import os
import sys
import numpy as np
import pandas as pd
from pandas import ExcelWriter
from datetime import datetime


def write_excel(outpath, data):
    writer = ExcelWriter(outpath)
    data.to_excel(writer, 'Sheet1')
    writer.save()


def load_currents_csv(file_path):
    data_frame = pd.read_csv(file_path, parse_dates=[0], infer_datetime_format=True)
    data_frame.sort_values(by='time', ascending=True, inplace=True)
    data_frame = data_frame.fillna(method='ffill')
    data_frame = data_frame.fillna(method='bfill')
    data_frame.index = pd.to_datetime(data_frame.pop('time'))
    return data_frame


def load_orbit_csv(file_path):
    data_frame = pd.read_csv(file_path, parse_dates=[0], infer_datetime_format=True)
    data_frame.columns = map(str.upper, data_frame.columns)
    data_frame.sort_values(by='TIME', ascending=True, inplace=True)
    data_frame.index = pd.to_datetime(data_frame.pop('TIME'))
    return data_frame


def load_platteaus_csv(file_path):
    data_frame = pd.read_csv(file_path, parse_dates=[0,1,2], infer_datetime_format=True)
    data_frame.columns = map(str.upper, data_frame.columns)
    data_frame.sort_values(by='KNOB_PLAT_START', ascending=True, inplace=True)
    return data_frame


def load_mcb_csv(file_path):
    data_frame = pd.read_csv(file_path, parse_dates=[0], infer_datetime_format=True)
    data_frame = data_frame.replace(to_replace='RUNNING', value=True)
    data_frame = data_frame.replace(to_replace='ARMED', value=False)
    data_frame = data_frame.replace(to_replace='IDLE', value=False)
    data_frame = data_frame.set_index(['circuit', 'time'])
    data_frame = data_frame.stack().unstack(0)
    data_frame.index = data_frame.index.droplevel(level=1)
    data_frame = data_frame.any(axis=1)
    return data_frame


def load_csv(file_path, filetype=None):
    if filetype is None:
        raise RuntimeError('The filetype is not defined for loading this csv file: ', file_path)
    elif filetype is 'currents':
        return load_currents_csv(file_path)
    elif filetype is 'orbit':
        return load_orbit_csv(file_path)
    elif filetype is 'platteaus':
        return load_platteaus_csv(file_path)
    elif filetype is 'mcb':
        return load_mcb_csv(file_path)


if __name__ == "__main__":
    sys.exit(-1)



