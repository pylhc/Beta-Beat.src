
from __future__ import print_function
import sys, os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use('ggplot')

from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
from matplotlib.axes import Axes
from matplotlib.patches import Polygon
import matplotlib.dates as mdates
    
from data_loader import load_csv 
from datetime import datetime
from datetime import timedelta 
from time import mktime
from polygon_interacter import PolygonInteractor 
import pandas as pd


KEYS_DICT = {
'BOSFU_B1':['LHC.BOFSU:EIGEN_FREQ_1_B1',
            'LHC.BOFSU:EIGEN_FREQ_2_B1',
            'LHC.BOFSU:COUPLING_ABS_B1'],
'BOSFU_B2':['LHC.BOFSU:EIGEN_FREQ_1_B2',
            'LHC.BOFSU:EIGEN_FREQ_2_B2',
            'LHC.BOFSU:COUPLING_ABS_B2'],
'BBQ_B1':[  'LHC.BQBBQ.CONTINUOUS.B1:EIGEN_FREQ_1',
            'LHC.BQBBQ.CONTINUOUS.B1:EIGEN_FREQ_2',
            'LHC.BQBBQ.CONTINUOUS.B1:COUPLING_ABS'],
'BBQ_B2':[  'LHC.BQBBQ.CONTINUOUS.B2:EIGEN_FREQ_1',
            'LHC.BQBBQ.CONTINUOUS.B2:EIGEN_FREQ_2',
            'LHC.BQBBQ.CONTINUOUS.B2:COUPLING_ABS'],
'BBQ_HS_B1':['LHC.BQBBQ.CONTINUOUS_HS.B1:EIGEN_FREQ_1',
             'LHC.BQBBQ.CONTINUOUS_HS.B1:EIGEN_FREQ_2',
             'LHC.BQBBQ.CONTINUOUS_HS.B1:COUPLING_ABS'],
'BBQ_HS_B2':['LHC.BQBBQ.CONTINUOUS_HS.B2:EIGEN_FREQ_1',
             'LHC.BQBBQ.CONTINUOUS_HS.B2:EIGEN_FREQ_2',
             'LHC.BQBBQ.CONTINUOUS_HS.B2:COUPLING_ABS']
             }


def load_tune_coupling_csv(file_path):
    data_frame = pd.read_csv(file_path, parse_dates=[0], infer_datetime_format=True)
    data_frame.set_index(['time'], inplace=True)
    return data_frame


if __name__ == '__main__':
    input_dir = '/afs/cern.ch/work/f/fcarlier/public/data/NL_test_data/'
    tune_filename  = os.path.join(input_dir,'data.BBQ.csv')
    platteaus_filename = os.path.join(input_dir, 'accepted_platteaus.csv')
    
    df = load_tune_coupling_csv(tune_filename)
    print (df)
    
    keys_B1 = [key for key in df.columns if 'B1' in key]
    keys_B2 = [key for key in df.columns if 'B2' in key]
   
    coupl_keys_B1 = [key for key in keys_B1 if 'COUPL' in key]
    coupl_keys_B2 = [key for key in keys_B2 if 'COUPL' in key]

    B1_keys = KEYS_DICT['BOSFU_B1']
    B2_keys = KEYS_DICT['BOSFU_B2']
    print (coupl_keys_B2)
    df[[B1_keys[0]]].dropna(how='all').plot.hist(bins=80)
    #df[['LHC.BOFSU:EIGEN_FREQ_1_B1']].dropna(how='all').plot.hist(bins=80)
    plt.show()


