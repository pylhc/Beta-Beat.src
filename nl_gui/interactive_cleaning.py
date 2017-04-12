#!/afs/cern.ch/work/o/omc/anaconda/bin/python

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


KEYS_DICT = {
'BOSFU':['LHC.BOFSU:EIGEN_FREQ_1_B1',
         'LHC.BOFSU:EIGEN_FREQ_2_B1',
         'LHC.BOFSU:COUPLING_ABS_B1',
         'LHC.BOFSU:EIGEN_FREQ_1_B2',
         'LHC.BOFSU:EIGEN_FREQ_2_B2',
         'LHC.BOFSU:COUPLING_ABS_B2'],
'BBQ':['LHC.BQBBQ.CONTINUOUS.B1:EIGEN_FREQ_1',
       'LHC.BQBBQ.CONTINUOUS.B1:EIGEN_FREQ_2',
       'LHC.BQBBQ.CONTINUOUS.B1:COUPLING_ABS',
       'LHC.BQBBQ.CONTINUOUS.B2:EIGEN_FREQ_1',
       'LHC.BQBBQ.CONTINUOUS.B2:EIGEN_FREQ_2',
       'LHC.BQBBQ.CONTINUOUS.B2:COUPLING_ABS'],
'BBQ_HS':['LHC.BQBBQ.CONTINUOUS_HS.B1:EIGEN_FREQ_1',
          'LHC.BQBBQ.CONTINUOUS_HS.B1:EIGEN_FREQ_2',
          'LHC.BQBBQ.CONTINUOUS_HS.B1:COUPLING_ABS',
          'LHC.BQBBQ.CONTINUOUS_HS.B2:EIGEN_FREQ_1',
          'LHC.BQBBQ.CONTINUOUS_HS.B2:EIGEN_FREQ_2',
          'LHC.BQBBQ.CONTINUOUS_HS.B2:COUPLING_ABS']
             }


class IterateCleaning(object):
    '''
    Author: Felix Carlier
    IterateCleaning
    '''
    def __init__(self, filenames):
        source = 'BBQ_HS'
        self.tune_df = load_csv(filenames[0], filetype='tune')
        self.platteaus_df = load_csv(filenames[1], filetype='platteaus')
        
        keys = KEYS_DICT[source]
        self.data_keys = {  'top_left':     keys[0],
                            'top_right':    keys[3], 
                            'middle_left':  keys[1], 
                            'middle_right': keys[4], 
                            'bottom_left':  keys[2], 
                            'bottom_right': keys[5]}
        self.cropped_data = {}
        self.hist_limits = {}
        self.span_limits = {}
        
        self.idx = 1 

        self.fig = plt.figure(figsize=(18,18))
        self.fig.patch.set_facecolor('white')
        self.fig.canvas.mpl_connect('key_press_event', self._next_plat_key)

        ax1 = self.fig.add_subplot(321)
        ax3 = self.fig.add_subplot(323)
        ax5 = self.fig.add_subplot(325)
        ax2 = self.fig.add_subplot(322)
        ax4 = self.fig.add_subplot(324)
        ax6 = self.fig.add_subplot(326)
        self.axes = {'top_left': ax1, 'top_right':ax2, 'middle_left':ax3, 
                     'middle_right':ax4, 'bottom_left':ax5, 'bottom_right':ax6}
        
        self._get_next_platteau()
        self._get_data_frames()
        self._get_hist_limits()
        self._make_plot()
        self.all_poly = self._make_all_poly()
        plt.show()
    
    def _get_next_platteau(self):
        self.idx += 1
        self.time_knob = self.platteaus_df['KNOB_PLAT_START'][self.idx]-timedelta(seconds=50)
        self.time_min_pl = self.platteaus_df['DATA_PLAT_START'][self.idx]
        self.time_max_pl = self.platteaus_df['DATA_PLAT_END'][self.idx]
        self.time_max = self.time_max_pl+timedelta(seconds=50)
        self.new_limits_B1 = [self.time_min_pl, self.time_max_pl]
        self.new_limits_B2 = [self.time_min_pl, self.time_max_pl]

    def _get_data_frames(self):
        for key in self.axes:
            cr_data = self.tune_df[[self.data_keys[key]]].loc[self.time_knob:self.time_max].dropna(how='all')
            self.cropped_data[key] = cr_data

    def _get_hist_limits(self):
        for key in self.axes:
            dfmin = self.cropped_data[key].min()[0]
            dfmax = self.cropped_data[key].max()[0]
            dfmean = self.cropped_data[key].mean()[0]
            dfstd = self.cropped_data[key].std()[0]
            
            self.hist_limits[key] = (dfmin, dfmax, dfmean, dfstd)            
            self.span_limits[key] = (dfmean-dfstd, dfmean+dfstd)

    def _make_plot(self):
        self.axes['top_left'].set_title('Beam 1')
        self.axes['top_right'].set_title('Beam 2')
        for key in self.axes:
            self.cropped_data[key].plot.hist(xlim=[self.hist_limits[key][0],self.hist_limits[key][1]],
                                             bins=200, 
                                             ax=self.axes[key])
        self.fig.canvas.draw()

    def _make_all_poly(self):
        all_poly = {} 
        for key in self.axes:
            span_temp = self.axes[key].axvspan(self.span_limits[key][0],self.span_limits[key][1], facecolor='g', alpha=0.2, animated=True)
            all_poly[key] = PolygonInteractor(self.axes[key], span_temp)
        return all_poly

    def _next_plat_key(self, event):
        if event.key == 'n':
            try:
                self._get_next_platteau()
                self._clear_plots()
                self._make_plot()
            except KeyError:
                wurstel

    def _clear_plots(self):
        self._get_data_frames()
        self._get_hist_limits()
        self._update_poly()
        for key in self.axes:
            self.axes[key].clear()

    def _update_poly(self):
        xmin_mask = np.array([1,1,0,0,1])
        xmax_mask = np.array([0,0,1,1,0])
        ylim = np.array([0,1,1,0,0])
        
        for key in self.span_limits:
            xlim = xmin_mask*self.span_limits[key][0] + xmax_mask*self.span_limits[key][1]
            self.all_poly[key].poly.xy = zip(xlim, ylim)

if __name__ == '__main__':
    input_dir = '/afs/cern.ch/work/f/fcarlier/public/data/NL_test_data/'
    tune_filename  = os.path.join(input_dir,'data.BBQ.csv')
    platteaus_filename = os.path.join(input_dir, 'accepted_platteaus.csv')
    filenames = [tune_filename, platteaus_filename]
    IterateCleaning(filenames)
