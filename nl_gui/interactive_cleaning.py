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


class IterateCleaning(object):
    '''
    Author: Felix Carlier
    
    IterateCleaning
    '''
    def __init__(self, filenames):
        self.source = ['BBQ_HS_B1','BBQ_HS_B2']
        self.tune_df = load_csv(filenames[0], filetype='tune')
        self.platteaus_df = load_csv(filenames[1], filetype='platteaus')
        
        self.B1_keys = KEYS_DICT[self.source[0]]
        self.B2_keys = KEYS_DICT[self.source[1]]
        
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
        self._make_plot()
        self.all_poly = self._make_all_poly()
        #self._update_span()
        plt.show()
    
    def _make_all_poly(self):
        span_d = {'top_left': (0.26, 0.31), 'top_right':(0.26, 0.31), 'middle_left':(0.28,0.33), 
                     'middle_right':(0.28, 0.33), 'bottom_left':(0, 0.015), 'bottom_right':(0,0.015)}

        ''' Create a dictionary of axvspan Polygons to be plotted on all axes'''
        all_poly = {} 
        for key in self.axes:
            span_temp = self.axes[key].axvspan(span_d[key][0],span_d[key][1], facecolor='g', alpha=0.2, animated=True)
            all_poly[key] = PolygonInteractor(self.axes[key], span_temp)
        return all_poly

    def _get_next_platteau(self):
        ''' Iterate through the guessed platteau times to define new plotting xlimits and platteau limits. '''
        self.idx += 1
        self.time_knob = self.platteaus_df['KNOB_PLAT_START'][self.idx]-timedelta(seconds=50)
        self.time_min_pl = self.platteaus_df['DATA_PLAT_START'][self.idx]
        self.time_max_pl = self.platteaus_df['DATA_PLAT_END'][self.idx]
        self.time_max = self.time_max_pl+timedelta(seconds=50)
        self.new_limits_B1 = [self.time_min_pl, self.time_max_pl]
        self.new_limits_B2 = [self.time_min_pl, self.time_max_pl]

    def _next_plat_key(self, event):
        '''
        Defines actions on (n)-key press. 
        1) Append accepted platteau values to platteau lists for B1 and B2
        2) Try: Get next platteau times, update spans and clear plots, generate data plots
        '''
        if event.key == 'n':
            try:
    #            self._get_next_platteau()
                self._clear_plots()
                self._make_plot()
            except KeyError:
                wurstel

    def _clear_plots(self):
        '''Updates vertical spans and clears the axes for new plot'''
        #self._update_span()
        for key in self.axes:
            self.axes[key].clear()

    def _get_data_frames(self):
        self.b1x = self.tune_df[[B1_keys[0]]].loc[self.time_knob:self.time_max].dropna(how='all')
        self.b1y = self.tune_df[[B1_keys[1]]].loc[self.time_knob:self.time_max].dropna(how='all')
        self.b1c = self.tune_df[[B1_keys[2]]].loc[self.time_knob:self.time_max].dropna(how='all')
        self.b2x = self.tune_df[[B2_keys[0]]].loc[self.time_knob:self.time_max].dropna(how='all')
        self.b2y = self.tune_df[[B2_keys[1]]].loc[self.time_knob:self.time_max].dropna(how='all')
        self.b2c = self.tune_df[[B2_keys[2]]].loc[self.time_knob:self.time_max].dropna(how='all')

    def _update_poly(self):
        '''
        Updates vertical spans. Creates two arrays for the x values of Polygon.xy using the mask. 
        Note, there are 5 values in Polygon.xy as the rectangle has to be closed i.e. xy[0] =x xy[4]. Ylimits to 0 and 1 will span the whole vertical space.
        The datetime objects (new_limits) have to be converted to matplotlib numbers, which uses a different reference than Unix time -> 2001/01/01
        '''
        xmin_mask = np.array([1,1,0,0,1])
        xmax_mask = np.array([0,0,1,1,0])
        
        self.b1x_lim = [self.b1x.min(), self.b1x.max(), self.b1x.mean(), self.b1x.std()]
        self.b1y_lim = [self.b1y.min(), self.b1y.max(), self.b1y.mean(), self.b1y.std()]
        self.b1c_lim = [self.b1c.min(), self.b1c.max(), self.b1c.mean(), self.b1c.std()]
        self.b2x_lim = [self.b2x.min(), self.b2x.max(), self.b2x.mean(), self.b2x.std()]
        self.b2y_lim = [self.b2y.min(), self.b2y.max(), self.b2y.mean(), self.b2y.std()]
        self.b2c_lim = [self.b2c.min(), self.b2c.max(), self.b2c.mean(), self.b2c.std()]

        xlim_B1 = xmin_mask*mdates.date2num(self.new_limits_B1[0]) + xmax_mask*mdates.date2num(self.new_limits_B1[1])
        xlim_B2 = xmin_mask*mdates.date2num(self.new_limits_B2[0]) + xmax_mask*mdates.date2num(self.new_limits_B2[1])
        ylim = np.array([0,1,1,0,0])
        for key in ['top_left', 'middle_left', 'bottom_left']:
            self.all_poly[key].poly.xy = zip(xlim_B1, ylim)
        for key in ['top_right', 'middle_right', 'bottom_right']:
            self.all_poly[key].poly.xy = zip(xlim_B2, ylim)

    def _make_plot(self):
        '''
        Plot the dataframes for xing angles, circuit currents, and orbits using the given platteau timestamps.
        The time ranged used is larger than the actual guessed platteau to get a better plot.
        '''
        
        self.axes['top_left'].set_title('Beam 1')
        self.axes['top_right'].set_title('Beam 2')
        
        self.b1x.plot.hist( xlim=[self.b1x_lim[0],self.b1x_lim[1]],
                            bins=50, 
                            ax=self.axes['top_left'])
        self.b1y.plot.hist( xlim=[self.b1y_lim[0],self.b1y_lim[1]],
                            bins=50, 
                            ax=self.axes['middle_left'])
        self.b1c.plot.hist( xlim=[self.b1c_lim[0],self.b1c_lim[1]],
                            bins=50, 
                            ax=self.axes['bottom_left'])
        self.b2x.plot.hist( xlim=[self.b2x_lim[0],self.b2x_lim[1]],
                            bins=50, 
                            ax=self.axes['top_right'])
        self.b2y.plot.hist( xlim=[self.b2y_lim[0],self.b2y_lim[1]],
                            bins=50, 
                            ax=self.axes['middle_right'])
        self.b2c.plot.hist( xlim=[self.b2c_lim[0],self.b2c_lim[1]],
                            bins=50, 
                            ax=self.axes['bottom_right'])
        self.fig.canvas.draw()

    def _normalize_data(self, df):
        '''Removes first value bias, and normalizes by range of values'''
        df = df-df.iloc[0]
        df = df/(df.max()-df.min())
        return df


if __name__ == '__main__':
    input_dir = '/afs/cern.ch/work/f/fcarlier/public/data/NL_test_data/'
    tune_filename  = os.path.join(input_dir,'data.BBQ.csv')
    platteaus_filename = os.path.join(input_dir, 'accepted_platteaus.csv')
    filenames = [tune_filename, platteaus_filename]
    IterateCleaning(filenames)
