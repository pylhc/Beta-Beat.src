#!/afs/cern.ch/work/o/omc/anaconda/bin/python

from __future__ import print_function
import sys, os
import numpy as np
import argparse
import matplotlib.pyplot as plt

import matplotlib
matplotlib.style.use('ggplot')

from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
from matplotlib.axes import Axes
from matplotlib.patches import Polygon
import matplotlib.dates as mdates
import pandas as pd
from data_loader import load_csv 
from datetime import datetime
from datetime import timedelta 
from time import mktime
from polygon_interacter import PolygonInteractor 


KEYS_DICT_B1 = {
'BOSFU':['LHC.BOFSU:EIGEN_FREQ_1_B1',
         'LHC.BOFSU:EIGEN_FREQ_2_B1',
         'LHC.BOFSU:COUPLING_ABS_B1'],
'BBQ':['LHC.BQBBQ.CONTINUOUS.B1:EIGEN_FREQ_1',
       'LHC.BQBBQ.CONTINUOUS.B1:EIGEN_FREQ_2',
       'LHC.BQBBQ.CONTINUOUS.B1:COUPLING_ABS'],
'BBQ_HS':['LHC.BQBBQ.CONTINUOUS_HS.B1:EIGEN_FREQ_1',
          'LHC.BQBBQ.CONTINUOUS_HS.B1:EIGEN_FREQ_2',
          'LHC.BQBBQ.CONTINUOUS_HS.B1:COUPLING_ABS']
             }


KEYS_DICT_B2 = {
'BOSFU':['LHC.BOFSU:EIGEN_FREQ_1_B2',
         'LHC.BOFSU:EIGEN_FREQ_2_B2',
         'LHC.BOFSU:COUPLING_ABS_B2'],
'BBQ':['LHC.BQBBQ.CONTINUOUS.B2:EIGEN_FREQ_1',
       'LHC.BQBBQ.CONTINUOUS.B2:EIGEN_FREQ_2',
       'LHC.BQBBQ.CONTINUOUS.B2:COUPLING_ABS'],
'BBQ_HS':['LHC.BQBBQ.CONTINUOUS_HS.B2:EIGEN_FREQ_1',
          'LHC.BQBBQ.CONTINUOUS_HS.B2:EIGEN_FREQ_2',
          'LHC.BQBBQ.CONTINUOUS_HS.B2:COUPLING_ABS']
             }


class IterateCleaning(object):
    '''
    Author: Felix Carlier
    IterateCleaning
    '''
    def __init__(self, filenames, output_cleaning_file, beam, source):
        self.tune_df = load_csv(filenames[0], filetype='tune')
        self.platteaus_df = load_csv(filenames[1], filetype='accepted_platteaus')
        
        self.data_summary = pd.DataFrame(index=np.arange(len(self.platteaus_df['B1_min'])), columns=['Qx_ave', 'Qx_std', 'Qy_ave', 'Qy_std', 'Coupl_ave', 'Coupl_std'])
        self.output_cleaning_file = output_cleaning_file
        self.beam = beam
        if self.beam == 1:
            keys = KEYS_DICT_B1[source]
        elif self.beam == 2:
            keys = KEYS_DICT_B2[source]
        
        self.data_keys = {  'top_left':     keys[0],
                            'middle_left':  keys[1], 
                            'bottom_left':  keys[2]}
        self.cropped_data = {}
        self.hist_limits = {}
        self.span_limits = {}
        
        self.idx = 0 

        self.fig = plt.figure(figsize=(18,18))
        self.fig.patch.set_facecolor('white')
        self.fig.canvas.mpl_connect('key_press_event', self._next_plat_key)
        self.fig.canvas.mpl_connect('button_press_event', self._determine_poly)

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
   
    def _determine_poly(self, event):
        if event.inaxes is self.axes['top_left']:
            self.all_poly['top_left'].button_press_callback(event)
        if event.inaxes is self.axes['middle_left']:
            self.all_poly['middle_left'].button_press_callback(event)
        if event.inaxes is self.axes['bottom_left']:
            self.all_poly['bottom_left'].button_press_callback(event)

    def _get_next_platteau(self):
        if self.beam == 1:
            self.start = self.platteaus_df['B1_min'][self.idx]
            self.end = self.platteaus_df['B1_max'][self.idx]
        elif self.beam == 2:
            self.start = self.platteaus_df['B2_min'][self.idx]
            self.end = self.platteaus_df['B2_max'][self.idx]
        
    def _get_data_frames(self):
        for key in ['top_left', 'middle_left', 'bottom_left']:
            self.cropped_data[key] = self.tune_df[[self.data_keys[key]]].loc[self.start:self.end].dropna(how='all')
        self.cropped_data['top_right'] = self.tune_df[[self.data_keys['top_left'], self.data_keys['middle_left']]].loc[self.start:self.end].dropna(how='all')
        self.cropped_data['middle_right'] = self.tune_df[[self.data_keys['bottom_left']]].loc[self.start:self.end].dropna(how='all')
        self.cropped_data['bottom_right'] = self.tune_df[[self.data_keys['bottom_left']]].loc[self.start:self.end].dropna(how='all')

    def _get_hist_limits(self):
        for key in self.axes:
            dfmin = self.cropped_data[key].min()[0]
            dfmax = self.cropped_data[key].max()[0]
            dfmean = self.cropped_data[key].mean()[0]
            dfstd = self.cropped_data[key].std()[0]
            
            self.hist_limits[key] = (dfmin, dfmax, dfmean, dfstd)            
            self.span_limits[key] = (dfmean-dfstd, dfmean+dfstd)
    
    def _get_clean_limits(self):
        self.clean_lim_Qx = [min(self.all_poly['top_left'].poly.xy[:,0]), max(self.all_poly['top_left'].poly.xy[:,0])]
        self.clean_lim_Qy = [min(self.all_poly['middle_left'].poly.xy[:,0]), max(self.all_poly['middle_left'].poly.xy[:,0])]
        self.clean_lim_Coupl = [min(self.all_poly['bottom_left'].poly.xy[:,0]), max(self.all_poly['bottom_left'].poly.xy[:,0])]
    
    def _summarize_cleaned_data(self):
        self._get_clean_limits()
        qx_data = self.cropped_data['top_left'].clip(lower=self.clean_lim_Qx[0], upper=self.clean_lim_Qx[1])
        qy_data = self.cropped_data['middle_left'].clip(lower=self.clean_lim_Qy[0], upper=self.clean_lim_Qy[1])
        coupl_data = self.cropped_data['bottom_left'].clip(lower=self.clean_lim_Coupl[0], upper=self.clean_lim_Coupl[1])
       
        self.data_summary.loc[self.idx] = qx_data.mean()[0], qx_data.std()[0], qy_data.mean()[0], qy_data.std()[0], coupl_data.mean()[0], coupl_data.std()[0], 

    def _make_plot(self):
        self.axes['top_left'].set_title('Beam 1')
        self.axes['top_right'].set_title('Beam 2')
        for key in ['top_left', 'middle_left', 'bottom_left']:
            self.cropped_data[key].plot.hist(xlim=[self.hist_limits[key][0],self.hist_limits[key][1]],
                                             bins=200, 
                                             ax=self.axes[key])
        for key in ['top_right', 'middle_right', 'bottom_right']:
            self.cropped_data[key].plot(ax=self.axes[key])
        self.fig.canvas.draw()

    def _make_all_poly(self):
        all_poly = {} 
        for key in ['top_left', 'middle_left', 'bottom_left']:
            span_temp = self.axes[key].axvspan(self.span_limits[key][0],self.span_limits[key][1], facecolor='g', alpha=0.2, animated=True)
            all_poly[key] = PolygonInteractor(self.axes[key], span_temp)
        return all_poly

    def _next_plat_key(self, event):
        if event.key == 'n':
            self._summarize_cleaned_data()
            #try:
            self.idx += 1
            try:
                self._get_next_platteau()
                self._clear_plots()
                self._make_plot()
            except KeyError:
                print('Last platteau cleaned, no more platteaus.') 
                self.data_summary.to_csv(self.output_cleaning_file)
                print('Cleaned data summary written to: ', self.output_cleaning_file)
                plt.close()

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
        
        for key in ['top_left', 'middle_left', 'bottom_left']:
            xlim = xmin_mask*self.span_limits[key][0] + xmax_mask*self.span_limits[key][1]
            self.all_poly[key].poly.xy = zip(xlim, ylim)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir", help="Input directory.",
        dest='input_dir',type=str,
    )
    parser.add_argument(
        "--beam", help="Beam to be analysed.",
        dest='beam',type=int,
    )
    parser.add_argument(
            "--source", help="Source for tune and coupling measurements: BOFSU, BBQ, BBQ_HS",
        dest='source',type=str,
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    input_dir = args.input_dir
    beam = args.beam
    source = args.source

    tune_filename  = os.path.join(input_dir,'data.BBQ.csv')
    platteaus_filename = os.path.join(input_dir,'accepted_platteaus.dat')
    
    filenames = [tune_filename, platteaus_filename]
    output_cleaning_file = os.path.join(input_dir, 'analysed_data_'+source+'.dat')
    IterateCleaning(filenames, output_cleaning_file, beam, source)
