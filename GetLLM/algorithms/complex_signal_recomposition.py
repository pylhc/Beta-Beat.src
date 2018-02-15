from __future__ import print_function
import __init__ 
import os, sys
import numpy as np
from Utilities import tfs_pandas


RDT_LIST = {'H':['f1001H', 'f1010H',  
                 'f3000H', 'f1200H', 'f1020H', 'f1002H',  
                 'f1101H', 'f2010H', 'f1110H', 'f2001H', 
                 'f4000H', 'f1300H', 'f2002H', 'f1120H',
                 'f1102H', 'f2020H', 'f3001H', 'f1210H', 
                 'f1220H', 'f3002H', 'f1130H', 'f2003H'
            ],
            'V':['f0110V', 'f1010V',
                 'f0111V', 'f1020V', 'f0120V', 'f1011V', 
                 'f0030V', 'f0012V', 'f0210V', 'f2010V',  
                 'f0220V', 'f0211V', 'f0040V', 'f0013V',
                 'f2020V', 'f2011V', 'f0130V', 'f1012V',
            ]}


def get_bpm_pairs(phase_d):
    model_phase_advances = phase_d['MODEL']
    df_idx = model_phase_advances.index()
    for bpm in df_idx:
        last_bpm_of_range = df_idx[df_idx.get_loc(bpm)+4] 
        three_bpm_array = model_phase_advances.loc[bpm:last_bpm_of_range]
        bpm2 = (np.abs(three_bpm_array-0.25)).idxmin()
        phase_d['MEAS'].loc[bpm, 'BPM2'] = bpm2 
    return phase_d


def get_line_label(rdt):
    _, j, k, l, m, plane = list(rdt)
    if plane == 'H':
        p, q = k-j+1, m-l
        label = (str(p)+str(q)).replace('-', '_')
    elif plane == 'V':
        p, q = k-j, m-l+1
        label = (str(p)+str(q)).replace('-', '_')
    return label


def get_lines(rdts):
    lines_info = []
    for rdt in rdts:
        label = get_line_label(rdt)
        lines_info.append((label, p, q, rdt))
    return lines_info


def get_signal_df(spectral_df, phase_d, lines):
    tp = 2.0*np.pi
    for line in lines:
        label = line[0] 
        spectral_df['SIG'+label] =  spectral_df['AMP'+label]*complex(np.cos(tp*spectral_df['PHASE'+label]), 
                                                                    np.sin(tp*spectral_df['PHASE'+label]))
    signals_df = get_complex_signal(spectral_df, phase_d, lines)
    return signals_df


def get_complex_signal(spectral_df, phase_d, lines):
    signals_df = phase_d['BPM2']

    for line in lines:
        label = line[0]
        sig_bpm1 = spectral_df.loc[:, 'SIG'+label]*complex(1,  1/np.tan(tp*spectral_df['delta']))
        sig_bpm2 = spectral_df.loc[df.loc[:, 'BPM2'], 'SIG'+label]*complex(0, -1/np.sin(tp*spectral_df['delta'])) 
        signals_df['CSIG'+label] = sig_bpm1 + sig_bpm2
    return signals_df


def concatenate_dataframes_to_multiIndex(dfs):
    keys_dfs = np.arange(len(dfs))
    data_set = pd.concat(dfs, keys=keys_dfs).swaplevel(i=-2, j=-1, axis=0)
    return data_set


def rdt_function_factory(rdt):
    _, j, k, l, m, plane = list(rdt)
    j, k, l, m = int(j), int(k), int(l), int(m)
    if plane == 'H':
        def rdt_function(x, f):
            return 2 * j * f * x[0]**((j+k-2)/2.) * x[1]**((l+m)/2.)
    elif plane == 'V':
        def rdt_function(x, f):
            return 2 * l * f * x[0]**((j+k)/2.) * x[1]**((l+m-2)/2.)
    return rdt_function


def get_rdt_amp():
    pass


def get_rdt_phase():
    pass


def get_resonance_driving_terms(data_set, actions, lines_info):
    for idx, line in enumerate(lines_info):
        fit_func = rdt_function_factory(line[3])
        for bpm in data_set.index[0]:
            line_data = data_set['CSIG'+line[0]].loc[bpm]
            popt, perr = fit_rdt_vs_actions(func, actions, line_data)
            data_set[line[3]].loc[bpm] = popt[0]
            data_set[line[3]+'_err'].loc[bpm] = perr[0]
    return data_set
    

def fit_rdt_vs_actions(func, actions, line_data):
    popt, pcov = curve_fit(func, actions, line_data)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr


def process_measurements(twiss_d_zero_dpp, phase_d, lines):
    dfs = []
    for spectral_df in twiss_d_zero_dpp:
        complex_signals = get_signal_df(spectral_df, phase_d, lines)    
        dfs.append(complex_signals)
    data_set = concatenate_dataframes_to_multiIndex(dfs)
    return data_set


def do_main(phase_d, twiss_d):
    actions = np.vstack((np.transpose(kick_x)[0]**2, np.transpose(kick_y)[0]**2))
    phase_d = get_bpm_pairs(phase_d)
    
    for plane in ['H', 'V']:
        lines_info = get_lines(RDT_LIST[plane])
        if plane == 'H':
            twiss_d_zero_dpp = twiss_d.zero_dpp_x
        elif plane == 'V':
            twiss_d_zero_dpp = twiss_d.zero_dpp_y
        data_set = process_measurements(twiss_d_zero_dpp, phase_d, lines_info)
        get_resonance_driving_terms(data_set, actions, lines_info)


if __name__ == '__main__':
    do_main()

 
