'''
Description: Do statitstic works on shunyi X band Radar
Author: Hejun Xie
Date: 2021-10-06 15:08:44
LastEditors: Hejun Xie
LastEditTime: 2021-10-07 17:07:39
'''

import pycwr
import numpy as np
from numpy import ma
from scipy.ndimage import gaussian_filter
import datetime as dt
import glob
import pyproj as pp
import pickle

from pycwr.io import read_auto

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

from ZJU_AERO.beam import compute_trajectory_radial

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'

# define the ground based radar geo. coords here
shunyi  = [116.6, 40.13]
haidian = [116.2833333, 39.98333333]
yanqing = [115.9666667, 40.45]
pinggu  = [117.1166667, 40.16666667]
nanjiao = [116.4666667, 39.8]

# define the vicinity range (half range)
vicinity = (2, 10)  # (tangential, radial)

def data_smooth(data, mask, sigma, truncate=30.0):
    data[mask] = 0.0
    if sigma is not None:
        data = gaussian_filter(data, sigma=sigma, truncate=truncate)
        
        # get WW as valid guassian weight
        W = data.copy() * 0.0 + 1.0
        W[mask] = 0.0
        WW = gaussian_filter(W, sigma=sigma, truncate=truncate) + 1e-10 # to avoild divide by zero error

        data = data / WW
    data[mask] = np.nan

    return data

def sweep_smooth(pyart_radar, start_radial, end_radial, field, sigma=None):
    sweep_data = pyart_radar.fields[field]['data'][start_radial:end_radial, :].data
    sweep_mask = pyart_radar.fields[field]['data'][start_radial:end_radial, :].mask
    sweep_data = data_smooth(sweep_data, sweep_mask, sigma)
    pyart_radar.fields[field]['data'][start_radial:end_radial, :] = sweep_data

def get_sweep_field(pyart_radar, start_radial, end_radial, field):
    return pyart_radar.fields[field]['data'][start_radial:end_radial, :]

def sweep_filter(pyart_radar, field_filtered, field_mask, threshold, start_radial, end_radial, lessthan=True):
    sweep_mask = pyart_radar.fields[field_mask]['data'][start_radial:end_radial, :].data
    if lessthan:
        mask = sweep_mask < threshold
    else:
        mask = sweep_mask >= threshold

    sweep_filtered = pyart_radar.fields[field_filtered]['data'][start_radial:end_radial, :]
    sweep_filtered[mask] = np.nan
    pyart_radar.fields[field_filtered]['data'][start_radial:end_radial, :] = \
        ma.array(sweep_filtered, mask=mask | sweep_filtered.mask)

def preprocess(radar_file, sweep):
    PRD = read_auto(radar_file)
    r = PRD.ToPyartRadar()

    radial_start_index = r.sweep_start_ray_index['data'][sweep]
    radial_end_index = r.sweep_end_ray_index['data'][sweep] + 1

    # print(PRD.fields[sweep])
    # print(PRD.fields[sweep].coords['azimuth'].data)
    # print(PRD.fields[sweep].coords['range'].data)
    # print(r.azimuth['data'][radial_start_index:radial_end_index])
    # print(r.range['data'])

    range_vec = r.range['data']
    elevation_angle = r.elevation['data'][radial_start_index]
    azimuth_angles = r.azimuth['data'][radial_start_index:radial_end_index]
    shunyi_coords = [shunyi[1], shunyi[0], 40.0]
    s, h, e = compute_trajectory_radial(range_vec, elevation_angle, shunyi_coords, 1)

    sweep_smooth(r, radial_start_index, radial_end_index, 'reflectivity', sigma=[2.0, 15.0])
    sweep_smooth(r, radial_start_index, radial_end_index, 'cross_correlation_ratio', sigma=[2.0, 15.0])
    zh      = get_sweep_field(r, radial_start_index, radial_end_index, 'reflectivity')
    rhohv   = get_sweep_field(r, radial_start_index, radial_end_index, 'cross_correlation_ratio')

    sweep_filter(r, 'differential_reflectivity', 'reflectivity', 20.0, radial_start_index, radial_end_index)
    sweep_filter(r, 'differential_reflectivity', 'cross_correlation_ratio', 0.95, radial_start_index, radial_end_index)

    sweep_smooth(r, radial_start_index, radial_end_index, 'differential_reflectivity', sigma=[2.0, 15.0])
    zdr     = get_sweep_field(r, radial_start_index, radial_end_index, 'differential_reflectivity')

    del r
    return zdr, zh, s, azimuth_angles

def get_ad(geoid, coords1, coords2):
    '''
    get azimuth and distance from coords1 to coords2 
    '''
    az12, az21, d = geoid.inv(coords1[0], coords1[1], coords2[0], coords2[1])
    if az12 < 0.0:
        az12 = az12 + 360.0
    return az12, d

def get_vicinity(idx_site, field):

    tempfield = np.zeros((field.shape[0]*3, field.shape[1]), dtype='float')
    tempfield[0:field.shape[0], :] = field
    tempfield[field.shape[0]:2*field.shape[0], :] = field
    tempfield[2*field.shape[0]:3*field.shape[0], :] = field

    vicinity_field = tempfield[field.shape[0]+idx_site[0]-vicinity[0]:field.shape[0]+idx_site[0]+vicinity[0]+1, \
                               idx_site[1]-vicinity[1]:idx_site[1]+vicinity[1]+1]

    return vicinity_field

def stats_radar_file(radar_file, sweep):
    
    zdr, zh, s_coords, azi_coords = preprocess(radar_file, sweep)

    geoid = pp.Geod(ellps = 'WGS84')
    a_yanqing, s_yanqing = get_ad(geoid, shunyi, yanqing)
    a_pinggu, s_pinggu   = get_ad(geoid, shunyi, pinggu)
    a_haidian, s_haidian = get_ad(geoid, shunyi, haidian)
    a_nanjiao, s_nanjiao = get_ad(geoid, shunyi, nanjiao)

    idx_yanqing = np.argmin(np.abs(azi_coords-a_yanqing)), np.argmin(np.abs(s_coords-s_yanqing))
    idx_pinggu = np.argmin(np.abs(azi_coords-a_pinggu)), np.argmin(np.abs(s_coords-s_pinggu))
    idx_haidian = np.argmin(np.abs(azi_coords-a_haidian)), np.argmin(np.abs(s_coords-s_haidian))
    idx_nanjiao = np.argmin(np.abs(azi_coords-a_nanjiao)), np.argmin(np.abs(s_coords-s_nanjiao))

    YQ = get_vicinity(idx_yanqing, zh), get_vicinity(idx_yanqing, zdr)
    PG = get_vicinity(idx_pinggu, zh), get_vicinity(idx_pinggu, zdr)
    HD = get_vicinity(idx_haidian, zh), get_vicinity(idx_haidian, zdr)
    NJ = get_vicinity(idx_nanjiao, zh), get_vicinity(idx_nanjiao, zdr)

    result = dict()
    result['54406-Yanqing'] = YQ
    result['54399-Haidian'] = HD
    result['54424-Pinggu'] = PG
    result['54511-Nanjiao'] = NJ
    return result

def stats_time_window(result_box):
    nwindows = len(result_box)
    result = dict()
    result['54406-Yanqing']     = np.zeros((nwindows, 4), dtype='float') + np.nan
    result['54399-Haidian']     = np.zeros((nwindows, 4), dtype='float') + np.nan
    result['54424-Pinggu']      = np.zeros((nwindows, 4), dtype='float') + np.nan
    result['54511-Nanjiao']     = np.zeros((nwindows, 4), dtype='float') + np.nan

    for iwindow in range(nwindows):
        for site in result.keys():
            ZHs = [result_box[iwindow][i][site][0] for i in range(len(result_box[iwindow]))]
            ZDRs = [result_box[iwindow][i][site][1] for i in range(len(result_box[iwindow]))]
            if len(ZHs) != 0:
                ZH, ZDR = np.stack(ZHs), np.stack(ZDRs)
                result[site][iwindow,0] = np.nanmean(ZH)
                result[site][iwindow,1] = np.nanstd(ZH)
                result[site][iwindow,2] = np.nanmean(ZDR)
                result[site][iwindow,3] = np.nanstd(ZDR)

    return result

def plot_with_nan(ax, x, y, yerr, **kwargs):
    mask = ~np.isnan(y)
    # ax.plot(x[mask], y[mask], **kwargs)

    kwargs['capsize'] = 5
    kwargs['elinewidth'] = 0.8
    eb = ax.errorbar(x[mask], y[mask], yerr=yerr[mask], **kwargs)
    eb[-1][0].set_linestyle(kwargs['ls'])


LOAD_RESULT = False

if __name__ == "__main__":

    # initialize time windows
    time_windows = np.arange(dt.datetime(2019,11,29,11), dt.datetime(2019,11,29,17), \
                                 dt.timedelta(seconds=1800)).astype(dt.datetime)

    if not LOAD_RESULT:
        radar_files = glob.glob('../../pathos/RADAR/19.11.29顺义X波段/29/BJXSY.20191129.??????.AR2.bz2')
        radar_files = radar_files[200:300]
        # radar_files = [radar_files[199]]
        # radar_files = radar_files[200:210]
        
        # initialize the result_box
        result_box = []
        for time_window in time_windows:
            result_box.append([])
        
        for ifile, radar_file in enumerate(radar_files):
            yyyymmdd = radar_file.split('.')[-4]
            hhmmss = radar_file.split('.')[-3]
            current_time = dt.datetime(int(yyyymmdd[0:4]),int(yyyymmdd[4:6]),int(yyyymmdd[6:8]),
                                       int(hhmmss[0:2]), int(hhmmss[2:4]), int(hhmmss[4:6]))
            print(current_time)
            window_idx = np.argmin(np.abs(current_time - time_windows))
            
            result_1time = stats_radar_file(radar_file, 1)
            result_box[window_idx].append(result_1time)

        with open("./result_box.pkl", "wb") as f:
            pickle.dump(result_box, f)
    else:
        with open("./result_box.pkl", "rb") as f:
            result_box = pickle.load(f)
    
    result = stats_time_window(result_box)

    fig = plt.figure(figsize=(10,12))
    time = np.arange(result['54406-Yanqing'].shape[0])

    xticklabels = [time_window.strftime('%Y-%m-%d %H') for time_window in time_windows[::2]]

    ax = plt.subplot(211)

    plot_with_nan(ax, time, result['54406-Yanqing'][:,0], result['54406-Yanqing'][:,1], lw=1.5, ls='-', c='k', label='54406-Yanqing', marker='x')
    plot_with_nan(ax, time, result['54424-Pinggu'][:,0], result['54424-Pinggu'][:,1],   lw=1.5, ls='--', c='k', label='54424-Pinggu', marker='x')
    plot_with_nan(ax, time, result['54399-Haidian'][:,0], result['54399-Haidian'][:,1], lw=1.5, ls='-', c='r', label='54399-Haidian', marker='x')
    plot_with_nan(ax, time, result['54511-Nanjiao'][:,0], result['54511-Nanjiao'][:,1], lw=1.5, ls='--', c='r', label='54511-Nanjiao', marker='x')

    ax.set_xticks(time[::2])
    ax.set_xticklabels(xticklabels, fontsize=9, rotation=30)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Reflectivity ZH dBZ]', fontsize=12)

    ax = plt.subplot(212)

    plot_with_nan(ax, time, result['54406-Yanqing'][:,2], result['54406-Yanqing'][:,3], lw=1.5, ls='-', c='k', label='54406-Yanqing', marker='x')
    plot_with_nan(ax, time, result['54424-Pinggu'][:,2], result['54424-Pinggu'][:,3],   lw=1.5, ls='--', c='k', label='54424-Pinggu', marker='x')
    plot_with_nan(ax, time, result['54399-Haidian'][:,2], result['54399-Haidian'][:,3], lw=1.5, ls='-', c='r', label='54399-Haidian', marker='x')
    plot_with_nan(ax, time, result['54511-Nanjiao'][:,2], result['54511-Nanjiao'][:,3], lw=1.5, ls='--', c='r', label='54511-Nanjiao', marker='x')

    ax.set_xticks(time[::2])
    ax.set_xticklabels(xticklabels, fontsize=9, rotation=30)

    ax.legend(loc='best', frameon=False, fontsize=12)
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Differential Reflectivity ZDR [dBZ]', fontsize=12)

    fig.tight_layout()
    plt.savefig('CloudRadarSiteStats.png', dpi=300)    
