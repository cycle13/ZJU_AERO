'''
Description: Do statitstic works on shunyi X band Radar
Author: Hejun Xie
Date: 2021-10-06 15:08:44
LastEditors: Hejun Xie
LastEditTime: 2021-10-06 19:16:37
'''

import pycwr
import numpy as np
from numpy import ma
from scipy.ndimage import gaussian_filter
import datetime as dt
import glob

from pycwr.io import read_auto

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'


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

    sweep_smooth(r, radial_start_index, radial_end_index, 'reflectivity', sigma=[2.0, 15.0])
    sweep_smooth(r, radial_start_index, radial_end_index, 'cross_correlation_ratio', sigma=[2.0, 15.0])
    zh      = get_sweep_field(r, radial_start_index, radial_end_index, 'reflectivity').flatten()
    rhohv   = get_sweep_field(r, radial_start_index, radial_end_index, 'cross_correlation_ratio').flatten()

    sweep_filter(r, 'differential_reflectivity', 'reflectivity', 20.0, radial_start_index, radial_end_index)
    sweep_filter(r, 'differential_reflectivity', 'cross_correlation_ratio', 0.95, radial_start_index, radial_end_index)

    sweep_smooth(r, radial_start_index, radial_end_index, 'differential_reflectivity', sigma=[2.0, 15.0])
    zdr     = get_sweep_field(r, radial_start_index, radial_end_index, 'differential_reflectivity').flatten()

    del r
    return zdr
    
def stats_radar_file(radar_file, sweep):

    yyyymmdd = radar_file.split('.')[-4]
    hhmmss = radar_file.split('.')[-3]

    current_time = dt.datetime(int(yyyymmdd[0:4]),int(yyyymmdd[4:6]),int(yyyymmdd[6:8]),
        int(hhmmss[0:2]), int(hhmmss[2:4]), int(hhmmss[4:6]))

    print(current_time)
    zdr = preprocess(radar_file, sweep)

if __name__ == "__main__":

    radar_files = glob.glob('../../pathos/RADAR/19.11.29顺义X波段/29/BJXSY.20191129.??????.AR2.bz2')
    # radar_files = radar_files[200:300]
    radar_files = [radar_files[199]]
    for radar_file in radar_files:
        stats_radar_file(radar_file, 1)
