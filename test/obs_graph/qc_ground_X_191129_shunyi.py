'''
Description: Do statitstic works on shunyi X band Radar
Author: Hejun Xie
Date: 2021-10-06 15:08:44
LastEditors: Hejun Xie
LastEditTime: 2021-10-07 18:43:57
'''

import pycwr
import pyart
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


# define some constants
fields  = ['reflectivity', 'differential_reflectivity', 'cross_correlation_ratio']
vrange  = {'reflectivity':  (0, 40),
        'differential_reflectivity': (-1.5, 1.5),
        'velocity': (-15, 15),
        'specific_differential_phase': (-1, 6),
        'cross_correlation_ratio': (0.85, 1.0)}
cmap    = {'reflectivity':  'pyart_Carbone11',
        'differential_reflectivity': 'pyart_RefDiff',
        'velocity': 'pyart_BuOr8',
        'specific_differential_phase': 'pyart_RefDiff',
        'cross_correlation_ratio': 'pyart_RefDiff'}
latex_name = {'reflectivity': r'$Z_{H}$',
        'differential_reflectivity': r'$Z_{DR}$',
        'velocity': r'$V_r$',
        'specific_differential_phase': r'$K_{DP}$',
        'cross_correlation_ratio': r'$\rho_{hv}$'}
short_name = {'reflectivity':  'ZH',
        'differential_reflectivity': 'ZDR',
        'velocity': 'RVEL',
        'specific_differential_phase': 'KDP',
        'cross_correlation_ratio': 'RHOHV'}
DEG = r'$^\circ$'

def plot_radar_file(pyart_radar, sweep, casename):
    radial_start_index = pyart_radar.sweep_start_ray_index['data'][sweep]
    radial_end_index = pyart_radar.sweep_end_ray_index['data'][sweep] + 1

    display = pyart.graph.RadarMapDisplayBasemap(pyart_radar)

    for field in fields:
        plt.figure()
        display.plot_ppi_map(field, sweep, vmin=vrange[field][0], vmax=vrange[field][1],
                            min_lon=114, max_lon=120, min_lat=38.0, max_lat=42.5,
                            lon_lines=np.arange(38, 42, 1), projection='lcc',
                            lat_lines=np.arange(113, 119, 1), resolution='h',
                            lat_0=pyart_radar.latitude['data'],
                            lon_0=pyart_radar.longitude['data'],
                            shapefile='./ChinaProvince/ChinaProvince',
                            cmap=cmap[field],
                            title= 'Elevation: {:.1f}'.format(pyart_radar.elevation['data'][radial_start_index]) + DEG + '\n' + \
                                latex_name[field])
        # plot range rings at 50, 100, 200km
        display.plot_range_ring(50., line_style='k-', lw=1.0)
        display.plot_range_ring(100., line_style='k--', lw=1.0)
        display.plot_range_ring(200., line_style='k-', lw=1.0)

        # plots cross hairs
        display.plot_line_xy(np.array([-300000.0, 300000.0]), np.array([0.0, 0.0]),
                            line_style='k-', lw=1.2)
        display.plot_line_xy(np.array([0.0, 0.0]), np.array([-300000.0, 300000.0]),
                            line_style='k-', lw=1.2)

        plt.rc('font', size=6)
        kwargs = {'marker':'+', 'color':'r', 'ms':3}
        display.plot_point(116.2833333, 39.98333333, label_text='54399 Haidian', **kwargs)
        display.plot_point(115.9666667, 40.45, label_text='54406 Yanqing', **kwargs)
        display.plot_point(117.1166667, 40.16666667, label_text='54424 Pinggu', **kwargs)
        display.plot_point(116.4666667, 39.8, label_text='54511 Nanjiao', **kwargs)
        plt.rc('font', size=8)

        plt.savefig('QCplot/{}/{}_ppi_obs.png'.format(casename, short_name[field]), dpi=300, bbox_inches='tight')

        plt.close()


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

def preprocess_ctl1(radar_file, sweep):
    PRD = read_auto(radar_file)
    r = PRD.ToPyartRadar()
    radial_start_index = r.sweep_start_ray_index['data'][sweep]
    radial_end_index = r.sweep_end_ray_index['data'][sweep] + 1
    zdr     = get_sweep_field(r, radial_start_index, radial_end_index, 'differential_reflectivity').flatten()

    plot_radar_file(r, sweep, 'ctl1')
    del r
    return zdr

def preprocess_ctl3(radar_file, sweep):
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

    plot_radar_file(r, sweep, 'ctl3')

    zdr     = get_sweep_field(r, radial_start_index, radial_end_index, 'differential_reflectivity').flatten()
    del r
    return zdr

def preprocess_ctl2(radar_file, sweep):
    PRD = read_auto(radar_file)
    r = PRD.ToPyartRadar()
    radial_start_index = r.sweep_start_ray_index['data'][sweep]
    radial_end_index = r.sweep_end_ray_index['data'][sweep] + 1
    sweep_smooth(r, radial_start_index, radial_end_index, 'differential_reflectivity', sigma=[2.0, 15.0])
    zdr     = get_sweep_field(r, radial_start_index, radial_end_index, 'differential_reflectivity').flatten()

    plot_radar_file(r, sweep, 'ctl2')
    del r
    return zdr

def preprocess_1(radar_file, sweep):
    PRD = read_auto(radar_file)
    r = PRD.ToPyartRadar()
    radial_start_index = r.sweep_start_ray_index['data'][sweep]
    radial_end_index = r.sweep_end_ray_index['data'][sweep] + 1

    sweep_smooth(r, radial_start_index, radial_end_index, 'reflectivity', sigma=[2.0, 15.0])
    sweep_smooth(r, radial_start_index, radial_end_index, 'cross_correlation_ratio', sigma=[2.0, 15.0])
    zh      = get_sweep_field(r, radial_start_index, radial_end_index, 'reflectivity').flatten()
    rhohv   = get_sweep_field(r, radial_start_index, radial_end_index, 'cross_correlation_ratio').flatten()

    sweep_smooth(r, radial_start_index, radial_end_index, 'differential_reflectivity', sigma=[2.0, 15.0])
    sweep_filter(r, 'differential_reflectivity', 'reflectivity', 20.0, radial_start_index, radial_end_index)
    sweep_filter(r, 'differential_reflectivity', 'cross_correlation_ratio', 0.95, radial_start_index, radial_end_index)

    zdr     = get_sweep_field(r, radial_start_index, radial_end_index, 'differential_reflectivity').flatten()

    plot_radar_file(r, sweep, 'prep1')
    del r
    return zdr

def preprocess_2(radar_file, sweep):
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

    plot_radar_file(r, sweep, 'prep2')
    del r
    return zdr
    
def stats_radar_file(radar_file, sweep):

    yyyymmdd = radar_file.split('.')[-4]
    hhmmss = radar_file.split('.')[-3]

    current_time = dt.datetime(int(yyyymmdd[0:4]),int(yyyymmdd[4:6]),int(yyyymmdd[6:8]),
        int(hhmmss[0:2]), int(hhmmss[2:4]), int(hhmmss[4:6]))

    print(current_time)

    zdr_c1 = preprocess_ctl1(radar_file, sweep)
    zdr_c2 = preprocess_ctl2(radar_file, sweep)
    zdr_c3 = preprocess_ctl3(radar_file, sweep)
    zdr_1 = preprocess_1(radar_file, sweep)
    zdr_2 = preprocess_2(radar_file, sweep)
    ZDR_c1 = zdr_c1[~np.isnan(zdr_c1)]
    ZDR_c2 = zdr_c2[~np.isnan(zdr_c2)]
    ZDR_c3 = zdr_c3[~np.isnan(zdr_c3)]
    ZDR_1 = zdr_1[~np.isnan(zdr_1)]
    ZDR_2 = zdr_2[~np.isnan(zdr_2)]

    zdr_bins = np.arange(-2.1, 2.3, 0.2)
    print(zdr_bins)
    
    fig = plt.figure(figsize=(6,8))
    ax = plt.gca()
    
    ax.hist(ZDR_c1, bins=zdr_bins, ec="k", fc="None", alpha=1.0, lw=1, histtype='step', ls='-')
    ax.hist(ZDR_c2, bins=zdr_bins, ec="k", fc="None", alpha=1.0, lw=1, histtype='step', ls='--')
    ax.hist(ZDR_c3, bins=zdr_bins, ec="k", fc="blue", alpha=0.5, lw=0.5)
    ax.hist(ZDR_2, bins=zdr_bins, ec="k", fc="green", alpha=0.5, lw=0.5)

    ax.set_xlabel('Differential Reflectivity ZDR [dBZ]', fontsize=14)
    ax.set_ylabel('Radar Gate Number', fontsize=14)
    ax.set_title('Shunyi X-Band Radar ZDR Histogram Distribution', fontsize=16)
    plt.savefig('ZDR_histogram.png', dpi=300)
    plt.close()

if __name__ == "__main__":

    radar_files = glob.glob('../../pathos/RADAR/19.11.29顺义X波段/29/BJXSY.20191129.??????.AR2.bz2')
    # radar_files = radar_files[200:300]
    radar_files = [radar_files[199]]
    for radar_file in radar_files:
        stats_radar_file(radar_file, 1)
