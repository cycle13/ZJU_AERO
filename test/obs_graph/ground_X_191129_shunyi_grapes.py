'''
Description: plot beijing 20191129 case
Author: Hejun Xie
Date: 2020-11-03 16:17:53
LastEditors: Hejun Xie
LastEditTime: 2021-10-04 17:52:39
'''

import pycwr
import pyart
import numpy as np
from scipy.ndimage import gaussian_filter
import datetime as dt
import glob

from pycwr.io import read_auto

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import matplotlib as mpl
mpl.use('Agg')
from pyart.graph import RadarMapDisplayBasemap
import matplotlib.pyplot as plt

# define some constants
fields  = ['reflectivity', 'differential_reflectivity', 'velocity', 'specific_differential_phase', 'cross_correlation_ratio']
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

def plot_radar_file(radar_file, sweep, directory):

    yyyymmdd = radar_file.split('.')[-4]
    hhmmss = radar_file.split('.')[-3]

    current_time = dt.datetime(int(yyyymmdd[0:4]),int(yyyymmdd[4:6]),int(yyyymmdd[6:8]),
        int(hhmmss[0:2]), int(hhmmss[2:4]), int(hhmmss[4:6]))

    print(current_time)
    
    PRD = read_auto(radar_file)
    r = PRD.ToPyartRadar()

    # print(r.fields)

    # print(PRD.scan_info)
    radial_start_index = r.sweep_start_ray_index['data'][sweep]
    radial_end_index = r.sweep_end_ray_index['data'][sweep] + 1
    sweep_smooth(r, radial_start_index, radial_end_index, 'cross_correlation_ratio', sigma=[2.0, 15.0])
    sweep_smooth(r, radial_start_index, radial_end_index, 'differential_reflectivity', sigma=[2.0, 15.0])
    sweep_smooth(r, radial_start_index, radial_end_index, 'reflectivity', sigma=[2.0, 15.0])
    sweep_smooth(r, radial_start_index, radial_end_index, 'specific_differential_phase', sigma=[2.0, 15.0])

    display = pyart.graph.RadarMapDisplayBasemap(r)
    plt.figure()

    for field in fields:
        display.plot_ppi_map(field, sweep, vmin=vrange[field][0], vmax=vrange[field][1],
                            min_lon=114, max_lon=120, min_lat=38.0, max_lat=42.5,
                            lon_lines=np.arange(38, 42, 1), projection='lcc',
                            lat_lines=np.arange(113, 119, 1), resolution='h',
                            lat_0=r.latitude['data'],
                            lon_0=r.longitude['data'],
                            shapefile='./ChinaProvince/ChinaProvince',
                            cmap=cmap[field],
                            title= 'Time: {}'.format(str(current_time)) + '\n' + \
                                'Elevation: {:.1f}'.format(r.elevation['data'][radial_start_index]) + DEG + '\n' + \
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

        plt.savefig('{}/{}/{}_ppi_obs.png'.format(directory, short_name[field],
            current_time.strftime('%Y-%m-%d_%H:%M_UTC')), dpi=300, bbox_inches='tight')

        plt.close()
    
    del r

if __name__ == "__main__":

    radar_files = glob.glob('../../pathos/RADAR/19.11.29顺义X波段/29/BJXSY.20191129.??????.AR2.bz2')
    radar_files = radar_files[200:300]
    # radar_files = [radar_files[199]]
    for radar_file in radar_files:
        plot_radar_file(radar_file, 1, './191129/')
        