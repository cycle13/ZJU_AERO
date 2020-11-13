'''
Description: plot beijing 20191129 case
Author: Hejun Xie
Date: 2020-11-03 16:17:53
LastEditors: Hejun Xie
LastEditTime: 2020-11-13 11:26:35
'''

import pycwr
import pyart
import numpy as np
import datetime as dt

from pycwr.io import read_auto

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# UTC 2019-11-29 12UTC
current_time = dt.datetime(2019,11,29,15)
file = r"../pathos/RADAR/19.11.29顺义X波段/29/BJXSY.20191129.150003.AR2.bz2"

DEG = r'$^\circ$'

PRD = read_auto(file)
r = PRD.ToPyartRadar()

print(PRD.scan_info)
# print(PRD.fields[0])
# print(r.fields.keys())
# exit()

import matplotlib as mpl
mpl.use('Agg')
from pyart.graph import RadarMapDisplayBasemap
display = pyart.graph.RadarMapDisplayBasemap(r)
import matplotlib.pyplot as plt
plt.figure()

fields  = ['reflectivity', 'differential_reflectivity', 'velocity']
vrange  = {'reflectivity':  (0, 40),
        'differential_reflectivity': (0, 2),
        'velocity': (-15, 15)}
cmap    = {'reflectivity':  'pyart_Carbone11',
        'differential_reflectivity': 'pyart_Carbone11',
        'velocity': 'pyart_BuOr8'}
latex_name = {'reflectivity': r'$Z_{H}$',
        'differential_reflectivity': r'$Z_{DR}$',
        'velocity': r'$V_r$'}
short_name = {'reflectivity':  'ZH',
        'differential_reflectivity': 'ZDR',
        'velocity': 'RVEL'}    

sweep = 0

for field in fields:
    display.plot_ppi_map(field, sweep, vmin=vrange[field][0], vmax=vrange[field][1],
                        min_lon=114, max_lon=120, min_lat=38.0, max_lat=42.5,
                        lon_lines=np.arange(38, 42, 1), projection='lcc',
                        lat_lines=np.arange(113, 119, 1), resolution='h',
                        lat_0=r.latitude['data'],
                        lon_0=r.longitude['data'],
                        cmap=cmap[field],
                        title= 'Time: {}'.format(str(current_time)) + '\n' + \
                            'Elevation: {:.1f}'.format(r.elevation['data'][sweep]) + DEG + '\n' + \
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

    plt.savefig('{}_ppi_obs_{}.png'.format(current_time.strftime('%Y-%m-%d-%HUTC'), short_name[field]), dpi=300, bbox_inches='tight')

    plt.close()
