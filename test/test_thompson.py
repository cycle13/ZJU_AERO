'''
Description: test for Thompson microsphysics scheme
Author: Hejun Xie
Date: 2020-08-22 12:36:55
LastEditors: Hejun Xie
LastEditTime: 2021-10-18 15:24:48
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# Global imports
import numpy as np
import pickle
import os
import glob
import datetime as dt

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Local imports
import ZJU_AERO
import pyart

LOAD_MODEL = False
LOAD_RADAR = False
DEG = r'$^\circ$'


fields  = ['ZH', 'RVEL', 'ZDR']
cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17'}
vrange  = {'ZH':  (0, 50),
        'ZDR': (-5.0, 5.0),
        'RVEL': (-15, 15)}
cmap    = {'ZH':  'pyart_Carbone11',
        'ZDR': 'pyart_RefDiff',
        'RVEL': 'pyart_BuOr8'}
latex_name = {'ZH': r'$Z_{H}$',
        'ZDR': r'$Z_{DR}$',
        'RVEL': r'$V_r$'}

if __name__ == "__main__":
    
    FILENAME = '../pathos/WRF/thompson/wrfout_d02_2021-08-08_00_00_00'
    a = ZJU_AERO.RadarOperator(options_file='./option_files/thompson_test.yml')
    a.load_model_file([FILENAME], load_datetime=dt.datetime(2021, 8, 8, 0), load_from_file=LOAD_MODEL, load_file='mdl.nc')

    if not LOAD_RADAR:
        r = a.get_PPI(elevations = 1.494)
        # r = a.get_PPI_test(elevations = 1.0)
        with open("./ppi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./ppi.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()

    from pyart.graph import RadarMapDisplayBasemap
    display = pyart.graph.RadarMapDisplayBasemap(r)
    plt.figure()

    for field in fields:
        display.plot_ppi_map(field, 0, vmin=vrange[field][0], vmax=vrange[field][1],
                            min_lon=104, max_lon=109, min_lat=28, max_lat=32,
                            lon_lines=np.arange(104, 109, 1), projection='lcc',
                            lat_lines=np.arange(28, 32, 1), resolution='h',
                            lat_0=r.latitude['data'],
                            lon_0=r.longitude['data'],
                            shapefile='./obs_graph/ChinaProvince/ChinaProvince',
                            cmap=cmap[field],
                            title= 'Time: {}'.format(a.get_pos_and_time()['time']) + '\n' + \
                                'Elevation: {:.1f}'.format(r.elevation['data'][0]) + DEG + '\n' + \
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

        plt.savefig('test_wsm6_{}.png'.format(field), dpi=300, bbox_inches='tight')
        plt.close()
