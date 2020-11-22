'''
Description: test GRAPES interface for radar operator
Author: Hejun Xie
Date: 2020-11-02 16:17:47
LastEditors: Hejun Xie
LastEditTime: 2020-11-22 15:02:49
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

# Local imports
import ZJU_AERO
import pyart

LOAD_MODEL = False
LOAD_RADAR = False
DEG = r'$^\circ$'


fields  = ['ZH', 'RVEL', 'ZDR']
cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17'}
vrange  = {'ZH':  (0, 40),
        'ZDR': (0, 2),
        'RVEL': (-15, 15)}
cmap    = {'ZH':  'pyart_Carbone11',
        'ZDR': 'pyart_Carbone11',
        'RVEL': 'pyart_BuOr8'}
latex_name = {'ZH': r'$Z_{H}$',
        'ZDR': r'$Z_{DR}$',
        'RVEL': r'$V_r$'}

if __name__ == "__main__":

    FOLDER = '../pathos/GRAPES/north_china_snowfall_20191112'
    data_file_list = glob.glob(FOLDER+os.sep+'*.nc')
    load_datetime = dt.datetime(2019,11,29,18)
    
    a = ZJU_AERO.RadarOperator(options_file='./option_files/grapes_interface.yml')
    a.load_model_file(data_file_list, load_datetime=load_datetime, load_from_file=LOAD_MODEL, load_file='mdl.nc')

    if not LOAD_RADAR:
        r = a.get_PPI(elevations = 0.5)
        # r = a.get_PPI_test(elevations = 1)
        with open("./ppi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./ppi.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()

    # exit()
    # test PyartRadop
    import matplotlib as mpl
    mpl.use('Agg')
    from pyart.graph import RadarMapDisplayBasemap
    display = pyart.graph.RadarMapDisplayBasemap(r)
    import matplotlib.pyplot as plt
    plt.figure()

    for field in fields:
        display.plot_ppi_map(field, 0, vmin=vrange[field][0], vmax=vrange[field][1],
                            min_lon=114, max_lon=120, min_lat=38.0, max_lat=42.5,
                            lon_lines=np.arange(38, 42, 1), projection='lcc',
                            lat_lines=np.arange(113, 119, 1), resolution='h',
                            lat_0=r.latitude['data'],
                            lon_0=r.longitude['data'],
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

        plt.savefig('{}_ppi_grapes_iitm_{}.png'.format(load_datetime.strftime('%Y-%m-%d-%HUTC'), field), dpi=300, bbox_inches='tight')

        plt.close()
