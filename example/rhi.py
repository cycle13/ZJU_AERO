'''
Description: Example to use RHI scan simulation
Author: Hejun Xie
Date: 2020-08-22 12:36:55
LastEditors: Hejun Xie
LastEditTime: 2021-10-21 19:27:44
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

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
    a = ZJU_AERO.RadarOperator(options_file='./option_files/example.yml')
    a.load_model_file([FILENAME], load_datetime=dt.datetime(2021, 8, 8, 0), load_from_file=LOAD_MODEL, load_file='mdl.nc')

    if not LOAD_RADAR:
        r = a.get_RHI(azimuths=40.0, elev_start=0.5, elev_stop=10.5, elev_step = 1.0)
        # r = a.get_PPI_test(elevations = 1.0)
        with open("./rhi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./rhi.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()

    from pyart.graph import RadarMapDisplayBasemap
    display = pyart.graph.RadarMapDisplayBasemap(r)
    plt.figure()

    for field in fields:
        fig = plt.figure(figsize=[10, 4])
        ax = fig.add_subplot(111)
        display.plot(field, 0, vmin=vrange[field][0], vmax=vrange[field][1],
                        cmap=cmap[field],
                        title=  'Time: {}'.format(a.get_pos_and_time()['time']) + '\n' + \
                                'Azimuth: {}'.format(r.azimuth['data'][0]) + DEG + '\n' + \
                                latex_name[field],
                        ax=ax)
        display.set_limits(ylim=[0, 15], xlim=[0, 150])
        plt.savefig('test_rhi_{}.png'.format(field), dpi=300, bbox_inches='tight')
        plt.close()
