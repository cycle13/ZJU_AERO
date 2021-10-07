'''
Description: plot W profile at cloudradar cite of Beijing
Author: Hejun Xie
Date: 2021-10-04 17:57:30
LastEditors: Hejun Xie
LastEditTime: 2021-10-07 18:28:21
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
plt.rcParams['font.family'] = 'serif'  

from ZJU_AERO.nwp.grapes import check_if_variables_in_file
from ZJU_AERO.nwp.grapes import get_grapes_variables as get_model_variables

FOLDER = '../../pathos/GRAPES/north_china_snowfall_20191112'
# var = 'W'
# var = 't'
var = 'RH'

if __name__ == "__main__":

    data_file_list = glob.glob(FOLDER+os.sep+'*.nc')
    load_datetime = dt.datetime(2019,11,29,18)

    ds = get_model_variables(data_file_list, ['W', 'T', 'QV', 'P'], load_datetime)
    ds['t'] = ds['T'] - 273.15
    ds['Es'] = 6.11 * 10 ** (7.5 * ds['t'] / (237.3 + ds['t']))
    ds['E'] = ds['QV'] * ds['P'] / (0.622 + ds['QV']) / 100
    ds['RH'] = ds['E'] / ds['Es'] * 100
 
    haidian = ds.interp(latitude=39.98333333, longitude=116.2833333)
    yanqing = ds.interp(latitude=40.45, longitude=115.9666667)
    pinggu = ds.interp(latitude=40.16666667, longitude=117.1166667)
    nanjiao = ds.interp(latitude=39.8, longitude=116.4666667)


    var = 'RH'
    fig = plt.figure(figsize=(10,8))
    ax = plt.gca()
    ax.plot(yanqing[var].data, yanqing.coords['z-levels']/1000., ls='--', color='k', label='54406 Yanqing')
    ax.plot(pinggu[var].data, pinggu.coords['z-levels']/1000., ls='-', color='k', label='54424 Pinggu')
    ax.plot(haidian[var].data, haidian.coords['z-levels']/1000., ls='-', color='r', label='54399 Haidian')
    ax.plot(nanjiao[var].data, nanjiao.coords['z-levels']/1000., ls='--', color='r', label='54511 Nanjiao')

    ax.legend(loc='best', frameon=False, fontsize=8)
    ax.set_xlim([40, 110])
    ax.set_ylim([0, 6])
    ax.set_xlabel('Relative Humidity [%]', fontsize=14)
    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_title('CloudRadar Site Vertical Profile Produced by GRAPES', fontsize=16)
    plt.savefig('cloudradarsite_profile_{}.png'.format(var), dpi=300)
    plt.close()

    var = 't'
    fig = plt.figure(figsize=(10,8))
    ax = plt.gca()
    ax.plot(yanqing[var].data, yanqing.coords['z-levels']/1000., ls='--', color='k', label='54406 Yanqing')
    ax.plot(pinggu[var].data, pinggu.coords['z-levels']/1000., ls='-', color='k', label='54424 Pinggu')
    ax.plot(haidian[var].data, haidian.coords['z-levels']/1000., ls='-', color='r', label='54399 Haidian')
    ax.plot(nanjiao[var].data, nanjiao.coords['z-levels']/1000., ls='--', color='r', label='54511 Nanjiao')

    ax.legend(loc='best', frameon=False, fontsize=8)
    ax.set_xlim([-20, 5])
    ax.set_ylim([0, 6])
    ax.set_xlabel(r'Temperature [$^\circ C$]', fontsize=14)
    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_title('CloudRadar Site Vertical Profile Produced by GRAPES', fontsize=16)
    plt.savefig('cloudradarsite_profile_{}.png'.format(var), dpi=300)
    plt.close()

    var = 'W'
    fig = plt.figure(figsize=(10,8))
    ax = plt.gca()
    ax.plot(yanqing[var].data, yanqing.coords['z-levels']/1000., ls='--', color='k', label='54406 Yanqing')
    ax.plot(pinggu[var].data, pinggu.coords['z-levels']/1000., ls='-', color='k', label='54424 Pinggu')
    ax.plot(haidian[var].data, haidian.coords['z-levels']/1000., ls='-', color='r', label='54399 Haidian')
    ax.plot(nanjiao[var].data, nanjiao.coords['z-levels']/1000., ls='--', color='r', label='54511 Nanjiao')

    ax.legend(loc='best', frameon=False, fontsize=8)
    ax.set_xlim([-1, 1])
    ax.set_ylim([0, 15])
    ax.set_xlabel('Vertical Wind [kg/kg]', fontsize=14)
    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_title('CloudRadar Site Vertical Profile Produced by GRAPES', fontsize=16)
    plt.savefig('cloudradarsite_profile_{}.png'.format(var), dpi=300)
    plt.close()
