'''
Description: simulate gpm dpr with grapes modelvar input
Author: Hejun Xie
Date: 2021-04-14 09:49:46
LastEditors: Hejun Xie
LastEditTime: 2021-04-16 20:54:01
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np
import pickle
# import pyart

# Local imports
import ZJU_AERO
import pyart
import datetime as dt
import os
import glob

LOAD_MODEL = True
LOAD_RADAR = True
DEG = r'$^\circ$'

cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17',
'KDP': 'pyart_EWilson17', 'PHIDP': 'pyart_Carbone42', 'RHOHV': 'pyart_GrMg16'}
clevels = {'ZH':np.arange(0,55,5), 'ZDR':np.linspace(0., 0.1, 50)}
units = {'ZH':'[dBZ]', 'ZDR':'[dBZ]'}
longname = {'ZH':'Horizontal Refelctivity', 'ZDR':'Differential Reflectivty'}

if __name__ == "__main__":

    FOLDER = '/mnt/f/modelvar/2020090500/'
    data_file_list = glob.glob(FOLDER+os.sep+'*.nc')
    load_datetime = dt.datetime(2020, 9, 5, 9)

    OBS_FILE = '/mnt/e/GPM/L2A/2A.GPM.Ku.V8-20180723.20200905-S083837-E101112.037051.V06A.HDF5'
    OBS_SLICE = (slice(2450, 2750), slice(0, 50))
    
    a = ZJU_AERO.RadarOperator(options_file='./option_files/gpm_dpr_sim.yml')
    a.load_model_file(data_file_list, load_datetime=load_datetime, load_from_file=LOAD_MODEL, load_file='mdl.nc')

    # print(a.dic_vars['T'])
    # a.close()
    # exit()

    if not LOAD_RADAR:
        # r = a.get_spaceborne_swath_test(OBS_FILE, slice=OBS_SLICE)
        r = a.get_spaceborne_swath(OBS_FILE, slice=OBS_SLICE)
        with open("./swath.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./swath.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()
    # exit()

    var = 'ZH'

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'
    from mpl_toolkits.basemap import Basemap

    fig, ax = plt.subplots(figsize=(10,8))

    slat = 18.0
    elat = 32.0
    slon = 125.0
    elon = 135.0

    zlevel = 5

    x = r.lons
    y = r.lats
    h = r.heights
    z = 10*np.log10(r.data[var])

    z[z<-10] = np.nan

    map = Basemap(projection='cyl', llcrnrlat=slat, urcrnrlat=elat, llcrnrlon=slon, urcrnrlon=elon, resolution='h', ax=ax)
    map.drawparallels(np.arange(slat, elat+1, 3.0), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
    map.drawmeridians(np.arange(slon, elon+1, 3.0), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])
    map.drawcoastlines()
    
    im = map.pcolormesh(x[...,zlevel], y[...,zlevel], z[...,zlevel], cmap='pyart_Carbone11', shading='auto')

    ax.set_title('Spaceborne Radar ' + var, fontsize=16)
    
    cb = fig.colorbar(im, ax=ax)
    cb.ax.set_ylabel('Reflectivity Factor [dBZ]', fontsize=14)

    plt.savefig('./gpm_dpr_swath', dpi=300, bbox_inches='tight')

    plt.close()
    del fig, ax
