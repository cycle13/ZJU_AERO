'''
Description: test vertical profile radial
and Full doppler scheme
Author: Hejun Xie
Date: 2021-02-27 19:52:15
LastEditors: Hejun Xie
LastEditTime: 2021-02-28 16:02:18
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

LOAD_MODEL = True
LOAD_RADAR = True
DEG = r'$^\circ$'


fields  = ['ZH', 'RVEL', 'ZDR']
cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17'}
vrange  = {'ZH':  (0, 40),
        'ZDR': (0, 0.5),
        'RVEL': (-15, 15)}
cmap    = {'ZH':  'pyart_Carbone11',
        'ZDR': 'pyart_Carbone11',
        'RVEL': 'pyart_BuOr8'}
latex_name = {'ZH': r'$Z_{H}$',
        'ZDR': r'$Z_{DR}$',
        'RVEL': r'$V_r$'}

def get_timelines(time_range, time_step):

    load_datetimes = list()
    load_datetime = time_range[0]
    while load_datetime <= time_range[1]:
        load_datetimes.append(load_datetime)
        load_datetime += time_step
    
    return load_datetimes
    

if __name__ == "__main__":

    FOLDER = '../pathos/GRAPES/north_china_snowfall_20191112'
    data_file_list = glob.glob(FOLDER+os.sep+'*.nc')
    
    # get load_datetimes
    time_step = dt.timedelta(hours=1)
    time_range = (dt.datetime(2019,11,29,0), dt.datetime(2019,11,30,0))

    load_datetimes = get_timelines(time_range, time_step)

    rs = list()
    for load_datetime in load_datetimes:
        load_timestr = load_datetime.strftime("%Y%m%d%H")
        mpkl = './vprof/mdl{}.nc'.format(load_timestr)
        rpkl = "./vprof/vprof{}.pkl".format(load_timestr)
        
        if not LOAD_RADAR:
            a = ZJU_AERO.RadarOperator(options_file='./option_files/vprof.yml')
            a.load_model_file(data_file_list, load_datetime=load_datetime, load_from_file=LOAD_MODEL, load_file=mpkl)
            
            r = a.get_VPROF()
            with open(rpkl, "wb") as f:
                pickle.dump(r, f)
            a.close()
        else:
            with open(rpkl, "rb") as f:
                r = pickle.load(f)
        
        rs.append(r)
    
    units = {'dBZ':'[dBZ]', 'V':'[m/s]', 'W':'[m/s]'}
    longname = {'dBZ':'Refelctivity', 'V':'Radial Velocity', 'W':'Spectrum Width'}

    from mpl_toolkits.basemap import Basemap
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'

    
    heights_raw = rs[0].heights_profile
    idx = np.where(heights_raw<=15000)[0] # 15km
    heights = heights_raw[idx] / 1000. # km
    
    dBZ = np.zeros((len(load_datetimes), len(heights)), dtype='float32')
    V = np.zeros((len(load_datetimes), len(heights)), dtype='float32')

    for ir,r in enumerate(rs):
        dBZ[ir,:] = 10 * np.log10(r.values['ZH'][idx])
        V[ir,:] = r.values['RVEL'][idx]

    time = range(len(load_datetimes))

    fig = plt.figure(figsize=(10,8))
    xticklabels = [load_datetime.strftime('%Y-%m-%d %H') for load_datetime in load_datetimes]

    ax = plt.subplot(211)
    pm = ax.pcolormesh(time, heights, dBZ.T, cmap='pyart_Carbone11', 
    vmin=-40, vmax=30, shading='auto')
    ax.set_xticks(time)
    ax.set_xticklabels([])
    ax.set_ylabel('Height [km]', fontsize=14)
    cb = fig.colorbar(pm, ax=ax)
    cb.ax.set_ylabel('{} {}'.format(longname['dBZ'], units['dBZ']), fontsize=14)

    ax = plt.subplot(212)
    pm = ax.pcolormesh(time, heights, V.T, cmap='pyart_BuOr8', 
    vmin=-2, vmax=0., shading='auto')
    ax.set_xticks(time)
    ax.set_xticklabels(xticklabels, fontsize=9, rotation=30)
    ax.set_ylabel('Height [km]', fontsize=14)
    cb = fig.colorbar(pm, ax=ax)
    cb.ax.set_ylabel('{} {}'.format(longname['V'], units['V']), fontsize=14)

    fig.tight_layout()

    plt.savefig('test_vprof.png', dpi=300)
