'''
Description: test vertical profile radial
and Full doppler scheme
Author: Hejun Xie
Date: 2021-02-27 19:52:15
LastEditors: Hejun Xie
LastEditTime: 2021-03-08 09:50:38
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
LOAD_RADAR = False
DEG = r'$^\circ$'

# np.set_printoptions(threshold=np.inf)

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
    a = ZJU_AERO.RadarOperator(options_file='./option_files/vprof.yml')

    for load_datetime in load_datetimes:

        # load_datetime = dt.datetime(2019,11,29,14)
        print(load_datetime)

        load_timestr = load_datetime.strftime("%Y%m%d%H")
        mpkl = './vprof/mdl{}.nc'.format(load_timestr)
        rpkl = "./vprof/vprof{}.pkl".format(load_timestr)
        
        a.load_model_file(data_file_list, load_datetime=load_datetime, load_from_file=LOAD_MODEL, load_file=mpkl)
        if not LOAD_RADAR:
            r = a.get_VPROF()
            with open(rpkl, "wb") as f:
                pickle.dump(r, f)
        else:
            with open(rpkl, "rb") as f:
                r = pickle.load(f)
        rs.append(r)

        # a.close()
        # exit()
    
    a.close()
    # exit()

    units = {'dBZ':'[dBZ]', 'V':'[m/s]', 'W':'[m/s]', 'DSP':'[dBZ]'}
    longname = {'dBZ':'Refelctivity', 'V':'Radial Velocity', 'W':'Spectrum Width', 'DSP':'Doppler Spectrum'}

    from mpl_toolkits.basemap import Basemap
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'

    
    heights_raw = rs[0].heights_profile
    idx = np.where(heights_raw<=15000)[0] # 15km
    heights = heights_raw[idx] / 1000. # km
    from ZJU_AERO.const import global_constants as constants
    varray = constants.VARRAY

    
    dBZ = np.zeros((len(load_datetimes), len(heights)), dtype='float32')
    V = np.zeros((len(load_datetimes), len(heights)), dtype='float32')
    W = np.zeros((len(load_datetimes), len(heights)), dtype='float32')
    DSP = np.zeros((len(load_datetimes), len(heights), len(varray)), dtype='float32')

    for ir,r in enumerate(rs):
        dBZ[ir,:] = 10 * np.log10(r.values['ZH'][idx])
        V[ir,:] = r.values['RVEL'][idx]
        W[ir,:] = r.values['SW'][idx]
        DSP[ir,:] = 10 * np.log10(r.values['DSPECTRUM'][idx,:])

    '''
    1. Plot Z, V, W
    '''
    time = range(len(load_datetimes))

    fig = plt.figure(figsize=(10,12))
    xticklabels = [load_datetime.strftime('%Y-%m-%d %H') for load_datetime in load_datetimes]

    ax = plt.subplot(311)
    pm = ax.pcolormesh(time, heights, dBZ.T, cmap='pyart_Carbone11', 
    vmin=-40, vmax=30, shading='auto')
    ax.set_xticks(time[::2])
    ax.set_xticklabels([])
    ax.set_ylabel('Height [km]', fontsize=14)
    cb = fig.colorbar(pm, ax=ax)
    cb.ax.set_ylabel('{} {}'.format(longname['dBZ'], units['dBZ']), fontsize=14)

    ax = plt.subplot(312)
    pm = ax.pcolormesh(time, heights, V.T, cmap='pyart_BuOr8', 
    vmin=-2, vmax=0., shading='auto')
    ax.set_xticks(time[::2])
    ax.set_xticklabels([])
    ax.set_ylabel('Height [km]', fontsize=14)
    cb = fig.colorbar(pm, ax=ax)
    cb.ax.set_ylabel('{} {}'.format(longname['V'], units['V']), fontsize=14)

    ax = plt.subplot(313)
    pm = ax.pcolormesh(time, heights, W.T, cmap='pyart_Carbone17', 
    vmin=0, vmax=2, shading='auto')
    ax.set_xticks(time[::2])
    ax.set_xticklabels(xticklabels[::2], fontsize=9, rotation=30)
    ax.set_ylabel('Height [km]', fontsize=14)
    ax.set_xlabel('Time')
    cb = fig.colorbar(pm, ax=ax)
    cb.ax.set_ylabel('{} {}'.format(longname['W'], units['W']), fontsize=14)

    fig.tight_layout()

    plt.savefig('test_vprof.png', dpi=300)
    plt.close()

    '''
    2. Plot DSP
    '''

    time_indices = [12, 14, 16, 18]

    for time_idx in time_indices: 
        
        time_file = load_datetimes[time_idx].strftime("%Y%m%d%H")
        time_title = load_datetimes[time_idx].strftime("%Y-%m-%d %H:00:00 UTC")

        fig = plt.figure(figsize=(5,8))
        ax = plt.subplot(111)

        idx_v = (varray<=0.) & (varray>=-5) # -10m/s < v < 0m/s
        v_valid = varray[idx_v]

        idx_h = (heights<=6) # h < 6km
        h_valid = heights[idx_h]

        plot_DSP =  DSP[time_idx,idx_h,:][:,idx_v]
        plot_DSP[plot_DSP<=-60] = np.nan

        plot_V = V[time_idx, idx_h]
        plot_W = W[time_idx, idx_h]

        pm = ax.pcolormesh(v_valid, h_valid, plot_DSP, # cmap='pyart_Carbone11', 
        vmin=-60, vmax=20, shading='auto')
        ax.set_xticks([-5, -4, -3, -2, -1, 0])
        ax.set_xticklabels(['-5','-4','-3','-2','-1','0'], fontsize=9)
        ax.set_ylabel('Height [km]', fontsize=14)
        ax.set_xlabel(r'Radial Velocity $V_{rad}$ [m/s]', fontsize=14)

        ax.plot(plot_V, h_valid, color='k')
        ax.plot(plot_V+plot_W, h_valid, color='k', ls='--')
        ax.plot(plot_V-plot_W, h_valid, color='k', ls='--')

        cb = fig.colorbar(pm, ax=ax)
        cb.ax.set_ylabel('{} {}'.format(longname['DSP'], units['DSP']), fontsize=14)

        ax.set_title(time_title, fontsize=16)

        plt.savefig('doppler_spectrum_{}.png'.format(time_file), dpi=300)
        plt.close()
