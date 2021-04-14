'''
Description: plot gpm-dpr observation
Author: Hejun Xie
Date: 2021-04-13 09:55:12
LastEditors: Hejun Xie
LastEditTime: 2021-04-13 11:05:47
'''

gpm_folder = '/mnt/e/GPM/L2A/'
# gpm_file = '/mnt/e/GPM/L2A/2A.GPM.Ku.V8-20180723.20200901-S085847-E103120.036989.V06A.HDF5'

import os
import h5py
import numpy as np
import pyart
# np.set_printoptions(threshold=np.inf)

if __name__ == "__main__":

    gpm_files = os.listdir(gpm_folder)

    for gpm_file in gpm_files:

        if gpm_file.split('.')[-1] != 'HDF5':
            continue

        print(gpm_file)

        gpm_f = h5py.File(os.path.join(gpm_folder, gpm_file), 'r')
        group = 'NS'
        
        binNoClutter = gpm_f[group]['PRE']['binClutterFreeBottom'][:]
        # print(binNoClutter)

        ZH_gpm = gpm_f[group]['PRE']['zFactorMeasured'][:]
        ZH_gpm[ZH_gpm<-1000] = float('nan') # Put Nan where data is missing

        [N,M,K] = ZH_gpm.shape
        j, k = np.meshgrid(np.arange(N), np.arange(M))

        ZH_gpm_grd = ZH_gpm[j.T, k.T, binNoClutter]
        ZH_gpm_grd[ZH_gpm_grd<-1000] = float('nan')

        mask = gpm_f[group]['PRE']['flagPrecip'][:]
        ZH_gpm_grd[mask] = np.nan

        # print(ZH_gpm.shape)
        # print(ZH_gpm_grd)


        GPMslice = (slice(0, N+1), slice(0, M+1))
        lat_2D = gpm_f[group]['Latitude'][GPMslice]
        lon_2D = gpm_f[group]['Longitude'][GPMslice]
        center_lat_sc = gpm_f[group]['navigation']['scLat'][GPMslice[0]]
        center_lon_sc = gpm_f[group]['navigation']['scLon'][GPMslice[0]]
        altitudes_sc = gpm_f[group]['navigation']['dprAlt'][GPMslice[0]]

        # print(lat_2D)
        # print(lon_2D)
        # print(center_lat_sc)
        # print(center_lon_sc)

        var = 'ZH'

        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams['font.family'] = 'serif'
        from mpl_toolkits.basemap import Basemap

        fig, ax = plt.subplots(figsize=(10,8))

        slat = 15.0
        elat = 35.0
        slon = 125.0
        elon = 150.0

        # slat = -60.0
        # elat = 60.0
        # slon = -180.0
        # elon = 180.0

        map = Basemap(projection='cyl', llcrnrlat=slat, urcrnrlat=elat, llcrnrlon=slon, urcrnrlon=elon, resolution='h', ax=ax)
        map.drawparallels(np.arange(slat, elat+1, 5.0), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
        map.drawmeridians(np.arange(slon, elon+1, 5.0), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])
        map.drawcoastlines()
        
        im = map.pcolormesh(lon_2D, lat_2D, ZH_gpm_grd, cmap='pyart_Carbone11', shading='auto')

        ax.set_title('Spaceborne Radar ' + var, fontsize=16)
        
        cb = fig.colorbar(im, ax=ax)
        cb.ax.set_ylabel('Reflectivity Factor [dBZ]', fontsize=14)

        plt.savefig('{}/swath_{}.png'.format(gpm_folder, gpm_file), dpi=300, bbox_inches='tight')

        plt.close()
        del fig, ax
    