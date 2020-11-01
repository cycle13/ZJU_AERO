'''
Description: test wrfout
Author: Hejun Xie
Date: 2020-08-28 17:03:50
LastEditors: Hejun Xie
LastEditTime: 2020-11-01 19:52:16
'''


import numpy as np
from numpy import ma
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import datetime as dt
from PIL import Image
import glob
import os


import pyWRF as pw

FILE_MDL = '../../../../cosmo_pol/pathos/WRF/wsm6_ERA5/wrfout_d01_2019-05-17_03_00_00'
model_ini_time = dt.datetime(2019, 5, 17, 3)
model_out_step = dt.timedelta(minutes=15)

def _add_title(ax, title):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])

    ax.text(0.5, 0.50, title, fontsize=20, ha='center', va='center')


if __name__ == '__main__':

    import netCDF4 as nc
    
    mdl = nc.Dataset(FILE_MDL, 'r')
    
    # (Time, bottom_top, south_north, west_east)
    p = mdl.variables['P'][0,...] + mdl.variables['PB'][0,...] # [Pa]
    lat = mdl.variables['XLAT'][0,...]
    lon = mdl.variables['XLONG'][0,...]

    times = np.arange(12, 60, 1)
    play_var    = 'Q_hydro'
    unit        = {'Q_hydro': 'kg*m^-2'}
    longname    = {'Q_hydro': 'Total Hydrometeor'}

    clevels_vars = {'Q_hydro': np.linspace(0., 30, 31)}

    # plot box
    slat, elat = 38, 42.
    slon, elon = 113, 119.

###############################################################################################################################

    # start plot model
    plt.rcParams['font.family'] = 'serif'

    for itime in times:
        
        file_h = pw.open_file(FILE_MDL)
        d = file_h.get_variable(['Zw','Q_v'], itime=itime, assign_heights=False)
        file_h.close()
        Z = d['Zw'].data
        Q = d['Q_v'].data
        dZ = Z[1:,...] - Z[:-1,...]

        var = np.einsum('jkl,jkl->kl', dZ, Q)

        plot_time = model_ini_time + model_out_step * itime
        plot_time_str = plot_time.strftime('%Y%m%d%H%Mz')

        fig = plt.figure(figsize=(10.0,7.0))

        ax_title = fig.add_axes([0.1, 0.90, 0.82, 0.08])
        ax_cf = fig.add_axes([0.16, 0.14, 0.76, 0.72])
        ax_cb = fig.add_axes([0.1, 0.05, 0.85, 0.03])

        _add_title(ax_title, 'WRF {} {}'.format(longname[play_var], plot_time_str))

        origin = 'lower'
        extend = 'max'
        cmap = 'jet'
        norm = None
        clevels = clevels_vars[play_var]
        ticks = clevels
        ticklabels = ['{:.1f}'.format(clevel) for clevel in clevels]

        map = Basemap(projection='cyl', llcrnrlat=slat, urcrnrlat=elat, llcrnrlon=slon, urcrnrlon=elon, resolution='l', ax=ax_cf)
        map.drawparallels(np.arange(slat, elat+1, 1), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
        map.drawmeridians(np.arange(slon, elon+1, 1), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])

        map.drawcoastlines()

        x, y = map(lon, lat)

        sliced_var = var
        sliced_var = ma.masked_where(sliced_var <= 0.1, sliced_var)

        CF = map.contourf(x, y, sliced_var, levels=clevels, cmap=cmap, origin=origin, extend=extend, norm=norm)
        CB = fig.colorbar(CF, cax=ax_cb, orientation='horizontal', ticks=ticks)
        CB.ax.set_xticklabels(ticklabels)
        for tick in CB.ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
            tick.label.set_rotation(45)

        CB.set_label('{} [{}]'.format(play_var, unit[play_var]), fontsize=14)

        plt.savefig('{}_{}.png'.format(play_var, plot_time_str), bbox_inches='tight', dpi=300)
        plt.close(fig)
    
    pic_files = glob.glob("{}_*.png".format(play_var))
    gif_file = "{}_ERA5_late.gif".format(play_var)
    imgs = [Image.open(ipic) for ipic in pic_files]
    imgs[0].save(gif_file, save_all=True, append_images=imgs, duration=2)

    for pic_file in pic_files:
        os.system("rm {}".format(pic_file))
    
    