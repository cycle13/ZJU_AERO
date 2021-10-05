'''
Description: plot grapes modelvar
Author: Hejun Xie
Date: 2021-10-04 19:03:31
LastEditors: Hejun Xie
LastEditTime: 2021-10-04 21:13:09
'''

# Global imports
import numpy as np
from numpy import ma
import glob
import os

from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Local import
from ..grapes import get_grapes_variables
from .graph_const import units, long_names, clevels, cmap


def _add_title(ax, title):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])

    ax.text(0.5, 0.50, title, fontsize=20, ha='center', va='center')


def graph_grapes_modelvar(model_file_list, plot_var, plot_time, plot_box, coords):
    '''
    Plot GRAPES model variables
    Params:
        model_file_list: A list of model file.
        plot_var: The variable to be ploted.
        plot_time: The datetime to load the variables.
        plot_box: [llc_lat, urc_lat, llc_lon, urc_lon]
        coords: [lat, lon]
    '''
    plt.rcParams['font.family'] = 'serif'

    llc_lat, urc_lat, llc_lon, urc_lon = plot_box
    lat, lon = coords

    plot_time_str = plot_time.strftime('%Y%m%d%H%Mz')
    
    # the plot_var as var
    ds = get_grapes_variables(model_file_list, ['W', 'U', 'V'], plot_time)
    W = ds.data_vars['W'] # [m/s]
    U = ds.data_vars['U'] # [m/s]
    V = ds.data_vars['V'] # [m/s]
    SPEED = np.sqrt(U*U + V*V) # [m/s]

    layer = 20
    u = U[layer, ...]
    v = V[layer, ...]
    w = W[layer, ...]
    speed = SPEED[layer, ...]

    print(W.coords)

    # var = W.data[20,...]
    var = W.coords['topograph']

    fig = plt.figure(figsize=(10.0,10.0))

    ax_title = fig.add_axes([0.1, 0.90, 0.82, 0.08])
    ax_cf = fig.add_axes([0.16, 0.14, 0.76, 0.72])
    ax_cb = fig.add_axes([0.1, 0.05, 0.85, 0.03])

    _add_title(ax_title, 'GRAPES {} {}'.format(long_names[plot_var], plot_time_str))

    origin = 'lower'
    extend = 'max'
    norm = None
    clevel = clevels[plot_var]
    ticks = clevel
    ticklabels = ['{:.1f}'.format(iclevel) for iclevel in clevel]

    map = Basemap(projection='lcc', llcrnrlat=llc_lat, urcrnrlat=urc_lat, llcrnrlon=llc_lon, urcrnrlon=urc_lon, 
    resolution='i', lat_1=30., lat_2=60., lat_0=34.83158, lon_0=115., ax=ax_cf)
    map.drawparallels(np.arange(llc_lat, urc_lat+1, 1), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
    map.drawmeridians(np.arange(llc_lon, urc_lon+1, 1), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])

    map.drawcoastlines()

    map.readshapefile('./ChinaProvince/ChinaProvince', 'shapefile', ax=ax_cf)

    x, y = map(lon, lat)

    CF = map.contourf(x, y, var, levels=clevel, cmap=cmap[plot_var], origin=origin, extend=extend, norm=norm)

    xx, yy = map.makegrid(v.shape[1], v.shape[0], returnxy=True)[2:4]
    STREAM = map.streamplot(xx, yy, u, v, color='orange', linewidth=1)

    CB = fig.colorbar(CF, cax=ax_cb, orientation='horizontal', ticks=ticks)
    CB.ax.set_xticklabels(ticklabels)
    for tick in CB.ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
        tick.label.set_rotation(45)

    CB.set_label('{} [{}]'.format(plot_var, units[plot_var]), fontsize=14)

    plt.savefig('{}_{}.png'.format(plot_var, plot_time_str), bbox_inches='tight', dpi=300)
    plt.close(fig)
