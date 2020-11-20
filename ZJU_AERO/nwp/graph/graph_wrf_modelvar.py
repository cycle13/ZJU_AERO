'''
Description: plot wrf modelvar
Author: Hejun Xie
Date: 2020-08-28 17:03:50
LastEditors: Hejun Xie
LastEditTime: 2020-11-07 12:58:31
'''

# Global imports
import numpy as np
from numpy import ma
from PIL import Image
import glob
import os

from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Local import
from ..wrf import get_wrf_variables
from .graph_const import units, long_names, clevels

def _add_title(ax, title):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])

    ax.text(0.5, 0.50, title, fontsize=20, ha='center', va='center')

def _gen_timeline(start_dt, end_dt, step_dt):
    timeline = []
    cur_dt = start_dt
    while cur_dt <= end_dt:
        timeline.append(cur_dt)
        cur_dt += step_dt

    return timeline

def graph_wrf_modelvar_timeline(model_file_list, play_var, plot_box, coords, 
    start_dt, step_dt, end_dt, gif_tag):
    timeline = _gen_timeline(start_dt, end_dt, step_dt)
    # for time in timeline:
    #     graph_wrf_modelvar(model_file_list, play_var, time, plot_box, coords)
    
    pic_files = glob.glob("{}_*.png".format(play_var))
    gif_file = "{}_{}.gif".format(play_var, gif_tag)
    imgs = [Image.open(ipic) for ipic in pic_files]
    imgs[0].save(gif_file, save_all=True, append_images=imgs, duration=2)

    for pic_file in pic_files:
        os.system("rm {}".format(pic_file))


def graph_wrf_modelvar(model_file_list, plot_var, plot_time, plot_box, coords):
    '''
    Plot WRF model variables
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
    ds = get_wrf_variables(model_file_list, ['Z', 'Qhydro', 'RHO'], plot_time)
    Z = ds.data_vars['Z'].data # [m]
    Q = ds.data_vars['Qhydro'] * ds.data_vars['RHO']  # [kg * m-3]
    dZ = np.pad(np.diff(Z, axis=0), ((0, 1), (0, 0), (0, 0)), 'edge')
    var = np.einsum('jkl,jkl->kl', dZ, Q)

    fig = plt.figure(figsize=(10.0,7.0))

    ax_title = fig.add_axes([0.1, 0.90, 0.82, 0.08])
    ax_cf = fig.add_axes([0.16, 0.14, 0.76, 0.72])
    ax_cb = fig.add_axes([0.1, 0.05, 0.85, 0.03])

    _add_title(ax_title, 'WRF {} {}'.format(long_names[plot_var], plot_time_str))

    origin = 'lower'
    extend = 'max'
    cmap = 'jet'
    norm = None
    clevel = clevels[plot_var]
    ticks = clevel
    ticklabels = ['{:.1f}'.format(iclevel) for iclevel in clevel]

    map = Basemap(projection='cyl', llcrnrlat=llc_lat, urcrnrlat=urc_lat, llcrnrlon=llc_lon, urcrnrlon=urc_lon, 
    resolution='l', ax=ax_cf)
    map.drawparallels(np.arange(llc_lat, urc_lat+1, 1), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
    map.drawmeridians(np.arange(llc_lon, urc_lon+1, 1), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])

    map.drawcoastlines()

    x, y = map(lon, lat)

    var = ma.masked_where(var<=0.1, var)

    CF = map.contourf(x, y, var, levels=clevel, cmap=cmap, origin=origin, extend=extend, norm=norm)
    CB = fig.colorbar(CF, cax=ax_cb, orientation='horizontal', ticks=ticks)
    CB.ax.set_xticklabels(ticklabels)
    for tick in CB.ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
        tick.label.set_rotation(45)

    CB.set_label('{} [{}]'.format(plot_var, units[plot_var]), fontsize=14)

    plt.savefig('{}_{}.png'.format(plot_var, plot_time_str), bbox_inches='tight', dpi=300)
    plt.close(fig)
