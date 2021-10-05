'''
Description: test ploting grapes modelvar
Author: Hejun Xie
Date: 2021-10-04 19:15:15
LastEditors: Hejun Xie
LastEditTime: 2021-10-04 21:12:25
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np
import netCDF4 as nc
import datetime as dt
import glob
import os

# Local imports
from ZJU_AERO.nwp.graph import graph_grapes_modelvar

FILE_MDL = '../../pathos/GRAPES/north_china_snowfall_20191112/grapes2019112900_00.nc'
FOLDER_MDL = '../../pathos/GRAPES/north_china_snowfall_20191112/'
plot_dt = dt.datetime(2019, 11, 29, 14)

if __name__ == '__main__':

    mdl = nc.Dataset(FILE_MDL, 'r')
    
    # XLAT, XLONG (Time, south_north, west_east)
    lat = mdl.variables['latitudes']
    lon = mdl.variables['longitudes']

    data_file_list = glob.glob(FOLDER_MDL+os.sep+'*.nc')

    TLAT, TLON = np.meshgrid(lat, lon)

    # plot box
    llc_lat, urc_lat = 38, 42.
    llc_lon, urc_lon = 113, 119.
    plot_box = tuple([llc_lat, urc_lat, llc_lon, urc_lon])
    coords = tuple([TLAT.T, TLON.T])

    graph_grapes_modelvar(data_file_list, 'topo', plot_dt, plot_box, coords)
