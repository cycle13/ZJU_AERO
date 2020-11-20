'''
Description: test graph_wrf_modelvar.py
Author: Hejun Xie
Date: 2020-11-07 11:05:47
LastEditors: Hejun Xie
LastEditTime: 2020-11-20 17:06:27
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np
import netCDF4 as nc
import datetime as dt


# Local imports
from ZJU_AERO.nwp.graph import graph_wrf_modelvar_timeline

FILE_MDL = '../pathos/WRF/wsm6_test/ERA5/wrfout_d01_2019-05-17_03_00_00'
start_dt = dt.datetime(2019, 5, 17, 6)
end_dt = dt.datetime(2019, 5, 17, 18)
step_dt = dt.timedelta(minutes=15)
play_var = 'Q_hydro'

if __name__ == '__main__':

    mdl = nc.Dataset(FILE_MDL, 'r')
    
    # XLAT, XLONG (Time, south_north, west_east)
    lat = mdl.variables['XLAT'][0,...]
    lon = mdl.variables['XLONG'][0,...]

    times = np.arange(12, 60, 1)

    # plot box
    llc_lat, urc_lat = 38, 42.
    llc_lon, urc_lon = 113, 119.
    plot_box = tuple([llc_lat, urc_lat, llc_lon, urc_lon])
    coords = tuple([lat, lon])

###############################################################################################################################

    graph_wrf_modelvar_timeline([FILE_MDL], play_var, plot_box, coords,
    start_dt, step_dt, end_dt, 'test')

