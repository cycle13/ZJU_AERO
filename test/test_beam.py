'''
@Description: test and plot exhausted beam solver
@Author: Hejun Xie
@Date: 2020-08-02 12:46:24
LastEditors: Hejun Xie
LastEditTime: 2020-11-11 10:13:20
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# global import
import numpy as np
# import pyWRF as pw
import datetime as dt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# Local imports
from pyCRO.beam import fixed_radius_KE, ODEZeng2014, ODEZeng2014_exhaustive
from pyCRO.utilities import DATAdecorator
from pyCRO.nwp.wrf import get_wrf_variables
from pyCRO.config import cfg
from pyCRO.constants import global_constants as constants

@DATAdecorator('./', False, './she.pkl')
def get_she():
    FILENAME = '../pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'
    ds = get_wrf_variables([FILENAME], ['N'], dt.datetime(2013,10,6,10))

    s1, h1, e1 = fixed_radius_KE(range_vec, elevation_angles, coords_radar)
    s2, h2, e2 = ODEZeng2014(range_vec, elevation_angles, coords_radar, ds.data_vars['N'])
    s3_pre, h3_pre, e3_pre = ODEZeng2014_exhaustive(range_vec, elevation_angles, 
    azimuth_angle, coords_radar, ds.data_vars['N'])
    s3, h3, e3 = s3_pre[0], h3_pre[0], e3_pre[0]
    return s1,s2,s3, h1,h2,h3, e1,e2,e3 

if __name__ == "__main__":

    # get global constants
    cfg.init('./option_files/wsm6_test.yml')
    cfg.CONFIG = cfg.sanity_check(cfg.CONFIG)
    constants.update()
    
    # range_vec = np.arange(0, 100*1000, 500)
    # elevation_angles = [0.2]
    # azimuth_angle = 120
    # coords_radar = [27.9, 120.8, 200]

    # range_vec = np.arange(75, 150*1000, 150)
    range_vec = np.arange(1, 150*1000, 150)
    elevation_angles = [1.0]
    azimuth_angle = 0
    coords_radar = [27.9, 120.8, 200]
    
    s1,s2,s3, h1,h2,h3, e1,e2,e3 = get_she()

    print(h1)
    print(h2)
    print(h3)

    plt.rcParams['font.family'] = 'serif'
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.set_title('Beam Track e={:4.1f}'.format(elevation_angles[0])+r'$^\circ$')

    ax.plot(s1, h1, color='k', label='43ERM')
    ax.plot(s2, h2, color='r', label='(Zeng, 2014) Simplified solver')
    ax.plot(s3, h3, color='b', label='(Zeng, 2014) Exhaustive solver')
    
    ax.set_xlabel('Distance [m]', fontsize=12)
    ax.set_ylabel('Height [m]', fontsize=12)

    ax.legend(frameon=False)
    
    plt.savefig('./beam_track.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    