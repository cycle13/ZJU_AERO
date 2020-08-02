'''
@Description: test and plot exhausted beam solver
@Author: Hejun Xie
@Date: 2020-08-02 12:46:24
@LastEditors: Hejun Xie
@LastEditTime: 2020-08-02 16:05:31
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# global import
import numpy as np
import pyWRF as pw
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# Local imports
from pyCRO.beam_propogation import fixed_radius_KE, ODEZeng2014, ODEZeng2014_exhausted
from pyCRO.utilities import DATAdecorator

@DATAdecorator('./', True, './she.pkl')
def get_she():
    FILENAME = '../../../cosmo_pol/pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'
    file_h = pw.open_file(FILENAME)
    d = file_h.get_variable(['N'], itime=10, assign_heights=True)

    s1, h1, e1 = fixed_radius_KE(range_vec, elevation_angles, coords_radar)
    s2, h2, e2 = ODEZeng2014(range_vec, elevation_angles, coords_radar, d['N'])
    s3_pre, h3_pre, e3_pre = ODEZeng2014_exhausted(range_vec, elevation_angles, 
    azimuth_angle, coords_radar, d['N'])
    s3, h3, e3 = s3_pre[0], h3_pre[0], e3_pre[0]
    return s1,s2,s3, h1,h2,h3, e1,e2,e3 

if __name__ == "__main__":
    
    range_vec = np.arange(0, 50*1000, 500)
    elevation_angles = [1.0]
    azimuth_angle = 120
    coords_radar = [27.9, 120.8, 200]
    
    s1,s2,s3, h1,h2,h3, e1,e2,e3 = get_she()

    print(h1)
    print(h2)
    print(h3)
