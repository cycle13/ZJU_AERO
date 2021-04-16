'''
Description: convert grapes file to nc for gpm test case
Author: Hejun Xie
Date: 2021-04-14 10:52:53
LastEditors: Hejun Xie
LastEditTime: 2021-04-16 20:07:29
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

from ZJU_AERO.nwp.grapes import convert_to_nc_general_ctl
from ZJU_AERO.nwp.grapes import convert_to_nc_specific_ctl

CTLFILE = '/mnt/f/modelvar/2020090500/model.ctl_202009050000900'
NC = '/mnt/f/modelvar/2020090500/modelvar202009050000900.nc'
convert_to_nc_specific_ctl(CTLFILE, NC, var=['zz', 'pi', 'th', 
    'u', 'v', 'w', 'Qv', 'Qc', 'Qr', 'Qi', 'Qs', 'Qg'], 
    geobox=[125, 135, 18, 32], v_grid_type='ml')

# CTLFILE = '../pathos/GRAPES/north_china_snowfall_20191112/grapes2019112900ctl'
# convert_to_nc_general_ctl(CTLFILE, var='all', v_grid_type='ml')
