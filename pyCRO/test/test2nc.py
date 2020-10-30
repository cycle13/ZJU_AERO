'''
Description: test transfer from GRAPES modelvar to nc
Author: Hejun Xie
Date: 2020-10-28 16:33:40
LastEditors: Hejun Xie
LastEditTime: 2020-10-30 18:55:03
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

from pyCRO.nwp.grapes import transf2nc_specific_ctl
from pyCRO.nwp.grapes import transf2nc_general_ctl

# CTLFILE = '/mnt/e/GRAPES/north_china_snowfall_20191112/grapes2019112900_10.ctl'
# NC = '/mnt/e/GRAPES/north_china_snowfall_20191112/grapes2019112900_10.nc'
# transf2nc_specific_ctl(CTLFILE, NC, var='all', v_grid_type='ml')

CTLFILE = '/mnt/e/GRAPES/north_china_snowfall_20191112/grapes2019112900ctl'
transf2nc_general_ctl(CTLFILE, var='all', v_grid_type='ml')
