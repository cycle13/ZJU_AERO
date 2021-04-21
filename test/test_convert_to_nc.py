'''
Description: test transfer from GRAPES modelvar to nc
Author: Hejun Xie
Date: 2020-10-28 16:33:40
LastEditors: Hejun Xie
LastEditTime: 2021-04-21 14:18:02
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

from ZJU_AERO.nwp.grapes import convert_to_nc_general_ctl
from ZJU_AERO.nwp.grapes import convert_to_nc_specific_ctl

# [A]. SPECIFIC
# CTLFILE = '/mnt/f/modelvar/2020090500/model.ctl_202009050000900'
# NC = '/mnt/f/modelvar/2020090500/modelvar202009050000900.nc'
# convert_to_nc_specific_ctl(CTLFILE, NC, var=['zz', 'pi', 'th', 
#     'u', 'v', 'w', 'Qv', 'Qc', 'Qr', 'Qi', 'Qs', 'Qg'], 
#     geobox=[125, 135, 18, 32], v_grid_type='ml')

# [B]. GENERAL
CTLFILE = '../pathos/GRAPES/north_china_snowfall_20191112/grapes2019112900ctl'
convert_to_nc_general_ctl(CTLFILE, var='all', v_grid_type='ml')
