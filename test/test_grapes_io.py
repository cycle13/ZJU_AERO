'''
Description: test nwp.grapes_io.py
Author: Hejun Xie
Date: 2020-11-01 11:12:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-11 10:15:02
'''


# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global import
import glob
import os
import datetime

FOLDER = '../pathos/GRAPES/north_china_snowfall_20191112'
data_file_list = glob.glob(FOLDER+os.sep+'*.nc')

from pyCRO.nwp.grapes import get_grapes_variables
from pyCRO.nwp.grapes import WGS_to_GRAPES

if __name__ == "__main__":
    
    # 1. U, V, W, T, P, Pw, RHO
    # ds = get_grapes_variables(data_file_list, 
    # ['U', 'V', 'W', 'T', 'P', 'Pw', 'RHO'], datetime.datetime(2019,11,29,3))

    # 2. QV_v, QR_v, QS_v, QG_v, QC_v, QI_v
    # ds = get_grapes_variables(data_file_list, 
    # ['QV_v', 'QR_v', 'QS_v', 'QG_v', 'QC_v', 'QI_v'], datetime.datetime(2019,11,29,3))

    # 3. N
    ds = get_grapes_variables(data_file_list, 
    ['N'], datetime.datetime(2019,11,29,3))

    print(ds.attrs['proj_info'])
    print(ds.data_vars['N'].attrs['proj_info'])

    a = WGS_to_GRAPES((39.69, 116.2), ds.attrs['proj_info'])
    b = WGS_to_GRAPES(((40.00, 117.0), (39.69, 116.2)), ds.attrs['proj_info'])

    print(a)
    print(b)
    