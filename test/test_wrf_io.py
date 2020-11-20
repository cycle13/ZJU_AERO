'''
Description: test nwp.wrf_io.py
Author: Hejun Xie
Date: 2020-11-01 11:12:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-11 10:17:32
'''


# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global import
import glob
import os
import datetime

FOLDER = '../pathos/WRF/wsm6_test/ERA_interim'
data_file_list = glob.glob(FOLDER+os.sep+'wrfout_d01_*')

from ZJU_AERO.nwp.wrf import get_wrf_variables
from ZJU_AERO.nwp.wrf import WGS_to_WRF

if __name__ == "__main__":
    
    # 1. U, V, W, T, P, Pw, RHO
    # ds = get_wrf_variables(data_file_list, 
    # ['U', 'V', 'W', 'T', 'P', 'Pw', 'RHO'], datetime.datetime(2019,5,17,10))

    # 2. QV_v, QR_v, QS_v, QG_v, QC_v, QI_v
    # ds = get_wrf_variables(data_file_list, 
    # ['QV_v', 'QR_v', 'QS_v', 'QG_v', 'QC_v', 'QI_v'], datetime.datetime(2019,5,17,10))

    # 3. N
    ds = get_wrf_variables(data_file_list, 
    ['N'], datetime.datetime(2019,5,17,10))

    print(ds)

    a = WGS_to_WRF((39.69, 116.2), ds.attrs)
    b = WGS_to_WRF(((40.00, 39.69), (117.0, 116.2)), ds.attrs)

    print(a)
    print(b)

