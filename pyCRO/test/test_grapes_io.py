'''
Description: test nwp.grapes_io.py
Author: Hejun Xie
Date: 2020-11-01 11:12:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-02 12:41:53
'''


# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global import
import glob
import os
import datetime

FOLDER = '/mnt/e/GRAPES/north_china_snowfall_20191112'
data_file_list = glob.glob(FOLDER+os.sep+'*.nc')

from pyCRO.nwp.grapes import get_grapes_variables

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

    print(ds)
    