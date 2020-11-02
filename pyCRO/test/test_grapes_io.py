'''
Description: test nwp.grapes_io.py
Author: Hejun Xie
Date: 2020-11-01 11:12:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-02 10:16:47
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
    ds = get_grapes_variables(data_file_list, ['U', 'V', 'W', 'T'], datetime.datetime(2019,11,29,3))
