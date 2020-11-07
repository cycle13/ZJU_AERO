'''
Description: test interpolation
Author: Hejun Xie
Date: 2020-08-15 20:59:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 23:15:28
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# global import
import numpy as np
# import pyWRF as pw
import copy
import datetime as dt

# Local imports
from pyCRO.config import cfg
from pyCRO.constants import global_constants as constants
from pyCRO.interpolation import get_interpolated_radial
from pyCRO.utilities import DATAdecorator
from pyCRO.nwp.wrf import get_wrf_variables, check_if_variables_in_file

BASE_VARIABLES = ['U','V','W','QR_v','QS_v','QG_v','QI_v','RHO','T']
MODEL_FILE = '../../../cosmo_pol/pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'

@DATAdecorator('./', False, './dic_variables.pkl')
def get_dic_variables(filename):
    
    vars_to_load = copy.deepcopy(BASE_VARIABLES)
    vars_ok = check_if_variables_in_file(['P','T','QV','QR','QC','QI','QS','QG','U','V','W'])
    vars_to_load.extend('N')

    print('Reading variables ' + str(vars_to_load) + ' from file')
    loaded_vars = get_wrf_variables([MODEL_FILE], vars_to_load, dt.datetime(2013,10,6,10))

    dic_vars = loaded_vars
    
    if 'N' in list(loaded_vars.data_vars.keys()):
            N=loaded_vars['N']
            dic_vars = dic_vars.drop_vars('N') # Remove N from the xarray dataset (we won't need it there)

    return dic_vars, N

# some unit tests
if __name__ == "__main__":

    # get global constants
    cfg.init('./option_files/interpolation.yml')
    cfg.CONFIG = cfg.sanity_check(cfg.CONFIG)
    constants.update()

    dic_variables, N = get_dic_variables(MODEL_FILE)

    azimuth = 120
    elevation = 1.0
     
    list_subradials = get_interpolated_radial(
        dic_variables, 
        azimuth, elevation, 
        N=N, list_refraction=None
    )
    
    # exit()

    print(len(list_subradials))
    print(list_subradials[4].dist_profile)
    print(list_subradials[4].heights_profile)
    print(list_subradials[4].elev_profile)
    # print(list_subradials[5].values)
