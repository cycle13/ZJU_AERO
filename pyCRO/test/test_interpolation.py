'''
Description: test interpolation
Author: Hejun Xie
Date: 2020-08-15 20:59:00
LastEditors: Hejun Xie
LastEditTime: 2020-08-15 22:05:58
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# global import
import numpy as np
import pyWRF as pw
import copy

# Local imports
from pyCRO.config import cfg
cfg.init('./option_files/interpolation.yml')
from pyCRO.interpolation import get_interpolated_radial
from pyCRO.utilities import DATAdecorator


BASE_VARIABLES = ['U','V','W','QR_v','QS_v','QG_v','QI_v','RHO','T']
MODEL_FILE = '../../../cosmo_pol/pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'

@DATAdecorator('./', True, './dic_variables.pkl')
def get_dic_variables(filename):
    
    file_h = pw.open_file(filename)
    vars_to_load = copy.deepcopy(BASE_VARIABLES)
    vars_ok = file_h.check_if_variables_in_file(['P','T','QV','QR','QC','QI','QS','QG','U','V','W'])
    vars_to_load.extend('N')


    print('Reading variables ' + str(vars_to_load) + ' from file')
    loaded_vars = file_h.get_variable(vars_to_load, itime=10, get_proj_info=True,
                                            shared_heights=False, assign_heights=True)
    
    # To deal with annoying issues with pickle
    for var in loaded_vars.values():
        var.file = None
    
    dic_vars = loaded_vars #  Assign to class
    if 'N' in loaded_vars.keys():
        N=loaded_vars['N']
        dic_vars.pop('N',None) # Remove N from the variable dictionnary (we won't need it there)
    file_h.close()

    print('-------done------')
    return dic_vars, N

# some unit tests
if __name__ == "__main__":

    dic_variables, N = get_dic_variables(MODEL_FILE)

    azimuth = 120
    elevation = 1.0
     
    list_subradials = get_interpolated_radial(
        dic_variables, 
        azimuth, elevation, 
        N=N, list_refraction=None
    )

    print(len(list_subradials))
    print(list_subradials[4].dist_profile)
    print(list_subradials[4].heights_profile)
    print(list_subradials[4].elev_profile)
    # print(list_subradials[5].values)
