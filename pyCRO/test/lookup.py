'''
Description: test lookup table 
Author: Hejun Xie
Date: 2020-08-19 22:17:23
LastEditors: Hejun Xie
LastEditTime: 2020-08-19 22:29:46
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global import
import numpy as np

# Local import
from pyCRO.lookup import load_all_lut

if __name__ == "__main__":

    lut_sz = load_all_lut('1mom', ['R', 'S', 'G'], 35.6, 'tmatrix_masc', folder_lut='../../../cosmo_pol/pathos/lut/')
    print(lut_sz['S'])

