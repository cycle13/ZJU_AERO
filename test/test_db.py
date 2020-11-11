'''
Description: test lookup table 
Author: Hejun Xie
Date: 2020-08-19 22:17:23
LastEditors: Hejun Xie
LastEditTime: 2020-11-11 10:14:23
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global import
import numpy as np

# Local import
from pyCRO.db import load_all_lut

if __name__ == "__main__":

    scattering_method = {'S':'tmatrix_masc', 'R':'tmatrix_masc', 'G': 'tmatrix_masc'}

    lut_sz = load_all_lut('1mom', ['R', 'S', 'G'], 35.6, scattering_method, folder_lut='../pathos/lut/')
    print(lut_sz['S'])

