'''
Description: hydrometeor interpolation
Author: Hejun Xie
Date: 2020-08-19 12:55:00
LastEditors: Hejun Xie
LastEditTime: 2020-12-11 10:36:01
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np

# Local imports
from ZJU_AERO.config.cfg import createConfig
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.hydro.hydrometeor import Snow, IceParticle, Rain, Graupel

if __name__ == "__main__":
    createConfig('./option_files/interpolation.yml')
    constants.update()

    s = Snow('1mom')
    i = IceParticle('1mom')
    r = Rain('1mom')

    print( r.get_aspect_ratio(np.array([9.0, 10.0, 12.0])) )
