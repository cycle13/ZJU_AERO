'''
Description: hydrometeor interpolation
Author: Hejun Xie
Date: 2020-08-19 12:55:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-20 17:06:13
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')


# Local imports
from ZJU_AERO.config import cfg
cfg.init('./option_files/interpolation.yml')

from ZJU_AERO.hydro.hydrometeor import Snow, IceParticle 

if __name__ == "__main__":
    s = Snow('1mom')
    i = IceParticle('1mom')
