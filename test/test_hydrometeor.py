'''
Description: hydrometeor interpolation
Author: Hejun Xie
Date: 2020-08-19 12:55:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-13 12:37:07
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')


# Local imports
from pyCRO.config import cfg
cfg.init('./option_files/interpolation.yml')

from pyCRO.hydrometeors.hydrometeor import Snow, IceParticle 

if __name__ == "__main__":
    s = Snow('1mom')
    i = IceParticle('1mom')
