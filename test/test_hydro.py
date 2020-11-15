'''
Description: hydrometeor interpolation
Author: Hejun Xie
Date: 2020-08-19 12:55:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 12:00:43
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')


# Local imports
from pyCRO.config import cfg
cfg.init('./option_files/interpolation.yml')

from pyCRO.hydro.hydrometeor import Snow, IceParticle 

if __name__ == "__main__":
    s = Snow('1mom')
    i = IceParticle('1mom')