'''
Description: test config subpackage
Author: Hejun Xie
Date: 2020-11-21 22:41:58
LastEditors: Hejun Xie
LastEditTime: 2020-11-22 11:19:14
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')


# Local import 
from ZJU_AERO.config.cfg import createConfig

if __name__ == "__main__":
    # cc = ConfigClass('./option_files/configtest.yml')
    # cc.check_sanity()
    # print(cc['radar'])

    createConfig('./option_files/configtest.yml')
    from ZJU_AERO.config.cfg import CONFIG
    print(CONFIG)

