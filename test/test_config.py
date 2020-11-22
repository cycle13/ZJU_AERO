'''
Description: test config subpackage
Author: Hejun Xie
Date: 2020-11-21 22:41:58
LastEditors: Hejun Xie
LastEditTime: 2020-11-22 10:35:02
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')


# Local import 
from ZJU_AERO.config.config_proc_new import createConfig

if __name__ == "__main__":
    # cc = ConfigClass('./option_files/configtest.yml')
    # cc.check_sanity()
    # print(cc['radar'])

    createConfig('./option_files/configtest.yml')
    from ZJU_AERO.config.config_proc_new import CONFIG
    print(CONFIG)

