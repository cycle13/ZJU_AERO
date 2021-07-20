'''
Description: test nonspherical hydrometeor
Author: Hejun Xie
Date: 2021-06-13 18:23:40
LastEditors: Hejun Xie
LastEditTime: 2021-07-20 20:53:01
'''


# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np
import time

# Local imports
from ZJU_AERO.config.cfg import createConfig
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.hydro.hydrometeor import Snow, NonsphericalSnow

if __name__ == "__main__":
    createConfig('./option_files/interpolation.yml')
    constants.update()

    s = Snow('1mom')
    ns = NonsphericalSnow('1mom', 'hexcol')

    ntest = 2666
    
    # T = np.array([250, 250, 250, 250])
    # QM = np.array([1e-5, 1e-4, 1e-3, 1e-2])

    T = np.ones((ntest), dtype='float32') * 253. # [K]
    QM = np.logspace(-5, -2, ntest) # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1)

    t1 = time.time()
    s.set_psd(T, QM)
    t2 = time.time()
    # print(s.lambda_)
    print(t2-t1)

    # exit()

    t1 = time.time()
    ns.set_psd(T, QM)
    t2 = time.time()
    # print(ns.lambda_)
    print(t2-t1)
