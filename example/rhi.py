'''
Description: test for scatter
Author: Hejun Xie
Date: 2020-08-22 12:36:55
LastEditors: Hejun Xie
LastEditTime: 2020-11-11 10:20:05
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np
import pickle
import datetime as dt
# import pyart

# Local imports
import ZJU_AERO
import pyart

LOAD_MODEL = False
LOAD_RADAR = False
DEG = r'$^\circ$'

cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17',
'KDP': 'pyart_EWilson17', 'PHIDP': 'pyart_Carbone42', 'RHOHV': 'pyart_GrMg16'}

if __name__ == "__main__":
    FILENAME = '../pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'
    a = ZJU_AERO.RadarOperator(options_file='./option_files/simulate_db_acess.yml')
    a.load_model_file([FILENAME], load_datetime=dt.datetime(2013, 10, 6, 10), load_from_file=LOAD_MODEL, load_file='mdl.nc')

    if not LOAD_RADAR:
        r = a.get_RHI(azimuths = 120, elev_step = 0.2)
        with open("./rhi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./rhi.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()
    # exit()
    
    from pyart.graph import RadarDisplay
    display = RadarDisplay(r)
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('Agg')
    
    fig = plt.figure(figsize=[10, 4])
    ax = fig.add_subplot(111)

    field = 'ZH'
    vrange = (0, 60)
    display.plot(field, 0, vmin=vrange[0], vmax=vrange[1],
                    cmap=cmap[field],
                    title= 'Time: {}'.format(a.get_pos_and_time()['time']) + '\n' + \
                            'Azimuth: {}'.format(r.azimuth['data'][0]) + DEG + '\n' + \
                            r'$Z_{H}$',
                    ax=ax)
    
    display.set_limits(ylim=[0, 17])
    
    plt.savefig('ZH_rhi_noatt.png',dpi=300, bbox_inches='tight')
