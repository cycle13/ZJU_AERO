'''
Description: test for scatter
Author: Hejun Xie
Date: 2020-08-22 12:36:55
LastEditors: Hejun Xie
LastEditTime: 2020-09-19 15:55:31
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global imports
import numpy as np
import pickle
# import pyart

# Local imports
import pyCRO
import pyart

LOAD_MODEL = True
LOAD_RADAR = True
DEG = r'$^\circ$'

cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17',
'KDP': 'pyart_EWilson17', 'PHIDP': 'pyart_Carbone42', 'RHOHV': 'pyart_GrMg16'}

if __name__ == "__main__":
    FILENAME = '../../../cosmo_pol/pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'
    a = pyCRO.RadarOperator(options_file='./option_files/simulate.yml')
    a.load_model_file(FILENAME, itime=10, load_pickle=LOAD_MODEL, pickle_file='mdl.pkl')

    if not LOAD_RADAR:
        r = a.get_PPI(elevations = 1)
        with open("./ppi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./ppi.pkl", "rb") as f:
            r = pickle.load(f)
    

    # exit()
    import matplotlib as mpl
    mpl.use('Agg')
    from pyart.graph import RadarMapDisplayBasemap
    display = pyart.graph.RadarMapDisplayBasemap(r)
    import matplotlib.pyplot as plt
    plt.figure()

    field = 'ZH'
    vrange = (0, 60)
    display.plot_ppi_map(field, 0, vmin=vrange[0], vmax=vrange[1],
                     min_lon=119, max_lon=122.5, min_lat=26.3, max_lat=29.5,
                     lon_lines=np.arange(119, 122.7, 1), projection='lcc',
                     lat_lines=np.arange(26.3, 29.5, 1), resolution='h',
                     lat_0=r.latitude['data'],
                     lon_0=r.longitude['data'],
                     cmap=cmap[field],
                     title= 'Time: {}'.format(a.get_pos_and_time()['time']) + '\n' + \
                            'Elevation: {}'.format(r.elevation['data'][0]) + DEG + '\n' + \
                            r'$Z_{H}$')
    # plot range rings at 10, 20, 30 and 40km
    display.plot_range_ring(50., line_style='k-', lw=1.0)
    display.plot_range_ring(100., line_style='k--', lw=1.0)
    display.plot_range_ring(150., line_style='k-', lw=1.0)

    # plots cross hairs
    display.plot_line_xy(np.array([-200000.0, 200000.0]), np.array([0.0, 0.0]),
                        line_style='k-', lw=1.2)
    display.plot_line_xy(np.array([0.0, 0.0]), np.array([-200000.0, 200000.0]),
                        line_style='k-', lw=1.2)

    # display.plot_point(r.longitude['data'], r.latitude['data'])
    
    plt.savefig('ZH_ppi.png',dpi=300, bbox_inches='tight')
    