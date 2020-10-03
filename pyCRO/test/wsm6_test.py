'''
Description: test for scatter
Author: Hejun Xie
Date: 2020-08-22 12:36:55
LastEditors: Hejun Xie
LastEditTime: 2020-10-03 19:37:04
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

LOAD_MODEL = False
LOAD_RADAR = False
DEG = r'$^\circ$'

cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17',
'KDP': 'pyart_EWilson17', 'PHIDP': 'pyart_Carbone42', 'RHOHV': 'pyart_GrMg16'}

if __name__ == "__main__":
    FILENAME = '../../../cosmo_pol/pathos/WRF/wsm6_test/wrfout_d01_2019-05-17_00_00_00'
    a = pyCRO.RadarOperator(options_file='./option_files/wsm6_test.yml')
    a.load_model_file(FILENAME, itime=40, load_pickle=LOAD_MODEL, pickle_file='mdl.pkl')

    # print(a.dic_vars['T'])

    if not LOAD_RADAR:
        r = a.get_PPI(elevations = 1.494, plot_engine='pycwr')
        with open("./ppi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./ppi.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()

    # test PycwrRadop
    # PRD = PycwrRadop('ppi',r)
    PRD = r

    print(PRD.scan_info)

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from pycwr.draw.RadarPlot import Graph, GraphMap

    import matplotlib as mpl
    mpl.use('Agg')

    isweep = 0

    elevation = PRD.scan_info['fixed_angle'][isweep]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    graph = GraphMap(PRD, ccrs.PlateCarree())
    graph.plot_ppi_map(ax, isweep, "dBZ", cmap="pyart_NWSRef")
    ax.set_title("PPI with map e={}".format(float(elevation)), fontsize=16)
    plt.savefig('pycwr_{}.png'.format(isweep), dpi=300)

    exit()
    
    # test PyartRadop
    import matplotlib as mpl
    mpl.use('Agg')
    from pyart.graph import RadarMapDisplayBasemap
    display = pyart.graph.RadarMapDisplayBasemap(r)
    import matplotlib.pyplot as plt
    plt.figure()

    field = 'ZH'
    vrange = (0, 60)
    display.plot_ppi_map(field, 0, vmin=vrange[0], vmax=vrange[1],
                     min_lon=113.5, max_lon=119, min_lat=37.5, max_lat=42.0,
                     lon_lines=np.arange(38, 42, 1), projection='lcc',
                     lat_lines=np.arange(113, 119, 1), resolution='h',
                     lat_0=r.latitude['data'],
                     lon_0=r.longitude['data'],
                     cmap=cmap[field],
                     title= 'Time: {}'.format(a.get_pos_and_time()['time']) + '\n' + \
                            'Elevation: {}'.format(r.elevation['data'][0]) + DEG + '\n' + \
                            r'$Z_{H}$')
    # plot range rings at 50, 100, 200km
    display.plot_range_ring(50., line_style='k-', lw=1.0)
    display.plot_range_ring(100., line_style='k--', lw=1.0)
    display.plot_range_ring(200., line_style='k-', lw=1.0)

    # plots cross hairs
    display.plot_line_xy(np.array([-300000.0, 300000.0]), np.array([0.0, 0.0]),
                        line_style='k-', lw=1.2)
    display.plot_line_xy(np.array([0.0, 0.0]), np.array([-300000.0, 300000.0]),
                        line_style='k-', lw=1.2)

    # display.plot_point(r.longitude['data'], r.latitude['data'])
    
    plt.savefig('ZH_ppi.png',dpi=300, bbox_inches='tight')
    