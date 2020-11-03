'''
Description: test GRAPES interface for radar operator
Author: Hejun Xie
Date: 2020-11-02 16:17:47
LastEditors: Hejun Xie
LastEditTime: 2020-11-03 11:58:38
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global imports
import numpy as np
import pickle
import os
import glob
import datetime as dt

# Local imports
import pyCRO
import pyart

LOAD_MODEL = False
LOAD_RADAR = False
DEG = r'$^\circ$'

cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17',
'KDP': 'pyart_EWilson17', 'PHIDP': 'pyart_Carbone42', 'RHOHV': 'pyart_GrMg16'}

if __name__ == "__main__":

    FOLDER = '/mnt/e/GRAPES/north_china_snowfall_20191112'
    data_file_list = glob.glob(FOLDER+os.sep+'*.nc')
    load_datetime = dt.datetime(2019,11,29,12)
    
    a = pyCRO.RadarOperator(options_file='./option_files/grapes_interface.yml')
    a.load_model_file(data_file_list, load_datetime=load_datetime, load_from_file=LOAD_MODEL, load_file='mdl.nc')

    if not LOAD_RADAR:
        r = a.get_PPI(elevations = 1, plot_engine='pycwr')
        # r = a.get_PPI_test(elevations = 1)
        with open("./ppi.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./ppi.pkl", "rb") as f:
            r = pickle.load(f)
    
    a.close()

    # exit()
    PRD = r

    print(PRD.scan_info)

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from pycwr.draw.RadarPlot import Graph, GraphMap

    import matplotlib as mpl
    mpl.use('Agg')

    fields = ['dBZ', 'ZDR', 'KDP', 'PhiDP', 'V', 'CC']
    PYCRO_name = {'dBZ': 'ZH',
            'ZDR':'ZDR',
            'KDP':'KDP',
            'PhiDP':'PHIDP',
            'CC': 'RHOHV',
            'V': 'RVEL'
            }

    units = {'dBZ': 'dBZ',
            'ZDR': 'dBZ',
            'KDP': 'degrees/km',
            'PhiDP': 'degrees',
            'V': 'm/s',
            'CC': '-'
    }

    isweep = 0
    elevation = PRD.scan_info['fixed_angle'][isweep]

    for field in fields:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        graph = GraphMap(PRD, ccrs.PlateCarree())
        graph.plot_ppi_map(ax, isweep, field)
        ax.set_title("Simulation {} [{}] \n e={:.3f} UTC: {}".format(PYCRO_name[field],
        units[field], float(elevation), str(load_datetime)), fontsize=15)
        plt.savefig('grapes_simulation_{}.png'.format(PYCRO_name[field]), dpi=300)

        del fig, ax, graph

