'''
Description: test pycwr
Author: Hejun Xie
Date: 2020-08-30 22:20:45
LastEditors: Hejun Xie
LastEditTime: 2020-10-08 11:12:41
'''

import pycwr

from pycwr.io import read_auto

file = r"/mnt/e/RADAR/19.5.17房山-2/BJXFS.20190517.100000.AR2.bz2"

PRD = read_auto(file)
PyartRadar = PRD.ToPyartRadar()

print(PRD.scan_info)
print(PRD.fields[0]['KDP'])

# exit()

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

isweep = 1
elevation = PRD.scan_info['fixed_angle'][isweep]

for field in fields:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    graph = GraphMap(PRD, ccrs.PlateCarree())
    graph.plot_ppi_map(ax, isweep, field)
    ax.set_title("Observation {} [{}] \n e={:.3f} UTC: 2019-05-17 10:00".format(PYCRO_name[field],
    units[field],float(elevation)), fontsize=15)
    plt.savefig('observation_{}.png'.format(PYCRO_name[field]), dpi=300)

    del ax, graph, fig
