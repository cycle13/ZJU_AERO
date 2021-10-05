'''
Description: Some graph constant for the graph package
Author: Hejun Xie
Date: 2020-11-07 12:18:44
LastEditors: Hejun Xie
LastEditTime: 2021-10-04 20:05:04
'''

import numpy as np

units           = {
    'Q_hydro': 'kg*m^-2',
    'W': 'm*m^-1',
    'topo': 'm'}

long_names      = {
    'Q_hydro': 'Total Hydrometeor',
    'W': 'W component wind',
    'topo': 'Topograph'}

clevels         = {
    'Q_hydro': np.linspace(0., 30, 31),
    'W': np.linspace(-2.0, 2.0, 41),
    'topo': np.linspace(-100, 2000, 22)}

cmap            = {
    'Q_hydro': 'jet',
    'W': 'bwr',
    'topo': 'terrain'}
