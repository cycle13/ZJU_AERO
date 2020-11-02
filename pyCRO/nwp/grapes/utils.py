'''
Description: Some utilities for grapes radar-operator interface
Author: Hejun Xie
Date: 2020-11-02 17:45:19
LastEditors: Hejun Xie
LastEditTime: 2020-11-02 18:58:26
'''

# Global import
import pyproj as pp
import numpy as np

def WGS_to_GRAPES(coords_WGS, proj_info):
    '''
    Description:
        Convert from WGS coordinates (lat, lon) to WRF coordinates (I, J)
    Params:
        coords_WGS: WGS coordinates, can be:
            (1). a single tuple (lat, lon)
            (2). g nested tuple or nested list ((lat, lon), (lat, lon), ...) or ([lat, lon], [lat, lon], ...)
            (3). a numpy ndarray [[lat, lon], [lat, lon], ...]
    Returns:
        GRAPES coordinates, for coords_WGS type
        (1) return: 1-D array [I, J]
        (2) and (3) return: 2-D array [[I, J], [I, J], ...]
    '''

    # convert into numpy.ndarray (ncoords, 2) or (lat, lon)
    if isinstance(coords_WGS, tuple) or isinstance(coords_WGS[0], tuple):
        coords_WGS = np.array(coords_WGS)

    x0, y0 = proj_info['LLC_LON'], proj_info['LLC_LAT']
    dx, dy = proj_info['DLON'], proj_info['DLAT']

    if len(coords_WGS.shape) == 1:
        x = coords_WGS[1]
        y = coords_WGS[0]
    else:
        x = coords_WGS[:, 1]
        y = coords_WGS[:, 0]

    I = (x-x0) / dx
    J = (y-y0) / dy

    if len(coords_WGS.shape) == 1:
        return np.array([I, J]).astype('float32')
    else:
        return np.array([I, J]).T.astype('float32')
