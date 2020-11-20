'''
Description: Some utilities for wrf radar-operator interface
Author: Hejun Xie
Date: 2020-11-05 19:31:44
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 22:04:15
'''


# Global import
import pyproj as pp
import numpy as np
import sys

def WGS_to_WRF(coords_WGS, proj_info):
    '''
    Description:
        Convert from WGS coordinates (lat, lon) to WRF coordinates (I, J)
        currently only for latlong projection
    Params:
        coords_WGS: WGS coordinates, can be:
            (1). a single tuple (lat, lon)
            (2). g nested tuple or nested list (lat_list, lon_list) 
            (3). a numpy ndarray [[lat1, lat2, ...], [lon1, lon2, ...]]
    Returns:
        GRAPES coordinates, for coords_WGS type
        (1) return: 1-D array [I, J]
        (2) and (3) return: 2-D array [[I1, J1], [I2, J2], ...]
    '''

    # convert into numpy.ndarray (ncoords, 2) or (lat, lon)
    if isinstance(coords_WGS, tuple) or isinstance(coords_WGS[0], tuple) or \
        isinstance(coords_WGS, list) or isinstance(coords_WGS[0], list):
        coords_WGS = np.array(coords_WGS)
    
    lon = coords_WGS[1]
    lat = coords_WGS[0]

    nx, ny = proj_info['nI'], proj_info['nJ']
    dx, dy = proj_info['DX'], proj_info['DY']

    # define the datum
    wrf_proj = pp.Proj(proj='lcc', # projection type: Lambert Conformal Conic
                       lat_1=proj_info['TRUELAT1'], lat_2=proj_info['TRUELAT2'], # Cone intersects with the sphere
                       lat_0=proj_info['MOAD_CEN_LAT'], lon_0=proj_info['STAND_LON'], # Center point
                       a=6370000, b=6370000) # This is it! The Earth is a perfect sphere
    
    # Easting and Northings of the domains center point
    wgs_proj = pp.Proj(proj='latlong', datum='WGS84')

    # to resolve python versions deprecation of pyproj
    python_version = sys.version_info
    if sys.version_info[0] >= 3:
        transformer = pp.Transformer.from_proj(wgs_proj, wrf_proj)
    if sys.version_info[0] >= 3:
        cen_lambert_x, cen_lambert_y = transformer.transform(proj_info['CEN_LON'], proj_info['CEN_LAT'])
    else:
        cen_lambert_x, cen_lambert_y = pp.transform(wgs_proj, wrf_proj, proj_info['CEN_LON'], proj_info['CEN_LAT'])
    
    if sys.version_info[0] >= 3:
        proj_lambert_x, proj_lambert_y = transformer.transform(lon, lat)
    else:
        proj_lambert_x, proj_lambert_y = pyproj.transform(wgs_proj, wrf_proj, lon, lat)


    x = proj_lambert_x - (cen_lambert_x - dx * (nx - 1) / 2.)
    y = proj_lambert_y - (cen_lambert_y - dy * (ny - 1) / 2.)

    I, J = x / dx, y / dy

    if len(coords_WGS.shape) == 1:
        return np.array([I, J]).astype('float32')
    else:
        return np.array([I, J]).T.astype('float32')

if __name__ == "__main__":
    # unit test
    dic_proj = {'TRUELAT1':30., 'TRUELAT2':60.,
                'MOAD_CEN_LAT':30., 'STAND_LON':125.,
                'CEN_LAT':28.8117, 'CEN_LON':123.2079,
                'DX':3000., 'DY':3000.,
                'nI':360, 'nJ':222}
    
    coords_WGS = np.array([[25.610306, 31.735355], [118.02133, 128.86807]])

    coords_WRF = WGS_to_WRF(coords_WGS, dic_proj)

    print(coords_WRF)
