'''
@Description: 3/4 earth radius model (see Doviak and Zrnic, p.21)
@Author: Hejun Xie
@Date: 2020-07-16 09:47:37
@LastEditors: Hejun Xie
@LastEditTime: 2020-07-16 11:00:34
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# global import
import numpy as np

# local import
from pyCRO.constants import global_constants as constants
from pyCRO.utilities import get_earth_radius

def fixed_radius_(range_vec, elevation_angle, coords_radar):
    '''
    Computes the trajectory of a radar beam along a specified radial wwith
    the simple 4/3 earth radius model (see Doviak and Zrnic, p.21)
    Args:
        range_vec: vector of all ranges along the radial [m]
        elevation_angle: elevation angle of the beam [degrees]
        coord_radar: radar coordinates as 3D tuple [lat, lon, alt], ex
            [47.35, 12.3, 682]

    Returns:
        s: vector of distance at the ground along the radial [m]
        h: vector of heights above ground along the radial [m]
        e: vector of incident elevation angles along the radial [degrees]
    '''
    # elevation_angle must be in radians in the formula
    elevation_angle = np.deg2rad(elevation_angle)
    KE = constants.KE

    altitude_radar=coords_radar[2]
    latitude_radar=coords_radar[1]

    # Compute earth radius at radar latitude
    RE = get_earth_radius(latitude_radar)
    # Compute height over radar of every range_bin
    temp = np.sqrt(range_vec ** 2 + (KE * RE) ** 2 + 2 * range_vec *
                  KE * RE * np.sin(elevation_angle))

    h = temp - KE * RE + altitude_radar
    # Compute arc distance of every range bin
    s = KE * RE * np.arcsin((range_vec * np.cos(elevation_angle)) /
                                   (KE * RE + h))
    e = elevation_angle + np.arctan(range_vec * np.cos(elevation_angle) /
                                    (range_vec*np.sin(elevation_angle) + KE *
                                     RE + altitude_radar))
    s = s.astype('float32')
    h = h.astype('float32')
    e = np.rad2deg(e.astype('float32'))

    return s,h,e

# unit tests
if __name__ == "__main__":
    range_vec = np.array([0, 500, 1000, 1500, 2000])
    elevation_angle = [1.0]
    coords_radar = [120, 30, 100]
    s, h, e =  fixed_radius_(range_vec, elevation_angle, coords_radar)
    print(s)
    print(h)
    print(e)
