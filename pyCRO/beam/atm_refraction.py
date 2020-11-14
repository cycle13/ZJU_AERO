'''
@Description: defines a set of functions to compute the path of a
radar beam while taking into account atmospheric refraction
@Author: Hejun Xie
@Date: 2020-07-16 17:44:29
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 11:55:51
'''

# Global imports
import numpy as np

# Local imports
from ..constants import global_constants as constants
from ..utilities import get_earth_radius
from .effective_earth_radius import effective_earth_radius
from .online_ode_exhaustive import Zeng2014_exhaustive
from .online_ode import Zeng2014


def compute_trajectory_radial(range_vec, elevation_angle, coords_radar,
                         refraction_method, N=None, azimuth_angle=None):
    '''
    Computes the trajectory of a radar beam along a specified radial by
    calling the appropriate subfunction depending on the desired method
    Args:
        range_vec: vector of all ranges along the radial [m]
        elevation_angle: elevation angle of the beam [degrees], be a list if refraction_method==3
        coord_radar: radar coordinates as 3D tuple [lat, lon, alt], ex
            [47.35, 12.3, 682]
        refraction_method: the method to be used to compute the trajectory
            can be either 1, for the 4/3 method, or 2 for the Zeng and Blahak
            (2014) ODE method
        N: atmospheric refractivity as a NWP variable, needs to be provided
            only for refraction_method == 2 and 3
        azimuth_angle: azimuth angle of the beam [degrees], needs to be provided
            only for refraction_method == 3

    Returns:
        s: vector of distance at the ground along the radial [m], be a list of s if refraction_method==3
        h: vector of heights above ground along the radial [m], be a list if of h if refraction_method==3
        e: vector of incident elevation angles along the radial [degrees], be a list of e if refraction_method==3
    '''

    if refraction_method==1:
        s, h, e = effective_earth_radius(range_vec, elevation_angle, coords_radar)
    elif refraction_method==2:
        s, h, e = Zeng2014(range_vec, elevation_angle, coords_radar, N)
    elif refraction_method==3:
        s, h, e = Zeng2014_exhaustive(range_vec, elevation_angle, azimuth_angle, coords_radar, N)

    return s, h, e


def compute_trajectory_spaceborne(elevation):
    '''
    Computes the trajectory of a spaceborne radar beam along a specified radial,
    currently atmospheric refraction is not taken into account
    Args:
        elevation: elevation angle of the radial in degrees

    Returns:
        s: vector of distance at the ground along the radial [m]
        h: vector of heights above ground along the radial [m]
        e: vector of incident elevation angles along the radial [degrees]
    '''
    from ..config.cfg import CONFIG

    # Get info about spaceborne radar position
    latitude = CONFIG['radar']['coords'][0]
    altitude_radar = CONFIG['radar']['coords'][2]
    max_range = CONFIG['radar']['range']
    radial_resolution = CONFIG['radar']['radial_resolution']

    elev_rad = np.deg2rad(elevation)

    # For spaceborne refraction is simply ignored...
    maxHeightMODEL = constants.MAX_MODEL_HEIGHT
    RE = get_earth_radius(latitude)
    # Compute maximum range to target (using cosinus law in the triangle
    # earth center-radar-target)

    range_vec=np.arange(radial_resolution/2.,max_range,radial_resolution)

    h = np.sqrt((altitude_radar + RE)**2 + range_vec**2 - 2*(altitude_radar + RE)*range_vec*np.sin(elev_rad)) - RE
    alpha = np.arcsin((range_vec * np.cos(elev_rad)) / (RE + h))
    s = RE * alpha
    e = elevation - np.rad2deg(alpha)

    in_lower_atm = [h < maxHeightMODEL]

    h = h[in_lower_atm]
    s = s[in_lower_atm]
    e = e[in_lower_atm]

    s = s.astype('float32')
    h = h.astype('float32')
    e = e.astype('float32')

    return s,h,e
