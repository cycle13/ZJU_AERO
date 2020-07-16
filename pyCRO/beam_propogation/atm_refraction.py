'''
@Description: defines a set of functions to compute the path of a
radar beam while taking into account atmospheric refraction
@Author: Hejun Xie
@Date: 2020-07-16 17:44:29
@LastEditors: Hejun Xie
@LastEditTime: 2020-07-16 17:47:55
'''

# local import
from .fixed_radius import fixed_radius_KE
from .ODE_solver import ODEZeng2014


def compute_trajectory_radial(range_vec, elevation_angle, coords_radar,
                         refraction_method, N = 0):
    '''
    Computes the trajectory of a radar beam along a specified radial by
    calling the appropriate subfunction depending on the desired method
    Args:
        range_vec: vector of all ranges along the radial [m]
        elevation_angle: elevation angle of the beam [degrees]
        coord_radar: radar coordinates as 3D tuple [lat, lon, alt], ex
            [47.35, 12.3, 682]
        refraction_method: the method to be used to compute the trajectory
            can be either 1, for the 4/3 method, or 2 for the Zeng and Blahak
            (2014) ODE method
        N: atmospheric refractivity as a COSMO variable, needs to be provided
            only for refraction_method == 2

    Returns:
        s: vector of distance at the ground along the radial [m]
        h: vector of heights above ground along the radial [m]
        e: vector of incident elevation angles along the radial [degrees]
    '''

    if refraction_method==1:
        s, h, e = fixed_radius_KE(range_vec, elevation_angle, coords_radar)
    elif refraction_method==2:
        s, h, e = ODEZeng2014(range_vec, elevation_angle, coords_radar, N)

    return s, h, e
