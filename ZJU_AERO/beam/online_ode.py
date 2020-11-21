'''
@Description: ODE (Odinary Differential Equation) solver 
for the trajectory of the radar beam
@Author: Hejun Xie
@Date: 2020-07-16 11:48:54
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 15:42:27
'''

# Global imports
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error
# import pyWRF as pw
from scipy.interpolate import interp1d
from scipy.integrate import odeint

# Local imports
from ..utils import get_earth_radius
from ..const import global_constants as constants

def _deriv_z(z, r, n_h_int, dn_dh_int, RE):
    '''
    Updates the state vector of the system of ODE used in the ODE refraction
    by Blahak and Zeng
    (2014)
    Args:
        z: state vector in the form of a tuple (height, sin(theta))
        r: range vector in m (not actually used in the state vector)
        n_h_spline: piecewise linear interpolator for the refractive index as
            a function of the altitude
        dn_dh_int: piecewise linear interpolator for the derivative of the
            refractive index as a function of the altitude
        RE: earth radius [m]
    Returns:
        An updated state vector
    '''
    # Computes the derivatives (RHS) of the system of ODEs
    h, u = z
    n = n_h_int(h)
    dn_dh = dn_dh_int(h)
    return [u, (-u ** 2 * ((1./n) * dn_dh + 1./ (RE + h)) +
                ((1. / n) * dn_dh + 1. / (RE + h)))]

def _piecewise_linear(x,y):
    '''
    Defines a piecewise linear interpolator, used to interpolate refractivity
    values between COSMO vertical coordinates
    Args:
        x: vector of independent variable
        y: vector of dependent variable

    Returns:
        A piecewise linear interpolator
    '''
    interpolator=interp1d(x,y)
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        if np.isscalar(xs):
            xs=[xs]
        return np.array([pointwise(xi) for xi in xs])

    return ufunclike


def Zeng2014(range_vec, elevation_angle, coords_radar, N):
    '''
    Computes the trajectory of a radar beam along a specified radial with the
    ODE method of Zeng and Blahak (2014)
    TODO: adapt for other NWP models
    Args:
        range_vec: vector of all ranges along the radial [m]
        elevation_angle: elevation angle of the beam [degrees]
        coord_radar: radar coordinates as 3D tuple [lat, lon, alt], ex
            [47.35, 12.3, 682]
        N: atmospheric refractivity as a COSMO variable

    Returns:
        s: vector of distance at the ground along the radial [m]
        h: vector of heights above ground along the radial [m]
        e: vector of incident elevation angles along the radial [degrees]
    '''

    from ..config.config_proc import CONFIG
    if CONFIG['nwp']['name'] == 'grapes':
        from ..nwp.grapes import WGS_to_GRAPES as WGS_to_MODEL
    elif CONFIG['nwp']['name'] == 'wrf':
        from ..nwp.wrf import WGS_to_WRF as WGS_to_MODEL

    # Get info about NWP coordinate system
    proj_MODEL = N.attrs
    
    # Convert WGS84 coordinates to MODEL coordinates --> model grid [I, J]
    coords_rad_in_MODEL = WGS_to_MODEL(coords_radar, proj_MODEL)

    # TODO: For model grid only, model grid (I, J)
    llc_MODEL = (0., 0.)
    res_MODEL = [1., 1.]

    # Get index of radar in NWP coordinates (I, J)
    pos_radar_bin = [(coords_rad_in_MODEL[0]-llc_MODEL[0]) / res_MODEL[0],
                    (coords_rad_in_MODEL[1]-llc_MODEL[1]) / res_MODEL[1]]

    # Get refractive index profile from refractivity estimated from NWP variables
    # data (bottom_top, south_north, west_east)
    n_vert_profile = 1 + (N.data[:,int(np.round(pos_radar_bin[1])),
                             int(np.round(pos_radar_bin[0]))]) * 1E-6
    # Get corresponding altitudes
    h = N.coords['z-levels'][:,int(np.round(pos_radar_bin[1])),
                                int(np.round(pos_radar_bin[0]))]

    # Get earth radius at radar latitude
    RE = get_earth_radius(coords_radar[0])

    # Create piecewise linear interpolation for n as a function of height
    n_h_int = _piecewise_linear(h, n_vert_profile)
    dn_dh_int = _piecewise_linear(h[0:-1], np.diff(n_vert_profile) / np.diff(h))

    z_0 = [coords_radar[2], np.sin(np.deg2rad(elevation_angle))]
    # Solve second-order ODE
    Z = odeint(_deriv_z, z_0, range_vec, args = (n_h_int, dn_dh_int, RE))
    h = Z[:,0] # Heights above ground
    e = np.arcsin(Z[:,1]) # Elevations
    s = np.zeros(h.shape) # Arc distances
    dR = range_vec[1]-range_vec[0]
    s[0] = 0

    for i in range(1,len(s)): # Solve for arc distances
        s[i] = s[i-1] + RE * np.arcsin((np.cos(e[i-1]) * dR) / (RE + h[i]))

    s = s.astype('float32')
    h = h.astype('float32')
    e = np.rad2deg(e.astype('float32'))
    return s, h, e

# unit tests
if __name__ == "__main__":

    FILENAME = '../../../cosmo_pol/pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'
    file_h = pw.open_file(FILENAME)
    d = file_h.get_variable(['N'], itime=10, assign_heights=True)
    
    range_vec = np.array([0, 500, 1000, 1500, 2000])
    elevation_angle = [1.0]
    coords_radar = [30, 120, 100]
    s, h, e =  Zeng2014(range_vec, elevation_angle, coords_radar, d['N'])
    print(s)
    print(h)
    print(e)
