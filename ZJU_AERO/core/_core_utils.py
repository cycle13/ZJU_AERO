'''
Description: Some utilities used in core submodule
Author: Hejun Xie
Date: 2021-03-01 09:57:09
LastEditors: Hejun Xie
LastEditTime: 2021-03-01 21:55:45
'''

import numpy as np
from ..const import global_constants as constants

rad2deg = 180. / np.pi
deg2rad = np.pi / 180.

def _get_varray():
    '''
    Get the radial velocity array 
    Unit: [m s-1]
    '''
    return constants.VARRAY

def _get_wavelength():
    '''
    Get wavelength Initialized in global constants
    Unit: [mm]
    '''
    return constants.WAVELENGTH

def _get_K_squared():
    '''
    Get Kw Initialized in global constants
    Unit: [-]
    '''
    return constants.KW

def _get_rho_0():
    '''
    Get the air density at MSL (rho_0) in global constants
    Unit: [kg m-3]
    '''
    return constants.RHO_0

def _get_refl_coeff():
    '''
    Get the coeffcient of reflectivity 
    Unit: [mm4]
    '''
    return _get_wavelength()**4 / (np.pi**5 * _get_K_squared())

def get_rcs_h(sz):
    '''
    Get the radar (back-scattering) cross section of horizontal polarization
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        rcs_h: radar cross section of horizontal polarization [mm2] (single) / [mm2 m-3] (bulk)
    '''
    return 2 * np.pi * (sz[...,0] - sz[...,1] - sz[...,2] + sz[...,3])
    
def get_rcs_v(sz):
    '''
    Get the radar (back-scattering) cross section of vertical polarization
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        rcs_h: radar cross section of vertical polarization [mm2] (single) / [mm2 m-3] (bulk)
    '''
    return 2 * np.pi * (sz[...,0] + sz[...,1] + sz[...,2] + sz[...,3])

def get_chv(sz):
    '''
    Get the back-scattering amplitude covariance
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        chv: back-scattering chv = shh * svv* (complex) [mm2] (single) / [mm2 m-3] (bulk)
    '''

    return (sz[..., 4] + sz[..., 7]) * 0.5 + (sz[..., 6] - sz[..., 5]) * 0.5 * 1j

def get_ecs_h(sz):
    '''
    Get the extinction cross section of horizontal polarization
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        ecs_h: extinction cross section of horizontal polarization [mm2] (single) / [mm2 m-3] (bulk)
    '''
    return _get_wavelength() * 2 * sz[..., 11]

def get_ecs_v(sz):
    '''
    Get the extinction cross section of vertical polarization
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        ecs_v: extinction cross section of vertical polarization [mm2] (single) / [mm2 m-3] (bulk)
    '''
    return _get_wavelength() * 2 * sz[..., 9]

def get_z_h(sz):
    '''
    Get the horizontal reflectivity factor [mm6 m-4]
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        zh: horizontal reflectivity factor [mm6 m-4]
    '''
    return get_rcs_h(sz) * _get_refl_coeff()

def get_z_v(sz):
    '''
    Get the vertical reflectivity factor [mm6 m-4]
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        zh: vertical reflectivity factor [mm6 m-4]
    '''
    return get_rcs_v(sz) * _get_refl_coeff()

def get_zdr(sz):
    '''
    Get the differential reflectivity factor [-]
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        zdr: differential reflectivity factor [-]
    '''
    return get_rcs_h(sz) / get_rcs_v(sz)

def get_kdp(sz):
    '''
    Get the differential phase shift [deg km-1]
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        kdp: differential phase shift [deg km-1]
    '''
    return 1e-3 * rad2deg * _get_wavelength() * (sz[...,10] - sz[..., 8])

def get_att_h(sz):
    '''
    Get the attenuation factor in horizontal polarization [dB km-1]
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        att_h: horizontal attenuation factor [dB km-1]
    '''

    '''
    4.343e-3 = 1e-3 * 10 * np.log10(np.e)
    a. factor 1e-3:  convert from [dB mm2 m-3] to [dB km-1]
    c. factor 10 * np.log10(np.e) : convert from exp-based to dB-based
    '''
    return 4.343e-3 * get_ecs_h(sz)

def get_att_v(sz):
    '''
    Get the attenuation factor in vertical polarization [dB km-1]
    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        att_v: vertical attenuation factor [dB km-1]
    '''

    '''
    4.343e-3 = 1e-3 * 10 * np.log10(np.e)
    a. factor 1e-3:  convert from [dB mm2 m-3] to [dB km-1]
    c. factor 10 * np.log10(np.e) : convert from exp-based to dB-based
    '''
    return 4.343e-3 * get_ecs_v(sz)

def get_delta_hv(sz):
    '''
    Get total differential phase shift upon backscattering [deg]

    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        delta_hv: total differential phase shift upon backscattering [deg]
    '''

    # Convert from FSA (Forward scattering alignment) to BSA (Backward scattering alignment)
    chv_fsa = get_chv(sz)
    chv_bsa = - chv_fsa.real - chv_fsa.imag * 1j
    return rad2deg * np.angle(chv_bsa)

def get_rho_hv(sz):
    '''
    Get backscattering coherence coefficient [-]

    Params:
        sz: packed Backward phase matrix Z [mm2] (single) / [mm2 m-3] (bulk) and
            Forward amplitude matrix S [mm] (single) / [mm m-3] (bulk)
            (..., sz_idx)
            sz_idx: [z11, z12, z21, z22, 
                     z33, z34, z43, z44,
                     svr, svi, shr, shi]
    Returns:
        rho_hv: backscattering coherence coefficient [-]
    '''
    shh2 = get_rcs_h(sz) / (4 * np.pi)
    svv2 = get_rcs_v(sz) / (4 * np.pi)
    chv = get_chv(sz)
    chv2 = chv.real**2 + chv.imag**2

    return np.sqrt( chv2 / (shh2 * svv2) ) 

def get_rho_corr(rho):
    '''
    Get the rho correction coefficient for terminal velocity
    Params:
        rho: the air density [kg m-3]
    Returns:
        rho_corr: rho correction coefficient for terminal velocity [-]
    '''
    return ( _get_rho_0() / rho )**(0.5)

def proj_vel(U, V, W, vf, theta, phi):
    """
    Gets the radial velocity from the 3D wind field and hydrometeor
    fall velocity
    Args:
        U: eastward wind component [m/s]
        V: northward wind component [m/s]
        W: vertical wind component [m/s]
        vf: terminal fall velocity [m/s]
        theta: elevation angle [degree]
        phi: azimuth angle [degree]

    Returns:
        The radial velocity, with reference to the radar beam
        positive values represent flow away from the radar [m/s]
    """

    theta = np.deg2rad(theta) # elevation
    phi   = np.deg2rad(phi) # azimuth

    return ((U*np.sin(phi) + V * np.cos(phi)) * np.cos(theta)
            + (W - vf) * np.sin(theta))

def proj_vel_back(U, V, W, vrad, theta, phi):
    """
    Gets the hydrometeor fall velocity from the 3D wind field and 
    radial velocity
    Args:
        U: eastward wind component [m/s]
        V: northward wind component [m/s]
        W: vertical wind component [m/s]
        vrad: radial velocity [m/s]
        theta: elevation angle [degree]
        phi: azimuth angle [degree]

    Returns:
        The terminal fall velocity of hydrometeors 
    """

    theta = np.deg2rad(theta) # elevation
    phi   = np.deg2rad(phi) # azimuth

    return W + (U * np.sin(phi) + V * np.cos(phi)) / np.tan(theta) - \
        vrad / np.sin(theta)
