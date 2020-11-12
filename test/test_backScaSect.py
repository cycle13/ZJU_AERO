'''
Description: test and plot the back scattering section of snow,
of IITM, TM, and Rayleigh scattering approximation.
Author: Hejun Xie
Date: 2020-11-12 17:03:03
LastEditors: Hejun Xie
LastEditTime: 2020-11-12 18:07:34
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global imports
import numpy as np
# np.set_printoptions(threshold=np.inf)
import xarray as xr

# Local imports
from pyCRO.db import load_lut
from pyCRO.hydrometeors.dielectric import dielectric_ice


def Ki_squared(f, T):
    """
    Computes the value of |Ki|^2 used in the definition of the refl. factor
    Args:
        f: the frequency [GHz]
        T: the Temperature [K]
    Returns:
        The value of |Ki|^2 at specified temperature and frequency.
    """

    m = dielectric_ice(T, f)
    k = (m**2-1)/(m**2 + 2)
    return np.abs(k)**2

def get_IITM_backScaSect(T, e):
    lut = load_lut('../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc', engine='xarray')
    sz = lut.lookup_line(t=[T], e=[e]).squeeze()

    radar_xsect_h = 2*np.pi*(sz[:,0]-sz[:,1]-sz[:,2]+sz[:,3]) # [mm2]
    print(radar_xsect_h)

    return radar_xsect_h
    
def get_TM_backScaSect(T, e):
    lut = load_lut('../pathos/lut/tmatrix_masc/lut_SZ_S_9_41_1mom.lut', engine='numpy')
    sz = lut.lookup_line(t=[T], e=[e]).squeeze()

    radar_xsect_h = 2*np.pi*(sz[:,0]-sz[:,1]-sz[:,2]+sz[:,3]) # [mm2]
    print(radar_xsect_h)

    return radar_xsect_h

def get_Rayleigh_backScaSect(T, f):
    Ki = Ki_squared(f, T)
    radar_xsect_h = np.pi**5 * D**6 / WAVELENGTH**4 * Ki * (RHO_S / RHO_I)**2 # [mm2]

    print(radar_xsect_h)
    return radar_xsect_h

if __name__ == '__main__':

    T = 253. # [K]
    e = 1.0  # [deg]
    f = 9.41 # [GHz]
    RHO_I = 916 # [kg m-3]
    RHO_S = 100 # [kg m-3]
    D = np.linspace(0.2, 20, 64) # [mm]
    C = 299792458. # [m s-1]
    WAVELENGTH = C / (f * 1.0E9) * 1000 # [mm]

    bss_IITM = get_IITM_backScaSect(T, e)
    # bss_TM = get_TM_backScaSect(T, e)
    bss_Rayleigh = get_Rayleigh_backScaSect(T, f)

    # start plot