'''
Description: test and plot the back scattering section of snow,
of IITM, TM, and Rayleigh scattering approximation.
Author: Hejun Xie
Date: 2020-11-12 17:03:03
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 11:12:54
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
from pyCRO.utilities import dielectric_ice


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
    radar_xsect_v = 2*np.pi*(sz[:,0]+sz[:,1]+sz[:,2]+sz[:,3]) # [mm2]

    return radar_xsect_h, radar_xsect_v
    
def get_TM_backScaSect(T, e):
    lut = load_lut('../pathos/lut/tmatrix_masc/lut_SZ_S_9_41_1mom.lut', engine='numpy')
    sz = lut.lookup_line(t=[T], e=[e]).squeeze()

    radar_xsect_h = 2*np.pi*(sz[:,0]-sz[:,1]-sz[:,2]+sz[:,3]) # [mm2]
    print(radar_xsect_h)

    return radar_xsect_h

def get_Rayleigh_backScaSect(T, f):
    Ki = Ki_squared(f, T)
    radar_xsect_h = np.pi**5 * D**6 / WAVELENGTH**4 * Ki * (RHO_S / RHO_I)**2 # [mm2]

    return radar_xsect_h

if __name__ == '__main__':

    T = 253. # [K]
    e = 1.0  # [deg]
    f = 9.41 # [GHz]
    RHO_I = 916 # [kg m-3]
    RHO_S = 100 # [kg m-3]
    D = np.linspace(0.2, 20, 64) # [mm]
    De = D * (RHO_S/RHO_I)**(1./3.) # [mm]
    C = 299792458. # [m s-1]
    WAVELENGTH = C / (f * 1.0E9) * 1000 # [mm]

    bss_h_IITM, bss_v_IITM = get_IITM_backScaSect(T, e)
    # bss_TM = get_TM_backScaSect(T, e)
    bss_Rayleigh = get_Rayleigh_backScaSect(T, f)

    # start plot
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.use('Agg')

    plt.rcParams['font.family'] = 'serif'
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.set_title('Horizontal Back scattering section')

    ax.plot(D, bss_h_IITM, color='r', ls='-', label='IITM orientation with preference Horizontal')
    ax.plot(D, bss_v_IITM, color='r', ls='--', label='IITM orientation with preference Vertical')
    ax.plot(D, bss_Rayleigh, color='k', ls='-', label='Rayleigh scattering Dmax')
    ax.plot(De, bss_Rayleigh, color='k', ls='--', label='Rayleigh scattering De')
    
    ax.set_xlabel(r'Dmax [$mm$]', fontsize=12)
    ax.set_ylabel(r'Back scattering section [$mm^{2}$]', fontsize=12)

    ax.set_yscale('log')

    ax.legend(frameon=False)
    
    plt.savefig('./backScaSect.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
