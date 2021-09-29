'''
Description: plot level A database sample
Author: Hejun Xie
Date: 2021-01-05 10:27:35
LastEditors: Hejun Xie
LastEditTime: 2021-09-12 16:43:30
'''

# Global imports
import xarray as xr
import numpy as np

# plot settings
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'

# FILE = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
# FILE = '../pathos/lut/tm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
FILE = '../pathos/lut/iitm_masc_snowflake/lut_SZ_S_9_41_1mom_LevelA.nc'

AR = 1.5 # [-]
Dmax = 20.0 # [mm]
T = 253.0 # [K]
Frequency = 9.41 # [GHz]
C = 299792458
lambda_ = C / (Frequency * 1E09) * 1000 # [mm]
# print(lambda_)

if __name__ == "__main__":
    
    with xr.open_dataset(FILE, engine="h5netcdf") as lut:

        lut_slice = lut.sel(aspect_ratio=AR, Dmax=Dmax, temperature=T, method='nearest')

        # (beta, elevation) [mm2]
        sigmav = (2 * np.pi) * (lut_slice['p11_bw'] + lut_slice['p22_bw'] + lut_slice['p12_bw'] + lut_slice['p21_bw']).data \
             * (lambda_/(2*np.pi))**2
        sigmah = (2 * np.pi) * (lut_slice['p11_bw'] + lut_slice['p22_bw'] - lut_slice['p12_bw'] - lut_slice['p21_bw']).data \
             * (lambda_/(2*np.pi))**2
        
        zdr = 10 * np.log10(sigmah / sigmav)
            
        beta = lut_slice.coords['beta'].data
        elevation = lut_slice.coords['elevation'].data


        TBETA,TELEVATION = np.meshgrid(beta,elevation)


        # start ploting
        clevels = np.linspace(min(sigmah.min(),sigmav.min()), max(sigmah.min(), sigmah.max()), 100)

        print(clevels)
        
        # zdr_abs_max = max(np.abs(zdr.min()), np.abs(zdr.max()))
        # clevels_zdr = np.linspace(-zdr_abs_max, zdr_abs_max, 100)

        clevels = np.arange(0, 2500, 5)
        # clevels = clevels * (Dmax / 20.0)**3
        clevels_zdr = np.arange(-2.0, 2.0, 0.05)

        # plot amplitude matrix output
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 6.5),constrained_layout=True)
        fig.suptitle('Oriented Particle Level-A Database sample', fontsize=18)

        CF1 = axs[0].contourf(TBETA, TELEVATION, sigmav.T, cmap='jet', levels=clevels, extend='both')
        CB1 = fig.colorbar(CF1, ax=axs[0], orientation='horizontal')
        CB1.set_label(r'Vertical Scattering Cross Section $\sigma_v$ [$mm^2$]', fontsize=14)
        axs[0].set_title(r'Vertical $\sigma_v$', fontsize=18)
        axs[0].set_xlabel(r'Euler Angle $\beta$ [degree]', fontsize=14)
        axs[0].set_ylabel(r'Elevation $e$ [degree]', fontsize=14)

        CF2 = axs[1].contourf(TBETA, TELEVATION, sigmah.T, cmap='jet', levels=clevels, extend='both')
        CB2 = fig.colorbar(CF2, ax=axs[1], orientation='horizontal')
        CB2.set_label(r'Horizontal Scattering Cross Section $\sigma_h$ [$mm^2$]', fontsize=14)
        axs[1].set_title(r'Horizontal $\sigma_h$', fontsize=18)
        axs[1].set_xlabel(r'Euler Angle $\beta$ [degree]', fontsize=14)
        axs[1].set_ylabel(r'Elevation $e$ [degree]', fontsize=14)

        CF3 = axs[2].contourf(TBETA, TELEVATION, zdr.T, cmap='RdBu_r', levels=clevels_zdr, extend='both')
        CB3 = fig.colorbar(CF3, ax=axs[2], orientation='horizontal')
        CB3.set_label(r'Differenrtial Reflectivity $Z_{DR}$ [dBZ]', fontsize=14)
        axs[2].set_title(r'Differenrtial Reflectivity $Z_{DR}$', fontsize=18)
        axs[2].set_xlabel(r'Euler Angle $\beta$ [degree]', fontsize=14)
        axs[2].set_ylabel(r'Elevation $e$ [degree]', fontsize=14)

        # plt.tight_layout()
        plt.savefig('levelA_test.png', dpi=300)
        plt.close()
