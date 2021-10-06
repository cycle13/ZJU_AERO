'''
Description: plot snowflake radar optical database
Author: Hejun Xie
Date: 2021-10-05 21:08:46
LastEditors: Hejun Xie
LastEditTime: 2021-10-05 21:43:48
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

FILE = '../pathos/lut/iitm_masc_snowflake_20.0/lut_SZ_S_9_41_1mom_LevelB.nc'

T = 253.0 # [K]
Frequency = 9.41 # [GHz]
C = 299792458 # [m/s]
ele = 1.0 # [deg] 
lambda_ = C / (Frequency * 1E09) * 1000 # [mm]

if __name__ == "__main__":
    
    with xr.open_dataset(FILE, engine="h5netcdf") as lut:
        lut_slice = lut.sel(elevation=ele, temperature=T, method='nearest')

        # (Dmax) [mm]
        sigmav = (2 * np.pi) * (lut_slice['p11_bw'] + lut_slice['p22_bw'] + lut_slice['p12_bw'] + lut_slice['p21_bw']).data \
             * (lambda_/(2*np.pi))**2
        sigmah = (2 * np.pi) * (lut_slice['p11_bw'] + lut_slice['p22_bw'] - lut_slice['p12_bw'] - lut_slice['p21_bw']).data \
             * (lambda_/(2*np.pi))**2

        kdp = lut_slice['s22_fw'].real - lut_slice['s11_fw'].real
        zdr = 10 * np.log10(sigmah / sigmav)

        Dmax = lut_slice.coords['Dmax'].data

        x = Dmax * np.pi / lambda_

        '''
        Choose the X axis
        '''
        # xaxis = x 
        # xaxisname = r'Size Parameter $x$'
        xaxis = Dmax
        xaxisname = r'Maximum Diameter $D_{max}$ [mm]'

        '''
        Start plot level B database
        '''        
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(xaxis, sigmah, label=r'$\sigma_h$ [$mm^2$]', ls='-', color='k')
        ax.plot(xaxis, sigmav, label=r'$\sigma_v$ [$mm^2$]', ls='--', color='k')
        ax.set_yscale('log')
        ax.set_title('Oriented Particle Level-B Database sample', fontsize=18)
        ax.set_xlabel(r'Maximum Diameter $D_{max}$ [mm]', fontsize=14)
        ax.set_ylabel(r'Back Scattering Sector [$mm^2$]', fontsize=14)

        ax2 = ax.twinx()
        ax2.plot(xaxis, zdr, label='ZDR [dBZ]', ls='-', color='r')
        ax2.set_ylabel(r'Differential Reflectivity Factor [dBZ]', fontsize=14)
        # ax2.plot(xaxis, kdp, label=r'KDP [$mm$]', ls='-', color='r')
        # ax2.set_ylabel(r'Real Part of Forward Scattering Amplitudes [$mm$]', fontsize=14)
        
        ax.legend(frameon=False, bbox_to_anchor=(0.2, 0.95))
        ax2.legend(frameon=False, bbox_to_anchor=(0.2, 0.85))
    
        plt.tight_layout()
        plt.savefig('Snowflake_database_levelB.png', dpi=300)
        plt.close()
