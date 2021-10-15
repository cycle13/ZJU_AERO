'''
Description: plot snowflake radar optical database
Author: Hejun Xie
Date: 2021-10-05 21:08:46
LastEditors: Hejun Xie
LastEditTime: 2021-10-15 19:57:34
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


# PARAMS = ['8.2', '9.5', '11.0', '14.0', '17.0', '20.0']
PARAMS = ['8.2', '20.0']
# COLORS = ['red', 'orange', 'yellow', 'green', 'blue', 'purple']
COLORS = ['red', 'purple']

# FILES = ['../pathos/lut/iitm_masc_snowflake_'+PARAM+'/lut_SZ_S_9_41_1mom_LevelB.nc' for PARAM in PARAMS]
FILES = ['../pathos/lut/iitm_masc_snowflake_'+PARAM+'/lut_SZ_S_35_0_1mom_LevelB.nc' for PARAM in PARAMS]
# HEX = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc'
HEX = '../pathos/lut/iitm_masc/lut_SZ_S_35_0_1mom_LevelB.nc'
# SPH = '../pathos/lut/tm_masc_release/lut_SZ_S_9_41_1mom_LevelB.nc'
SPH = '../pathos/lut/tm_masc_release/lut_SZ_S_35_0_1mom_LevelB.nc'

T = 253.0 # [K]
# Frequency = 9.41 # [GHz]
Frequency = 35.0 # [GHz]
C = 299792458 # [m/s]
ele = 1.0 # [deg] 
lambda_ = C / (Frequency * 1E09) * 1000 # [mm]
nD = 64


def get_zdr(FILE):
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
     
     return zdr, Dmax


if __name__ == "__main__":


     ZDR = np.zeros((nD, len(PARAMS)+2), dtype='float')
     for ifile, FILE in enumerate(FILES):
          zdr, Dmax = get_zdr(FILE)
          ZDR[:,ifile] = zdr
     
     zdr, Dmax = get_zdr(HEX)
     ZDR[:,-2] = zdr

     zdr, Dmax = get_zdr(SPH)
     ZDR[:,-1] = zdr

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

     for iparam, param in enumerate(PARAMS):
          ax.plot(xaxis, ZDR[:,iparam], label='n2/n1={}'.format(param), ls='-', color=COLORS[iparam], marker='x', markersize=3)
     
     ax.plot(xaxis, ZDR[:,-2], label='hexagon', ls='--', color='k', marker='x', markersize=3)
     ax.plot(xaxis, ZDR[:,-1], label='spheroid', ls='-', color='k', marker='x', markersize=3)
     
     ax.set_xlabel(r'Maximum Diameter $D_{max}$ [mm]', fontsize=14)
     ax.set_ylabel(r'Differential Reflectivity [$dBZ$]', fontsize=14)

     ax.legend(frameon=False, bbox_to_anchor=(0.2, 0.95))

     plt.tight_layout()
     plt.savefig('Snowflake_database_levelB.png', dpi=300)
     plt.close()
