'''
Description: test level A and Level B database
Author: Hejun Xie
Date: 2020-12-30 15:11:58
LastEditors: Hejun Xie
LastEditTime: 2021-01-02 19:07:54
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import xarray as xr
import numpy as np

# Local imports
from ZJU_AERO.db import load_lut

# plot settings
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'


IITM_LEVELA = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
IITM_LEVELB = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc'
IITM_LEVELB_EXP = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB_exp/t253.0e1.0.nc'
TM_LEVELA = '../pathos/lut/tm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
TM_LEVELB = '../pathos/lut/tm_masc/lut_SZ_S_9_41_1mom_LevelB.nc'

def plot(ptype, iitm_data, tm_data, Dmax, picname):

    fig, ax = plt.subplots(figsize=(12,6))
    
    ax.plot(Dmax,   tm_data, color='r', label='IITM hexagon'  )
    ax.plot(Dmax, iitm_data, color='b', label='EBCM spheroid' )

    ax.set_xlabel(r'Dmax [$mm$]', fontsize=14)

    if ptype == 'ZDR':
        ax.set_ylabel(r'$Z_{DR}$ [dBZ]', fontsize=14)
    elif ptype == 'KDP':
        ax.set_ylabel(r'$K_{DP}$ [-]', fontsize=14)
        # ax.set_yscale('log')

    ax.legend(loc='best', frameon=False, fontsize=10)

    plt.savefig(picname, dpi=300)
    plt.close()

if __name__ == "__main__":

    '''
    0. LEVELB API
    '''
    # db = load_lut(IITM_LEVELB, engine='xarray')
    # sz = db.lookup_line(e=[1.], t=[253.]).squeeze()
    # db.close()

    # zdr = (sz[:,0]-sz[:,1]-sz[:,2]+sz[:,3]) / (sz[:,0]+sz[:,1]+sz[:,2]+sz[:,3])
    # kdp = sz[:,10]-sz[:,8]
    # ZDR = 10 * np.log10(zdr)

    # print(ZDR)
    # print(kdp)
    # exit()
    

    '''
    1. LEVELB ZDR
    '''
    with xr.open_dataset(IITM_LEVELB, engine="h5netcdf") as iitm_levelb:
        Dmax = iitm_levelb.coords['Dmax'].data
        # iitm_slice = iitm_levelb.sel(temperature=253, elevation=1.0, Dmax=20.0, method='nearest')
        iitm_slice = iitm_levelb.sel(temperature=253, elevation=1.0, method='nearest')
        iitm_zdr = (iitm_slice['p11_bw'] - iitm_slice['p12_bw'] - \
           iitm_slice['p21_bw'] + iitm_slice['p22_bw']) / \
          (iitm_slice['p11_bw'] + iitm_slice['p12_bw'] + \
           iitm_slice['p21_bw'] + iitm_slice['p22_bw'])
        iitm_ZDR = 10 * np.log10(iitm_zdr.data)
        iitm_kdp = iitm_slice['s22_fw'].data.real - iitm_slice['s11_fw'].data.real
    
    with xr.open_dataset(TM_LEVELB, engine="h5netcdf") as tm_levelb:
        tm_slice = tm_levelb.sel(temperature=253, elevation=1.0, method='nearest')
        tm_zdr = (tm_slice['p11_bw'] - tm_slice['p12_bw'] - \
           tm_slice['p21_bw'] + tm_slice['p22_bw']) / \
          (tm_slice['p11_bw'] + tm_slice['p12_bw'] + \
           tm_slice['p21_bw'] + tm_slice['p22_bw'])
        tm_ZDR = 10 * np.log10(tm_zdr.data)
        tm_kdp = tm_slice['s22_fw'].data.real - tm_slice['s11_fw'].data.real
    
    plot('ZDR', iitm_ZDR, tm_ZDR, Dmax, 'ZDR_LEVELB.png')
    plot('KDP', iitm_kdp, tm_kdp, Dmax, 'KDP_LEVELB.png')
    # exit()

    '''
    2. LEVELB EXP
    '''
    
    # with xr.open_dataset(IITM_LEVELB, engine="h5netcdf") as iitm_levelb:
    #     iitm_slice = iitm_levelb.sel(temperature=253, elevation=1.0, method='nearest')
    #     print(iitm_slice['p12_bw'].data)
    
    # with xr.open_dataset(IITM_LEVELB_EXP, engine="h5netcdf") as iitm_levelb_exp:
    #     iitm_slice_exp = iitm_levelb_exp.interp(std_a=40.0, std_b=-0.077, lambda_a=8.42, 
    #         lambda_b=-0.57, mu_a=0.053, mu_b=0.79)
    #     # iitm_slice_exp = iitm_levelb_exp.sel(std_a=40.0, std_b=-0.08, lambda_a=6.0, 
    #     #     lambda_b=-0.6, mu_a=0.055, mu_b=0.8, method='nearest')
    #     print(iitm_slice_exp['p12_bw'].data)
    # exit()

    '''
    3. LEVELA ZDR
    '''
    AR = 1.5
    BETA = 30.0

    wavelength = 299792458. / (9.41E09) * 1000 # [mm]
    
    with xr.open_dataset(IITM_LEVELA, engine="h5netcdf") as iitm_levela:

        iitm_slice = iitm_levela.sel(temperature=253, aspect_ratio=AR,
        beta=BETA, elevation=1.0, method='nearest')

        iitm_zh = 2 * np.pi * (iitm_slice['p11_bw'] - iitm_slice['p12_bw'] - \
            iitm_slice['p21_bw'] + iitm_slice['p22_bw'])
        
        iitm_zv = 2 * np.pi * (iitm_slice['p11_bw'] + iitm_slice['p12_bw'] + \
            iitm_slice['p21_bw'] + iitm_slice['p22_bw'])
        
        iitm_zdr = iitm_zh / iitm_zv

        iitm_kdp = iitm_slice['s22_fw'].data.real - iitm_slice['s11_fw'].data.real
    
    with xr.open_dataset(TM_LEVELA, engine="h5netcdf") as tm_levela:

        tm_slice = tm_levela.sel(temperature=253, aspect_ratio=AR,
        beta=BETA, elevation=1.0, method='nearest')

        tm_zh = 2 * np.pi * (tm_slice['p11_bw'] - tm_slice['p12_bw'] - \
            tm_slice['p21_bw'] + tm_slice['p22_bw'])
        
        tm_zv = 2 * np.pi * (tm_slice['p11_bw'] + tm_slice['p12_bw'] + \
            tm_slice['p21_bw'] + tm_slice['p22_bw'])
        
        tm_zdr = tm_zh / tm_zv

        tm_kdp = tm_slice['s22_fw'].data.real - tm_slice['s11_fw'].data.real

    iitm_ZDR = 10 * np.log10(iitm_zdr.data)
    tm_ZDR = 10 * np.log10(tm_zdr.data)
    
    plot('ZDR', iitm_ZDR, tm_ZDR, Dmax, 'ZDR_LEVELA.png')
    plot('KDP', iitm_kdp, tm_kdp, Dmax, 'KDP_LEVELA.png')
