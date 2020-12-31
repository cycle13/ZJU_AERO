'''
Description: test level A and Level B database
Author: Hejun Xie
Date: 2020-12-30 15:11:58
LastEditors: Hejun Xie
LastEditTime: 2020-12-31 19:25:49
'''

import xarray as xr
import numpy as np

IITM_LEVELA = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
IITM_LEVELB = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc'
IITM_LEVELB_EXP = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB_exp/t253.0e1.0.nc'
TM_LEVELA = '../pathos/lut/tm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
TM_LEVELB = '../pathos/lut/tm_masc/lut_SZ_S_9_41_1mom_LevelB.nc'

if __name__ == "__main__":
    
    with xr.open_dataset(IITM_LEVELB, engine="h5netcdf") as iitm_levelb:
        iitm_slice = iitm_levelb.sel(temperature=253, elevation=1.0, method='nearest')
        print(iitm_slice['p12_bw'].data)
    
    with xr.open_dataset(IITM_LEVELB_EXP, engine="h5netcdf") as iitm_levelb_exp:
        iitm_slice_exp = iitm_levelb_exp.interp(std_a=40.0, std_b=-0.077, lambda_a=8.42, 
            lambda_b=-0.57, mu_a=0.053, mu_b=0.79)
        # iitm_slice_exp = iitm_levelb_exp.sel(std_a=40.0, std_b=-0.08, lambda_a=6.0, 
        #     lambda_b=-0.6, mu_a=0.055, mu_b=0.8, method='nearest')
        print(iitm_slice_exp['p12_bw'].data)
    
    exit()

    AR = 3.0
    BETA = 30.0
    
    with xr.open_dataset(IITM_LEVELA, engine="h5netcdf") as iitm_levela:

        iitm_slice = iitm_levela.sel(temperature=253, aspect_ratio=AR,
        beta=BETA, elevation=1.0, method='nearest')

        iitm_zh = 2 * np.pi * (iitm_slice['p11_bw'] - iitm_slice['p12_bw'] - \
            iitm_slice['p21_bw'] + iitm_slice['p22_bw'])
        
        iitm_zv = 2 * np.pi * (iitm_slice['p11_bw'] + iitm_slice['p12_bw'] + \
            iitm_slice['p21_bw'] + iitm_slice['p22_bw'])
        
        iitm_zdr = iitm_zh / iitm_zv
    
    with xr.open_dataset(TM_LEVELA, engine="h5netcdf") as tm_levela:

        tm_slice = tm_levela.sel(temperature=253, aspect_ratio=AR,
        beta=BETA, elevation=1.0, method='nearest')

        tm_zh = 2 * np.pi * (tm_slice['p11_bw'] - tm_slice['p12_bw'] - \
            tm_slice['p21_bw'] + tm_slice['p22_bw'])
        
        tm_zv = 2 * np.pi * (tm_slice['p11_bw'] + tm_slice['p12_bw'] + \
            tm_slice['p21_bw'] + tm_slice['p22_bw'])
        
        tm_zdr = tm_zh / tm_zv

    print(np.log10(iitm_zdr.data))
    print(np.log10(tm_zdr.data))
