'''
Description: test the radar operator core by
ploting intrinsic variables ZH, ZDR, KDP against
hydrometeor concentration: SNOW, GRAUPEL, RAIN, etc.
Author: Hejun Xie
Date: 2020-11-14 18:03:56
LastEditors: Hejun Xie
LastEditTime: 2021-07-20 20:12:39
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np

# Local imports
from ZJU_AERO.core._rdop_scatter import one_rad_one_hydro, get_pol_from_sz 
from ZJU_AERO.interp import Radial
from ZJU_AERO.db import load_lut
from ZJU_AERO.hydro.hydrometeor import Snow, NonsphericalSnow
from ZJU_AERO.config import createConfig
from ZJU_AERO.const import global_constants as constants

if __name__ == "__main__":
    # get global constants
    createConfig('./option_files/wsm6_test.yml')
    from ZJU_AERO.config.cfg import CONFIG
    constants.update()

    hydro_name = 'S'
    
    db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_9_41_1mom_LevelB.nc', engine='xarray')
    hydro_istc = Snow('1mom')

    # db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_35_0_1mom_LevelB.nc', engine='xarray')
    # hydro_istc = NonsphericalSnow('1mom', 'spheroid')
    
    # db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_35_0_1mom_LevelB.nc', engine='xarray')
    # hydro_istc = NonsphericalSnow('1mom', 'spheroid')

    # db = load_lut('../pathos/lut/iitm_masc/lut_SZ_S_35_0_1mom_LevelB.nc', engine='xarray')
    # hydro_istc = NonsphericalSnow('1mom', 'hexcol')

    # db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_9_41_1mom_LevelB.nc', engine='xarray')
    # hydro_istc = NonsphericalSnow('1mom', 'spheroid')

    # db = load_lut('../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc', engine='xarray')
    # hydro_istc = NonsphericalSnow('1mom', 'hexcol')
    
    # set profile
    ntest = 1000
    e = np.ones((ntest), dtype='float32') * 1.0 # [deg]
    dic_values = dict()
    dic_values['T'] = np.ones((ntest), dtype='float32') * 253. # [K]
    dic_values['QS_v'] = np.logspace(-5, -2, ntest) # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1)
    
    test_rad = Radial(dic_values, elev_profile=e,
    mask=None, lats_profile=None, lons_profile=None,
    dist_ground_profile=None, heights_profile=None)

    # print(test_rad.values)
    # print(test_rad.elev_profile)

    '''
    sz_psd_integ: 
            dimension: (ntest, 12) 
            unit: Z[mm2 m-3]; S[mm m-3]
    '''
    valid_data, sz_psd_integ = one_rad_one_hydro(test_rad, hydro_name, hydro_istc,
    db, simulate_doppler=False, ngates=ntest)

    '''
    ZH: radar refl. factor at hor. pol. in linear units [mm6 m-3]
    ZV: radar refl. factor at vert. pol. in linear units [mm6 m-3]
    ZDR: diff. refl. = z_h / z_v [-]
    RHOHV: copolar. corr. coeff [-]
    KDP: spec. diff. phase shift upon propagation [deg km-1]
    AH: spec. att. at hor. pol. [dB km-1]
    AV: spec. att. at vert. pol. [dB km-1]
    DELTA_HV: total phase shift upon backscattering [deg]
    '''
    ZH, ZV, ZDR, RHOHV, KDP, AH, AV, DELTA_HV = get_pol_from_sz(sz_psd_integ)
    dbZH = 10 * np.log10(ZH) # [dBZ]
    dbZDR = 10 * np.log10(ZDR) # [dB]
    
    print(dbZH)
    print(dbZDR)
    print(KDP)

    db.close()
