'''
Description: plot snowflake radar optical database
Author: Hejun Xie
Date: 2021-10-05 21:08:46
LastEditors: Hejun Xie
LastEditTime: 2021-10-15 20:05:30
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
from ZJU_AERO.hydro.hydrometeor import Snow, NonsphericalSnow, Graupel, NonsphericalGraupel
from ZJU_AERO.config import createConfig
from ZJU_AERO.const import global_constants as constants

ntest = 100
e = np.ones((ntest), dtype='float32') * 1.0 # [deg]

def test_core(hydro_istc, hydro_name, db):

	dic_values = dict()
	dic_values['T'] = np.ones((ntest), dtype='float32') * 253. # [K]
	if hydro_name == 'S':
		dic_values['QS_v'] = np.logspace(-5, -2, ntest) # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1)
	elif hydro_name == 'G':
		dic_values['QG_v'] = np.logspace(-5, -2, ntest) # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1)
	
	test_rad = Radial(dic_values, elev_profile=e,
	mask=None, lats_profile=None, lons_profile=None,
	dist_ground_profile=None, heights_profile=None)

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

	# print(dbZH)
	
	db.close()
	
	return dbZH, dbZDR, KDP


if __name__ == "__main__":
    # get global constants
    createConfig('../test/option_files/wsm6_test.yml')
    from ZJU_AERO.config.cfg import CONFIG
    constants.update()

    freq_str = '9_41'

    hydro_istc = NonsphericalSnow('1mom', 'snowflake', 8.2)
    db = load_lut('../pathos/lut/iitm_masc_snowflake_8.2/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN_8_2, ZDR_SN_8_2, KDP_SN_8_2 = test_core(hydro_istc, 'S', db)

    hydro_istc = NonsphericalSnow('1mom', 'snowflake', 9.5)
    db = load_lut('../pathos/lut/iitm_masc_snowflake_9.5/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN_9_5, ZDR_SN_9_5, KDP_SN_9_5 = test_core(hydro_istc, 'S', db)

    hydro_istc = NonsphericalSnow('1mom', 'snowflake', 11.0)
    db = load_lut('../pathos/lut/iitm_masc_snowflake_11.0/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN_11, ZDR_SN_11, KDP_SN_11 = test_core(hydro_istc, 'S', db)

    hydro_istc = NonsphericalSnow('1mom', 'snowflake', 14.0)
    db = load_lut('../pathos/lut/iitm_masc_snowflake_14.0/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN_14, ZDR_SN_14, KDP_SN_14 = test_core(hydro_istc, 'S', db)

    hydro_istc = NonsphericalSnow('1mom', 'snowflake', 17.0)
    db = load_lut('../pathos/lut/iitm_masc_snowflake_17.0/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN_17, ZDR_SN_17, KDP_SN_17 = test_core(hydro_istc, 'S', db)

    hydro_istc = NonsphericalSnow('1mom', 'snowflake', 20.0)
    db = load_lut('../pathos/lut/iitm_masc_snowflake_20.0/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN_20, ZDR_SN_20, KDP_SN_20 = test_core(hydro_istc, 'S', db)
	
    hydro_istc = NonsphericalSnow('1mom', 'hexcol')
    db = load_lut('../pathos/lut/iitm_masc/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_HN, ZDR_HN, KDP_HN = test_core(hydro_istc, 'S', db)

    hydro_istc = NonsphericalSnow('1mom', 'spheroid')
    db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
    ZH_SN, ZDR_SN, KDP_SN = test_core(hydro_istc, 'S', db)

    # Plot imports
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'    

    fig, axes = plt.subplots(3, 1, figsize=(10,12), sharex=True)
    # fig.subplots_adjust(hspace=0)

    # ZH
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN_8_2, color='r', label='Snow Snowflake n2/n1=8.2', ls='-')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN_9_5, color='orange', label='Snow Snowflake n2/n1=9.5', ls='-')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN_11, color='yellow', label='Snow Snowflake n2/n1=11.0', ls='-')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN_14, color='green', label='Snow Snowflake n2/n1=14.0', ls='-')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN_17, color='blue', label='Snow Snowflake n2/n1=17.0', ls='-')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN_20, color='purple', label='Snow Snowflake n2/n1=20.0', ls='-')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_HN, color='k', label='Snow Hexcol', ls='--')
    axes[0].plot(np.logspace(-5, -2, ntest), ZH_SN, color='k', label='Snow Spheroid', ls='-')
    axes[0].set_ylabel(r'$Z_{H}$ [dBZ]', fontsize=14)

    # ZDR
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN_8_2, color='r', label='Snow Snowflake n2/n1=8.2', ls='-')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN_9_5, color='orange', label='Snow Snowflake n2/n1=9.5', ls='-')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN_11, color='yellow', label='Snow Snowflake n2/n1=11.0', ls='-')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN_14, color='green', label='Snow Snowflake n2/n1=14.0', ls='-')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN_17, color='blue', label='Snow Snowflake n2/n1=17.0', ls='-')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN_20, color='purple', label='Snow Snowflake n2/n1=20.0', ls='-')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_HN, color='k', label='Snow Hexcol', ls='--')
    axes[1].plot(np.logspace(-5, -2, ntest), ZDR_SN, color='k', label='Snow Spheroid', ls='-')
    axes[1].set_ylabel(r'$Z_{DR}$ [dBZ]', fontsize=14)

    # KDP
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN_8_2, color='r', label='Snow Snowflake n2/n1=8.2', ls='-')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN_9_5, color='orange', label='Snow Snowflake n2/n1=9.5', ls='-')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN_11, color='yellow', label='Snow Snowflake n2/n1=11.0', ls='-')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN_14, color='green', label='Snow Snowflake n2/n1=14.0', ls='-')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN_17, color='blue', label='Snow Snowflake n2/n1=17.0', ls='-')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN_20, color='purple', label='Snow Snowflake n2/n1=20.0', ls='-')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_HN, color='k', label='Snow Hexcol', ls='--')
    axes[2].plot(np.logspace(-5, -2, ntest), KDP_SN, color='k', label='Snow Spheroid', ls='-')
    axes[2].set_ylabel(r'$K_{DP}$ [$deg \cdot km^{-1}$]', fontsize=14)

    axes[2].set_xscale('log')
    axes[2].set_xlabel(r'Ice water content [$kg m-3$]', fontsize=14)

    axes[0].legend(loc='best', frameon=False, fontsize=8)
    plt.savefig('Snowflake_database_levelC.png', dpi=300)
    plt.close()
