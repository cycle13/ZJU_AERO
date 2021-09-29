'''
Description: test graupel and snow
Author: Hejun Xie
Date: 2021-06-23 20:47:23
LastEditors: Hejun Xie
LastEditTime: 2021-09-29 17:14:32
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

	print(dbZH)
	
	db.close()
	
	return dbZH, dbZDR, KDP


if __name__ == "__main__":
    # get global constants
	createConfig('./option_files/wsm6_test.yml')
	from ZJU_AERO.config.cfg import CONFIG
	constants.update()

	# freq_str = '9_41'
	freq_str = '35_0'
	
	hydro_istc = Snow('1mom')
	db = load_lut('../pathos/lut/iitm_masc/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_S_HO, ZDR_S_HO, KDP_S_HO = test_core(hydro_istc, 'S', db)

	hydro_istc = NonsphericalSnow('1mom', 'snowflake', 20.0)
	db = load_lut('../pathos/lut/iitm_masc_snowflake_20.0/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_S_SN_20, ZDR_S_SN_20, KDP_S_SN_20 = test_core(hydro_istc, 'S', db)

	hydro_istc = NonsphericalSnow('1mom', 'snowflake', 8.2)
	db = load_lut('../pathos/lut/iitm_masc_snowflake_8.2/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_S_SN_8_2, ZDR_S_SN_8_2, KDP_S_SN_8_2 = test_core(hydro_istc, 'S', db)

	hydro_istc = NonsphericalSnow('1mom', 'hexcol')
	db = load_lut('../pathos/lut/iitm_masc/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_S_HN, ZDR_S_HN, KDP_S_HN = test_core(hydro_istc, 'S', db)

	hydro_istc = Snow('1mom')
	db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_S_SO, ZDR_S_SO, KDP_S_SO = test_core(hydro_istc, 'S', db)

	hydro_istc = NonsphericalSnow('1mom', 'spheroid')
	db = load_lut('../pathos/lut/tm_masc_release/lut_SZ_S_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_S_SN, ZDR_S_SN, KDP_S_SN = test_core(hydro_istc, 'S', db)

	hydro_istc = NonsphericalGraupel('1mom', 'spheroid')
	db = load_lut('../pathos/lut/tm_masc/lut_SZ_G_'+freq_str+'_1mom_LevelB.nc', engine='xarray')
	ZH_G, ZDR_G, KDP_G = test_core(hydro_istc, 'G', db)

	# Plot imports
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	plt.rcParams['font.family'] = 'serif'    

	fig, axes = plt.subplots(3, 1, figsize=(10,12), sharex=True)
	# fig.subplots_adjust(hspace=0)

    # ZH
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_G, color='r', label='Graupel', ls='-')
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_S_SO, color='k', label='Snow Spheroid old m-D scheme', ls='-')
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_S_SN, color='b', label='Snow Spheroid new m-D scheme', ls='-')
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_S_SN_20, color='g', label='Snow Snowflake new m-D scheme n2/n1=20.0', ls='--')
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_S_SN_8_2, color='g', label='Snow Snowflake new m-D scheme n2/n1=8.2', ls='-')
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_S_HO, color='k', label='Snow HexCol old m-D scheme', ls='--')
	axes[0].plot(np.logspace(-5, -2, ntest), ZH_S_HN, color='b', label='Snow HexCol new m-D scheme', ls='--')
	
	axes[0].set_ylabel(r'$Z_{H}$ [dBZ]', fontsize=14)

	# ZDR
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_G, color='r', label='Graupel', ls='-')
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_S_SO, color='k', label='Snow Spheroid old m-D scheme', ls='-')
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_S_SN, color='b', label='Snow Spheroid new m-D scheme', ls='-')
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_S_SN_20, color='g', label='Snow Snowflake new m-D scheme n2/n1=20.0', ls='--')
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_S_SN_8_2, color='g', label='Snow Snowflake new m-D scheme n2/n1=8.2', ls='-')
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_S_HO, color='k', label='Snow HexCol old m-D scheme', ls='--')
	axes[1].plot(np.logspace(-5, -2, ntest), ZDR_S_HN, color='b', label='Snow HexCol new m-D scheme', ls='--')

	axes[1].set_ylabel(r'$Z_{DR}$ [dBZ]', fontsize=14)

	# KDP
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_G, color='r', label='Graupel', ls='-')
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_S_SO, color='k', label='Snow Spheroid old m-D scheme', ls='-')
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_S_SN, color='b', label='Snow Spheroid new m-D scheme', ls='-')
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_S_SN_20, color='g', label='Snow Snowflake new m-D scheme n2/n1=20.0', ls='--')
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_S_SN_8_2, color='g', label='Snow Snowflake new m-D scheme n2/n1=8.2', ls='-')
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_S_HO, color='k', label='Snow HexCol old m-D scheme', ls='--')
	axes[2].plot(np.logspace(-5, -2, ntest), KDP_S_HN, color='b', label='Snow HexCol new m-D scheme', ls='--')

	axes[2].set_ylabel(r'$K_{DP}$ [$deg \cdot km^{-1}$]', fontsize=14)
	axes[2].set_yscale('log')
	
	axes[2].set_xscale('log')
	axes[2].set_xlabel(r'Ice water content [$kg m-3$]', fontsize=14)
	
	axes[0].legend(loc='best', frameon=False, fontsize=8)

	plt.savefig('S_and_G', dpi=300)
	plt.close()
