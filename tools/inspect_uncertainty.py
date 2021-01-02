'''
Description: Inspect the uncertainty of polarmitric radar variables (Zdr, Kdp, etc.)
of a given mass mixing ratio, temperature and elevation angle.
The uncertainty is caused by particle orientation and aspect ratio parameters.
Author: Hejun Xie
Date: 2020-12-31 10:33:16
LastEditors: Hejun Xie
LastEditTime: 2021-01-02 11:05:16
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import xarray as xr
import numpy as np
import itertools

# plot settings
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'

# colormap interpolation
from matplotlib import cm
from scipy import interpolate

# ZJU AERO
from ZJU_AERO.core._rdop_scatter import one_rad_one_hydro, get_pol_from_sz 
from ZJU_AERO.interp import Radial
from ZJU_AERO.db._lut_ssp_xarray import Lookup_table_xarray
from ZJU_AERO.hydro.hydrometeor import Snow
from ZJU_AERO.config import createConfig
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.utils import DATAdecorator

# EXP_FILE = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB_exp/t253.0e1.0.nc'
EXP_FILE = '../pathos/lut/tm_masc/lut_SZ_S_9_41_1mom_LevelB_exp/t253.0e1.0.nc'
nQ = 100
Qs = np.logspace(-5, -2, nQ) # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1) 


def get_colormap(cmap_name):
    cm_org = cm.__getattribute__(cmap_name)
    colormap = np.empty((256, 3), dtype='float64')

    for i in range(256):
        colormap[i,:] = np.array(cm_org(i))[0:3]
    return colormap

def interp_colormap(colormap, npts):
    x = np.arange(256)
    xnew = np.arange(npts) / (npts-1) * 255.
    
    interped_colormap = np.empty((npts, 3), dtype='float64')

    for iband in range(3):
        y = colormap[:,iband]
        f = interpolate.interp1d(x, y)
        ynew = f(xnew)
        interped_colormap[:,iband] = ynew

    return interped_colormap

class exp_lookup_table_xarray(Lookup_table_xarray):
    '''
    The Lookup_table class used to store scattering properties of hydrometeors
    and perform queries of the scattering properties, this class assumes all
    stored data to be defined on regular grids, actually a wrapper for xarray
    '''

    def __init__(self, filename):
        super(exp_lookup_table_xarray, self).__init__(filename)
        # print(self._ds)
    
    def set_params(self, **kwargs):
        self.params = kwargs
    
    def lookup_line(self, **kwargs):
        '''
        Query the xarray levelB database with interp method.
        get the interpolated and flattened data table
        '''
        interped_data = {}
        # Adapt the kwargs for iitm xarray database
        iitm_kwargs = self.params
        for key, value in kwargs.items():
            if key == 'e' or key == 't':
                nQ = len(value)
            else:
                iitm_kwargs[key] = xr.DataArray(value, dims='line')
        
        # interp the xarray data into queried hyperplane
        for var in self._ds.data_vars.keys():
            interped_data[var] = self._ds[var].interp(**iitm_kwargs).values
            # iitm_kwargs['method'] = 'nearest'
            # interped_data[var] = self._ds[var].sel(**iitm_kwargs).values

        flattened_data = self._flatten_SZ_matrix(interped_data)

        flattened_data = np.expand_dims(flattened_data, axis=0).repeat(nQ, axis=0)
        return flattened_data


# @DATAdecorator('./', True, './uctt_iitm.pkl')
@DATAdecorator('./', True, './uctt_tm.pkl')
def get_uctt():
    # get global constants
    createConfig('../test/option_files/wsm6_test.yml')
    from ZJU_AERO.config.cfg import CONFIG
    constants.update()
    
    db = exp_lookup_table_xarray(EXP_FILE)
    hydro_istc = Snow('1mom')
    hydro_name = 'S'

    # set profile
    e = np.ones((nQ), dtype='float32') * 1.0 # [deg]
    dic_values = dict()
    dic_values['T'] = np.ones((nQ), dtype='float32') * 253. # [K]
    dic_values['QS_v'] = Qs # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1) 

    test_rad = Radial(dic_values, elev_profile=e,
    mask=None, lats_profile=None, lons_profile=None,
    dist_ground_profile=None, heights_profile=None)

    # formulate the iteration
    # params_dict = {
    #     'std_a': db.get_axis_value('std_a')[0:2],
    #     'std_b': db.get_axis_value('std_b')[0:2],
    #     'lambda_a': db.get_axis_value('lambda_a')[2:4],
    #     'lambda_b': db.get_axis_value('lambda_b')[2:4],
    #     'mu_a': db.get_axis_value('mu_a')[0:2],
    #     'mu_b': db.get_axis_value('mu_b')[0:2],
    # }

    params_dict = {
        'std_a': db.get_axis_value('std_a'),
        'std_b': db.get_axis_value('std_b'),
        'lambda_a': db.get_axis_value('lambda_a'),
        'lambda_b': db.get_axis_value('lambda_b'),
        'mu_a': db.get_axis_value('mu_a'),
        'mu_b': db.get_axis_value('mu_b'),
    }

    params_shape = np.array([len(value) for value in params_dict.values()])

    iter_idx    = tuple([range(len(value)) for value in params_dict.values()])
    iter_keys   = params_dict.keys()
    key2ikey = dict()
    for ikey, key in enumerate(iter_keys):
        key2ikey[key] = ikey
    
    # formulate the numpy array
    
    out_shape = [nQ]
    out_shape.extend(params_shape)
    nsample = np.product(params_shape)
    out_zdr = np.empty(out_shape, dtype='float64')
    out_kdp = np.empty(out_shape, dtype='float64')
    out_zh = np.empty(out_shape, dtype='float64')

    # formulate the xarray
    coords = [("nQ", dic_values['QS_v'])]
    coords.extend(params_dict.items())

    out_ZDR = xr.DataArray(out_zdr, coords=coords)
    out_KDP = xr.DataArray(out_kdp, coords=coords)
    out_ZH = xr.DataArray(out_zh, coords=coords)

    # print(out_ZDR) 

    # start iteration
    for idx_list in itertools.product(*iter_idx):
        lambda_a = params_dict['lambda_a'][idx_list[key2ikey['lambda_a']]]
        lambda_b = params_dict['lambda_b'][idx_list[key2ikey['lambda_b']]]
        mu_a = params_dict['mu_a'][idx_list[key2ikey['mu_a']]]
        mu_b = params_dict['mu_b'][idx_list[key2ikey['mu_b']]]
        std_a = params_dict['std_a'][idx_list[key2ikey['std_a']]]
        std_b = params_dict['std_b'][idx_list[key2ikey['std_b']]]

        params = dict(std_a=std_a, std_b=std_b, 
        lambda_a=lambda_a, lambda_b=lambda_b, 
        mu_a=mu_a, mu_b=mu_b)

        message = ['{}={:>.3f}'.format(item[0], item[1]) for item in params.items()]
        print(' '.join(message))

        db.set_params(**params)

        '''
        sz_psd_integ: 
                dimension: (nQ, 12) 
                unit: Z[mm2 m-3]; S[mm m-3]
        '''
        valid_data, sz_psd_integ = one_rad_one_hydro(test_rad, hydro_name, hydro_istc,
        db, simulate_doppler=False, ngates=nQ)

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
        dbZH    = 10 * np.log10(ZH) # [dB]
        dbZDR   = 10 * np.log10(ZDR) # [dB]

        if np.isnan(dbZDR).all():
            print('NaN')

        out_ZDR.loc[params] = dbZDR
        out_ZH.loc[params] = dbZH
        out_KDP.loc[params] = KDP        
    
    db.close()
    return out_ZDR, out_KDP, out_ZH, nsample

if __name__ == "__main__":

    uctt_ZDR, uctt_KDP, uctt_ZH, nsample = get_uctt()
    uctt_ZDR = uctt_ZDR.data.reshape((nQ, nsample))
    uctt_KDP = uctt_KDP.data.reshape((nQ, nsample)) 
    uctt_ZH = uctt_ZH.data.reshape((nQ, nsample)) 
    print(uctt_ZDR.shape)

    plot_Q = Qs[0:100:9] * 1e3 # [kg/m^3] --> [g/m^3]

    '''
    1. ZDR
    '''
    fig, ax = plt.subplots(figsize=(12,6))
    
    ZDR_list = [uctt_ZDR[i,:] for i in range(0,100,9)]
    ZDR_list = [ZDR[~np.isnan(ZDR)] for ZDR in ZDR_list]
    inds = [ x+1 for x in range(len(plot_Q))]

    parts = ax.violinplot(ZDR_list, showmeans=False, showmedians=False, showextrema=False, 
    widths=0.7)
    
    colormap = get_colormap('viridis')
    mycolors = interp_colormap(colormap, len(plot_Q))
    for ipc, pc in enumerate(parts['bodies']):
        pc.set_facecolor(mycolors[ipc])
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    ax.set_title('Uncertainty of ZDR', fontsize=16)

    ax.set_xlabel(r'Snow concentration Qs [$g \cdot m^{-3}$]', fontsize=14)
    ax.set_ylabel(r'$Z_{DR}$ [dBZ]', fontsize=14)

    ax.set_xticks(inds)
    ax.set_xticklabels(['{:>.3f}'.format(Q) for Q in plot_Q])
    # ax.set_yscale('log')

    ZDR = np.array(ZDR_list, dtype='float32')
    quartile1, medians, quartile3 = np.percentile(ZDR, [25, 50, 75], axis=1)
    whiskers_min, whiskers_max = np.min(ZDR, axis=1), np.max(ZDR, axis=1)
    ax.scatter(inds, medians, marker='o', color='white', s=15, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=3.5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    plt.savefig('uncertainty_ZDR.png', dpi=300)
    plt.close()

    '''
    2. ZH
    '''
    fig, ax = plt.subplots(figsize=(12,6))
    
    ZH_list = [uctt_ZH[i,:] for i in range(0,100,9)]
    ZH_list = [ZH[~np.isnan(ZH)] for ZH in ZH_list]
    inds = [ x+1 for x in range(len(plot_Q))]

    parts = ax.violinplot(ZH_list, showmeans=False, showmedians=False, showextrema=False, 
    widths=0.7)
    
    colormap = get_colormap('viridis')
    mycolors = interp_colormap(colormap, len(plot_Q))
    for ipc, pc in enumerate(parts['bodies']):
        pc.set_facecolor(mycolors[ipc])
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    ax.set_title('Uncertainty of ZH', fontsize=16)

    ax.set_xlabel(r'Snow concentration Qs [$g \cdot m^{-3}$]', fontsize=14)
    ax.set_ylabel(r'$Z_{H}$ [dBZ]', fontsize=14)

    ax.set_xticks(inds)
    ax.set_xticklabels(['{:>.3f}'.format(Q) for Q in plot_Q])
    # ax.set_yscale('log')

    ZH = np.array(ZH_list, dtype='float32')
    quartile1, medians, quartile3 = np.percentile(ZH, [25, 50, 75], axis=1)
    whiskers_min, whiskers_max = np.min(ZH, axis=1), np.max(ZH, axis=1)
    ax.scatter(inds, medians, marker='o', color='white', s=5, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=3.5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    plt.savefig('uncertainty_ZH.png', dpi=300)
    plt.close()

    '''
    3. KDP
    '''
    fig, ax = plt.subplots(figsize=(12,6))
    
    KDP_list = [uctt_KDP[i,:] for i in range(0,100,9)]
    KDP_list = [KDP[~np.isnan(KDP)] for KDP in KDP_list]
    inds = [ x+1 for x in range(len(plot_Q))]

    parts = ax.violinplot(KDP_list, showmeans=False, showmedians=False, showextrema=False, 
    widths=0.7)
    
    colormap = get_colormap('viridis')
    mycolors = interp_colormap(colormap, len(plot_Q))
    for ipc, pc in enumerate(parts['bodies']):
        pc.set_facecolor(mycolors[ipc])
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    ax.set_title('Uncertainty of KDP', fontsize=16)

    ax.set_xlabel(r'Snow concentration Qs [$g \cdot m^{-3}$]', fontsize=14)
    ax.set_ylabel(r'$K_{DP}$ [$deg \cdot km^{-1}$]', fontsize=14)

    ax.set_xticks(inds)
    ax.set_xticklabels(['{:>.3f}'.format(Q) for Q in plot_Q])
    ax.set_yscale('log')
    ax.set_ylim([5e-4, 1e2])

    KDP = np.array(KDP_list, dtype='float32')
    quartile1, medians, quartile3 = np.percentile(KDP, [25, 50, 75], axis=1)
    whiskers_min, whiskers_max = np.min(KDP, axis=1), np.max(KDP, axis=1)
    ax.scatter(inds, medians, marker='o', color='white', s=15, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=3.5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    plt.savefig('uncertainty_KDP.png', dpi=300)
    plt.close()
