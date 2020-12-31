'''
Description: Inspect the uncertainty of polarmitric radar variables (Zdr, Kdp, etc.)
of a given mass mixing ratio, temperature and elevation angle.
The uncertainty is caused by particle orientation and aspect ratio parameters.
Author: Hejun Xie
Date: 2020-12-31 10:33:16
LastEditors: Hejun Xie
LastEditTime: 2020-12-31 19:33:16
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import xarray as xr
import numpy as np

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

EXP_FILE = '../pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB_exp/t253.0e1.0.nc'


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
        print(self._ds)
    
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

if __name__ == "__main__":

    # get global constants
    createConfig('../test/option_files/wsm6_test.yml')
    from ZJU_AERO.config.cfg import CONFIG
    constants.update()
    
    db = exp_lookup_table_xarray(EXP_FILE)
    hydro_istc = Snow('1mom')
    hydro_name = 'S'

    # set profile
    ntest = 100
    e = np.ones((ntest), dtype='float32') * 1.0 # [deg]
    dic_values = dict()
    dic_values['T'] = np.ones((ntest), dtype='float32') * 253. # [K]
    dic_values['QS_v'] = np.logspace(-5, -2, ntest) # [kg m-3] (1E-5, 1E-2) [g m-3] (1E-2, 1E+1) 

    test_rad = Radial(dic_values, elev_profile=e,
    mask=None, lats_profile=None, lons_profile=None,
    dist_ground_profile=None, heights_profile=None)


    std_a_list, std_b_list = db.get_axis_value('std_a'), db.get_axis_value('std_b')
    lambda_a_list, lambda_b_list = db.get_axis_value('lambda_a'), db.get_axis_value('lambda_b') 
    mu_a_list, mu_b_list = db.get_axis_value('mu_a'), db.get_axis_value('mu_b') 

    # params = dict(std_a=40.0, std_b=-0.077, lambda_a=8.42, lambda_b=-0.57, mu_a=0.053, mu_b=0.79)
    params = dict(std_a=40.0, std_b=-0.08, lambda_a=8.0, lambda_b=-0.5, mu_a=0.055, mu_b=0.8)
    db.set_params(**params)

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
    


    # test
    # Dmax = db.get_axis_value('Dmax')
    # print(Dmax)
    # interped = db.lookup_line(e=[1.5,2.5], t=[213, 223])
    # print(interped.shape)
    # print(interped[0,:,0])
    # db.close()
    
    # with xr.open_dataset(EXP_FILE, engine='h5netcdf') as exp:
    #     params['method'] = 'nearest'
    #     exp_selcted = exp.sel(**params)
    #     p11 = exp_selcted['p11_bw'].data
    #     print(p11)

