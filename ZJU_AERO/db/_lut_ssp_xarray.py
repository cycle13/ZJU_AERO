'''
Description: submodule to deal with the new xarray 
ssp lookup tables, iitm_masc.
Author: Hejun Xie
Date: 2020-10-12 11:30:13
LastEditors: Hejun Xie
LastEditTime: 2021-03-06 19:42:50
'''

# Global import
import xarray as xr
import numpy as np


def load_lut_xarray(filename):
    '''
    Loads an instance of the Loookup_table class by loading the LevelB xarray database,
    previously saved to the drive with the sz_lut function in compute_lut_sz.py
    Args:
        filename: the complete filename (with path), indicating where the
            lookup table is stored

    Returns:
        lut: the lookup table as an instance of the class Lookup_table Class
        (see below)
    '''

    lut = Lookup_table_xarray(filename)
    return lut

class Lookup_table_xarray:
    '''
    The Lookup_table class used to store scattering properties of hydrometeors
    and perform queries of the scattering properties, this class assumes all
    stored data to be defined on regular grids, actually a wrapper for xarray
    '''
    
    def __init__(self, filename):
        self.type = 'xarray'
        self._ds = xr.open_dataset(filename, engine='h5netcdf')
    
    def close(self):
        self._ds.close()
    
    def get_axis_value(self, axis_name):
        return self._ds.coords[axis_name].values
    
    def lookup_line(self, **kwargs):
        '''
        Query the xarray levelB database with interp method.
        get the interpolated and flattened data table
        '''
        interped_data = {}
        # Adapt the kwargs for iitm xarray database
        iitm_kwargs = {}
        for key, value in kwargs.items():
            if key == 'e':
                iitm_kwargs['elevation'] = xr.DataArray(value, dims='line')
            elif key == 't':
                iitm_kwargs['temperature'] = xr.DataArray(value, dims='line')
            else:
                iitm_kwargs[key] = xr.DataArray(value, dims='line')
        # interp the xarray data into queried hyperplane
        for var in self._ds.data_vars.keys():
            interped_data[var] = self._ds[var].interp(**iitm_kwargs).values

        flattened_data = self._flatten_SZ_matrix(interped_data)
        return flattened_data
    
    def _flatten_SZ_matrix(self, interped_data):
        '''
        Flatten the SZ matrix interpolated from xarray levelB matrix
        return a table with dimension (n_valid_gates, nbins_D, 12)
        '''
        shape_single = list(interped_data['p11_bw'].shape)
        shape_single.append(12)
        flattened_shape = tuple(shape_single)
        flattened_data = np.empty(flattened_shape, dtype='float32')
        
        flattened_map_real = {  'p11_bw':0,     'p12_bw':1,     'p13_bw':-1,    'p14_bw':-1,
                                'p21_bw':2,     'p22_bw':3,     'p23_bw':-1,    'p24_bw':-1,
                                'p31_bw':-1,    'p32_bw':-1,    'p33_bw':4,     'p34_bw':5,
                                'p41_bw':-1,    'p42_bw':-1,    'p43_bw':6,     'p44_bw':7}
        flattened_map_cplx = {  's11_fw':8,     's12_fw':-1,    's21_fw':-1,    's22_fw':10}

        # flatten the data
        for key,value in flattened_map_real.items():
            if value != -1:
                flattened_data[...,value] = interped_data[key]

        for key,value in flattened_map_cplx.items():
            if value != -1:
                flattened_data[...,value] = interped_data[key].real
                flattened_data[...,value+1] = interped_data[key].imag
        
        return flattened_data

if __name__ == "__main__":
    # lut = load_lut_xarray('/home/xhj/wkspcs/Radar-Operator/cosmo_pol/pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc')
    lut = load_lut_xarray('/home/xhj/wkspcs/Radar-Operator/cosmo_pol/pathos/lut/tm_masc/lut_SZ_R_35_0_1mom_LevelB.nc')
    
    Dmax = lut.get_axis_value('Dmax')
    print(Dmax)

    # interped = lut.lookup_line(e=[1.5,2.5], t=[213, 223])
    interped = lut.lookup_line(e=[1.5,2.5], t=[283, 293])
    print(interped.shape)
    print(interped)

    lut.close()
