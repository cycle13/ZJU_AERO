'''
Description: define some derived variables.
The derived variables are useful for radar operator applications
Author: Hejun Xie
Date: 2020-11-01 10:41:10
LastEditors: Hejun Xie
LastEditTime: 2020-11-02 10:40:41
'''

# Global imports
import numpy as np
import xarray as xr

# Local imports
from .. import DerivedVarWorkstation
from .interp_to_regular_grid import to_reg_grid
from .grapes_constant import grapes_constant, derived_var_unit

class GRAPESDerivedVar(DerivedVarWorkstation):
    def __init__(self, fhandle, regular_grids, raw_varnames, time_idx):
        '''
        Params:
            fhandle: The file handle offered by xarray dataset
                interface to netCDF4 file.
            regular_grids: The regular model grids.
            raw_varnames: The list variables names that can be accessed
                from netCDF4 format GRAPES file straightforwardly.
            time_idx: The time index we want to extract.
        Returns:
            A GRAPES derived variable work station.
        '''
        super(GRAPESDerivedVar, self).__init__(fhandle)
        self.regular_grids = regular_grids
        self.raw_varnames = raw_varnames
        self.time_idx = time_idx
    
    def _get_raw_var(self, varname):
        '''
        Description:
            raw.data: an numpy array.
            raw.dims: something like '('level', 'lat', 'lon')'.
        '''
        print("Raw var: {}".format(varname))
        raw = self.fhandle.data_vars[varname].isel(time=self.time_idx)
        return to_reg_grid(raw.data, varname, raw.dims, self.regular_grids)

    def _get_derived_var(self, varname):
        print("Derived var: {}".format(varname))
        raw_map = {'U':'u', 'V':'v', 'W':'w'}
        if varname in raw_map.keys():
            tmp_var = self.get_var( raw_map[varname] )
        elif varname == 'T':
            # unit: [K] = [K] * [-]
            tmp_var = self.get_var('th') * self.get_var('pi')
        else:
            return self._no_derived_var(varname)
        
        # assign the name and unit
        tmp_var.name = varname
        tmp_var.attrs['unit'] = derived_var_unit[varname]
        return tmp_var
