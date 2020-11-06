'''
Description: define some derived variables.
The derived variables are useful for radar operator applications
Author: Hejun Xie
Date: 2020-11-05 19:29:35
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 21:33:38
'''

# Global imports
import numpy as np
import xarray as xr

# Local imports
from .. import DerivedVarWorkstation
from .interp_to_regular_grid import to_reg_grid
from .wrf_constant import wrf_constant as const
from .wrf_constant import derived_var_unit as unit
from .wrf_constant import derived_var_long_name as long_name
from .wrf_constant import raw_var_map

class WRFDerivedVar(DerivedVarWorkstation):
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
            A WRF derived variable work station.
        '''
        super(WRFDerivedVar, self).__init__(fhandle)
        self.regular_grids = regular_grids
        self.raw_varnames = raw_varnames
        self.time_idx = time_idx

    def _get_raw_var(self, varname):
        '''
        Description:
            raw.data: An numpy array.
            raw.dims: Something like '(bottom_top, south_north, west_east)'.
        '''
        print("Raw var: {}".format(varname))
        raw = self.fhandle.data_vars[raw_var_map[varname]].isel(Time=self.time_idx)
        return to_reg_grid(raw.data, varname, raw.dims, self.regular_grids)

    def _get_derived_var(self, varname):
        print("Derived var: {}".format(varname))

        mass_density_map = {'QV_v':'QV', 'QR_v':'QR', 'QS_v':'QS', 
                            'QG_v':'QG', 'QC_v':'QC', 'QI_v':'QI'}
        
        if varname in mass_density_map.keys():
            # unit: [kg*m-3]
            tmp_var = self.get_var(mass_density_map[varname]) * self.RHO
        elif varname == 'P':
            # [Pa]
            tmp_var = self.PP + self.PB
        elif varname == 'T':
            # [K]
            tmp_var = (self.PTP + self.PTB) * (self.P / self.P0) ** const['kappa']
        elif varname == 'Z':
            # [m] = [m2*s-2] / [m*s-2]
            tmp_var = (self.GPP + self.GPB) / const['g']
        elif varname == 'Pw':
            # unit: [Pa] = [Pa] * [kg*kg-1] / [kg*kg-1]
            tmp_var = self.P * self.QV / (self.QV*(1.-const['Rdv'])+const['Rdv'])
        elif varname == 'Qhydro':
            # unit: [kg*kg-1]
            tmp_var = self.QR + self.QC + self.QS + self.QG + self.QI
        elif varname == 'RHO':
            # unit: [kg*m-3]
            denominator = self.T * const['Rd'] * \
            (1. + self.QV * (const['Rvd']-1.-self.Qhydro))
            tmp_var = self.P / denominator
        elif varname == 'N':
            # unit: [-] * 1E-6
            tmp_var = (77.6/self.T)*(0.01*self.P+4810*(0.01*self.Pw)/self.T)
        else:
            return self._no_derived_var(varname)

        # assign the name and unit
        tmp_var.name = varname
        tmp_var.attrs['unit'] = unit[varname]
        tmp_var.attrs['long_name'] = long_name[varname]
        return tmp_var


def check_if_variables_in_file(varname_list):
    '''
    Description:
        Check if all the variable in varname_list
        are avaiable in this model file.
    Params:
        varname_list: A list of variable names.
    Returns:
        A bool flag indicating if all the variables 
        are available in the model file.
    '''
    pass
