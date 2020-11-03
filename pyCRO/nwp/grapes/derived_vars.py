'''
Description: define some derived variables.
The derived variables are useful for radar operator applications
Author: Hejun Xie
Date: 2020-11-01 10:41:10
LastEditors: Hejun Xie
LastEditTime: 2020-11-02 22:40:54
'''

# Global imports
import numpy as np
import xarray as xr

# Local imports
from .. import DerivedVarWorkstation
from .interp_to_regular_grid import to_reg_grid
from .grapes_constant import grapes_constant as const
from .grapes_constant import derived_var_unit as unit
from .grapes_constant import derived_var_long_name as long_name

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
        raw_map = { 'U':'u', 'V':'v', 'W':'w', 
                    'QV':'Qv', 'QR':'Qr', 'QS':'Qs', 
                    'QG':'Qg', 'QC':'Qc', 'QI':'Qi'}
        mass_density_map = {'QV_v':'Qv', 'QR_v':'Qr', 'QS_v':'Qs', 
                            'QG_v':'Qg', 'QC_v':'Qc', 'QI_v':'Qi'}
        if varname in raw_map.keys():
            tmp_var = self.get_var(raw_map[varname])
        elif varname in mass_density_map.keys():
            # unit: [kg*m-3]
            tmp_var = self.get_var(mass_density_map[varname]) * self.RHO
        elif varname == 'T':
            # unit: [K] = [K] * [-]
            tmp_var = self.th * self.pi
        elif varname == 'P':
            # unit: [Pa] = [Pa] * [-]
            tmp_var = const['P0'] * (self.pi ** (1./const['kappa']))
        elif varname == 'Pw':
            # unit: [Pa] = [Pa] * [kg*kg-1] / [kg*kg-1]
            tmp_var = self.P * self.Qv / (self.Qv*(1.-const['Rdv'])+const['Rdv'])
        elif varname == 'Qhydro':
            # unit: [kg*kg-1]
            tmp_var = self.Qr + self.Qc + self.Qs + self.Qg + self.Qi
        elif varname == 'RHO':
            # unit: [kg*m-3]
            denominator = self.T * const['Rd'] * \
            (1. + self.Qv * (const['Rvd']-1.-self.Qhydro))
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

    vars_ok = True
    for varname in varname_list:
        if varname not in long_name.keys():
            vars_ok = False
            break
    
    return vars_ok
