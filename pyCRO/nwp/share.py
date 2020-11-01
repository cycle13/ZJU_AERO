'''
Description: define the some shared utilities for GRAPES
and WRF radar operator interface.
Author: Hejun Xie
Date: 2020-11-01 19:51:00
LastEditors: Hejun Xie
LastEditTime: 2020-11-01 23:51:38
'''


class DerivedVarWorkstation(object):
    '''
    Description:
         A base class defining how to get the derived var 
         for radar-operator-nwp-model nterfaces.
    '''
    def __init__(self, fhandle):
        self.fhandle = fhandle
        # A forced attribute, must be specified in specfic classes 
        self.raw_varnames = list() 
        self.known_vars = dict()

    def get_var(self, varname):
        
        var = self._check_register(varname)

        if var is not None:
            return var
        
        if varname in self.raw_varnames:
            var = self._get_raw_var(varname)
        else:
            var = self._get_derived_var(varname)
        
        self._register_var(varname, var)

        return var

    def close(self):
        for var in self.known_vars.values():
            del var
        del self.known_vars

    def _register_var(self, varname, data):
        if varname not in self.known_vars.keys():
            self.known_vars[varname] = data

    def _check_register(self, varname):
        if varname in self.known_vars.keys():
            return self.known_vars[varname]

    def _no_derived_var(self, varname):
        raise ValueError('Could not compute derived variable {}, please specify a valid variable name'.format(varname))

    def _get_raw_var(self, varname):
        pass

    def _get_derived_var(self, varname):
        pass

