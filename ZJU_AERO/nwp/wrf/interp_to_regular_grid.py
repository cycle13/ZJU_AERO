'''
Description: interp the wrf model data to regular grid.
By regular grid, we mean ...
Author: Hejun Xie
Date: 2020-11-05 19:30:49
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 21:03:42
'''

import numpy as np
import xarray as xr

def to_reg_grid(data, varname, dim_names, regular_coords):
    '''
    Description:
        All the netCDF4 data will be passed into this function to determine if
        the numpy data have to be interpolated into regular model grid.
        model level variables.
        Some irregular dimensions 
        (a). bottom_top_stag: bottom_top + 1 
        (b). south_north_stag: south_north + 1
        (c). west_east_stag: west_east + 1
    Params:
        data: The numpy array to be interpolated into regular grid,
            i.e., horizontal regular grid or vertical full grid.
        varname: The name of the variable, will be used to determine which
            interpolation rule to follow.
        dim_names: The name of dimensions of the numpy array, respectively.
            It will be used to determine which dimension will be interpolated.
            Something like '(bottom_top, south_north, west_east)'
        regular coords: The regular grid coordinates to prepare the new
            dataArray defined on regular grids.
    Returns:
        An interpolated xarray dataArray defined on regular grids.

    '''

    stag_dimensions = ['bottom_top_stag', 'south_north_stag', 'west_east_stag']
    avg_dimension_map = {
        'bottom_top_stag': 'bottom_top',
        'south_north_stag': 'south_north',
        'west_east_stag': 'west_east',
    } 

    dim_regular = [] 
    for dim_idx, dim_name in enumerate(dim_names):
        if dim_name in stag_dimensions:
            data = _get_stag_avg(dim_idx, data)
            dim_regular.append(avg_dimension_map[dim_name])
        else:
            dim_regular.append(dim_name)
    
    out_xr = xr.DataArray(data, dims=dim_regular)

    if 'south_north' in dim_regular and 'west_east' in dim_regular:
        out_xr.coords['latitude']   = (('south_north', 'west_east'), regular_coords['latitude'].data)
        out_xr.coords['longitude']  = (('south_north', 'west_east'), regular_coords['longitude'].data)

    return out_xr


def _get_stag_avg(I, data):
    '''
    Description: Make stag average for the dimension I of numpy-Array data
    Params:
        I: The dimension index to make stag slices.
        data: The numpy data to make stag average.
    Returns:
        data: The numpy data average after stag average
    '''

    def ___simple_avg(d1,d2):
        return 0.5*(d1+d2)

    N = len(data.shape)
    slice1, slice2 = _make_stag_slice(I, N)
    return ___simple_avg(data[slice1], data[slice2])

def _make_stag_slice(I, N):
    '''
    Description: Make stag slice for the dimension I of N-Dimension numpy-Array
    Params:
        I: The dimension index to make stag slices.
        N: Total dimensions numbers of the numpy array.
    Returns:
        stag_slice1: slice 1 , tuple of slice object.
        stag_slice2: slice 2 , tuple of slice object.
    Example:
        (bottom_top, south_north, west_east_stag) --> I = 2, N = 3
        stag_slice1 = [:,:,:-1]
        stag_slice2 = [:,:,1:]
    '''

    ele_slice1 = slice(None, -1)
    ele_slice2 = slice(1, None)
    ele_all = slice(None)
    
    stag_slice1, stag_slice2 = [], []
    
    for i in range(N):
        if i == I:
            stag_slice1.append(ele_slice1)
            stag_slice2.append(ele_slice2)
        else:
            stag_slice1.append(ele_all)
            stag_slice2.append(ele_all)
    
    return tuple(stag_slice1), tuple(stag_slice2)
