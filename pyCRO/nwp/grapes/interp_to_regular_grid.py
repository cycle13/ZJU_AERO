'''
Description: interp the grapes model data to regular grid.
By regular grid, we mean the full grid of vertical grid (mass grid),
and the regular grid in horizontal C-grid.
Author: Hejun Xie
Date: 2020-11-01 18:45:30
LastEditors: Hejun Xie
LastEditTime: 2020-11-01 23:43:37
'''

import numpy as np
import xarray as xr

def to_reg_grid(data, varname, dim_names, regular_coords):
    '''
    Description:
        All the netCDF4 data will be passed into this function to determine if
        the numpy data have to be interpolated into regular model grid.
        model level variables 
        (a). u, v: on half vertical model levels, but not available on top and bottom levels. (nz-1)
        (b). pi: on half vertical model levels, available on top and bottom levels. (nz+1)
        (c). other mass variables: available on full vertical model levels. (nz)
    Params:
        data: The numpy array to be interpolated into regular grid,
            i.e., horizontal regular grid or vertical full grid.
        varname: The name of the variable, will be used to determine which
            interpolation rule to follow.
        dim_names: The name of dimensions of the numpy array, respectively.
            It will be used to determine which dimension will be interpolated. 
        regular coords: The regular grid coordinates to prepare the new
            dataArray defined on regular grids.
    Returns:
        An interpolated xarray dataArray defined on regular grids

    TODO: implement the dim_names
    '''
    # 1. vertical interpolation
    if varname in ['u', 'v']:
        data = _to_reg_grid_ver_interp_uv(data, dim_names)
    elif varname in ['pi']:
        data = _to_reg_grid_ver_interp_pi(data, dim_names)
    # 2. horizontal interpolation
    if varname in ['u']:
        data = _to_reg_grid_hor_interp_u(data, dim_names)
    elif varname in ['v']:
        data = _to_reg_grid_hor_interp_v(data, dim_names)
    # 3. generate the xarray dataArray
    out_xr = xr.DataArray(data, 
    coords=[("levels", regular_coords['levels']),
            ("latitude", regular_coords['latitude']),
            ("longitude", regular_coords['longitude'])]
    )

    return out_xr
    

def _to_reg_grid_hor_interp_u(data, dim_names):
    '''
    Assume the data dimension as (levels, latitude, longitude)
    levels: bottom-top
    latitude: west-east
    longitude: south-north
    '''
    tmp_data = 0.5 * (data[:,:,1:] + data[:,:,:-1])
    return np.pad(tmp_data, ((0, 0), (0, 0), (0, 1)), 'edge')

def _to_reg_grid_hor_interp_v(data, dim_names):
    '''
    Assume the data dimension as (levels, latitude, longitude)
    levels: bottom-top
    latitude: west-east
    longitude: south-north
    '''
    tmp_data = 0.5 * (data[:,1:,:] + data[:,:-1,:])
    return np.pad(tmp_data, ((0, 0), (0, 1), (0, 0)), 'edge')

def _to_reg_grid_ver_interp_uv(data, dim_names):
    '''
    Assume the data dimension as (levels, latitude, longitude)
    levels: bottom-top
    latitude: west-east
    longitude: south-north
    '''
    # may be surface wind should be zero
    tmp_data = 0.5 * (data[1:,:,:] + data[:-1,:,:])
    return np.pad(tmp_data, ((1, 1), (0, 0), (0, 0)), 'edge')

def _to_reg_grid_ver_interp_pi(data, dim_names):
    '''
    Assume the data dimension as (levels, latitude, longitude)
    levels: bottom-top
    latitude: west-east
    longitude: south-north
    '''
    # zm = sqrt(z1*z2) 
    # 1. isothermal layer assumption
    # 2. Half layer is evenly spaced by full layer geometrically
    tmp_data = np.sqrt(data[1:,:,:] * data[:-1,:,:])
    return tmp_data
