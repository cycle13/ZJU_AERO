'''
Description: wrf_io.py read WRF model variables in netCDF4 format.
This module can read in derived variables xarray dataset, as defined in derived_vars.py.
Author: Hejun Xie
Date: 2020-11-05 19:28:25
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 21:34:18
'''

# Global imports
import xarray as xr
import numpy as np
import datetime as dt

# Local imports
from .derived_vars import WRFDerivedVar
from .wrf_constant import raw_var_map


def get_wrf_variables(data_file_list, varnames_list, ex_datetime):
    '''
    Params:
        data_file_list: A group of filename containing the datetime to be extracted.
        varname_list: A group of variable names we want to extract.
        ex_datetime: The datetime of the variable we want to extract.
    Returns:
        An xarray dataset containing all the variables in varname_list,
        With coordinates as geopotential height, latitude, longitude.
    '''

    for data_file in data_file_list:
        with xr.open_dataset(data_file) as ds:

            '''
            Dataset Example:
            1. ds.dims: {'Time': 73, 
            'south_north': 199, 'west_east': 199, 'bottom_top': 32, 
            'south_north_stag': 200, 'west_east_stag': 200, 'bottom_top_stag': 33, 
            'soil_layers_stag': 4}
            2. list(ds.coords.keys()): ['XLAT', 'XLONG', 'XTIME', 'XLAT_U', 'XLONG_U', 'XLAT_V', 'XLONG_V']
            3. list(ds.data_vars.keys()): ['Times', 'U', 'V', 'W', 'PH', 'PHB', 'T', 'P', 'PB' ...]
            4. ds.attrs: {'CEN_LAT': 40.000008, 'CEN_LON': 116.0, 'TRUELAT1': 30.0, 'TRUELAT2': 60.0, 
            'MOAD_CEN_LAT': 40.000008, 'STAND_LON': 116.0, 'POLE_LAT': 90.0, 'POLE_LON': 0.0, 
            'WEST-EAST_GRID_DIMENSION': 200, 'SOUTH-NORTH_GRID_DIMENSION': 200, 'BOTTOM-TOP_GRID_DIMENSION': 33, 
            'DX': 3000.0, 'DY': 3000.0, 'MAP_PROJ_CHAR': 'Lambert Conformal', ...}
            '''
            
            times_str = ds.data_vars['Times'].data[:]
            times_dt = [dt.datetime.strptime(time_str.decode('UTF-8'), '%Y-%m-%d_%H:%M:%S') for time_str in times_str]
            
            if ex_datetime in times_dt:
                found_file = data_file
                found_time_idx = times_dt.index(ex_datetime)
                print('Successfully found queried time in file:{}, time_index={}'.format(found_file, found_time_idx))
                break
    
    if 'found_file' not in locals().keys():
        raise ValueError("Queried time not available in this set of data files")

    output_ds = _get_wrf_variables(found_file, varnames_list, found_time_idx)
    output_ds.attrs['time'] = str(ex_datetime)

    return output_ds

def _get_wrf_variables(data_file, varname_list, time_idx):
    '''
    Params:
        data_file: A filename containing the queried data.
        varname_list: A group of variable names we want to extract.
        time_idx: The index of queried time in the datafile.
    Returns:
        An xarray dataset containing all the variables in varname_list,
        With coordinates as geopotential height, latitude, longitude.
    '''
    
    with xr.open_dataset(data_file) as ds: 

        # 1. Get regular coords and raw_varnames to generate the DerivedVarWorkStation
        regular_coords = dict()
        regular_coords['latitude'] = ds.coords['XLAT'].isel(Time=time_idx)
        regular_coords['longitude'] = ds.coords['XLONG'].isel(Time=time_idx)

        # here we use white namelist for raw variables as registered in raw_var_map
        raw_varnames = []
        file_wrf_names = list(ds.data_vars.keys())
        for raw_name, wrf_name in raw_var_map.items():
            if wrf_name in file_wrf_names:
                raw_varnames.append(raw_name)

        # 2.1 Get the derived vars
        dv = WRFDerivedVar(ds, regular_coords, raw_varnames, time_idx)
        prepare_ds = dict()
        for varname in varname_list:
            prepare_ds[varname] = dv.get_var(varname)

        # 2.2 Get the coords
        z = dv.get_var('Z') # get z-levels on regular grids 
        tp = dv.get_var('TP') # get topograph on regular grids
        dv.close()

    # 3. Make the xarray data set
    output_ds = xr.Dataset(prepare_ds)
    output_ds.coords['z-levels'] = (("bottom_top", "south_north", "west_east"), z.data)
    output_ds.coords['topograph'] = (("south_north", "west_east"), tp.data)

    # 4. Add the projection information (now only latlon)
    proj_info = {'TRUELAT1':ds.attrs['TRUELAT1'], 'TRUELAT2':ds.attrs['TRUELAT2'],
                    'MOAD_CEN_LAT':ds.attrs['MOAD_CEN_LAT'], 'STAND_LON':ds.attrs['STAND_LON'],
                    'CEN_LAT':ds.attrs['CEN_LAT'], 'CEN_LON':ds.attrs['CEN_LON'],
                    'DX':ds.attrs['DX'], 'DY':ds.attrs['DY'],
                    'nI':ds.dims['west_east'], 'nJ':ds.dims['south_north'],
                    'PROJ_TYPE':ds.attrs['MAP_PROJ_CHAR']}
    output_ds.attrs.update(proj_info)
    if 'N' in output_ds.data_vars.keys():
        output_ds.data_vars['N'].attrs.update(proj_info)

    return output_ds
