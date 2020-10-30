# -*- coding: utf-8 -*- 
'''
@Author: wanghao
@Date: 2019-01-22 15:31:50
LastEditors: Hejun Xie
LastEditTime: 2020-10-30 20:36:26
@Description: 转换数据在一个文件中的数据
'''

# Global import
import sys
import os
from netCDF4 import Dataset
import numpy as np
import time
import datetime
import pandas as pd
import re
import glob
import copy
import datetime as dt

# Local import
from .CTLReader import CTLReader

def transf2nc_specific_ctl(ctlname, nc_out, var='all', v_grid_type='ml'):
    '''
    Params:
        ctlname: The ctl filename, should be non-general ctl file.
        nc_out: The output netCDF4 file to be saved as.
        var: The variavbles you want to read, set it as 'all' if you want to read all variables.   
        v_grid_type: 'ml' or 'pl':
            1. 'ml', model level variables 
                (a). u, v: on half vertical model levels, but not available on top and bottom levels. (nz-1)
                (b). pi: on half vertical model levels, available on top and bottom levels. (nz+1)
                (c). other mass variables: available on full vertical model levels. (nz)
            2. 'pl', pressure level variables
                all variables on full vertical model pressure levels. (nz)
    Returns:
        No returns, but save the data in netCDF4 format.
    '''

    if v_grid_type not in ['ml', 'pl']:
        raise KeyError('Invalid v_grid_type:{}'.format(v_grid_type))
    
    data = CTLReader(ctlname)
    times_tmp = [itime.strftime("%Y%m%d%H") for itime in data.variables['time']]

    # 提取单位信息
    cur_path = os.path.dirname(os.path.realpath(__file__))
    vartable_path = cur_path + os.sep + 'GRAPES_VAR_TABLE.txt'
    df = pd.read_table(vartable_path, sep=',')
    var_unit = dict()
    for index, row in df.iterrows():
        ivar, unit = row
        var_unit[ivar] = unit

    # 提取维度信息
    nt = data.dimensions['time']
    nz = data.dimensions['levels']
    ny = data.dimensions['latitude']
    nx = data.dimensions['longitude']
    if v_grid_type == 'ml':
        nz -= 1 # set nz as the mass level

    ncfile = Dataset(nc_out,'w')

    ncfile.createDimension('ntime',nt)
    ncfile.createDimension('nlevel',nz)
    ncfile.createDimension('nlat',ny)
    ncfile.createDimension('nlon',nx)
    if v_grid_type == 'ml':
        ncfile.createDimension('nlevel_uv',nz-1)
        ncfile.createDimension('nlevel_pi',nz+1)

    times     = ncfile.createVariable('times',np.int,('ntime'))
    levels    = ncfile.createVariable('levels',np.float32,('nlevel'))
    latitude  = ncfile.createVariable('latitude',np.float32,('nlat'))
    longitude = ncfile.createVariable('longitude',np.float32,('nlon'))
    if v_grid_type == 'ml':
        levels_uv   = ncfile.createVariable('levels_uv',np.float32,('nlevel_uv'))
        levels_pi   = ncfile.createVariable('levels_pi',np.float32,('nlevel_pi'))
    
    times[:]     = [int(itime) for itime in times_tmp] 
    latitude[:]  = data.variables['latitude'][:]
    longitude[:] = data.variables['longitude'][:]
    
    if v_grid_type == 'ml':
        '''
        0.5 is the bottom levele
        0.5 + nz is the top level
        '''
        levels_uv[:]    = np.arange(nz-1) + 1.5
        levels_pi[:]    = np.arange(nz+1) + 0.5 
        levels[:]       = np.arange(nz) + 1.0
    elif v_grid_type == 'pl':
        levels[:]       = data.variables['levels'][:]

    # Read all variables added by Xie 2020/10/30
    if var == 'all':
        var = data.varname

    for ivar in var:
        if data.variables[ivar].dimensions['levels'] == 1:
            locals()[ivar] = ncfile.createVariable(ivar,np.float32,('ntime','nlat','nlon'))
            locals()[ivar][:] = data.variables[ivar][:]
            locals()[ivar].units = var_unit[ivar]
            locals()[ivar].long_name = data.variables[ivar].attributes['long_name']
        else:
            if ivar in ['u', 'v']:
                locals()[ivar] = ncfile.createVariable(ivar,np.float32,('ntime','nlevel_uv','nlat','nlon'))
            elif ivar in ['pi']:
                locals()[ivar] = ncfile.createVariable(ivar,np.float32,('ntime','nlevel_pi','nlat','nlon'))
            else:
                locals()[ivar] = ncfile.createVariable(ivar,np.float32,('ntime','nlevel','nlat','nlon'))
            locals()[ivar][:] = data.variables[ivar][:]
            locals()[ivar].units = var_unit[ivar]
            locals()[ivar].long_name = data.variables[ivar].attributes['long_name']

    # Variable Attributes
    latitude.units  = 'degree_north'
    longitude.units = 'degree_east'
    if v_grid_type == 'pl':
        levels.units = 'hPa'
    elif v_grid_type == 'ml':
        levels.units = 'level_number'
        levels_uv.units = 'level_number'
        levels_pi.units = 'level_number'
    times.units     = 'days since 2000-01-01 00:00:00'
    times.calendar  = 'gregorian'
    times.incr      = '{}'.format(data.crement['time'].seconds/3600)

    # Global Attributes
    ncfile.description = 'Transf postvar data to NC'
    ncfile.history = 'Created ' + time.ctime(time.time())
    ncfile.source = 'netCDF4 python module tutorial'

    ncfile.close()

def transf2nc_general_ctl(ctlname, var='all', v_grid_type='ml'):
    '''
    Params:
        ctlname: The general ctl filename
        var: The variavbles you want to read, set it as 'all' if you want to read all variables.
        v_grid_types: See transf2nc_specific_ctl.
    Returns:
        No returns, but save the data in netCDF4 format.
    '''
    specfic_ctl_filenames, nc_outs = _gen_specific_ctl(ctlname)
    for specfic_ctl_filename, nc_out in zip(specfic_ctl_filenames, nc_outs):
        print('{}-->{}'.format(specfic_ctl_filename, nc_out))
        transf2nc_specific_ctl(specfic_ctl_filename, nc_out, var, v_grid_type)
        os.system('rm {}'.format(specfic_ctl_filename))

def _gen_specific_ctl(general_ctl):
    '''
    Params:
        general_ctl: The general ctl filename.
    Returns:
        A group of specfic ctl filenames, and the output netCDF4 filenames 
        generated by this function.
    '''

    # 1. Read in the general ctl file
    with open(general_ctl, 'r') as fin:
       ctl_content = fin.read()

    # 2. Get some information from general ctl file
    p = re.compile("{}\s+(.*)".format('dset'))
    m = p.search(ctl_content)
    if m.group(1)[0] == '^':
        path = os.path.dirname(general_ctl)
        filename = path + os.sep + m.group(1)[1:]
    else:
        filename = m.group(1)[:]
    
    p = re.compile("{}\s+(\d+)\s+linear\s+(\w+)\s+(\w+)\s".format('tdef'))
    m = p.search(ctl_content)
    nfile = int(m.group(1))
    init_dt = dt.datetime.strptime(m.group(2), '%Hz%d%b%Y')

    # 3. Get all the specific data files 
    # TODO: A crude match, maybe we should fix it in the future
    match = filename.split('%f')[0] + '?'*int(filename.split('%f')[1])
    datfiles = glob.glob(match)
    if len(datfiles) != nfile:
        raise ValueError('Find too many or too little specfic data files {} as indicated by general ctl file {}'.format(len(datfiles), nfile))

    # 4. Do replaces and ouput the specfic ctl file 
    nc_outs = []
    specfic_ctl_filenames = []

    for datfile in datfiles:
        nc_outs.append(datfile+'.nc')
        specfic_ctl_filenames.append(datfile+'.ctl')

        # do some replace jobs
        temp_ctl_content = copy.copy(ctl_content)
        p = re.compile(r"(dset\s+\^*).*")
        m1 = p.sub('\g<1>{}'.format(datfile.split(os.sep)[-1]), temp_ctl_content)

        # TODO: A crude time string, maybe we should fix it in the future
        dat_dt = init_dt + dt.timedelta(hours=int(datfile.split('_')[-1]))
        timestr = dat_dt.strftime('%Hz%d%b%Y').upper()

        p = re.compile(r"(tdef\s+)\d+(\s+linear\s+)\w+(\s+\w+\s)")
        m2 = p.sub('\g<1>{}\g<2>{}\g<3>'.format(1, timestr), m1)

        with open(datfile+'.ctl', 'w') as fout:
            fout.write(m2)
    
    return specfic_ctl_filenames, nc_outs
