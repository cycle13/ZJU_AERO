'''
Description: computing lookup table using pytmatrix
Author: Hejun Xie
Date: 2020-12-06 10:21:24
LastEditors: Hejun Xie
LastEditTime: 2021-04-08 10:46:05
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import os
import numpy as np
import xarray as xr
from pytmatrix import orientation
from pytmatrix.tmatrix import Scatterer
import multiprocessing as mp 
from scipy import stats
import yaml

# Local imports
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.interp import quadrature
from ZJU_AERO.hydro import create_hydrometeor
from ZJU_AERO.utils import dielectric_ice, dielectric_water, dielectric_mixture


# Whether to regenerate  all lookup tables, even if already present
FORCE_REGENERATION_SCATTER_TABLES = True


with open('db.yaml', 'r') as yamlfile:
    db_cfg = yaml.load(yamlfile, Loader=yaml.SafeLoader)

# Load DBS settings
for key in db_cfg['DBS']:
    locals()[key] = db_cfg['DBS'][key]

# Load DIMENSIONS settings:
for key in db_cfg['DIMENSIONS']:
    if key in ['ELEVATIONS', 'BETAS']:
        params = db_cfg['DIMENSIONS'][key]
        locals()[key] = np.arange(params[0], params[1], params[2])
    else:
        locals()[key] = db_cfg['DIMENSIONS'][key]
    
    if isinstance(locals()[key], list):
        locals()[key] = np.array(locals()[key])

# the key in levelA and levelB dictionaries
real_variables      = [ 'p11_bw', 'p12_bw', 'p13_bw', 'p14_bw',
                        'p21_bw', 'p22_bw', 'p23_bw', 'p24_bw',
                        'p31_bw', 'p32_bw', 'p33_bw', 'p34_bw',
                        'p41_bw', 'p42_bw', 'p43_bw', 'p44_bw']

real_variables_map  =  { 'p11_bw':(0,0), 'p12_bw':(0,1), 'p13_bw':(0,2), 'p14_bw':(0,3),
                         'p21_bw':(1,0), 'p22_bw':(1,1), 'p23_bw':(1,2), 'p24_bw':(1,3),
                         'p31_bw':(2,0), 'p32_bw':(2,1), 'p33_bw':(2,2), 'p34_bw':(2,3),
                         'p41_bw':(3,0), 'p42_bw':(3,1), 'p43_bw':(3,2), 'p44_bw':(3,3)}

complex_variables   =  [ 's11_fw', 's12_fw', 
                         's21_fw', 's22_fw']

complex_variables_map   = { 's11_fw':(0,0), 's12_fw':(0,1), 
                            's21_fw':(1,0), 's22_fw':(1,1)}

def dielectric_solid(t, f, rho_solid):
    """
    Compute the complex dielectric constant of snow (a mixture of ice and air), based on
    the article of H.Liebe: "A model for the complex permittivity of
    water at frequencies below 1 THz"
    Args:
        rho_solid: density of solid hydrometeors in kg/m^3
        t: temperature in K
        f: frequency in GHz
    Returns:
        m: the complex dielectric constant m = x + iy
    """
    rho_ice = 918 # [kg/m^3]

    frac_ice_volume = rho_solid / rho_ice
    mix = [frac_ice_volume, 1-frac_ice_volume]

    return dielectric_mixture(mix, [dielectric_ice(t, f), constants.M_AIR])

def _gen_one_tm(hydrom_type, list_elevation, list_beta, D, T, AR, frequency, wavelength):
    '''
    Params:
        scatt: A scatter instance defined in pytmatrix.
        list_elevation: list of radar elevation
        list_beta: list of Euler angle beta
    '''
    # Get m
    if hydrom_type in ['S']:
        m = dielectric_solid(T, frequency, 100)
    elif hydrom_type in ['G']:
        m = dielectric_solid(T, frequency, 500)
    elif hydrom_type in ['I']:
        m = dielectric_solid(T, frequency, 300)
    elif hydrom_type in ['R']:
        m = dielectric_water(T, frequency)

    # Get configuration for one T-matrix
    scatt = Scatterer(radius=D/2., radius_type=Scatterer.RADIUS_MAXIMUM, 
            wavelength=wavelength, m=m, axis_ratio=AR, ndgs=10)

    temp_dic = dict()
    dims = ['elevation', 'beta']
    coords = {'elevation': list_elevation, 
              'beta': list_beta}
    size = tuple([len(coords[dim]) for dim in dims])

    for real_variable in real_variables:
        temp_dic[real_variable] = (dims, np.empty(size, dtype='float32'))
    for complex_variable in complex_variables:
        temp_dic[complex_variable] = (dims, np.empty(size, dtype='complex64'))
    temp_lut = xr.Dataset(temp_dic, coords)
    
    for elevation in list_elevation:
        # elevation = 1.0

        geom_back=(90-elevation, 180-(90-elevation), 0., 180, 0.0, 0.0)
        geom_forw=(90-elevation, 90-elevation, 0., 0.0, 0.0, 0.0)

        for beta in list_beta:
            # beta = 30.0
            
            scatt.beta_p = np.array([beta], dtype='float32')
            scatt.beta_w = np.array([1.], dtype='float32')
            scatt.orient = orientation.orient_averaged_fixed
            scatt.n_alpha = NALPHA

            loc_dict = dict(elevation=elevation, beta=beta)
            
            # Back Scattering Z
            scatt.set_geometry(geom_back)
            Z = scatt.get_Z()
            # zdr = (Z[0,0] - Z[0,1] - Z[1,0] + Z[1,1]) / (Z[0,0] + Z[0,1] + Z[1,0] + Z[1,1])
            # ZDR = 10 * np.log10(zdr)
            # print(Z)
            # print(ZDR)
            
            for real_variable in real_variables:
                real_index = real_variables_map[real_variable]
                temp_lut[real_variable].loc[loc_dict] = Z[real_index[0], real_index[1]]

            # Forward Scattering S
            scatt.set_geometry(geom_forw)
            S = scatt.get_S()
            # kdp = S[1,1].real - S[0,0].real 
            # print(S)
            # print(kdp)

            for complex_variable in complex_variables:
                complex_index = complex_variables_map[complex_variable]
                temp_lut[complex_variable].loc[loc_dict] = S[complex_index[0], complex_index[1]]
            
            # exit()

    # exit()

    del scatt

    if ARS[hydrom_type] != 'SINGLE':
        return (temp_lut, D, T, AR)
    else:
        return (temp_lut, D, T, None)

def assign_dataset_pieces(pack):
    ds, D, T, AR = pack[0], pack[1], pack[2], pack[3]
    
    if AR is not None:
        print('Dmax={:>.3f}, temperature={:>.3f}, aspect_ratio={:>.3f}'.format(D, T, AR))
        loc_dict = dict(Dmax=D, temperature=T, aspect_ratio=AR)
    else:
        print('Dmax={:>.3f}, temperature={:>.3f}'.format(D, T))
        loc_dict = dict(Dmax=D, temperature=T)
    
    '''
    There is an annoying issue with the value assignment
    for xarray DataArray using '.loc', we just transpose the two-dimension array
    to fix the issue.
    '''
    for real_variable in real_variables:
        levela_lut[real_variable].loc[loc_dict] = np.transpose(ds[real_variable].data)
    for complex_variable in complex_variables:
        levela_lut[complex_variable].loc[loc_dict] = np.transpose(ds[complex_variable].data)
    
    del ds
    return

def gen_levelA(hydrom_type, frequency, levela_name_lut):
    '''
        Generate levelA scattering property database with ECBM T-matrix solver
        see pytmatrix
        Params:
            hydrom_type: 'S', 'G', 'R', 'I'
            frequency: frequencies for which to obtain the lookup tables, in GHz.
            levela_name_lut: the filename of level A database (Output).
    '''

    global levela_lut

    scheme = '1mom'
    hydrom = create_hydrometeor(hydrom_type, scheme)
    
    list_D = np.linspace(hydrom.d_min, hydrom.d_max, NUM_DIAMETERS).astype('float32') # [mm]
    wavelength=constants.C / (frequency * 1E09) * 1000 # [mm]
    list_elevation = ELEVATIONS # [deg]
    list_beta = BETAS # [deg]
    # Get list_temperature [K]
    if hydrom_type in ['S','G','I']:
        list_temperature = TEMPERATURES_SOL 
    elif hydrom_type in ['R']:
        list_temperature = TEMPERATURES_LIQ
    # Get list_AR [ar<1. for snowplates]
    if ARS[hydrom_type] != 'SINGLE':
        params = ARS[hydrom_type]
        list_AR = np.arange(params[0], params[1], params[2])
    else:
        list_AR = [None]

    # start formulating the levelB database
    dims = ['aspect_ratio', 'temperature', 'Dmax', 'beta', 'elevation']
    coords = {'aspect_ratio': list_AR,
        'temperature': list_temperature,
        'Dmax': list_D, 
        'beta': list_beta,
        'elevation': list_elevation}
    if isinstance(list_AR, list): # == [None]
        coords.pop('aspect_ratio')
        dims.remove('aspect_ratio')
    # print(coords)
    # exit()
    size = tuple([len(coords[dim]) for dim in dims])
    datadic = {}
    for real_variable in real_variables:
        datadic[real_variable] = (dims, np.empty(size, dtype='float32'))
    for complex_variable in complex_variables:
        datadic[complex_variable] = (dims, np.empty(size, dtype='complex64'))
    levela_lut = xr.Dataset(datadic, coords)

    # pool = mp.Pool(processes=mp.cpu_count(), maxtasksperchild=1)
    pool = mp.Pool(processes=mp.cpu_count())
    for D in list_D:
        for T in list_temperature:
            for AR in list_AR:
                if AR is None:
                    AR = hydrom.get_aspect_ratio(D)
                # pack = _gen_one_tm(hydrom_type, list_elevation, list_beta, D, T, AR, frequency, wavelength)
                # assign_dataset_pieces(pack)
                # exit()
                '''
                1. The args should have no user defined class
                '''            
                args = (hydrom_type, list_elevation, list_beta, D, T, AR, frequency, wavelength)
                pool.apply_async(_gen_one_tm, args=args, callback=assign_dataset_pieces)
    
    # Gather processes of multiprocess
    pool.close()
    pool.join()

    # pack = _gen_one_tm(hydrom_type, list_elevation, list_beta, 20.0, 253.0, 1.5, frequency, wavelength)
    # assign_dataset_pieces(pack)

    # print(levela_lut)
    print(levela_name_lut)

    levela_lut.to_netcdf(levela_name_lut, engine="h5netcdf")

    del levela_lut    


if __name__ == "__main__":
    '''
       Create all lookup tables for the specified hydrometeor types and
        microphysical schemes (currently only 1mom scheme implemented...)
    '''

    scheme = '1mom'
    
    for frequency in FREQUENCIES:
        for hydrom_type in HYDROM_TYPES:
            
            # The name of the lookup table is lut_SZ_<hydro_name>_<freq>_<scheme>_<level_name>.nc
            levela_name_lut = (FOLDER_LUT+"lut_SZ_"+hydrom_type+'_'+
                    str(frequency).replace('.','_')+'_'+scheme+"_LevelA"+".nc")
            
            if (FORCE_REGENERATION_SCATTER_TABLES
                or not os.path.exists(levela_name_lut)):
                msg = '''
                Generating scatter table for 1 moment scheme,
                hydrometeor = {:s}
                freq = {:s}
                '''.format(hydrom_type, str(frequency))
                print(msg)

                gen_levelA(hydrom_type, frequency, levela_name_lut)

