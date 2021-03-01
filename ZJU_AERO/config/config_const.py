'''
Description: Define some config constants for 
config subpackage.
Author: Hejun Xie
Date: 2020-11-21 11:07:03
LastEditors: Hejun Xie
LastEditTime: 2021-03-01 18:58:40
'''

# Local imports, see config_types.py for the definition of Range and TypeList
from .config_types import Range, TypeList

'''
Defines the valid values in a dictionnary, if a field is present in
VALID_VALUES but absent from DEFAULTS, this means that it mandatory and an
error will be returned if no valid value is provided, ex. frequency
IF DEFAULTS valus is None, then no correction by sanity check will be made
'''

DEFAULTS={
    'radar':
        {'type': 'ground',
        'coords': None,
        'range': None,
        'radial_resolution': 500,
        'sensitivity': [-5, 10000],
        '3dB_beamwidth': 1.,
        'PRI': 700, # Pulse Repetition Interval [us]
        'FFT_length': 256,
        'nyquist_velocity': None},
    'refraction':
        {'scheme': 1},
    'integration':
        {'scheme': 1,
        'nv_GH': 9,
        'nh_GH': 3,
        'weight_threshold': 1.},
    'nwp':
        {'name': 'wrf',
         'modeltop': 'default'},
    'core':
        {'simulate_doppler': 1,
         'doppler_scheme': 1,
         'doppler_spectrum': 0,
        'engine': 'rdop'},
    'microphysics':
        {'scheme': '1mom',
         'with_melting': 0,
         'with_ice_crystals': 1,
         'with_attenuation': 1,
         'scattering': 'tmatrix_masc',
         'scattering_S': 'default',
         'scattering_R': 'default',
         'scattering_I': 'default',
         'scattering_G': 'default'}
    }

VALID_VALUES={
    'radar':
        {'type': ['ground', 'spaceborne'],
        'coords': [None, TypeList([float, int], dim=[3])],
        'frequency':[2.7, 5.6, 9.41, 9.8, 13.6, 35.6],
        'range': [None, Range(5000, 500000)],
        'radial_resolution': Range(25, 5000),
        'sensitivity': [TypeList([float, int], dim=[3]), 
                        TypeList([float, int], dim=[2]),
                        float],
        '3dB_beamwidth': Range(0.1, 10.),
        'PRI': [None, Range(10, 3000)],
        'FFT_length': [None, Range(16, 2048)],
        'nyquist_velocity': [None, float]}, # TODO: nyquist velocity may have to do with the radar elevation
    'refraction':
        {'scheme': [1, 2, 3]},
    'integration': # TODO: weight_threshold; scheme 2, 3, 4 melting layer: 'ml'
        {'scheme': [1, 2, 3, 4, 'ml'],
        'nv_GH': range(1, 31, 2),
        'nh_GH': range(1, 31, 2),
        'weight_threshold':  Range(0.0001, 1.)},
    'nwp':
        {'name': ['wrf', 'grapes'],
         'modeltop': [float, 'default']}, # default: wrf=20km, grapes=30km
    'core':
        {'simulate_doppler': [0, 1],
         'doppler_scheme': [1, 2],
         'doppler_spectrum': [0, 1],
        'engine': ['rdop', 'wrfda']},
    'microphysics':
        {'scheme': ['1mom', '2mom'],
        'with_ice_crystals': [0, 1],
        'with_melting': [0, 1],
        'with_attenuation': [0, 1],
        'scattering':    ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_S':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_R':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_I':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_G':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'folder_lut': str}
    }
