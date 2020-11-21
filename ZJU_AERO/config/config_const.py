'''
Description: Define some config constants for 
config subpackage.
Author: Hejun Xie
Date: 2020-11-21 11:07:03
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 17:24:27
'''

# Local imports, see utilities.py for the definition of Range and TypeList
from .config_types import Range, TypeList


DEFAULTS={
    'radar':
        {'type':'ground',
        'range':150000,
        'radial_resolution':500,
        'PRI':700,
        'FFT_length':256,
        'sensitivity':[-5,10000],
        '3dB_beamwidth':1.,
        'nyquist_velocity': None,
        'frequency': 9.41},
    'refraction':
        {'scheme':1},
    'integration':
        {'scheme':1,
        'nv_GH':9,
        'nh_GH':3,
        'weight_threshold':1.},
    'nwp': # modeltop: wrf=20km, grapes=30km
        {'name': 'wrf',
         'modeltop': 'default'},
    'core':
        {'simulate_doppler':1,
        'scheme': 1,
        'engine':'rdop'},
    'microphysics':
        {'scheme':'1mom',
         'with_melting': 0,
         'with_ice_crystals':1,
         'with_attenuation': 1,
         'scattering': 'tmatrix_masc',
         'scattering_S':'default',
         'scattering_R':'default',
         'scattering_I':'default',
         'scattering_G':'default'}
    }

'''
Defines the valid values in a dictionnary, if a field is present in
VALID_VALUES but absent from DEFAULTS, this means that it mandatory and an
error will be returned if no valid value is provided, ex. frequency
'''

VALID_VALUES={
    'radar':
        {'type': ['ground', 'spaceborne'],
        'coords': [TypeList([float, int],[3]), None],
        'frequency':[2.7,5.6,9.41,9.8,13.6,35.6],
        'range': [Range(5000,500000), None],
        'radial_resolution': Range(25,5000),
        'PRI': Range(10,3000),
        'FFT_length':[Range(16,2048), None],
        'sensitivity': [TypeList([float,int],[3]),TypeList([float,int],[2]),
                        float],
        '3dB_beamwidth': Range(0.1,10.),
        'nyquist_velocity': [float, None]},
    'refraction':
        {'scheme':[1,2,3]},
    'integration': # TODO: weight_threshold; scheme 2, 3, 4 melting layer: 'ml'
        {'scheme':[1,2,3,4,'ml'],
        'nv_GH': range(1,31,2),
        'nh_GH': range(1,31,2),
        'weight_threshold':  Range(0.0001,1.)},
    'nwp':
        {'name': ['wrf', 'grapes'],
         'modeltop': [float, 'default']},
    'core':
        {'simulate_doppler': [0,1],
        'scheme': 1,
        'engine':['rdop', 'wrfda']},
    'microphysics':
        {'scheme':['1mom','2mom'],
        'with_ice_crystals':[0,1],
        'with_melting':[0,1],
        'with_attenuation':[0,1],
        'scattering':    ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_S':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_R':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_I':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'scattering_G':  ['tmatrix_masc', 'iitm_masc', 'default'],
        'folder_lut': [str]}
    }
