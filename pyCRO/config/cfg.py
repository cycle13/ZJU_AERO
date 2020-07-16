'''
@Description: parses the user configuration from a .yml file, defines defaults
and checks if entered values are valid
@Author: Hejun Xie
@Date: 2020-07-16 10:10:37
@LastEditors: Hejun Xie
@LastEditTime: 2020-07-16 10:27:18
'''


# Global imports
import numpy as np
import yaml
import copy
import builtins
import re
from textwrap import dedent

# Local imports, see utilities.py for the definition of Range and TypeList
from cosmo_pol.utilities import (Range, TypeList)

'''
Initialize CONFIG, which is a global variable, because it needs to be
accesses everywhere in the radar operator (this might not be the most
pythonic way of doing it, but I couldn't find a better way...) '
'''

global CONFIG
CONFIG = None

'''
Defines the defaults values in a dictionnary, if a field is present in
VALID_VALUES but absent from DEFAULTS, this means that it mandatory and an
error will be returned if no valid value is provided, ex. frequency
'''

DEFAULTS={
    'radar':
        {'type':'ground',\
        'range':150000,\
        'radial_resolution':500,\
        'PRI':700,\
        'FFT_length':256,\
        'sensitivity':[-5,10000],\
        '3dB_beamwidth':1.,\
        'K_squared': 0.93,
        'antenna_speed': 0.2,
        'frequency': 9.41},\
     'refraction':
        {'scheme':1},\
    'integration':
        {'scheme':1,\
        'nv_GH':9,\
        'nh_GH':3,\
        'n_gaussians':7,\
        'weight_threshold':1.,\
        'nr_GH':7,\
        'na_GL':7},\
    'doppler':
        {'scheme':1,\
        'turbulence_correction':0,\
        'motion_correction':0},\
    'microphysics':
        {'scheme':'1mom',\
         'with_melting': 0,\
         'with_ice_crystals':1,\
         'with_attenuation': 1,\
         'scattering':'tmatrix_masc',\
         'folder_lut':None}\
    }

'''
Defines the valid values in a dictionnary, if a field is present in
VALID_VALUES but absent from DEFAULTS, this means that it mandatory and an
error will be returned if no valid value is provided, ex. frequency
'''

VALID_VALUES={
    'radar':
        {'type': ['ground', 'GPM'],\
        'coords': TypeList([float, int],[3]),\
        'frequency':[2.7,5.6,9.41,9.8,13.6,35.6],\
        'range': Range(5000,500000),\
        'radial_resolution': Range(25,5000),\
        'PRI': Range(10,3000),\
        'FFT_length':Range(16,2048),\
        'sensitivity': [TypeList([float,int],[3]),TypeList([float,int],[2]),
                        float],\
        '3dB_beamwidth': Range(0.1,10.),\
         'K_squared': [float, None],\
         'antenna_speed': Range(1E-6,10.)},\
    'refraction':
        {'scheme':[1,2]},\
    'integration':
        {'scheme':[1,2,3,4,'ml'],\
        'nv_GH': range(1,31,2),\
        'nh_GH': range(1,31,2),\
        'n_gaussians': range(1,13,2),\
        'weight_threshold':  Range(0.0001,1.),\
        'nr_GH': range(1,31,2),\
        'na_GL': range(1,31,2)},\
    'doppler':
        {'scheme':[1,2,3],\
        'turbulence_correction':[0,1],
        'motion_correction':[0,1]},\
    'microphysics':
        {'scheme':['1mom','2mom'],\
        'with_ice_crystals':[0,1],\
        'with_melting':[0,1],\
        'with_attenuation':[0,1],\
        'scattering':['tmatrix_masc','tmatrix','dda'],\
        'folder_lut': [None, str]}\
    }


def _check_validity(input_value, valid_value):
    '''
    Checks the validity a specific key in the user specified configuration
    Args:
        input_value: the value provided by the user
        valid_value: the valid value specified in the VALID_VALUES dictionary
            for this particular key

    Returns:
        flag_valid: a boolean that states if the provided value is valid
    '''
    
    flag_valid = False
    if type(input_value) == list:
        # If input is a list, check all elements in the list
        flag_valid = all([_check_validity(i,valid_value) for i in input_value])

    else:
        # Check if valid value is a type
        if type(valid_value) == builtins.type:
            flag_valid = type(input_value) == valid_value
        # Check if valid value is a list or a Range
        elif type(valid_value) == list:
            flag_valid = any([_check_validity(input_value,i) for i in valid_value])
        elif type(valid_value) == Range or type(valid_value) == range:
            flag_valid = input_value in valid_value
        # Check if valid value is a string with a regex
        elif type(valid_value) == str and valid_value[0:5] == '-reg-':
            # See if input matches regex (the \Z is used to match end of string)
            if re.match(valid_value[5:]+'\Z',input_value):
                flag_valid = True
        # Last possibility is TypeList
        elif type(valid_value) == TypeList:
            flag_valid = valid_value == input_value
        else:
            # Last possibility is that valid_value is a single value
            flag_valid = valid_value == input_value
    
    return flag_valid


def init(options_file):
    '''
    Initialites the user CONFIG, by reading a yml file and parsing all its
        values
    Args:
        options_file: name of the .yml user configuration to read

    Returns:
        No output but updates the CONFIG global
    '''
    global CONFIG
    CONFIG = None
    try:
        with open(options_file, 'r') as ymlfile:
            CONFIG = yaml.load(ymlfile)
    except Exception as e:
        CONFIG = copy.deepcopy(DEFAULTS)
        msg = '''
        'Could not find or read {:s}, using default options...
        The error was:
        '''
        print(dedent(msg))
        print(e)

        return

    return CONFIG

def sanity_check(config):
    '''
    Check all all keys in the provided config dictionary, if they qre
        invalid replace them with the default value, unless they are mandatory
        (no default)
    Args:
        config: the user specified configuration in the form of a dictionary

    Returns:
        The parsed user input in the form of a dictionary
    '''

    # Parsing values
    for section in VALID_VALUES:
        if section not in config.keys():
            config[section] = {}

        for key in VALID_VALUES[section]:
            if key not in config[section].keys():
                config[section][key] = DEFAULTS[section][key]

            if key in VALID_VALUES[section].keys():
                flag_ok = _check_validity(config[section][key],
                                         VALID_VALUES[section][key])

                if not flag_ok:
                    valid = VALID_VALUES[section][key]
                    if type(valid) == list:
                        valid_str = [str(l) for l in valid]
                    else:
                        valid_str = valid

                    msg = '''
                    Invalid value entered for key: {:s}/{:s}
                    The value must be: {:s}
                    '''.format(section,key, str(valid_str))
                    print(dedent(msg))

                    if key in DEFAULTS[section].keys():
                        print('The default value {:s} was assigned'
                                  .format(str(DEFAULTS[section][key])))

                    else:
                        msg = '''
                        This key is mandatory, please provide a
                        valid value, aborting...
                        '''
                        raise ValueError(dedent(msg))

    return config
