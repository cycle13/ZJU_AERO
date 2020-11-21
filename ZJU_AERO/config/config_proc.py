'''
Description: Define the global CONFIG, and define
how to initialize and check it.
Author: Hejun Xie
Date: 2020-11-21 11:06:42
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 21:15:47
'''


# Global imports
import re
import yaml
import copy
import numpy as np
import builtins
from textwrap import dedent

# Local imports, see utilities.py for the definition of Range and TypeList
from .config_const import DEFAULTS, VALID_VALUES
from .config_types import Range, TypeList

'''
Initialize CONFIG, which is a global variable, because it needs to be
accesses everywhere in the radar operator (this might not be the most
pythonic way of doing it, but I couldn't find a better way...) '
'''

global CONFIG
CONFIG = None

def config_init(options_file):
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
            print('Loading user config file: {}'.format(options_file))
            CONFIG = yaml.load(ymlfile, Loader=yaml.SafeLoader)
    except Exception as e:
        CONFIG = copy.deepcopy(DEFAULTS)
        msg = '''
        'Could not find or read {}, using default options...
        The error was
        '''.format(options_file)
        print(dedent(msg))
        print(e)

        return

    return CONFIG

def config_sanity_check(config):
    '''
    Check all all keys in the provided config dictionary, if they qre
        invalid replace them with the default value, unless they are mandatory
        (no default)
    Args:
        config: the user specified configuration in the form of a dictionary

    Returns:
        The parsed user input in the form of a dictionary
    '''

    print('Start checking user configuration')
    for section in VALID_VALUES:
        # Initialize this section if user wants default settings
        if section not in config.keys():
            config[section] = {}

        for key in VALID_VALUES[section]:
            # Assign this key-value if user wants default settings
            if key not in config[section].keys():
                config[section][key] = DEFAULTS[section][key]

            # Check validity
            if key in VALID_VALUES[section].keys():
                flag_ok = _check_validity(config[section][key], VALID_VALUES[section][key])

                if not flag_ok:
                    # Get the string specifying the correct input
                    valid = VALID_VALUES[section][key]
                    if type(valid) == list:
                        valid_str = ' or '.join(['<'+str(l)+'>' for l in valid])
                    else:
                        valid_str = valid
                    msg = '''
                    Invalid value entered for key: {}/{}
                    The value must be: {}
                    '''.format(section, key, valid_str)
                    print(dedent(msg))

                    # Assign the key-value if this key is not mandatory
                    if key in DEFAULTS[section].keys():
                        msg = '''
                        The default value {} was assigned
                        '''.format(DEFAULTS[section][key])
                        print(dedent(msg))
                        config[section][key] = DEFAULTS[section][key]
                    else:
                        msg = '''
                        This key is mandatory, please provide a
                        valid value, aborting...
                        '''
                        raise ValueError(dedent(msg))
    
    _print_config(config)
    return config

def _print_config(config):
    '''
    Print the config after sanity checking
    '''

    for section in VALID_VALUES:
        print('Section: {}'.format(section))
        for key in VALID_VALUES[section]:
            print(' '*4 + '{}: {}'.format(key, config[section][key]))


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

    # check if input value is None
    if input_value is None:
        return (None in valid_value)

    # 1. If valid value is list, then check recursively 
    if type(valid_value) == list:
        flag_valid = any([_check_validity(input_value, i) for i in valid_value])
    # 2. Check if valid value is a build type
    elif type(valid_value) == builtins.type:
        flag_valid = (type(input_value) == valid_value)
    # 3. Check if valid value is a Range (range)
    elif type(valid_value) == Range or type(valid_value) == range:
        flag_valid = input_value in valid_value
    # 4. Check if valid value is a string with a regex
    elif type(valid_value) == str and valid_value[0:5] == '-reg-':
        if re.match(valid_value[5:], input_value):
            flag_valid = True
    # 5. Last possibility is TypeList
    elif type(valid_value) == TypeList:
        flag_valid = (valid_value == input_value)
    # 6. Last possibility is that valid_value is a single value, ex: any string like 'iitm_masc'
    else:
        flag_valid = (valid_value == input_value)

    return flag_valid
