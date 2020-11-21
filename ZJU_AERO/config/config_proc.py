'''
Description: Define the global CONFIG, and define
how to initialize and check it.
Author: Hejun Xie
Date: 2020-11-21 11:06:42
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 17:24:05
'''


# Global imports
import numpy as np
import yaml
import copy
import builtins
import re
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
            # See if input matches regex (the $ is used to match end of string)
            if re.match(valid_value[5:],input_value):
                flag_valid = True
        # Last possibility is TypeList
        elif type(valid_value) == TypeList:
            flag_valid = valid_value == input_value
        else:
            # Last possibility is that valid_value is a single value
            flag_valid = valid_value == input_value
    
    return flag_valid


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
            CONFIG = yaml.load(ymlfile, Loader=yaml.SafeLoader)
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
