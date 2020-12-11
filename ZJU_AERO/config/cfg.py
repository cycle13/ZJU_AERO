'''
Description: Define the global CONFIG, and define
how to initialize and check it.
Author: Hejun Xie
Date: 2020-11-21 11:06:42
LastEditors: Hejun Xie
LastEditTime: 2020-12-11 10:00:04
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


def createConfig(options_file):
    '''
    Initialize CONFIG, which is a global variable, because it needs to be
    accesses everywhere in the radar operator (this might not be the most
    pythonic way of doing it, but I couldn't find a better way...) '
    '''
    global CONFIG
    CONFIG = ConfigClass(options_file)
    CONFIG.check_sanity()

class ConfigClass(object):
    def __init__(self, options_file):
        '''
        Initialites a user ConfigClass instance, by reading a yml file and parsing all its
            values
        Args:
            options_file: name of the .yml user configuration to read

        Returns:
            A ConfigClass instance.
        '''
        
        self.__config = None

        try:
            with open(options_file, 'r') as ymlfile:
                print('Loading user config file: {}'.format(options_file))
                self.__config = yaml.load(ymlfile, Loader=yaml.SafeLoader)
        except Exception as e:
            self.__config = copy.deepcopy(DEFAULTS)
            msg = '''
            'Could not find or read {}, using default options...
            The error was
            '''.format(options_file)
            print(dedent(msg))
            print(e)

    def __getitem__(self, key):
        '''
        Define how you can use a ConfigClass instance as if it was a dictionary.
        '''
        return self.__config[key]
    
    def __str__(self):
        '''
        Define how you can print a ConfigClass instance.
        '''

        str_list = []
        for section in VALID_VALUES:
            str_list.append('Section: {}'.format(section))
            for key in VALID_VALUES[section]:
                str_list.append(' '*4 + '{}: {}'.format(key, self.__config[section][key]))
        
        return '\n'.join(str_list)
    
    def check_sanity(self):
        '''
        Check all all keys in this config instance, if they are
            invalid replace them with the default value, unless they are mandatory
            (no default)
        Args:
            No args, but check self.__config.

        Returns:
            No returns, but update self.__config.
        '''

        print('Start checking user configuration')
        for section in VALID_VALUES:
            # Initialize this section if user wants default settings
            if section not in self.__config.keys():
                self.__config[section] = {}

            for key in VALID_VALUES[section]:
                # Deal with absent key
                if key not in self.__config[section].keys():
                    # NOT MANDATORY
                    if key in DEFAULTS[section].keys():
                        self.__config[section][key] = DEFAULTS[section][key]
                    # MANDATORY
                    else:
                        msg = '''
                            This key: {} is mandatory,
                            please provide a valid value, aborting...
                            '''.format(key)
                        raise ValueError(dedent(msg)) 

                # Check validity
                if key in VALID_VALUES[section].keys():
                    flag_ok = self._check_validity(self.__config[section][key],
                                        VALID_VALUES[section][key])

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
                        # IF DEFAULTS valus is None, then no correction by sanity check will be made
                        if (key in DEFAULTS[section].keys()) and (DEFAULTS[section][key] is not None):
                            msg = '''
                            The default value {} was assigned
                            '''.format(DEFAULTS[section][key])
                            print(dedent(msg))
                            self.__config[section][key] = DEFAULTS[section][key]
                        else:
                            msg = '''
                            This key is mandatory, or no proper defaults can be provided, 
                            please provide a valid value, aborting...
                            '''
                            raise ValueError(dedent(msg))    
        print(self)

    @staticmethod
    def _check_validity(input_value, valid_value):
        '''
        A static method checking the validity a specific key in the user specified configuration.
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
            flag_valid = any([ConfigClass._check_validity(input_value, i) for i in valid_value])
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

