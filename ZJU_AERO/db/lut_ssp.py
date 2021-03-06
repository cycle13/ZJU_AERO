'''
Description: defines the Lookup Class as well as a set of functions
used to load and save single scattering lookup tables
Author: Hejun Xie
Date: 2020-08-19 22:09:15
LastEditors: Hejun Xie
LastEditTime: 2021-09-29 15:56:50
'''


# Global imports
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error
import os

# Local imports
from ._lut_ssp_numpy import load_lut_numpy
from ._lut_ssp_xarray import load_lut_xarray


def load_all_lut(scheme, list_hydrom, frequency, scattering_method, folder_lut=None):
    '''
    Loads all scattering lookup tables for the user specified parameters
    Args:
        scheme: the microphysical scheme, either '1mom' or '2mom'(not implemented yet)
        list_hydrom: the list of hydrometeors for which scattering properties
            should be obtained: 'R': rain, 'S': snow, 'G': graupel. Ex: ['R','S','G']
        frequency: the frequency in GHz, make sure the lookup tables have
            been previously computed for the corresponding frequency!
        scattering_method: A dictionary describing the scattering method that is used by every
            hydrometeors in list_hydrom, can be either 'tmatrix_masc' or 'iitm_masc',
            which correspond to subfolders in the lookup folder. You could add more...
            Ex: {'S':'iitm_masc', 'R':'tmatrix_masc', ...}
        folder_lut: The folder that contains lookup table

    Returns:
        lut_sz: dictionary containing the lookup table for every hydrometeor
            type given in 'list_hydrom', the lookup tables are instances of
            the class Lookup_table (see below)
    '''

    # Get current directory
    if folder_lut is None:
        folder_lut = os.path.dirname(os.path.realpath(__file__)) + os.sep
    lut_sz = {}

    for h in list_hydrom:
        folder_lut_method = folder_lut + scattering_method[h] + os.sep

        freq_str = str(frequency).replace('.','_')
        if scattering_method[h] in ['tmatrix_masc']:
            name = 'lut_SZ_' + h + '_' + freq_str + '_' + scheme + '.lut'
        else:
            name = 'lut_SZ_' + h + '_' + freq_str + '_' + scheme + '_' + 'LevelB' + '.nc'
        print(folder_lut_method + name)
        try:
            engine = 'numpy' if scattering_method[h] in ['tmatrix_masc'] else 'xarray'
            lut_sz[h] = load_lut(folder_lut_method + name, engine=engine)
        except:
            msg = """
            Could not find lookup table for scheme = {:s}, hydrometeor =
            {:s}, frequency = {:f} and scattering method = {:s}
            """.format(scheme, h, frequency, scattering_method)
            raise IOError(dedent(msg))

    return lut_sz

def load_lut(filename, engine):
    if engine == 'numpy':
        return load_lut_numpy(filename)
    elif engine == 'xarray':
        return load_lut_xarray(filename)

