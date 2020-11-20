'''
Description:  computes all radar observables for a given radial
Author: Hejun Xie
Date: 2020-08-20 22:01:10
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 12:29:13
'''

# Global imports
import numpy as np
np.warnings.filterwarnings('ignore')

# Local imports
from ..const import global_constants as constants

from ._rdop_scatter import get_radar_observables_rdop
from ._wrfda_scatter import get_radar_observables_wrfda


def get_radar_observables(list_subradials, lut_sz):
    """
    Computes Doppler and polarimetric radar variables for all subradials
    over ensembles of hydrometeors and integrates them over all subradials at
    the end, by different radar operator core engines.
    Args:
        list_subradials: list of subradials (Radial claass instances) as
            returned by the interpolation.py code
        lut_sz: Lookup tables for all hydrometeor species as returned by the
            load_all_lut function in the lut submodule

    Returns:
        A radial class instance containing the integrated radar observables
    """

    from ..config.cfg import CONFIG
    core_engine = CONFIG['core']['engine']

    if core_engine == 'rdop':
        return get_radar_observables_rdop(list_subradials, lut_sz)
    elif core_engine == 'wrfda':
        return get_radar_observables_wrfda(list_subradials)
    else:
        raise ValueError('Invalid core engine: {}'.format(core_engine)) 


def cut_at_sensitivity(list_subradials):
    '''
    Censors simulated measurements where the reflectivity falls below the
    sensitivity specified by the user, see the wiki for how to define
    the sensitivity in the configuration files
    Args:
        list_subradials: a list of subradials containing the computed radar
            observables
    Returns:
         the list_subradials but censored with the radar sensitivity
    '''
    from ..config.cfg import CONFIG
    sens_config = CONFIG['radar']['sensitivity']
    if not isinstance(sens_config,list):
        sens_config = [sens_config]

    if len(sens_config) == 3: # Sensitivity - gain - snr
        threshold_func = lambda r: (sens_config[0] + constants.RADAR_CONSTANT_DB
                         + sens_config[2] + 20*np.log10(r/1000.))
    elif len(sens_config) == 2: # ZH - range
        threshold_func = lambda r: ((sens_config[0] -
                                     20 * np.log10(sens_config[1] / 1000.)) +
                                     20 * np.log10(r / 1000.))
    elif len(sens_config) == 1: # ZH
        threshold_func = lambda r: sens_config[0]

    else:
        print('Sensitivity parameters are invalid, cannot cut at specified sensitivity')
        print('Enter either a single value of refl (dBZ) or a list of',
              '[refl (dBZ), distance (m)] or a list of [sensitivity (dBm), gain (dBm) and snr (dB)]')

        return list_subradials

    if isinstance(list_subradials[0],list): # Loop on list of lists
        for i,sweep in enumerate(list_subradials):
            for j,b, in enumerate(sweep):
                rranges = (CONFIG['radar']['radial_resolution'] *
                           np.arange(len(b.dist_profile)))
                mask = 10*np.log10(b.values['ZH']) < threshold_func(rranges)
                for k in b.values.keys():
                    if k in constants.SIMULATED_VARIABLES:
                        if k == 'DSPECTRUM':
                            logspectrum = 10 * np.log10(list_subradials[i][j].values[k])
                            thresh = threshold_func(rranges)
                            thresh =  np.tile(thresh,
                                              (logspectrum.shape[1],1)).T
                            list_subradials[i][j].values[k][logspectrum < thresh] = np.nan
                        else:
                            list_subradials[i][j].values[k][mask] = np.nan

    else:
        for i, subradial in enumerate(list_subradials): # Loop on simple list
            rranges = (CONFIG['radar']['radial_resolution'] *
                       np.arange(len(subradial.dist_profile)))
            mask = 10 * np.log10(subradial.values['ZH']) < threshold_func(rranges)
            for k in subradial.values.keys():
                if k in constants.SIMULATED_VARIABLES:
                    list_subradials[i].values[k][mask] = np.nan
    return list_subradials
