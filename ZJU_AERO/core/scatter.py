'''
Description:  computes all radar observables for a given radial
Author: Hejun Xie
Date: 2020-08-20 22:01:10
LastEditors: Hejun Xie
LastEditTime: 2020-11-22 19:46:56
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

    # deal with single value case ZH(dBZ)
    if not isinstance(sens_config, list):
        sens_config = [sens_config]

    if len(sens_config) == 3: # [Sensitivity(dBm), Gain(dBm), SNR(dB)]
        threshold_func = lambda r: sens_config[0] + sens_config[1] + sens_config[2] \
                                    + 20 * np.log10(r / 1000.)
    elif len(sens_config) == 2: # [ZH(dBZ), range(m)]
        threshold_func = lambda r: sens_config[0] \
                                    - 20 * np.log10(sens_config[1] / 1000.) \
                                    + 20 * np.log10(r / 1000.)
    elif len(sens_config) == 1: # [ZH(dBZ)]
        threshold_func = lambda r: sens_config[0]
    else:
        print('Sensitivity parameters are invalid, cannot cut at specified sensitivity')
        print('Enter either a single value of refl (dBZ) or a list of',
              '[refl (dBZ), distance (m)] or a list of [sensitivity (dBm), gain (dBm) and snr (dB)]')
        return list_subradials
    
    rranges = constants.RANGE_RADAR

    # Case if list_subradial is a list of radar sweeps
    if isinstance(list_subradials[0], list):
        for isweep, sweep in enumerate(list_subradials):
            for isubradial, subradial in enumerate(sweep):
                subradial = list_subradials[isweep][isubradial]
                _cut_at_sensitivity(subradial, rranges, threshold_func)
    # Case if list_subradial is a single sweep
    else:
        for isubradial, subradial in enumerate(list_subradials):
            subradial = list_subradials[isubradial]
            _cut_at_sensitivity(subradial, rranges, threshold_func)
            
    return list_subradials

def _cut_at_sensitivity(subradial, rranges, threshold_func):
    '''
    Cut the the SIMULATED_VARIABLES by radar sensitivity
    threshold function.

    Args:
        subradial: A single Radial instance.
        rranges: Radar gates distance[m].
        threshold_func: A threshold function.
    Returns:
        A modified subradial instance setting all the 
        SIMULATED_VARIABLES to np.nan where ZH(dBZ) is below 
        the threshold function threshold_func(r):
        ZH(dBZ) < threshold_func(r).
    '''
    
    # Also get rid of radar gates where ZH == np.nan
    mask = np.array( 10 * np.log10(subradial.values['ZH']) < threshold_func(rranges) ) | \
            np.isnan(subradial.values['ZH'])
    for varname in subradial.values.keys():
        if varname in constants.SIMULATED_VARIABLES:
            subradial.values[varname][mask] = np.nan
