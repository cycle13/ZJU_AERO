'''
Description: The WRFDA module for computing radar observables
from the interpolated radials given by NWP models, see (Dowell et al. 2011)
'Ensemble Kalman Filter Assimilation of Radar Observations of the 8 May 2003 
Oklahoma City Supercell: Influences of Reflectivity Observations on Storm-Scale Analyses'
Used for comparision with radar operator
Author: Hejun Xie
Date: 2020-10-12 10:46:54
LastEditors: Hejun Xie
LastEditTime: 2020-10-12 17:22:17
'''

# Global imports
import numpy as np
np.warnings.filterwarnings('ignore')

# Local imports
from ..config.cfg import CONFIG
from ..interpolation import Radial
from ..constants import global_constants as constants

def get_radar_observables_wrfda(list_subradials):
    """
    Computes radar variables for all subradials, using WRFDA method, see (Dowell et al. 2011)
    Args:
        list_subradials: list of subradials (Radial claass instances) as
            returned by the interpolation.py code, only simulate Horizontal Reflextivity (ZH)
    
    Returns:
        A radial class instance containing the integrated radar observables
    """
    
    # Get dimensions of subradials
    num_beams = len(list_subradials) # Number of subradials (quad. pts)
    idx_0 = int(num_beams/2) # Index of central subradial

    central_radial = list_subradials[idx_0]

    '''
        Q<hydro> [kg/kg] * RHO [kg/m3] = Q<hydro>_v [kg/m3]
    '''

    qr = central_radial.values['QR_v']
    qs = central_radial.values['QS_v']
    qg = central_radial.values['QG_v']

    zr = 3.63 * 1e9 * qr**1.75
    zs = 9.80 * 1e8 * qs**1.75
    zg = 4.33 * 1e10 * qg**1.75

    ZH = zr + zs + zg

    # Create outputs
    rad_obs = {}
    rad_obs['ZH'] = ZH

    # Finally get vectors of distances, height and lat/lon at the central beam
    heights_radar = central_radial.heights_profile
    distances_radar = central_radial.dist_profile
    lats = central_radial.lats_profile
    lons = central_radial.lons_profile
    mask = central_radial.mask

    # Create final radial
    radar_radial = Radial(rad_obs, mask, lats, lons, distances_radar,
                          heights_radar)

    return radar_radial

