'''
Description: defines a certain number of constants, for hydrometeors
in the Thompson (used in WRF and GRAPES) microphysical scheme
Author: Hejun Xie
Date: 2020-08-19 08:52:23
LastEditors: Hejun Xie
LastEditTime: 2021-10-15 20:08:30
'''


import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error
import scipy.special as spe

# Local imports
from . import global_constants as constants

'''
The list of constants are the following
N0: Intercept parameter in the exponential or gamma PSD [mm-1 m-3]
MU: Shape parameter in the PSD [-]
BM: exponent in the mass diameter relation M = AM * D^BM [-]
AM: intercept in the mass diameter relation M = AM * D^BM [mm^-BM kg]
BV: exponent in the velocity diameter relation V = AV * D^BV [-]
AV: intercept in the velocity diameter relation V = AV * D^BV [mm^-BM m s-1]
D_MIN: minimum considered diameter in the PSD
D_MAX: maximum considered diameter in the PSD
LAMBDA_FACTOR: constant factor used in the integration of the mass to get
    lambda of the PSD, is computed from other parameters, should not be
    changed. Note that it is precomputed to save computation time
VEL_FACTOR: constant factor used in the integration of the fall velocity
    over the PSD, is computed from other parameters, should not be
    changed
NTOT_FACTOR: constant factor used in the integration of the PSD,
     is computed from other parameters, should not be
    changed

For Thompson Microphysics scheme, the density of hydrometeor sphere is fixed,
that means BM=3 for sure. 
'''

# Graupel
RHO_G = 5.0E2   # [kg m-3]
N0_G = 4.0E3    # [mm-1 m-3]
BM_G = 3.0      # [-]
AM_G = np.pi / 6.0 * RHO_G * (1000**-BM_G) # [kg mm-BM_G]
BV_G = 0.80     # [-]
AV_G = 330.0 * (1000**-BV_G) # [m/s mm-BV_G]
MU_G = 0.0      # [-]
D_MIN_G = 0.2   # [mm]
D_MAX_G = 15    # [mm]
LAMBDA_FACTOR_G     = spe.gamma(BM_G + MU_G + 1)
VEL_FACTOR_G        = spe.gamma(BV_G + MU_G + 1)
NTOT_FACTOR_G       = spe.gamma(       MU_G + 1)

# Snow
RHO_S = 1.0E2   # [kg m-3]
# but we use (Field, 2005) snow PSD
# N0_S = 2.0E3 * np.exp(1.2E-1*(T-constants.T0))  # [mm-1 m-3]
BM_S = 3.0      # [-]
AM_S = np.pi / 6.0 * RHO_S * (1000**-BM_S) # [kg mm-BM_S]
BV_S = 0.41     # [-]
AV_S = 11.72 * (1000**-BV_S) # [m/s mm-BV_S]
MU_S = 0.0      # [-]
D_MIN_S = 0.2   # [mm]
D_MAX_S = 20    # [mm]
LAMBDA_FACTOR_S     = spe.gamma(BM_S + MU_S + 1)
VEL_FACTOR_S        = spe.gamma(BV_S + MU_S + 1)
NTOT_FACTOR_S       = spe.gamma(       MU_S + 1)

# Rain
RHO_R = 1.0E3   # [kg m-3]
N0_R = 8.0E3    # [mm-1 m-3]
BM_R = 3.0      # [-]
AM_R = np.pi / 6.0 * RHO_R * (1000**-BM_R) # [kg mm-BM_R]
BV_R = 0.8      # [-]
AV_R = 841.9 * (1000**-BV_R) # [m/s mm-BV_R]
MU_R = 0.0      # [-]
D_MIN_R = 0.1   # [mm]
D_MAX_R = 15     # [mm]
LAMBDA_FACTOR_R     = spe.gamma(BM_R + MU_R + 1)
VEL_FACTOR_R        = spe.gamma(BV_R + MU_R + 1)
NTOT_FACTOR_R       = spe.gamma(       MU_R + 1)

# ICE crystals
# For crystal we use (Field, 2005) particle size distribution, 
# not (generalized) gamma distribution
RHO_I  = 9.16E2 # [kg m-3]
FRAC_I = 0.3 # Fracture that ice crystal occupies in its circumscribed sphere (TODO) 
BM_I = 3.0   # [-]
AM_I = np.pi / 6.0 * RHO_I * FRAC_I * (1000 ** -BM_I) # [kg mm-BM_I]
# Zero fall speed for ice crystals, maybe problematic for doppler spectrum simulation
BV_I = 0.0    # [-]
AV_I = 0.0    # [m/s mm-BV_I]
MU_I = 0.0      # [-]
D_MIN_I = 0.05  # [mm]
D_MAX_I = 2     # [mm]

# This is the phi function for the double-normalized PSD for moments 2,3
# Taken from Field et al. (2005)
PHI_23_I = lambda x : (490.6*np.exp(-20.78 * x) +
                       17.46 * x**(0.6357) * np.exp(-3.290 * x))
