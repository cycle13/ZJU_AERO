'''
@Description: global constants for pyCRO
@Author: Hejun Xie
@Date: 2020-07-16 09:53:33
LastEditors: Hejun Xie
LastEditTime: 2020-08-15 21:53:29
'''

# Global import
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error

# Local import
from pyCRO.config.cfg import CONFIG


'''
--Independent constants--
EPS: machine epsilon, numbers which differ by less than machine epsilon are
    numerically the same
SIMULATED_VARIABLES:  list of all variables simulated by the radar operator
C: speed of light [m/s]
RHO_W: density of liquid water [kg/mm3 ]
RHO_I: density of ice [kg/mm3 ]
T0 : water freezing temperature [K]
A : Constant used to compute the spectral width due to turbulence
    see Doviak and Zrnic (p.409) [-]
RHO_0: Density of air at the sea level [kg/m3]
KE: 4/3 parameter in the 4/3 refraction model
MAX_MODEL_HEIGHT: The maximum height above which simulated radar gates are
    immediately discarded, i.e. there is not chance the model simulates
    anything so high. This is used when interpolation COSMO variables to the
    GPM beam, because GPM is at 407 km from the Earth, so there is no
    need at all to simulate all radar gates (this would make more than
    3000 gates at Ka band...)
T0: freezing temperature of water in K
'''

class Constant_class(object):
    def __init__(self):

        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Numerical parameters
        self.EPS = np.finfo(np.float).eps
        self.SIMULATED_VARIABLES = ['ZH','DSPECTRUM','RVEL','ZV','PHIDP',
                                    'ZDR','RHOHV','KDP']
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        # Physical parameters
        self.C = 299792458.
        self.RHO_W = 1000./(1000**3)
        self.RHO_I = 916./(1000**3)
        self.A = 1.6
        self.RHO_0 = 1.225
        self.KE = 4./3.
        self.MAX_MODEL_HEIGHT = 35000
        self.T0 = 273.15

        # secondary parameters
        self.RANGE_RADAR=np.arange(
            CONFIG['radar']['radial_resolution']/2.,
            CONFIG['radar']['range'],
            CONFIG['radar']['radial_resolution'])


global_constants = Constant_class()
