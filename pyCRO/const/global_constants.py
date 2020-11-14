'''
@Description: global constants for pyCRO
@Author: Hejun Xie
@Date: 2020-07-16 09:53:33
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 13:32:16
'''

# Global import
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error

# Local import
from ..utils import K_squared



class ConstantClass(object):
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
        anything so high. This is used when interpolation MODEL variables to the
        spaceborne radar beam, because spaceborne is at 407 km from the Earth, 
        so there is no need at all to simulate all radar gates (this would make 
        more than 3000 gates at Ka band...)
    T0: freezing temperature of water in K
    A_AR_LAMBDA_AGG: intercept parameter a in the power-law relation defining the
        value of the Lambda in the gamma distribution for aggregate aspect-ratios
        as a function of diameter
        Lambda(D) = a * D^b
    B_AR_LAMBDA_AGG: exponent parameter b in the power-law relation defining the
        value of the Lambda in the gamma distribution for aggregate aspect-ratios
        as a function of diameter
        Lambda(D) = a * D^b
    A_AR_LAMBDA_GRAU: intercept parameter a in the power-law relation defining the
        value of the Lambda in the gamma distribution for graupel aspect-ratios
        as a function of diameter
        Lambda(D) = a * D^b
    B_AR_LAMBDA_GRAU: exponent parameter b in the power-law relation defining the
        value of the Lambda in the gamma distribution for graupel aspect-ratios
        as a function of diameter
        Lambda(D) = a * D^b
    A_AR_M_AGG: intercept parameter a in the power-law relation defining the
        value of the M in the gamma distribution for aggregate aspect-ratios
        as a function of diameter
        M(D) = a * D^b
    B_AR_M_AGG: exponent parameter b in the power-law relation defining the
        value of the M in the gamma distribution for aggregate aspect-ratios
        as a function of diameter
        M(D) = a * D^b
    A_AR_M_GRAU: intercept parameter a in the power-law relation defining the
        value of the M in the gamma distribution for graupel aspect-ratios
        as a function of diameter
        M(D) = a * D^b
    B_AR_M_GRAU: exponent parameter b in the power-law relation defining the
        value of the M in the gamma distribution for graupel aspect-ratios
        as a function of diameter
        M(D) = a * D^b
    A_CANT_STD_AGG: intercept parameter a in the power-law relation defining the
        value of the standard deviation of aggregates orientations as a function
        of the diameter
        sigma_o(D) = a * D^b
    B_CANT_STD_AGG: exponent parameter b in the power-law relation defining the
        value of the standard deviation of aggregates orientations as a function
        of the diameter
        sigma_o(D) = a * D^b
    A_CANT_STD_GRAU: intercept parameter a in the power-law relation defining the
        value of the standard deviation of graupels orientations as a function
        of the diameter
        sigma_o(D) = a * D^b
    B_CANT_STD_GRAU: exponent parameter b in the power-law relation defining the
        value of the standard deviation of graupels orientations as a function
        of the diameter
        sigma_o(D) = a * D^b
    '''
    def __init__(self):
        '''
        Initialize the constant class by assigning the raw constants
        '''
        self._init_raw_const()
    
    def update(self):
        '''
        Update the derived part of the constant using new CONFIG.
        '''
        global CONFIG
        from ..config.cfg import CONFIG
        if CONFIG is not None:
            self._init_derived_const()
    
    def _init_raw_const(self):
        '''
        initialize the raw constants
        '''

        # 1. Numerical parameters
        self.EPS = np.finfo(np.float).eps
        self.SIMULATED_VARIABLES = ['ZH','DSPECTRUM','RVEL','ZV','PHIDP',
                                    'ZDR','RHOHV','KDP']
        
        # 2. Physical parameters
        self.C = 299792458.
        self.RHO_W = 1000./(1000**3)
        self.RHO_I = 916./(1000**3)
        self.A = 1.6
        self.RHO_0 = 1.225
        self.KE = 4./3.
        self.T0 = 273.15

        # 3.Power laws based on MASC observations
        # 3.1 Axis-ratios:
        # Aggregates (Snow)
        self.A_AR_LAMBDA_AGG =  8.42003348664
        self.B_AR_LAMBDA_AGG =  -0.568465084269
        self.A_AR_M_AGG =  0.0527252217284
        self.B_AR_M_AGG = 0.792594862923

        # Graupel
        self.A_AR_LAMBDA_GRAU =  1.97869286543
        self.B_AR_LAMBDA_GRAU =  -0.426770312328
        self.A_AR_M_GRAU =  0.0743715480794
        self.B_AR_M_GRAU =  0.672627814141

        # 3.2 Canting angles std:
        self.A_CANT_STD_AGG = 30.2393875
        self.B_CANT_STD_AGG = -0.077397563
        self.A_CANT_STD_GRAU = 26.65795932
        self.B_CANT_STD_GRAU = -0.10082787

        # 4. some missing global constant
        self.T_K_SQUARED = 283.15       # [K] Tk = 10 Celsius degree
        self.M_AIR = 1 + 0j             # [-] Mair ~ 1

    def _init_derived_const(self):
        '''
        initialize some derived constants
        '''
        self.WAVELENGTH = self.C/(CONFIG['radar']['frequency']*1E09)*1000 # [mm]
        self.KW = K_squared(CONFIG['radar']['frequency'], self.T_K_SQUARED)

        # get model top 
        if CONFIG['nwp']['modeltop'] == 'default':
            if CONFIG['nwp']['name'] == 'wrf':
                self.MAX_MODEL_HEIGHT = 20000
            elif CONFIG['nwp']['name'] == 'grapes':
                self.MAX_MODEL_HEIGHT = 30000
        else:
            self.MAX_MODEL_HEIGHT = CONFIG['nwp']['modeltop']

        # get radar range 
        if CONFIG['radar']['type'] in ['ground']:
            self.RANGE_RADAR=np.arange(
                CONFIG['radar']['radial_resolution']/2.,
                CONFIG['radar']['range'],
                CONFIG['radar']['radial_resolution'])

global_constants = ConstantClass()
