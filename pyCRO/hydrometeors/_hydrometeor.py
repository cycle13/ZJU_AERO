'''
Description: Hydrometeor base class,
should not be initialized directly
Author: Hejun Xie
Date: 2020-11-13 13:04:15
LastEditors: Hejun Xie
LastEditTime: 2020-11-13 13:05:26
'''

# Global import
import numpy as np
np.seterr(divide='ignore')


class _Hydrometeor(object):
    '''
    Base class for rain, snow and graupel, should not be initialized
    directly
    '''
    def __init__(self, scheme):
        """
        Create a Hydrometeor Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)

        Returns:
            A Hydrometeor class instance (see below)
        """
        self.precipitating = True
        self.scheme = scheme
        self.d_max = None # Max diam in the integration over D
        self.d_min = None # Min diam in the integration over D
        self.nbins_D = 1024 # Number of diameter bins used in the numerical integrations

        # Power-law parameters
        self.a = None
        self.b = None
        self.alpha = None
        self.beta = None

        # PSD parameters
        self.lambda_ = None # [mm-1]
        self.N0 = None      # [mm-1 m-3] 
        self.mu = None      # [-]
        self.nu = None      # [-]

        # Scattering parameters
        self.canting_angle_std = None

        # Integration factors
        self.lambda_factor = None
        self.vel_factor = None
        self.ntot_factor = None

    def get_N(self, D):
        """
        Returns the PSD in mm-1 m-3
        Args:
            D: vector or matrix (melting particle , not implemented yet) of diameters in mm
                for 1mom solid / liquid hydrometeors, D has dimension as (nbins_D)
                lambda_ has dimension as (n_valid_gates)
                N0 is scalar since it is 1mom scheme and it is fixed 

        Returns:
            N: the number of particles for every diameter, has dimension of (ngates, nbins_D)
        """
        
        operator = lambda x,y: np.outer(x,y)
        # print('lambda {}'.format(self.lambda_.shape))
        # print('D {}'.format(D.shape))
        # print('N0 {}'.format(self.N0.shape))
        if np.isscalar(self.N0):
            return self.N0 * D**self.mu * \
                            np.exp(-operator(self.lambda_, D**self.nu))
        else:
            return operator(self.N0, D**self.mu) * \
                            np.exp(-operator(self.lambda_, D**self.nu))
        
    def get_V(self,D):
        """
        Returns the terminal fall velocity i m/s
        Args:
            D: vector or matrix of diameters in mm

        Returns:
            V: the terminal fall velocities, same dimensions as D
        """
        V = self.alpha * D**self.beta
        return V

    def get_D_from_V(self,V):
        """
        Returns the diameter for a specified terminal fall velocity by
        simply inverting the power-law relation
        Args:
            V: the terminal fall velocity in m/s

        Returns:
            the corresponding diameters in mm
        """
        return (V / self.alpha) ** (1. / self.beta)

    def integrate_V(self):
        """
        Integrates the terminal fall velocity over the PSD
        Args:
            None

        Returns:
            v: the integral of V(D) * N(D)
            n: the integratal of the PSD: N(D) with the same dimension as self.lambda_
        """

        v = (self.vel_factor * self.N0 * self.alpha / self.nu * \
             self.lambda_ ** (-(self.beta + self.mu + 1) / self.nu))
        
        n = self.ntot_factor * self.N0 / self.nu * \
             self.lambda_ ** (-(self.mu + 1) / self.nu)
        
        if np.isscalar(v):
            v = np.array([v])
            n = np.array([n])

        return v, n

    def get_M(self,D):
        """
        Returns the mass of a particle in kg
        Args:
            D: vector or matrix of diameters in mm

        Returns:
            m: the particle masses, same dimensions as D
        """
        return self.a * D**self.b

    def set_psd(self,*args):
        """
        Sets the particle size distribution parameters
        Args:
            *args: for the one-moment scheme, see function redefinition in
                the more specific classes Rain, Snow, Graupel and Hail
                for the two-moment scheme, a tuple (QN, QM), where
                QN is the number concentration in m-3 (not used currently...)
                TODO: Not implemented yet

        Returns:
            No return but initializes class attributes for psd estimation
        """
        if len(args) == 2 and self.scheme == '2mom':
            pass