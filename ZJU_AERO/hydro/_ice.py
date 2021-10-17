'''
Description: hydrometeor ice
Author: Hejun Xie
Date: 2020-11-13 12:14:03
LastEditors: Hejun Xie
LastEditTime: 2021-10-17 21:01:23
'''

# Global imports
import numpy as np
np.seterr(divide='ignore')

# Local imports
from ..const import global_constants as constants
from ._hydrometeor import _Hydrometeor


class IceParticle(_Hydrometeor):
    '''
    Class for ice crystals
    '''
    def __init__(self, scheme, scheme_name='wsm6'):
        """
        Create a IceParticle Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            scheme_name: microphysics scheme name. Ex: wsm6, wdm6, lin, thompson...
        Returns:
            An IceParticle class instance (see below)
        """

        if scheme_name == 'wsm6':
            from ..const import constants_wsm6 as constants_1mom
        elif scheme_name == 'thompson':
            from ..const import constants_thompson as constants_1mom

        self.scheme = scheme
        self.scheme_name = scheme_name
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_I
        self.d_min = constants_1mom.D_MIN_I
        self.list_D = np.linspace(self.d_min, self.d_max, self.nbins_D)

        # Power-law parameters
        self.a = constants_1mom.AM_I
        self.b = constants_1mom.BM_I
        self.alpha = constants_1mom.AV_I
        self.beta = constants_1mom.BV_I
        self.f = constants_1mom.FV_I

        # PSD parameters
        self.N0 = None      # Taken from Field et al (2005), not true N0
        self.lambda_ = None # Taken from Field et al (2005), not true lambda
        self.mu = constants_1mom.MU_I
        self.nu = 1

        # Scattering parameters
        # See Noel and Sassen 2004
        self.canting_angle_std = 5.

        # Others
        self.lambda_factor = None
        self.ntot_factor = None
        self.vel_factor = None

        # Some exceptional const. attributes or functions by constants module
        if scheme_name in ['wsm6', 'thompson']:
            self.PHI_23_I = constants_1mom.PHI_23_I

    def get_N(self,D):
        """
        Returns the PSD in mm-1 m-3 using the moment explained in the paper
        based on the double-normalized PSD for moments 2,3, from
        Field et al. (2005)
        Args:
            D: vector or matrix of diameters in mm
            self.lambda_ in dimension of (n_valid_gates)
            D in dimension of (nbins_D)

        Returns:
            N: the number of particles for every diameter, with dimension (n_valid_gates, nbins_D)
        """
        
        # Taken from Field et al (2005)
        x = np.outer(self.lambda_, D) / 1000. # [-]
        return self.N0[:,None] * self.PHI_23_I(x)
        
    def integrate_V(self):
        """
        Integrates the terminal fall velocity over the PSD
        Args:
            None

        Returns:
            v: the integral of V(D) * N(D)
            n: the integratal of the PSD: N(D)
        """
        # Again no analytical solution here
        D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        dD = D[1] - D[0]
        N = self.get_N(D)
        V = self.get_V(D)
        v = np.sum(N*V) * dD
        n = np.sum(N)   * dD
        if np.isscalar(v):
            v = np.array([v])
            n = np.array([n])
        return v, n

    def get_mom_2(self, T, Qn, n):
        """
        Get second moment from third moment using the best fit provided
        by Field et al (2005):
            Qn = a(n, Tc)Q2^{b(n, Tc)}, 
            so Q2 = {Qn/a(n, Tc)}^{1./b(n, Tc)}
        Args:
            T: temperature in K
            Qn: n-order moment, i.e., integration of D^n * N(D),
             in SI unit [m^(n-3)]
            n: order of Qn

        Returns:
            Q2: 2-order moment, i.e., integration of D^2 * N(D),
             in SI unit [m^-1]
        """
        T = T - constants.T0 # Convert to celcius
        a = 5.065339 - 0.062659 * T -3.032362 * n + 0.029469 * T * n \
            - 0.000285 * T**2  + 0.312550 * n**2 + 0.000204 * T**2*n \
            + 0.003199 * T* n**2 - 0.015952 * n**3

        a  = 10**(a)
        b = 0.476221 - 0.015896 * T + 0.165977 * n + 0.007468 * T * n \
            - 0.000141 * T**2 + 0.060366 * n**2 + 0.000079 * T**2 * n \
            + 0.000594 * T * n**2 -0.003577 * n**3

        return (Qn/a) ** (1./b)

    def set_psd(self, arg1, arg2):
        """
        Sets the particle size distribution parameters
        Args:
            *args: for the one-moment scheme, a tuple (T, QM) containing the
                temperature in K and the mass concentration QM,
                for the two-moment scheme, a tuple (QN, QM), where QN is the
                number concentration in m-3 (not used currently, not implemented yet...)

        Returns:
            No return but initializes class attributes for psd estimation
        """

        # Reference is Field et al. (2005)
        # We use the two-moment normalization to get a PSD
        # if one moment scheme, arg1 = T, arg2 = Q(M)I
        # if two moments scheme, arg1 =  QNI, arg2 = Q(M)I

        T  = arg1.astype(np.float64)
        QM = arg2.astype(np.float64) # [kg m-3]
        Qn = QM / self.a * (1.0E-3)**self.b # [m^(n-3)]
        
        # get first estimated N0 and lambda_
        Q2 = self.get_mom_2(T, Qn, self.b) # [m^(-1)]
        N0 = Q2**((self.b+1)/(self.b-2)) * Qn**((2+1)/(2-self.b)) # [m-4] 
        N0 *= 1.0E-3 # convert from [m-4] to [mm-1 m-3]
        lambda_ = (Q2/Qn)**(1/(self.b-2))

        # Apply correction factor to match third moment
        D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        dD = D[1]-D[0]
        x = lambda_[:,None] * D.T / 1000 # [-]
        N = N0[:,None] * self.PHI_23_I(x)
        QM_est = np.nansum(self.a * D**self.b * N, axis=1) * dD
        N0 = N0 * (QM / QM_est)

        # assign the attributes
        self.N0 = N0.T
        self.lambda_ = lambda_.T
        self.ntot = np.nansum(N, axis=1) * dD

        # vectorlize
        if np.isscalar(self.lambda_):
            self.lambda_ = np.array([self.lambda_])
        if np.isscalar(self.ntot):
            self.ntot = np.array([self.ntot])

    def get_aspect_ratio(self,D):
        """
        Return the aspect ratio for specified diameters, based on  Auer and
        Veal (1970)
        Args:
            D: the vector of diameters

        Returns:
            The aspect-ratios defined by the smaller dimension over the
            larger dimension
        """
        ar = 11.3 * D**0.414 * 1000**(-0.414)
        return 1/ar

