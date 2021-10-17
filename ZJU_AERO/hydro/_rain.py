'''
Description: hydrometeor rain
Author: Hejun Xie
Date: 2020-11-13 12:13:05
LastEditors: Hejun Xie
LastEditTime: 2021-10-17 16:20:20
'''

# Global imports
import numpy as np
np.seterr(divide='ignore')

# Local imports
from ..const import global_constants as constants
from ._hydrometeor import _Hydrometeor


class Rain(_Hydrometeor):
    '''
    Class for raindrops
    '''
    def __init__(self, scheme, scheme_name='wsm6'):
        """
        Create a Rain Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            scheme_name: microphysics scheme name. Ex: wsm6, wdm6, lin, thompson...
        Returns:
            A Rain class instance (see below)
        """

        if scheme_name == 'wsm6':
            from ..const import constants_wsm6 as constants_1mom
        elif scheme_name == 'thompson':
            from ..const import constants_thompson as constants_1mom
        
        self.scheme = scheme
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_R
        self.d_min = constants_1mom.D_MIN_R

        # Power-law parameters
        self.a = constants_1mom.AM_R
        self.b = constants_1mom.BM_R
        self.alpha = constants_1mom.AV_R
        self.beta = constants_1mom.BV_R

        # PSD parameters
        self.lambda_ = None
        self.N0 = constants_1mom.N0_R
        self.mu = constants_1mom.MU_R
        self.nu = 1.0

        # Canting angle stdev, taken from Bringi
        self.canting_angle_std = 7.

        # Others
        self.lambda_factor = constants_1mom.LAMBDA_FACTOR_R
        self.vel_factor = constants_1mom.VEL_FACTOR_R
        self.ntot_factor = constants_1mom.NTOT_FACTOR_R
        self.ntot = None # Total number of particles

    def set_psd(self, *args):
        """
        Sets the particle size distribution parameters
        Args:
            *args: for the one-moment scheme, a tuple (QM) containing the
                mass concentration, for the two-moment scheme, a tuple
                (QN, QM), where QN is the number concentration in m-3
                (not used currently and not implemented yet...)

        Returns:
            No return but initializes class attributes for psd estimation,
            specfically, it initialize self.lambda_ and self.ntot
        """
        if len(args) == 2 and self.scheme == '2mom':
            super(Rain,self).set_psd(*args)
        elif self.scheme == '1mom':
            QM = args[0]
            with np.errstate(divide='ignore'):
                # QM = N0 * a * lambda^-(b+mu+1) * lambda_factor * N0 * a
                _lambda = (self.N0 * self.a * self.lambda_factor / QM) \
                    ** (1. / (self.b + self.mu + 1))
                _lambda = np.array(_lambda)
                _lambda[args[0]==0] = np.nan
                self.lambda_ = _lambda
                self.ntot = self.ntot_factor * self.N0 * \
                    self.lambda_ ** (-(self.mu + 1))
                
                if np.isscalar(self.lambda_):
                    self.lambda_ = np.array([self.lambda_])
                if np.isscalar(self.ntot):
                    self.ntot = np.array([self.ntot])
        else:
            msg = '''
            Invalid call to function, if scheme == ''2mom'',
            input must be tuple of (QN,QM) if scheme == '1mom',
            input must be (QM)
            '''
            print(dedent(msg))

    def get_aspect_ratio(self, D):
        """
        Return the aspect ratio for specified diameters, based on Thurai et
        al (2007)
        Args:
            D: the vector of diameters

        Returns:
            The aspect-ratios defined by the smaller dimension over the
            larger dimension, in the same dimension as D, i.e. (nbins_D)
        """

        if np.isscalar(D):
            D = np.array([D])

        ar = np.zeros((len(D),))
        ar[D<0.7] = 1.0
        mid_diam = np.logical_and(D<1.5,D>=0.7)
        ar[mid_diam] = (1.173 - 0.5165 * D[mid_diam] + 0.4698*D[mid_diam]**2
            - 0.1317*D[mid_diam]**3 - 8.5e-3*D[mid_diam]**4)

        ar[D>=1.5] = (1.065 - 6.25e-2 * D[D>=1.5] - 3.99e-3 * D[D>=1.5] ** 2 +
            7.66e-4 * D[D>=1.5] ** 3 - 4.095e-5 * D[D>=1.5] ** 4)

        # This model tends to diverge for large drops so we threshold it to
        # a reasonable max drop size (10mm)
        ar[D>=10.] = (1.065 - 6.25e-2 * 10.- 3.99e-3 * 10. ** 2 +
            7.66e-4 * 10. ** 3 - 4.095e-5 * 10. ** 4)

        return 1./ar
