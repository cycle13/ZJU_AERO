'''
Description: hydrometeor snow
Author: Hejun Xie
Date: 2020-11-13 12:13:17
LastEditors: Hejun Xie
LastEditTime: 2021-10-17 20:59:40
'''

# Global imports
import numpy as np
np.seterr(divide='ignore')
from textwrap import dedent

# Local imports
from ..const import global_constants as constants
from ._hydrometeor import _Hydrometeor, _NonsphericalHydrometeor


class Snow(_Hydrometeor):
    '''
    Class for snow in the form of aggregates
    '''
    def __init__(self, scheme, scheme_name='wsm6'):
        """
        Create a Snow Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            scheme_name: microphysics scheme name. Ex: wsm6, wdm6, lin, thompson...
        Returns:
            A Snow class instance (see below)
        """

        if scheme_name == 'wsm6':
            from ..const import constants_wsm6 as constants_1mom
        elif scheme_name == 'thompson':
            from ..const import constants_thompson as constants_1mom

        self.scheme = scheme
        self.scheme_name = scheme_name
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_S
        self.d_min = constants_1mom.D_MIN_S
        self.list_D = np.linspace(self.d_min, self.d_max, self.nbins_D)

        # Power-law parameters
        self.a = constants_1mom.AM_S
        self.b = constants_1mom.BM_S
        self.alpha = constants_1mom.AV_S
        self.beta = constants_1mom.BV_S
        self.f = constants_1mom.FV_S

        # PSD parameters
        self.N0 = None # For snow, N0 use relation by Field et al. 2005 (QJRMS)
        self.lambda_ = None
        self.mu = constants_1mom.MU_S
        self.nu = 1

        # Scattering parameters
        self.canting_angle_std = 20.

        # Others
        self.lambda_factor = constants_1mom.LAMBDA_FACTOR_S
        self.vel_factor = constants_1mom.VEL_FACTOR_S
        self.ntot_factor = constants_1mom.NTOT_FACTOR_S
        self.ntot = None

    def set_psd(self,*args):
        """
        Sets the particle size distribution parameters
        Args:
            *args: for the one-moment scheme, a tuple (T,QM) containing the
                temperatue in K and the mass concentration QM,
                for the two-moment scheme, a tuple (QN, QM), where QN is the
                number concentration in m-3 (not used currently, not implemented yet...)

        Returns:
            No return but initializes class attributes for psd estimation
        """
        if len(args) == 2 and self.scheme == '2mom':
            super(Snow,self).set_psd(*args)

        elif len(args) == 2 and self.scheme == '1mom':
            T, QM = args[0], args[1]
            # For N0 use relation by Field et al. 2005 (QJRMS)
            N0_23 = 5.65E5 * np.exp(-0.107*(T-constants.T0)) # [m-4] 
            self.N0 = 13.5 * N0_23 / 1.0E3 # [mm^-1 m^-3]
            with np.errstate(divide='ignore'):
                _lambda = (self.a * self.N0 * self.lambda_factor / QM) ** \
                    (1. / (self.b + self.mu + 1))
                _lambda = np.array(_lambda) # in m-1
                _lambda[args[1] == 0] = np.nan
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
            input must be (T,QM)
            '''
            print(dedent(msg))

    def get_aspect_ratio(self, D):
        """
        Return the aspect ratio for specified diameters, based on Brandes
         et al (2007)
        Args:
            D: the vector of diameters

        Returns:
            The aspect-ratios defined by the smaller dimension over the
            larger dimension
        """
        ar = (0.01714 * D + 0.8467) # Brandes et al 2007 (Colorado snowstorms)
        return 1.0 / ar

    def get_aspect_ratio_pdf_masc(self,D):
        """
        Returns the parameters of the gamma distribution of aspect-ratios
        of snow for the specific diameter.
        the gamma pdf has the following form, where x is the aspect ratio
        p(x) = ((a - loc)^(lamb-1) exp(-(x - 1)/mu)) / mu^lamb* Gamma(lamb)

        Args:
            D: vector of diameters in mm

        Returns:
            lamb: lambda parameter of the gamma pdf
            loc: location parameter of the gamma pdf
            mu: shape parameter of the gamma pdf

        """
        lambd = constants.A_AR_LAMBDA_AGG * D**constants.B_AR_LAMBDA_AGG
        loc = np.ones(len(lambd))
        mu = constants.A_AR_M_AGG * D**constants.B_AR_M_AGG
        return lambd, loc, mu

    def get_canting_angle_std_masc(self, D):
        """
        Returns the standard deviation of the distribution of canting angles
        of snow hydrometeors for a given diameter
        Args:
            D: vector of diameters in mm

        Returns:
            the standard deviations of the distribution of canting angles
            for the specified diameters
        """
        cant_std = constants.A_CANT_STD_AGG * D**constants.B_CANT_STD_AGG
        return cant_std

class NonsphericalSnow(_NonsphericalHydrometeor, Snow):
    '''
    Class for snow in the form of aggregates,
    but of nonspherical hydrometeor type, i.e., free a and b
    '''
    def __init__(self, scheme, shape, param=None, scheme_name='wsm6'):
        """
            Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            scheme_name: microphysics scheme name. Ex: wsm6, wdm6, lin, thompson...
            
            shape: shape of a particle
                1. hexcol: hexagonal column 
                2. TODO
            
            param: param of that shape

            Returns:
                A nonspherical Hydrometeor class instance (see below)
        """
        super(NonsphericalSnow, self).__init__(scheme, scheme_name=scheme_name)

        if shape not in ['hexcol', 'spheroid', 'snowflake']:
            msg = """
            Invalid Nonspherical snow shape
            """
            raise ValueError(dedent(msg))
            
        self.shape = shape
        self.param = param
        self.list_D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        self.asp_wgt = self.get_asp_wgt(self.list_D)
        self.M = self.get_M(self.list_D)

    def get_list_asp(self):
        '''
        Aspect ratio list of snow
        '''
        return np.arange(1.1, 5.2, 0.2)
    
    def _get_mass(self, D, asp):
        '''
        compute the mass of a nonspherical particle

        Args:
            D: maximum dimension of a particle
            asp: aspect ratio (asp > 1.0)

        Returns:
            Mass of a particle
        '''
        
        if self.shape in ['spheroid']:
            rho = constants.RHO_S # [kg mm-3]
        else:
            rho = constants.RHO_I
        V = self._get_volumn(D, asp)
        
        return V * rho
    
    def set_psd(self, *args):
        """
        Sets the particle size distribution parameters
        Args:
            *args: for the one-moment scheme, a tuple (T,QM) containing the
                temperatue in K and the mass concentration QM,
                for the two-moment scheme, a tuple (QN, QM), where QN is the
                number concentration in m-3 (not used currently, not implemented yet...)

        Returns:
            No return but initializes class attributes for psd estimation
        """

        ''''''

        T, QM = args[0], args[1]
        # For N0 use relation by Field et al. 2005 (QJRMS)
        N0_23 = 5.65E5 * np.exp(-0.107*(T-constants.T0)) # [m-4] 
        self.N0 = 13.5 * N0_23 / 1.0E3 # [mm^-1 m^-3]

        if np.isscalar(QM):
            _lambda = self.solve_lambda(QM, self.N0)
        else:
            ngates = len(QM)
            _lambda = np.zeros((ngates), dtype='float32')
            for igate, n0, qm in zip(range(ngates), self.N0, QM):
                if qm < constants.QM_THRESHOLD:
                    continue
                _lambda[igate] = self.solve_lambda(qm, n0)
            
        _lambda = np.array(_lambda) # [m-1]
        _lambda[QM < constants.QM_THRESHOLD] = np.nan

        self.lambda_ = _lambda
        self.ntot = self.ntot_factor * self.N0 * \
            self.lambda_ ** (-(self.mu + 1))
