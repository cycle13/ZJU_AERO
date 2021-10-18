'''
Description: hydrometeor grauple
Author: Hejun Xie
Date: 2020-11-13 12:13:29
LastEditors: Hejun Xie
LastEditTime: 2021-10-18 17:07:23
'''

# Global imports
import numpy as np
np.seterr(divide='ignore')
from textwrap import dedent

# Local imports
from ..const import global_constants as constants
from ._hydrometeor import _Hydrometeor, _NonsphericalHydrometeor


class Graupel(_Hydrometeor):
    '''
    Class for graupel
    '''
    def __init__(self, scheme, scheme_name='wsm6'):
        """
        Create a Graupel Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            scheme_name: microphysics scheme name. Ex: wsm6, wdm6, lin, thompson...
        Returns:
            A Graupel class instance (see below)
        """

        if scheme_name == 'wsm6':
            from ..const import constants_wsm6 as constants_1mom
        elif scheme_name == 'thompson':
            from ..const import constants_thompson as constants_1mom

        self.scheme = scheme
        self.scheme_name = scheme_name
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_G
        self.d_min = constants_1mom.D_MIN_G
        self.list_D = np.linspace(self.d_min, self.d_max, self.nbins_D)

        # Power-law parameters
        self.a = constants_1mom.AM_G
        self.b = constants_1mom.BM_G
        self.alpha = constants_1mom.AV_G
        self.beta = constants_1mom.BV_G
        self.f = constants_1mom.FV_G

        # PSD parameters
        self.N0 = constants_1mom.N0_G
        self.lambda_ = None
        self.mu = constants_1mom.MU_G
        self.nu = 1

        # Scattering parameters
        self.canting_angle_std = 40.

        # Others
        self.lambda_factor = constants_1mom.LAMBDA_FACTOR_G
        self.vel_factor = constants_1mom.VEL_FACTOR_G
        self.ntot_factor = constants_1mom.NTOT_FACTOR_G

        # Exceptional constants defined for particular microphysics schemes:
        if self.scheme_name == 'thompson':
            self.N0_UL = constants_1mom.N0_UL_G
            self.N0_LL = constants_1mom.N0_LL_G
            self.Q0 = constants_1mom.Q0_G
    
    def set_N0(self, QM):
        """
        Sets the interception parameter N0.
        Args:
            QM: Mass concentration of hydrometeor graupel.
        
        a). For microphysics scheme wsm6, N0 is a constant for hydrometeor rain
        b). For microphysics scheme thompson, N0 is the function of graupel mixing ratio (see Reference [1].),
            or mass concentration, equivalently.
            N0_G = max(10, min(0.2/QG, 5.0E3)) # [mm-1 m-3]
        References:
        [1].Thompson, Gregory, et al. "Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. 
        Part II: Implementation of a new snow parameterization." Monthly Weather Review 136.12 (2008): 5095-5115.
        """
        if self.scheme_name == 'wsm6':
            pass # already set in __init__()
        elif self.scheme_name == 'thompson':
            self.N0 = np.array(self.Q0 / QM)
            self.N0[self.N0 >= self.N0_UL] = self.N0_UL
            self.N0[self.N0 <= self.N0_LL] = self.N0_LL

    def set_psd(self,*args):
        """
        Sets the particle size distribution parameters
        Args:
            *args: for the one-moment scheme, a tuple (QM) containing the mass
                concentration QM, for the two-moment scheme, a tuple (QN, QM),
                where QN is the number concentration in m-3
                (not used currently, not implemented yet...)

        Returns:
            No return but initializes class attributes for psd estimation
        """
        if len(args) == 2 and self.scheme == '2mom':
            super(Graupel,self).set_psd(*args)
        elif self.scheme == '1mom':
            QM = args[0]
            self.set_N0(QM)
            with np.errstate(divide='ignore'):
                _lambda = (self.N0 * self.a * self.lambda_factor / QM) ** \
                                    (1. / (self.b + self.mu + 1))
                _lambda = np.array(_lambda)
                _lambda[args[0] == 0] = np.nan
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

    def get_aspect_ratio(self,D):
        """
        Return the aspect ratio for specified diameters, based on Garrett (2015)
        Args:
            D: the vector of diameters

        Returns:
            The aspect-ratios defined by the smaller dimension over the
            larger dimension
        """
        if np.isscalar(D):
            D=np.array([D])
        ar = 0.9*np.ones(len(D),)
        return 1.0/ar

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
        lamb = constants.A_AR_LAMBDA_GRAU * D**constants.B_AR_LAMBDA_GRAU
        loc =  np.ones(len(lamb))
        mu = constants.A_AR_M_GRAU * D**constants.B_AR_M_GRAU
        return lamb, loc, mu

    def get_canting_angle_std_masc(self,D):
        """
        Returns the standard deviation of the distribution of canting angles
        of snow hydrometeors for a given diameter
        Args:
            D: vector of diameters in mm

        Returns:
            the standard deviations of the distribution of canting angles
            for the specified diameters
        """
        cant_std = constants.A_CANT_STD_GRAU * D**constants.B_CANT_STD_GRAU
        return cant_std

class NonsphericalGraupel(_NonsphericalHydrometeor, Graupel):
    '''
    Class for snow in the form of graupel,
    but of nonspherical hydrometeor type, i.e., free a and b
    '''
    def __init__(self, scheme, shape, scheme_name='wsm6'):
        """
            Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            scheme_name: microphysics scheme name. Ex: wsm6, wdm6, lin, thompson...
            
            shape: shape of a particle
                1. hexcol: hexagonal column 
                2. TODO

            Returns:
                A nonspherical Hydrometeor class instance (see below)
        """
        super(NonsphericalGraupel, self).__init__(scheme, scheme_name=scheme_name)

        if shape not in ['spheroid']:
            msg = """
            Invalid Nonspherical graupel shape
            """
            raise ValueError(dedent(msg))
            
        self.shape = shape
        self.list_D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        self.asp_wgt = self.get_asp_wgt(self.list_D)
        self.M = self.get_M(self.list_D)

    def get_list_asp(self):
        '''
        Aspect ratio list of snow
        '''
        return np.arange(1.1, 3.1, 0.1)
    
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
            rho = constants.RHO_G # [kg mm-3]
        else:
            rho = constants.RHO_I
        
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

        QM = args[0]
        
        if np.isscalar(QM):
            _lambda = self.solve_lambda(QM, self.N0)
        else:
            ngates = len(QM)
            _lambda = np.zeros((ngates), dtype='float32')
            for igate, qm in zip(range(ngates), QM):
                if qm < constants.QM_THRESHOLD:
                    continue
                _lambda[igate] = self.solve_lambda(qm, self.N0)
            
        _lambda = np.array(_lambda) # [m-1]
        _lambda[QM < constants.QM_THRESHOLD] = np.nan

        self.lambda_ = _lambda
        self.ntot = self.ntot_factor * self.N0 * \
            self.lambda_ ** (-(self.mu + 1))
