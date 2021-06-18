'''
Description: hydrometeor grauple
Author: Hejun Xie
Date: 2020-11-13 12:13:29
LastEditors: Hejun Xie
LastEditTime: 2021-06-18 11:21:07
'''

# Global imports
import numpy as np
np.seterr(divide='ignore')

# Local imports
from ..const import global_constants as constants
from ..const import constants_wsm6 as constants_1mom
from ._hydrometeor import _Hydrometeor, _NonsphericalHydrometeor


class Graupel(_Hydrometeor):
    '''
    Class for graupel
    '''
    def __init__(self, scheme):
        """
        Create a Graupel Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
        Returns:
            A Graupel class instance (see below)
        """

        self.scheme = scheme
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_G
        self.d_min = constants_1mom.D_MIN_G

        # Power-law parameters
        self.a = constants_1mom.AM_G
        self.b = constants_1mom.BM_G
        self.alpha = constants_1mom.AV_G
        self.beta = constants_1mom.BV_G

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

class NonsphericalGraupel(Graupel, _NonsphericalHydrometeor):
    '''
    Class for snow in the form of graupel,
    but of nonspherical hydrometeor type, i.e., free a and b
    '''
    def __init__(self, scheme, shape):
        """
            Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
            
            shape: shape of a particle
                1. hexcol: hexagonal column 
                2. TODO

            Returns:
                A nonspherical Hydrometeor class instance (see below)
        """
        super(NonsphericalGraupel, self).__init__(scheme)

        if shape not in ['spheroid']:
            msg = """
            Invalid Nonspherical graupel shape
            """
            return ValueError(dedent(msg))
            
        self.shape = shape
        self.list_D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        self.asp_wgt = self.get_asp_wgt(self.list_D)

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
        
        rho = constants.RHO_I # [kg mm-3]
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
