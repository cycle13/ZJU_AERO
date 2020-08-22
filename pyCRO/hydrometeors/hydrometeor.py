# -*- coding: utf-8 -*-

'''
Description: hydrometeors.py: Provides classes with relevant functions for all
hydrometeor types considered in the radar operator.
Computes all diameter dependent properties (orientation, aspect-ratio,
dielectric constants, velocity, mass... 
Author: Hejun Xie
Date: 2020-08-18 09:37:31
LastEditors: Hejun Xie
LastEditTime: 2020-08-22 22:18:19
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global import 
import numpy as np
np.seterr(divide='ignore')

from scipy.integrate import trapz
from scipy.optimize import fsolve, root, newton, brentq
from scipy.stats import gamma as gampdf
from scipy.special import gamma, erf
from scipy.interpolate import interp1d
from textwrap import dedent

# Local import
from pyCRO.constants import global_constants as constants
from pyCRO.constants import constants_1mom
from pyCRO.hydrometeors.dielectric import dielectric_ice, dielectric_water, dielectric_mixture
from pyCRO.utilities import vlinspace

def create_hydrometeor(hydrom_type, scheme = '1mom'):
    """
    Creates a hydrometeor class instance, for a specified microphysical
    scheme
    Args:
        hydrom_type: the hydrometeor types, can be either
            'R': rain, 'S': snow aggregates, 'G': graupel
        scheme: microphysical scheme to use, can be either '1mom' (operational
           one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
    Returns:
        A hydrometeor class instance (see below)
    """

    if  hydrom_type == 'R':
       return Rain(scheme)
    elif hydrom_type == 'S':
        return Snow(scheme)
    elif hydrom_type == 'G':
        return Graupel(scheme)
    elif hydrom_type == 'I':
        return IceParticle(scheme)
    else:
        msg = """
        Invalid hydrometeor type, must be R, S, G, I
        """
        return ValueError(dedent(msg))

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
        self.lambda_ = None
        self.N0 = None
        self.mu = None
        self.nu = None

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
            return self.N0*D**self.mu * np.exp(- operator(self.lambda_,D**self.nu))
        else:
            return operator(self.N0,D**self.mu) * \
                            np.exp(- operator(self.lambda_,D**self.nu))
        
    def get_V(self,D):
        """
        Returns the terminal fall velocity i m/s
        Args:
            D: vector or matrix of diameters in mm

        Returns:
            V: the terminal fall velocities, same dimensions as D
        """
        V = self.alpha * D ** self.beta
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
            v: the integral of V(D) * (D)
            n: the integratal of the PSD: N(D) with the same dimension as self.lambda_
        """

        v = (self.vel_factor * self.N0 * self.alpha / self.nu *
             self.lambda_ ** (-(self.beta + self.mu + 1) / self.nu))
        
        n = self.ntot_factor*self.N0/self.nu*self.lambda_**(-(self.mu+1)/self.nu)
        
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
        return self.a*D**self.b

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


class _Solid(_Hydrometeor):
    '''
    Base class for snow, graupel, should not be initialized
    directly
    '''
    
    def get_fractions(self,D):
        """
        Returns the volumic fractions of pure ice and air within
        the melting snow particles
        Args:
            D: vector of diameters in mm

        Returns:
            A nx2 matrtix of fractions, where n is the number of dimensions.
            The first column is the fraction of ice, the second the
            fraction of air
        """
        # Uses COSMO mass-diameter rule
        f_ice = 6*self.a/(np.pi*constants.RHO_I)*D**(self.b-3)
        f_air = 1-f_ice
        return [f_ice, f_air]

    def get_m_func(self,T,f):
        """
        Returns a function giving the dielectric constant as a function of
        the diameters, depending on the frequency and the temperature.
        Used for the computation of scattering properties (see lut submodule)
        Args:
            T: temperature in K
            f: frequency in GHz

        Returns:
            A lambda function f(D) giving the dielectric constant as a function
            of the particle diameter
        """
        def func(D, T, f):
            frac = self.get_fractions(D)
            m_ice = dielectric_ice(T,f)
            return dielectric_mixture(frac, [m_ice, constants.M_AIR])
        return lambda D: func(D, T, f)

class Rain(_Hydrometeor):
    '''
    Class for raindrops
    '''
    def __init__(self, scheme):
        """
        Create a Rain Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
        Returns:
            A Rain class instance (see below)
        """
        
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
        self.canting_angle_std = 10.

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
            with np.errstate(divide='ignore'):
                _lambda = np.array((self.lambda_factor/args[0])**
                                   (1./(4.+self.mu)))
                _lambda[args[0]==0] = np.nan
                self.lambda_ = _lambda
                self.ntot = (self.ntot_factor *
                             self.N0 * self.lambda_**(self.mu - 1))
                
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
        # a reasonable max drop size
        ar[D>=self.d_max] = ar[D<=self.d_max][-1]

        return 1./ar

    def get_m_func(self,T,f):
        return lambda D: dielectric_water(T, f)


class Snow(_Solid):
    '''
    Class for snow in the form of aggregates
    '''
    def __init__(self, scheme):
        """
        Create a Snow Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
        Returns:
            A Snow class instance (see below)
        """

        self.scheme = scheme
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_S
        self.d_min = constants_1mom.D_MIN_S

        # Power-law parameters
        self.a = constants_1mom.AM_S
        self.b = constants_1mom.BM_S
        self.alpha = constants_1mom.AV_S
        self.beta = constants_1mom.BV_S

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
            # For N0 use relation by Field et al. 2005 (QJRMS)
            self.N0 = 13.5*(5.65*10**5*np.exp(-0.107*(args[0]-273.15)))/1000 # mm^-1 m^-3
            with np.errstate(divide='ignore'):
                _lambda = np.array((self.a * self.N0 * self.lambda_factor
                                    / args[1]) ** (1. / (self.b + 1))) # in m-1
                _lambda[args[1] == 0] = np.nan
                self.lambda_ = _lambda
                self.ntot = (self.ntot_factor * self.N0 *
                             self.lambda_ ** (self.mu - 1))
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
        lambd = constants.A_AR_LAMBDA_AGG*D**constants.B_AR_LAMBDA_AGG
        loc = np.ones(len(lambd))
        mu = constants.A_AR_M_AGG*D**constants.B_AR_M_AGG
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
        cant_std = constants.A_CANT_STD_AGG * D ** constants.B_CANT_STD_AGG
        return cant_std

class Graupel(_Solid):
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
            with np.errstate(divide='ignore'):
                _lambda = np.array((self.lambda_factor/args[0]) **
                                    (1./(4.+self.mu)))
                _lambda[args[0] == 0] = np.nan
                self.lambda_ = _lambda
                self.ntot = self.ntot_factor * self.N0 * self.lambda_**(self.mu - 1)
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
        lamb = constants.A_AR_LAMBDA_GRAU * D ** constants.B_AR_LAMBDA_GRAU
        loc =  np.ones(len(lamb))
        mu = constants.A_AR_M_GRAU*D**constants.B_AR_M_GRAU
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
        cant_std = constants.A_CANT_STD_GRAU*D**constants.B_CANT_STD_GRAU
        return cant_std

class IceParticle(_Solid):
    '''
    Class for ice crystals
    '''
    def __init__(self,scheme):
        """
        Create a IceParticle Class instance
        Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
        Returns:
            An IceParticle class instance (see below)
        """
        self.scheme = scheme
        self.nbins_D = 1024

        self.d_max = constants_1mom.D_MAX_I
        self.d_min = constants_1mom.D_MIN_I

        # Power-law parameters
        self.a = constants_1mom.AM_I
        self.b = constants_1mom.BM_I
        self.alpha = constants_1mom.AV_I
        self.beta = constants_1mom.BV_I

        # PSD parameters
        self.N0 = None      # Taken from Field et al (2005)
        self.lambda_ = None
        self.mu = constants_1mom.MU_I
        self.nu = 1

        # Scattering parameters
        # See Noel and Sassen 2004
        self.canting_angle_std = 5.

        # Others
        self.lambda_factor = constants_1mom.LAMBDA_FACTOR_I
        self.ntot_factor = constants_1mom.NTOT_FACTOR_I
        self.vel_factor = constants_1mom.VEL_FACTOR_I

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
        x = np.outer(self.lambda_, D) / 1000.
        return self.N0[:,None] * constants_1mom.PHI_23_I(x)
        
    def integrate_V(self):
        """
        Integrates the terminal fall velocity over the PSD
        Args:
            None

        Returns:
            v: the integral of V(D) * (D)
            n: the integratal of the PSD: N(D)
        """
        # Again no analytical solution here
        D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        dD = D[1] - D[0]
        N = self.get_N(D)
        v =  np.sum(N * self.get_V(D)) * dD
        n = np.sum(N) * dD
        if np.isscalar(v):
            v = np.array([v])
            n = np.array([n])
        return v, n

    def get_mom_2(self, T, QM):
        """
        Get second moment from third moment using the best fit provided
        by Field et al (2005)
        Args:
            T: temperature in K
            QM: mass concentration in kg/m3

        Returns:
            The estimated moment of orf_{\text{vol}}^{\text{water}} & =  f_{\text{wet}}^m \frac{\rho^{\text{water}}}{\rho^{m}} \\der 2 of the PSD
        """
        n = 3
        T = T - constants.T0 # Convert to celcius
        a = 5.065339 - 0.062659 * T -3.032362 * n + 0.029469 * T * n \
            - 0.000285 * T**2  + 0.312550 * n**2 + 0.000204 * T**2*n \
            + 0.003199 * T* n**2 - 0.015952 * n**3

        a  = 10**(a)
        b = 0.476221 - 0.015896 * T + 0.165977 * n + 0.007468 * T * n \
            - 0.000141 * T**2 + 0.060366 * n**2 + 0.000079 * T**2 * n \
            + 0.000594 * T * n**2 -0.003577 * n**3

        return (QM/a) ** (1/b)

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

        QM = arg2.astype(np.float64)
        
        # get first estimated N0 and lambda_
        T = arg1
        Q2 = self.get_mom_2(T,QM/constants_1mom.BM_I)
        N0 = Q2**((self.b + 1)/(self.b - 2)) * QM**((2 + 1)/(2 - self.b))
        N0 /= 10**5 # From m^-1 to mm^-1
        lambda_ = (Q2/QM) ** (1/(self.b - 2))

        # Apply correction factor to match third moment
        D = np.linspace(self.d_min, self.d_max, self.nbins_D)
        x = lambda_[:,None] * D.T / 1000
        N = N0[:,None] * constants_1mom.PHI_23_I(x)
        QM_est = np.nansum(self.a * D ** self.b * N, axis = 1) * (D[1]-D[0])
        N0 = N0/QM_est * QM

        # assign the attributes
        self.N0 = N0.T
        self.lambda_ = lambda_.T
        self.ntot = np.nansum(N, axis = 1) *(D[1]-D[0])

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
        ar = 11.3 * D ** 0.414 * 1000 ** (-0.414)
        return 1/ar
