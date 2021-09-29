'''
Description: Hydrometeor base class,
should not be initialized directly
Author: Hejun Xie
Date: 2020-11-13 13:04:15
LastEditors: Hejun Xie
LastEditTime: 2021-09-29 16:06:00
'''

# Global import
import numpy as np
np.seterr(divide='ignore')
from scipy import stats
from scipy import optimize

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
        self.list_D = None

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

        if ~np.isnan(self.vel_factor) and ~np.isnan(self.ntot_factor):
            v = (self.vel_factor * self.N0 * self.alpha / self.nu * \
                self.lambda_ ** (-(self.beta + self.mu + 1) / self.nu))
            
            n = self.ntot_factor * self.N0 / self.nu * \
                self.lambda_ ** (-(self.mu + 1) / self.nu)
        else:
            # Particle size distribution [mm-1 m-3]
            N = self.get_N(self.list_D)
            # Terminal velocity [m s-1] 
            v_f = self.get_V(self.list_D)

            # perform the integration
            v = np.trapz(np.multiply(v_f, N), self.list_D, axis=1)
            n = np.trapz(N, self.list_D, axis=1)
        
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
    
class _NonsphericalHydrometeor(_Hydrometeor):
    '''
    Base class for nonspherical hydrometeor.
    Complicated Definition of the m-D, vt-D relation of nonspherical particles can be found here
    Here three features of the particles are preserved as compared with spherical hydrometeor:
    1. The particle size distribution (PSD)
    2. The aspect ratio distribution (ARD) observed by particle camera
    3. The total water content 
    '''

    def __init__(self, scheme):
        """
            Args:
            scheme: microphysical scheme to use, can be either '1mom' (operational
               one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)

            Returns:
                A nonspherical Hydrometeor class instance (see below)
        """
        super(_NonsphericalHydrometeor, self).__init__(scheme)
        self.a = np.nan
        self.b = np.nan
        self.ntot_factor = np.nan
        self.vel_factor = np.nan
        self.list_D = None
        self.asp_wgt = None
        self.shape = None

    def get_list_asp(self):
        """
        Return the aspect ratio grids of a particle
        Defined by a specific kind of hydrometeor

        Returns:
            list_asp: A list of aspect ratio (asp > 1.0)
        """
        pass

    def _get_volumn(self, D, asp):
        '''
        Compute the volumn of a nonspherical particle

        Args:
            D: maximum dimension of a particle
            asp: aspect ratio (asp > 1.0)

        Returns:
            Volumn of a particle [mm3]
        '''

        if self.shape == 'hexcol':
            L = D / np.sqrt(asp**2 + 1)
            a = (asp * L) / 2
            V = L * 3 / 2 * np.sqrt(3) * a**2
        elif self.shape == 'spheroid':
            V = (np.pi / 6) * D ** 3 / asp
        elif self.shape == 'snowflake':
            L = D / np.sqrt(asp**2 + 1)
            a = (asp * L) / 2
            n1 = 0.25
            n2 = n1 * self.param # n2 / n1 range from 8.2 to 20 (plate to snowflake)

            deg2rad = np.pi / 180.
            dph = 1.0     # in degree
            nsym = 6      # 6-fold symmetry
            ph = np.arange(0.0, 360.0 / nsym, dph) * deg2rad # in degree
            ratio = ( 2.0**(n2/2.0 - 1.0) * (np.abs(np.cos(1.5*ph))**n2 + np.abs(np.sin(1.5*ph))**n2) ) \
                        ** (-1.0/(2.0*n1))
            
            if np.isscalar(asp):
                rr = a * ratio
            else:
                rr = np.outer(a, ratio)

            ss = 0.5 * rr * rr * dph * deg2rad
            if np.isscalar(asp):
                s = np.sum(ss) * nsym
            else:
                s = np.sum(ss, axis=1) * nsym
            
            V = s * L
        
        return V

    def _get_mass(self, D, asp):
        '''
        Compute the mass of a nonspherical particle

        Args:
            D: maximum dimension of a particle
            asp: aspect ratio (asp > 1.0)

        Returns:
            Mass of a particle
        need to be initialized
        '''
        pass
    
    def get_asp_wgt(self, list_D):
        """
        Return the aspect ratio distribution of a particle
        Args:
            list_D: list of maximum dimension of a particle [mm]
        
        Returns:
            asp_wgt: Aspect ratio distribution weight 
                dimension: (nD, nAR)
        """

        if np.isscalar(list_D):
            list_D = np.array([list_D], dtype='float32')

        list_asp = self.get_list_asp()
    
        ar_lambda, ar_loc, ar_mu = self.get_aspect_ratio_pdf_masc(list_D)
        asp_wgt = np.empty((len(list_D), len(list_asp)), dtype='float32')

        for l in zip(range(len(list_D)), ar_lambda, ar_loc, ar_mu):
            gamm = stats.gamma(l[1],l[2],l[3])
            wei = gamm.pdf(list_asp)
            wei /= np.sum(wei) # renormalization
            asp_wgt[l[0], :] = wei

        return asp_wgt

    def get_M(self, list_D):
        """
        Returns the mass of a particle
        Args:
            list_D: list of maximum dimension of a particle [mm]

        Returns:
            m: the averaged particle masses of aspect ratio distribution, 
            same dimensions as D [kg]
        """

        if np.isscalar(list_D):
            list_D = np.array([list_D], dtype='float32')
        
        list_asp = self.get_list_asp() # (nAR)
        asp_wgt = self.asp_wgt # (nD, nAR)

        M = np.zeros((len(list_D)), dtype='float32') # [kg]

        for iD, D in enumerate(list_D):
            m = self._get_mass(D, list_asp)
            M[iD] = np.sum(m * asp_wgt[iD, :])
        
        return M

    def solve_lambda(self, QM, N0):
        '''
        Solve lambda of exponential particle size distribution, 
        N(D) = N0 * exp(- lambda * D)
        according to given QM and N0, and M-D relationship,
        by zero-finding method of Newton optimization

        Args:
            QM: mass content # [kg m^-3]
            N0: intercept of exponential particle size distribution # [mm^-1 m^-3]

        Returns:
            Lambda: the slope of exponential particle size ditribution [mm-1]
        '''

        list_D = self.list_D
        dD = list_D[1] - list_D[0]

        M = self.M # (nD)

        def QM_P0(x):
            return np.sum(np.exp(-x*list_D) * M) * N0 * dD - QM
        
        def QM_P1(x):
            return - np.sum(np.exp(-x*list_D) * M * list_D) * N0 * dD
        
        def QM_P2(x):
            return np.sum(np.exp(-x*list_D) * M * list_D**2) * N0 * dD
        try:
            _lambda = optimize.newton(QM_P0, 1.0, fprime=QM_P1, fprime2=QM_P2, tol=1e-4)
        except RuntimeError:
            print(QM)

        return _lambda
