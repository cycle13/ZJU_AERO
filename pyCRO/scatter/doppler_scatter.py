'''
Description:  computes all radar observables for a given radial
Author: Hejun Xie
Date: 2020-08-20 22:01:10
LastEditors: Hejun Xie
LastEditTime: 2020-08-22 22:24:44
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global imports
import numpy as np
np.warnings.filterwarnings('ignore')
from scipy.ndimage.filters import gaussian_filter
import copy
from scipy.optimize import fsolve

# Local imports
from pyCRO.config.cfg import CONFIG
from pyCRO.interpolation import Radial
from pyCRO.hydrometeors import create_hydrometeor
from pyCRO.constants import global_constants as constants
from pyCRO.utilities import nansum_arr, sum_arr, vlinspace, nan_cumprod, nan_cumsum, aliasing

def proj_vel(U, V, W, vf, theta,phi):
    """
    Gets the radial velocity from the 3D wind field and hydrometeor
    fall velocity
    Args:
        U: eastward wind component [m/s]
        V: northward wind component [m/s]
        W: vertical wind component [m/s]
        vf: terminal fall velocity averaged over all hydrometeors [m/s]
        theta: elevation angle in degrees
        phi: azimuth angle in degrees

    Returns:
        The radial velocity, with reference to the radar beam
        positive values represent flow away from the radar
    """
    return ((U*np.sin(phi) + V * np.cos(phi)) * np.cos(theta)
            + (W - vf) * np.sin(theta))

def get_radar_observables(list_subradials, lut_sz):
    """
    Computes Doppler and polarimetric radar variables for all subradials
    over ensembles of hydrometeors and integrates them over all subradials at
    the end
    Args:
        list_subradials: list of subradials (Radial claass instances) as
            returned by the interpolation.py code
        lut_sz: Lookup tables for all hydrometeor species as returned by the
            load_all_lut function in the lut submodule

    Returns:
        A radial class instance containing the integrated radar observables
    """

    # Get scheme setup
    global doppler_scheme
    global microphysics_scheme

    # Get info from user config
    # some schme are displayed here but not implemented yet...
    from pyCRO.config.cfg import CONFIG
    doppler_scheme = CONFIG['doppler']['scheme'] # TODO
    microphysics_scheme = CONFIG['microphysics']['scheme'] # TODO only '1mom'
    with_ice_crystals = CONFIG['microphysics']['with_ice_crystals']
    att_corr = CONFIG['microphysics']['with_attenuation']
    radial_res = CONFIG['radar']['radial_resolution']
    # TODO only Guass-Hermite intergration (scheme 1) in interpolation.py
    integration_scheme = CONFIG['integration']['scheme'] 
    KW = CONFIG['radar']['K_squared']
    nyquist_velocity = CONFIG['radar']['nyquist_velocity']

    # Get dimensions of subradials
    num_beams = len(list_subradials) # Number of subradials (quad. pts)
    idx_0 = int(num_beams/2) # Index of central subradial
    # Nb of gates in final integrated radial ( = max length of all subradials)
    n_gates = max([len(l.dist_profile) for l in list_subradials])
    # Get elevation and az of center subradial (for nyquist velocity)
    elev_0 = list_subradials[idx_0].quad_pt[1]
    az_0 = list_subradials[idx_0].quad_pt[0]

    # get hydrometeor types
    hydrom_types = ['R','S','G']
    if with_ice_crystals:
        hydrom_types.extend('I')

    # Initialize
    # Create dictionnary with all hydrometeor Class instances
    dic_hydro = {}
    for h in hydrom_types:
        dic_hydro[h] = create_hydrometeor(h, microphysics_scheme)
        # Add info on number of bins to use for numerical integrations
        # Needs to be the same as in the lookup tables
        _nbins_D = lut_sz[h].value_table.shape[-2]
        _dmin = lut_sz[h].axes[2][0]
        _dmax = lut_sz[h].axes[2][-1]

        dic_hydro[h].nbins_D = _nbins_D
        dic_hydro[h].d_max = _dmax
        dic_hydro[h].d_min = _dmin
    
    # Initialize integrated scattering matrix, see lut submodule for info
    # about the 12 columns
    sz_integ = np.zeros((n_gates,len(hydrom_types),12),
                        dtype = 'float32') + np.nan
    
    # Intialize Doppler variables
    # average terminal velocity
    rvel_avg = np.zeros(n_gates,) + np.nan
    total_weight_rvel = np.zeros(n_gates,)

    ###########################################################################
    for i, subrad in enumerate(list_subradials): # Loop on subradials (quad pts)
        v_integ = np.zeros(n_gates,) # Integrated fall velocity
        n_integ = np.zeros(n_gates,) # Integrated number of particles

        for j, h in enumerate(hydrom_types): # Loop on hydrometeors
            
            """
            Since lookup tables are defined for angles in [0,90], we have to
            check if elevations are larger than 90°, in that case we take
            180-elevation by symmetricity. Also check if angles are smaller
            than 0, in that case, flip sign
            """
            elev_lut = subrad.elev_profile
            elev_lut[elev_lut > 90] = 180 - elev_lut[elev_lut > 90]
            # Also check if angles are smaller than 0, in that case, flip sign
            elev_lut[elev_lut < 0] = - elev_lut[elev_lut < 0]
            
            T = subrad.values['T']

            '''
            Part 1 : Compute the PSD of the particles
            '''

            QM = subrad.values['Q'+h+'_v'] # Get mass densities
            valid_data = QM > 0
            # spped up
            if not np.any(valid_data):
                continue # Skip
            
            # 1 Moment case
            if microphysics_scheme == '1mom':
                if h in ['S','I'] :
                    # For snow and ice crystals, we need T and QM
                    dic_hydro[h].set_psd(T[valid_data], QM[valid_data])
                else: # Rain and graupel
                    dic_hydro[h].set_psd(QM[valid_data])
            
            list_D = lut_sz[h].axes[lut_sz[h].axes_names['d']]
            dD = list_D[1] - list_D[0]
            
            # Compute particle numbers N(D) for all diameters
            N = dic_hydro[h].get_N(list_D)
            
            '''
            Part 2: Query of the scattering Lookup table
            '''
            sz = lut_sz[h].lookup_line(e = elev_lut[valid_data],
                                        t = T[valid_data])

            '''
            Part 3 : Integrate the SZ coefficients over PSD
            '''
            # sz (n_valid_gates, nbins_D, 12)
            # N  (n_valid_gates, nbins_D)
            sz_psd_integ = np.einsum('ijk,ij->ik',sz,N) * dD
            
            # Check for special cases where the beam is truncated by model top or topo
            if len(valid_data) < n_gates:
                valid_data = np.pad(valid_data,(0,n_gates - len(valid_data)),
                                    mode = 'constant',
                                    constant_values = False)
            
            sz_integ[valid_data,j,:] = nansum_arr(sz_integ[valid_data,j,:],
                                                      sz_psd_integ *
                                                      subrad.quad_weight)
            
            '''
            Part 4 : Doppler
            '''
            # Get terminal velocity integrated over PSD
            vh,n = dic_hydro[h].integrate_V()

            v_integ[valid_data] = nansum_arr(v_integ[valid_data],vh)
            n_integ[valid_data] = nansum_arr(n_integ[valid_data],n)
            
        ########################################################################### (inner loop hydrometeor finished)
        
        """ For every beam, we get the average fall velocity for
        all hydrometeors and the resulting radial velocity """
        
        # Obtain hydrometeor average fall velocity
        v_hydro = v_integ/n_integ
        # Add density weighting
        v_hydro*(subrad.values['RHO']/subrad.values['RHO'][0])**(0.5)

        # Get radial velocity knowing hydrometeor fall speed and U,V,W from model
        theta_deg = subrad.elev_profile
        phi_deg = subrad.quad_pt[0]
        theta = np.deg2rad(theta_deg) # elevation
        phi = np.deg2rad(phi_deg)       # azimuth
        proj_wind = proj_vel(subrad.values['U'],subrad.values['V'],
                            subrad.values['W'], v_hydro, theta, phi)

        # Get mask of valid values
        total_weight_rvel = sum_arr(total_weight_rvel, ~np.isnan(proj_wind)*subrad.quad_weight)
        # Average radial velocity for all sub-beams
        rvel_avg = nansum_arr(rvel_avg, (proj_wind) * subrad.quad_weight)
    
    ########################################################################### (outer loop subbeam finished)
    
    '''
    Here we derive the final quadrature integrated radar observables after
    integrating all scattering properties over hydrometeor types and
    all subradials
    '''

    # some over all hydrometeors (n_valid_gates, n_hydrometeors, 12) --> (n_valid_gates, 12)
    sz_integ = np.nansum(sz_integ,axis=1)
    sz_integ[sz_integ == 0] = np.nan

    # Get radar observables
    ZH, ZV, ZDR, RHOHV, KDP, AH, AV, DELTA_HV = get_pol_from_sz(sz_integ, KW)

    PHIDP = nan_cumsum(2 * KDP) * radial_res/1000. + DELTA_HV

    if att_corr:
        # AH and AV are in dB so we need to convert them to linear
        ZV *= nan_cumprod(10**(-0.1*AV*(radial_res/1000.))) # divide to get dist in km
        ZH *= nan_cumprod(10**(-0.1*AH*(radial_res/1000.)))
        ZDR = ZH / ZV
    
    rvel_avg /= total_weight_rvel

    # Apply aliasing if wanted
    if nyquist_velocity != None:
        # Note that the real elevation angle might vary for a given
        # ray due to refraction, but the nyquist velocity is computed based on the first
        # elevation angle
        nyq = nyquist_velocity
        rvel_avg = aliasing(rvel_avg, nyq)

    ###########################################################################

    '''
    Create the final Radial class instance containing all radar observables
    '''

    # Create outputs
    rad_obs = {}
    rad_obs['ZH'] = ZH
    rad_obs['ZDR'] = ZDR
    rad_obs['ZV'] = ZV
    rad_obs['KDP'] = KDP
    rad_obs['DELTA_HV'] = DELTA_HV
    rad_obs['PHIDP'] = PHIDP
    rad_obs['RHOHV'] = RHOHV
    # Add attenuation at every gate
    rad_obs['ATT_H'] = AH
    rad_obs['ATT_V'] = AV
    # doppler
    rad_obs['RVEL'] = rvel_avg
    
    '''
    Once averaged , the meaning of the mask is the following
    mask == -1 : all beams are below topography
    mask == 1 : all beams are above COSMO top
    mask > 0 : at least one beam is above COSMO top
    We will keep only gates where no beam is above COSMO top and at least
    one beam is above topography
    '''

    # Sum the mask of all beams to get overall average mask
    mask = np.zeros(n_gates,)
    for i,beam in enumerate(list_subradials):
        mask = sum_arr(mask,beam.mask[0:n_gates], cst = 1) # Get mask of every Beam

    mask /= float(num_beams)
    mask[np.logical_and(mask > -1, mask <= 0)] = 0

    # Finally get vectors of distances, height and lat/lon at the central beam
    heights_radar = list_subradials[idx_0].heights_profile
    distances_radar = list_subradials[idx_0].dist_profile
    lats = list_subradials[idx_0].lats_profile
    lons = list_subradials[idx_0].lons_profile

    # Create final radial
    radar_radial = Radial(rad_obs, mask, lats, lons, distances_radar,
                          heights_radar)

    return radar_radial

def get_pol_from_sz(sz, KW):
    '''
    Computes polarimetric radar observables from integrated scattering properties
    Args:
        sz: integrated scattering matrix, with an arbitrary number of rows
            (gates) and 12 columns (seet lut submodule)
        KW: the refractive factor of water, usually 0.93 for radar applications

    Returns:
         z_h: radar refl. factor at hor. pol. in linear units [mm6 m-3]
         z_v: radar refl. factor at vert. pol. in linear units [mm6 m-3]
         zdr: diff. refl. = z_h / z_v [-]
         rhohv: copolar. corr. coeff [-]
         kdp: spec. diff. phase shift upon propagation [deg km-1]
         ah: spec. att. at hor. pol. [dB km-1]
         av: spec. att. at vert. pol. [dB km-1]
         delta_hv: total phase shift upon backscattering [deg]
    '''

    wavelength = constants.WAVELENGTH

    from pyCRO.config.cfg import CONFIG
    K_squared = CONFIG['radar']['K_squared']

    # Horizontal reflectivity
    radar_xsect_h = 2*np.pi*(sz[:,0]-sz[:,1]-sz[:,2]+sz[:,3])
    z_h = wavelength**4/(np.pi**5*K_squared)*radar_xsect_h

    # Vertical reflectivity
    radar_xsect_v = 2*np.pi*(sz[:,0]+sz[:,1]+sz[:,2]+sz[:,3])
    z_v = wavelength**4/(np.pi**5*K_squared)*radar_xsect_v

    # Differential reflectivity
    zdr = radar_xsect_h/radar_xsect_v

    # Differential phase shift
    kdp = 1e-3 * (180.0/np.pi) * wavelength * (sz[:,10]-sz[:,8])

    # Attenuation
    ext_xsect_h = 2 * wavelength * sz[:,11]
    ext_xsect_v = 2 * wavelength * sz[:,9]
    ah = 4.343e-3 * ext_xsect_h
    av = 4.343e-3 * ext_xsect_v

    # Copolar correlation coeff.
    a = (sz[:,4] + sz[:,7])**2 + (sz[:,6] - sz[:,5])**2
    b = (sz[:,0] - sz[:,1] - sz[:,2] + sz[:,3])
    c = (sz[:,0] + sz[:,1] + sz[:,2] + sz[:,3])
    rhohv = np.sqrt(a / (b*c))

    # Backscattering differential phase
    delta_hv = np.arctan2(sz[:,5] - sz[:,6], -sz[:,4] - sz[:,7])

    return z_h,z_v,zdr,rhohv,kdp,ah,av,delta_hv

def cut_at_sensitivity(list_subradials):
    '''
    Censors simulated measurements where the reflectivity falls below the
    sensitivity specified by the user, see the wiki for how to define
    the sensitivity in the configuration files
    Args:
        list_subradials: a list of subradials containing the computed radar
            observables
    Returns:
         the list_subradials but censored with the radar sensitivity
    '''
    from pyCRO.config.cfg import CONFIG
    sens_config = CONFIG['radar']['sensitivity']
    if not isinstance(sens_config,list):
        sens_config = [sens_config]

    if len(sens_config) == 3: # Sensitivity - gain - snr
        threshold_func = lambda r: (sens_config[0] + constants.RADAR_CONSTANT_DB
                         + sens_config[2] + 20*np.log10(r/1000.))
    elif len(sens_config) == 2: # ZH - range
        threshold_func = lambda r: ((sens_config[0] -
                                     20 * np.log10(sens_config[1] / 1000.)) +
                                     20 * np.log10(r / 1000.))
    elif len(sens_config) == 1: # ZH
        threshold_func = lambda r: sens_config[0]

    else:
        print('Sensitivity parameters are invalid, cannot cut at specified sensitivity')
        print('Enter either a single value of refl (dBZ) or a list of',
              '[refl (dBZ), distance (m)] or a list of [sensitivity (dBm), gain (dBm) and snr (dB)]')

        return list_subradials

    if isinstance(list_subradials[0],list): # Loop on list of lists
        for i,sweep in enumerate(list_subradials):
            for j,b, in enumerate(sweep):
                rranges = (CONFIG['radar']['radial_resolution'] *
                           np.arange(len(b.dist_profile)))
                mask = 10*np.log10(b.values['ZH']) < threshold_func(rranges)
                for k in b.values.keys():
                    if k in constants.SIMULATED_VARIABLES:
                        if k == 'DSPECTRUM':
                            logspectrum = 10 * np.log10(list_subradials[i][j].values[k])
                            thresh = threshold_func(rranges)
                            thresh =  np.tile(thresh,
                                              (logspectrum.shape[1],1)).T
                            list_subradials[i][j].values[k][logspectrum < thresh] = np.nan
                        else:
                            list_subradials[i][j].values[k][mask] = np.nan

    else:
        for i, subradial in enumerate(list_subradials): # Loop on simple list
            rranges = (CONFIG['radar']['radial_resolution'] *
                       np.arange(len(subradial.dist_profile)))
            mask = 10 * np.log10(subradial.values['ZH']) < threshold_func(rranges)
            for k in subradial.values.keys():
                if k in constants.SIMULATED_VARIABLES:
                    list_subradials[i].values[k][mask] = np.nan
    return list_subradials
