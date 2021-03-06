'''
Description: The radar operator module for computing radar observables
from the interpolated radials given by NWP models
Author: Hejun Xie
Date: 2020-10-12 10:45:48
LastEditors: Hejun Xie
LastEditTime: 2021-06-23 21:32:01
'''

# Global imports
import numpy as np
import time
np.warnings.filterwarnings('ignore')

# Local imports
from ..interp import Radial
from ..hydro import create_hydrometeor
from ..utils import nansum_arr, sum_arr, vlinspace, nan_cumprod, nan_cumsum, aliasing
from . import _core_utils as utils
from ._doppler_spectrum import get_doppler_spectrum

def get_radar_observables_rdop(list_subradials, lut_sz):
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
    global simulate_doppler, simulate_spectrum, doppler_scheme
    global microphysics_scheme

    # Get info from user config
    # Some schme are displayed here but not implemented yet...
    from ..config.cfg import CONFIG
    doppler_scheme = CONFIG['core']['doppler_scheme'] 
    microphysics_scheme = CONFIG['microphysics']['scheme'] # TODO 2-moment scheme
    with_ice_crystals = CONFIG['microphysics']['with_ice_crystals']
    att_corr = CONFIG['microphysics']['with_attenuation']
    radial_res = CONFIG['radar']['radial_resolution']
    integration_scheme = CONFIG['integration']['scheme'] # TODO only Guass-Hermite intergration (scheme 1) 
    nyquist_velocity = CONFIG['radar']['nyquist_velocity']
    simulate_doppler = CONFIG['core']['simulate_doppler']
    simulate_spectrum = CONFIG['core']['doppler_spectrum']

    # No doppler for spaceborne radar
    if CONFIG['radar']['type'] == 'spaceborne':
        simulate_doppler = False
        simulate_spectrum = False

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
    dic_hydro = init_hydro_istc(hydrom_types, lut_sz, microphysics_scheme)
    
    # Initialize integrated scattering matrix, see lut submodule for info about the 12 columns
    sz_integ = np.zeros((n_gates,len(hydrom_types),12), dtype='float32') + np.nan
    
    # Intialize Doppler variables average terminal velocity
    if simulate_doppler:
        rvel_avg = np.zeros((n_gates,), dtype='float32') + np.nan
        total_weight_rvel = np.zeros((n_gates,), dtype='float32')
    # Initialize Full Doppler spectrum
    if simulate_spectrum:
        varray = utils._get_varray()
        doppler_spectrum = np.zeros((n_gates, len(varray)), dtype='float32') + np.nan

    for i, subrad in enumerate(list_subradials): # Loop on subradials (quad pts)

        if simulate_doppler:
            vn_integ = np.zeros((n_gates,), dtype='float32') # Integrated V(D)*N(D)
            n_integ = np.zeros((n_gates,), dtype='float32') # Integrated N(D)
        if simulate_spectrum:
            ah_subrad = np.zeros((n_gates,), dtype='float32') # attenuation coefficient of a sub-radial [dB km-1]

        for j, h in enumerate(hydrom_types): # Loop on hydrometeors
            
            '''
            Part I
            Get the integrated SZ for one hydrometeor and one sub-radial
            '''
            return_pack = one_rad_one_hydro(subrad, h, dic_hydro[h], lut_sz[h], 
                simulate_doppler=simulate_doppler, ngates=n_gates)

            if return_pack is None:
                continue
            else:
                if simulate_doppler:
                    valid_data, sz_psd_integ, vn_psd_integ, n_psd_integ = return_pack
                else:
                    valid_data, sz_psd_integ = return_pack
            
            '''
            Part II
            Add integrated items to sz_integ, vn_integ, n_integ, which counts all hydrometeor types (and
            all subradials for sz_integ but not for vn_integ and n_integ).
            That means the hydrometeor (weighted) averaged falling speed are computed for every sub-radial,
            but the bulk scattering amplitude matrix S and bulk scattering phase matrix Z are computed a 
            set of radials.
            '''
            sz_integ[valid_data, j, :] = nansum_arr(sz_integ[valid_data,j,:], sz_psd_integ * subrad.quad_weight)
            
            if simulate_doppler:
                vn_integ[valid_data] = nansum_arr(vn_integ[valid_data], vn_psd_integ)
                n_integ[valid_data] = nansum_arr(n_integ[valid_data], n_psd_integ)
            
            '''
            Part III
            Get the attenuation coefficient of this sub-radial
            for Full doppler spectrum solver and add that up for all hydrometeor types 
            '''
            if simulate_spectrum and att_corr:
                ah_per_hydro = utils.get_att_h(sz_psd_integ)
                ah_subrad = nansum_arr(ah_subrad, ah_per_hydro)
    
        ########################################################################### (inner loop hydrometeor finished)
        
        """ For every subradial, we get the average fall velocity for
        all hydrometeors and the resulting radial velocity """
        
        if simulate_doppler:
            # Obtain hydrometeor average fall velocity
            v_hydro = vn_integ / n_integ
            # Add density adjustments
            # rho_air in unit [kg*m-3]
            v_hydro *= utils.get_rho_corr(subrad.values['RHO'])

            # Get radial velocity knowing hydrometeor fall speed and U,V,W from model
            theta   = subrad.elev_profile
            phi     = subrad.quad_pt[0]
            proj_wind = utils.proj_vel(subrad.values['U'], subrad.values['V'],
                                subrad.values['W'], v_hydro, theta, phi)

            # Get mask of valid values
            total_weight_rvel = nansum_arr(total_weight_rvel, ~np.isnan(proj_wind) * subrad.quad_weight)
            # Average radial velocity for all sub-beams
            rvel_avg = nansum_arr(rvel_avg, proj_wind * subrad.quad_weight)
        
        """ For every subradial we get the doppler spectrum """
        if simulate_spectrum:
            '''
            subrad_spectrum: doppler spectrum
                unit: [mm6 m-3]
                dimension: (n_gates, len_FFT)
            A_ACCU_H: Accumulated Attenuation in horizontal polarization (two-way)
                unit: [dB]
                dimension: (n_gates) --> (ngates, 1)
            '''
            # exit()
            subrad_spectrum = get_doppler_spectrum(subrad, dic_hydro, lut_sz)

            if att_corr:
                A_ACCU_H =  nan_cumsum( - 2 * ah_subrad * (radial_res / 1000.) )
                subrad_spectrum *= 10 ** (0.1 * A_ACCU_H[:,None])
            
            # two-dimension array is not compatible with nansum_arr 
            doppler_spectrum = np.nansum([doppler_spectrum, subrad_spectrum * subrad.quad_weight], axis=0)
    
    ########################################################################### (outer loop sub-radial finished)
    
    '''
    Here we derive the final quadrature integrated radar observables after
    integrating all scattering properties over hydrometeor types and
    all subradials
    '''

    # some over all hydrometeors (n_valid_gates, n_hydrometeors, 12) --> (n_valid_gates, 12)
    sz_integ = np.nansum(sz_integ, axis=1)
    sz_integ[sz_integ == 0] = np.nan

    # Get radar observables
    ZH, ZV, ZDR, RHOHV, KDP, AH, AV, DELTA_HV = get_pol_from_sz(sz_integ)
    PHIDP = nan_cumsum(2 * KDP) * (radial_res / 1000.) + DELTA_HV
    
    if att_corr:
        A_ACCU_V = nan_cumsum( - 2 * AV * (radial_res / 1000.) )
        A_ACCU_H = nan_cumsum( - 2 * AH * (radial_res / 1000.) )
        # AV and AH are in dB so we need to convert them to linear
        ZV *= 10 ** (0.1 * A_ACCU_V)
        ZH *= 10 ** (0.1 * A_ACCU_H)
        ZDR = ZH / ZV
    
    if simulate_doppler:
        rvel_avg /= total_weight_rvel

        # Apply aliasing if wanted
        if nyquist_velocity != None:
            # Note that the real elevation angle might vary for a given
            # ray due to refraction, but the nyquist velocity is computed based on the first
            # elevation angle
            nyq = nyquist_velocity
            rvel_avg = aliasing(rvel_avg, nyq)
    
    if simulate_spectrum:
        '''
        doppler_spectrum: doppler spectrum
                unit: [mm6 m-3]
                dimension: (n_gates, len_FFT)
        varray: discernible velocity of FFT scheme
                unit [m s-1]
                dimension: (len_FFT)
        rvel_avg, avel_std: radial velocity V and spectrum width W
                unit [m s-1]
                dimension: (n_gates)
        '''
        vd = np.trapz(np.multiply(doppler_spectrum, varray), varray, axis=1)
        d = np.trapz(doppler_spectrum, varray, axis=1)
        rvel_avg = vd / d

        rvel_avg_expand = np.tile(np.expand_dims(rvel_avg, axis=1), (1,len(varray)))
        # varray_expand = np.tile(np.expand_dims(varray, axis=0), (n_gates,1))
        delta_v_square = (varray - rvel_avg_expand)**2
        rvel_var = np.trapz(np.multiply(doppler_spectrum, delta_v_square), varray, axis=1) / d
        rvel_std = np.sqrt(rvel_var)

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
    if simulate_doppler or simulate_spectrum:
        rad_obs['RVEL'] = rvel_avg
    if simulate_spectrum:
        rad_obs['DSPECTRUM'] = doppler_spectrum
        rad_obs['SW'] = rvel_std
    
    '''
    Once averaged , the meaning of the mask is the following
    mask == -1 : all beams are below topography
    mask == 1 : all beams are above NWP model top
    mask > 0 : at least one beam is above NWP model top
    We will keep only gates where no beam is above NWP model top and at least
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

def one_rad_one_hydro(rad, hydro_name, hydro_istc, db, simulate_doppler=True, ngates=-1):
    '''
    Work unit, make the core less procedure-oriented. 
    Integrate the sz of one subradial and one hydrometeor over particle size distribution.
    Args:
        rad: A single radial to integrate sz over psd.
        hydro_name: A single hydrometeor short name like 'R', 'S', 'G' or 'I'. 
        hydro_istc: A single hydrometeor instance to integrate sz over psd.
        db: The database (lookup table of that hydrometeor type).
        simulate_doppler: Boolean flag indicating whether to integrate V(D)N(D) and N(D).
        ngates: Apply it to -1 if you want ngates = rad.dist_profile.
    Returns:
        if simulate_doppler is True, return valid_data, sz_psd_integ, vn_psd_integ, n_psd_integ;
        else if simulate_doppler is True only return valid_data, sz_psd_integ;
        else QM == 0., return None.
    '''

    """
    Part 1 . Get e and T, QM from Radial instances
    Since lookup tables are defined for angles in [0,90], we have to
    check if elevations are larger than 90??, in that case we take
    180-elevation by symmetricity. Also check if angles are smaller
    than 0, in that case, flip sign
    """
    # 1.1 e
    e = rad.elev_profile
    e[e > 90]   = 180 - e[e > 90]
    e[e < 0 ]   = - e[e < 0]

    # 1.2 T    
    T = rad.values['T']

    # 1.3 QM
    QM = rad.values['Q' + hydro_name + '_v'] # Get mass densities
    valid_data = QM > 0
    # speed up if QM == 0.
    if not np.any(valid_data):
        return None

    '''
    Part 2 : Compute the particle size distribution of the particles
    in unit [mm-1 m-3];
    in dimension (n_valid_gates, nbins_D).
    TODO: Only WSM6 (a type of one moment scheme implemented)
    '''

    # 2.1 set psd for hydrometeor instances
    if hydro_name in ['S','I']: # For snow and ice crystals, we need T and QM
        hydro_istc.set_psd(T[valid_data], QM[valid_data])
    else: # For Rain and graupels, we only need mass concentration QM
        hydro_istc.set_psd(QM[valid_data])

    # 2.2 get D axis, and increment dD in [mm], from scattering property database
    if db.type == 'numpy':
        list_D = db.axes[db.axes_names['d']]
    elif db.type == 'xarray':
        list_D = db.get_axis_value('Dmax')
    dD = list_D[1] - list_D[0]
    
    # 2.3 Compute particle size distribution N(D) for all diameters in [mm-1 m-3]
    N = hydro_istc.get_N(list_D)
    
    '''
    Part 2: Query of the scattering database,
    in unit [mm2] for Z and [mm] for S;
    in dimension (n_valid_gates, nbins_D, 12).
    '''
    sz = db.lookup_line(e=e[valid_data], t=T[valid_data])


    '''
    Part 3 : Integrate the SZ coefficients over PSD.
    sz:     
            dimension: (n_valid_gates, nbins_D, 12) 
            unit: Z[mm2]; S[mm]
    N:       
            dimension: (n_valid_gates, nbins_D) 
            unit: [mm-1 m-3]
    dD:      
            dimension: scalar 
            unit: [mm]
    sz_psd_integ: 
            dimension: (n_valid_gates, 12) 
            unit: Z[mm2 m-3]; S[mm m-3]
    '''
    sz_psd_integ = np.einsum('ijk,ij->ik', sz, N) * dD
    
    '''
    Part 4 :
    Pad valid data to make it in the same dimension as radar_range
    Check for special cases where the beam is truncated by model top or topo
    '''
    if ngates == -1:
        ngates = len(rad.dist_profile)
    if len(valid_data) < ngates:
        valid_data = np.pad(valid_data, (0, ngates - len(valid_data)), 
            mode = 'constant', constant_values = False)

    '''
    Part 5 :
    if simulate doppler then:
    doppler_scheme == 1:
        Integrate V(D)*N(D) and N(D) over particle size distribution
    doppler_scheme == 2:
        Integrate V(D)*N(D)*rcs_h(D) and N(D)*rcs_h(D) over particle size distribution

    Get vn_psd_integ, n_psd_integ.
    vn_psd_integ: 
            dimension: (n_valid_gates) 
            unit: [m/s m-3]
    n_psd_integ: 
            dimension: (n_valid_gates) 
            unit: Z[m-3]
    '''
    if simulate_doppler:
        if doppler_scheme == 1:
            vn_psd_integ, n_psd_integ = hydro_istc.integrate_V()
        elif doppler_scheme == 2:
            # Horizontal back-scattering cross section [mm2] 
            rcs = utils.get_rcs_h(sz) # (n_valid_gates, nbins_D)
            # Terminal velocity [m s-1] 
            v_f = hydro_istc.get_V(list_D) # (nbins_D)

            # perform the integration
            vn_psd_integ = np.trapz(np.multiply(v_f, N * rcs), list_D, axis=1)
            n_psd_integ = np.trapz(N * rcs, list_D, axis=1)
    
    '''
    Part 6 : return integrated data
    '''
    if simulate_doppler:
        return valid_data, sz_psd_integ, vn_psd_integ, n_psd_integ
    else:
        return valid_data, sz_psd_integ

def get_pol_from_sz(sz):
    '''
    Computes polarimetric radar observables from integrated scattering properties
    constants.KW: the refractive factor of water, usually 0.93 for radar applications
    Args:
        sz: integrated scattering matrix, with an arbitrary number of rows
            (gates) and 12 columns (see lut submodule)
        

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

    # Horizontal reflectivity
    z_h = utils.get_z_h(sz)

    # Vertical reflectivity
    z_v = utils.get_z_v(sz)

    # Differential reflectivity
    zdr = utils.get_zdr(sz)

    # Differential phase shift
    kdp = utils.get_kdp(sz)

    # Attenuation
    ah = utils.get_att_h(sz)
    av = utils.get_att_v(sz)
    
    # Copolar correlation coeff.
    rho_hv = utils.get_rho_hv(sz)

    # Backscattering differential phase
    delta_hv = utils.get_delta_hv(sz)

    return z_h, z_v, zdr, rho_hv, kdp, ah, av, delta_hv

def init_hydro_istc(hydrom_types, lut_sz, microphysics_scheme):
    '''
    initialize hydrometeor instance
    Args:
        hydrom_types: A list of hydrometeor names.
        lut_sz: A dictionary of database for hydrometeors scattering properties.
        microphysics_scheme: The microphysics scheme of NWP.
    
    Returns:
        A dictionary of hydrometeor instance.
    '''

    dic_hydro = dict()
    
    for h in hydrom_types:
        dic_hydro[h] = create_hydrometeor(h, microphysics_scheme)
        # Add info on number of bins to use for numerical integrations
        # Needs to be the same as in the lookup tables
        if lut_sz[h].type == 'numpy':
            _nbins_D = lut_sz[h].value_table.shape[-2]
            _dmin = lut_sz[h].axes[2][0]
            _dmax = lut_sz[h].axes[2][-1]
        elif lut_sz[h].type == 'xarray':
            _nbins_D = len(lut_sz[h].get_axis_value('Dmax'))
            _dmin = lut_sz[h].get_axis_value('Dmax')[0]
            _dmax = lut_sz[h].get_axis_value('Dmax')[-1]

        dic_hydro[h].nbins_D = _nbins_D
        dic_hydro[h].d_max = _dmax
        dic_hydro[h].d_min = _dmin
    
    return dic_hydro
