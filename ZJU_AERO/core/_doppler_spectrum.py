'''
Description: compute doppler spectrum
Author: Hejun Xie
Date: 2021-03-01 20:25:31
LastEditors: Hejun Xie
LastEditTime: 2021-03-06 21:54:29
'''

import numpy as np
import copy
from . import _core_utils as utils
from .spectrum_integ import get_refl 

def get_doppler_spectrum(subrad, dic_hydro, lut_sz):
    '''
    Computes the reflectivity within every bin of the Doppler spectrum
    Args:
        subrad: a subradial, containing all necessary info
        dic_hydro: a dictionary containing the hydrometeor Class instances
        lut_sz: dictionary containing the lookup tables of scattering
            properties for every hydrometeor type
    Returns:
         refl: array of size [n_gates, len_FFT] containing the reflectivities
             at every range gate and for every velocity bin 
             unit: [mm6 m-3]
    '''

    # Get dimensions
    n_gates = len(subrad.dist_profile) # Number of radar gates

    varray = utils._get_varray()

    # Initialize matrix of reflectivities
    refl = np.zeros((n_gates, len(varray)), dtype='float32')

    # Get all elevations
    elev = subrad.elev_profile
    # Get azimuth angle
    phi = subrad.quad_pt[0]

    # Since lookup tables are defined for angles >0, we have to check
    # if angles are larger than 90°, in that case we take 180-elevation
    # by symmetricity
    elev_lut = copy.deepcopy(elev)
    elev_lut[elev_lut>90] = 180 - elev_lut[elev_lut>90]
    # Also check if angles are smaller than 0, in that case, flip sign
    elev_lut[elev_lut<0] = - elev_lut[elev_lut<0]

    rho_corr = utils.get_rho_corr(subrad.values['RHO'])

    # np.set_printoptions(threshold=np.inf)
    # print(subrad.values['T'][:100])
    # print(subrad.values['QR_v'][:100])
    # print(subrad.values['QS_v'][:100])
    # print(subrad.values['QG_v'][:100])
    
    for i in range(n_gates):
        
        
        # check if the radial is below the topograph or above the model top
        if subrad.mask[i] != 0:
            continue
        
        dic_hydrom_gate = {}
        for j,h in enumerate(dic_hydro.keys()):    
            # set psd for valid hydrometeors at this radar gate
            Q = subrad.values['Q' + h + '_v'][i]
            T = subrad.values['T'][i]
            if Q > 0:
                dic_hydrom_gate[h] = dic_hydro[h]
                if h in ['S', 'I']:
                    dic_hydrom_gate[h].set_psd(np.array([T]), np.array([Q]))
                else:
                    dic_hydrom_gate[h].set_psd(np.array([Q]))

        for j,h in enumerate(dic_hydrom_gate.keys()):
            db = lut_sz[h]
            if db.type == 'numpy':
                list_D = db.axes[db.axes_names['d']]
            elif db.type == 'xarray':
                list_D = db.get_axis_value('Dmax')
            
            # Initialize matrix of radar cross sections, N and D
            n_d_bins = len(list_D)
            rcs     = np.zeros((n_d_bins,), dtype='float32') + np.nan
            N       = np.zeros((n_d_bins,), dtype='float32') + np.nan
            D       = np.zeros((n_d_bins,), dtype='float32') + np.nan

            # Get N and D for all hydrometeors that are present
            D = np.linspace(dic_hydrom_gate[h].d_min,
                            dic_hydrom_gate[h].d_max,
                            dic_hydrom_gate[h].nbins_D)

            D_min   = D[0]
            step_D  = D[1] - D[0]
            N[:]  = dic_hydrom_gate[h].get_N(D)
        
            # Compute RCS for all hydrometeors that are present
            # NOTE: lookup_line() only takes in a list and output a list
            sz = db.lookup_line(e = [elev_lut[i]],
                                t = [subrad.values['T'][i]])[0]
            # get RCS
            rcs[:] = utils.get_rcs_h(sz)

            # Important we use symetrical elevations only for lookup querying, not
            # for actual trigonometrical velocity estimation
            Da, Db, idx = _get_diameter_from_rad_vel(dic_hydrom_gate[h], phi, elev[i],
                            subrad.values['U'][i],
                            subrad.values['V'][i],
                            subrad.values['W'][i],
                            rho_corr[i])
        
            # print(h)
            # print(Da)
            # print(Db)
            # print(idx)
            # print(N)
            # print(rcs)

            try:
                arguments_c_code = (len(idx), Da, Db, rcs, N, step_D, D_min)
                refl_hydro = get_refl(*arguments_c_code)[1]
                # print(refl_hydro)
                refl[i,idx] += refl_hydro
                # print(refl[i,idx])
            except:
                print('An error occured in the Doppler spectrum calculation...')
                raise
            
        # exit()
        
        # Add reflectivity coefficient for that radar gate
        refl[i,:] *= utils._get_refl_coeff()

    return refl
    

def _get_diameter_from_rad_vel(hydro_istc, phi, theta, U, V, W, rho_corr):
    '''
    Retrieves the diameters corresponding to all radial velocity bins, by
    getting the corresponding terminal velocity and inverting the
    diameter-velocity relations
    Args:
        hydro_istc: a hydrometeor Class instances
        phi: the azimuth angle in degrees
        theta: the elevation angle in degrees
        rho_corr: the correction for density [v_true = v_formula * rho_corr]
        U: the eastward wind component of NWP model
        V: the northward wind component of NWP model
        W: the vertical wind component of NWP model
    Returns:
         Da: the diameters corresponding to the left edge of all velocity bins [len_valid_FFT]
         Db: the diameters corresponding to the right edge of all velocity bins [len_valid_FFT]
         idx: the indices of the radar gates where Da and Db are valid values [len_valid_FFT]
             (i.e. between dmin and dmax)
    '''

    varray = utils._get_varray()
    vf_true = utils.proj_vel_back(U, V, W, varray, theta, phi)
    vf_formula = (1. / rho_corr) * vf_true

    # We are only interested in positive fall speeds
    valid_idx = np.where(vf_formula >= 0)[0] 
    vf_formula = vf_formula[valid_idx]

    '''
    D, Da, Db: diamater bins edge of hydrometeors
        dimension: [len_FFT]
        Unit: [mm]
    '''

    D = np.zeros((len(valid_idx),), dtype='float32')

    # Get D bins from V bins
    D[:] = hydro_istc.get_D_from_V(vf_formula)
    # Threshold to valid diameters
    D[D >= hydro_istc.d_max] = hydro_istc.d_max
    D[D <= hydro_istc.d_min] = hydro_istc.d_min

    # we dont know D is ascending or descending on axis [len_FFT]
    # Array of left bin limits
    Da = np.minimum(D[0:-1], D[1:])
    # Array of right bin limits
    Db = np.maximum(D[0:-1], D[1:])

    # Get indice of Dbins where at least one hydrometeor bin width is larger than 0
    mask = np.where((Db-Da) > 0.0)[0]

    return Da[mask], Db[mask], valid_idx[mask]
