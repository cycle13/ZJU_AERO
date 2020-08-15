'''
Description: interpolatiom.py: Provides routines for the trilinear interpolation of
COSMO variables to the radar gates
Author: Hejun Xie
Date: 2020-08-15 11:07:01
LastEditors: Hejun Xie
LastEditTime: 2020-08-15 16:56:31
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global imports
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error
import pyproj
import pickle
import pyWRF as pw
import scipy.interpolate as interp
from scipy.ndimage import gaussian_filter
from textwrap import dedent

# Local imports
from pyCRO.interpolation import Radial, get_all_radar_pts
from pyCRO.constants import global_constants as constants
from pyCRO.utilities import nansum_arr, sum_arr
from pyCRO.beam_propogation import compute_trajectory_radial

def integrate_radials(list_subradials):
    '''
    Integrates a set of radials corresponding to different quadrature points
    Args:
        list_subradials: list of Radial class instances corresponding to
            different subradials

    Returns:
        integrated_radial: an integrated Radial class instance

    '''

    num_subradials = len(list_subradials)
    list_variables = list_subradials[0].values.keys()

    integrated_variables = {}
    for k in list_variables:
        integrated_variables[k] = np.array([np.nan])
        sum_weights=0
        for i in list_subradials:
            sum_weights += i.quad_weight
        for i in list_subradials:
            integrated_variables[k] = nansum_arr(integrated_variables[k],
                                i.values[k] * i.quad_weight / sum_weights)

    # Get index of central beam
    idx_0 = int(num_subradials / 2.)

    # Sum the mask of all beams to get overall average mask
    mask = np.zeros(num_subradials,)
    for i,p in enumerate(list_subradials):
        mask = sum_arr(mask, p.mask)

    '''
    Once averaged , the meaning of the mask is the following
    mask == -1 : all subradials are below topography
    mask == 1 : all subradials are above COSMO top
    mask > 0 : at least one subradials is above COSMO top
    We will keep only gates where no subradials is above COSMO top and at least
    one subradials is above topography
    '''
    mask /= float(num_subradials)
    mask[np.logical_and(mask > -1, mask <= 0)] = 0

    heights_radar = list_subradials[idx_0].heights_profile
    distances_radar = list_subradials[idx_0].dist_profile
    lats = list_subradials[idx_0].lats_profile
    lons = list_subradials[idx_0].lons_profile

    # Create new beam with summed up variables
    integrated_radial = Radial(integrated_variables, mask, lats, lons,
                             distances_radar, heights_radar)

    return integrated_radial


def get_interpolated_radial(dic_variables, azimuth, elevation, N = None,
                            list_refraction = None):
    '''
    Interpolates a radar radial using a specified quadrature and outputs
    a list of subradials
    Args:
        dic_variables: dictionary containing the COSMO variables to be
            interpolated
        azimuth: the azimuth angle in degrees (phi) of the radial
        elevation: the elevation angle in degrees (theta) of the radial
        N : if the differential refraction scheme by Zeng and Blahak (2014) is
            used, the refractivity of the atmosphere must be provided as an
            additional COSMO variable
        list_refraction : To save time, a list of (s,h,e) tuples corresponding
            to the dist at ground, height above ground and incident elev. ang.
            for all quadrature points (output of atmospheric refraction)
            can be provided, in which case the atmospheric refraction will
            not be recalculated. This should be done only if the elevation
            angle is the same from one interpolated radial to the other
            (PPI). Also this is not possible for quadrature schemes 2 and 6
            which have irregular grids.
    Returns:
        list_subradials: a list of Radial class instances containing all
            subradials corresponding to all quadrature points
                         defined along the specified radial
    '''

    list_variables = dic_variables.values()
    keys = dic_variables.keys()

    # Get options
    from cosmo_pol.config.cfg import CONFIG
    bandwidth_3dB = CONFIG['radar']['3dB_beamwidth']
    integration_scheme = CONFIG['integration']['scheme']
    refraction_method = CONFIG['refraction']['scheme']
    
    # Calculate quadrature weights
    if integration_scheme == 1: # Classical single gaussian scheme
        nh_GH = int(CONFIG['integration']['nh_GH'])
        nv_GH = int(CONFIG['integration']['nv_GH'])

        # Get GH points and weights
        sigma = bandwidth_3dB/(2*np.sqrt(2*np.log(2)))

        pts_hor, weights_hor=np.polynomial.hermite.hermgauss(nh_GH)
        pts_hor = pts_hor*sigma

        pts_ver, weights_ver=np.polynomial.hermite.hermgauss(nv_GH)
        pts_ver = pts_ver*sigma

        weights = np.outer(weights_hor*sigma,weights_ver*sigma)
        weights *= np.abs(np.cos(np.deg2rad(pts_ver)))
        sum_weights = np.sum(weights.ravel())
        weights /= sum_weights # Normalize weights

        beam_broadening=nh_GH>1 or nv_GH>1 # Boolean for beam-broadening (if only one GH point : No beam-broadening)
    
    
    # Initialize list of subradials
    list_subradials = []

    # Get rrange of the radar
    rranges = constants.RANGE_RADAR

    # Get the subradial trajaxtory of the radar
    list_refraction = []

    # Get coordinates of virtual radar
    radar_pos = CONFIG['radar']['coords']

    for pt in pts_ver:
        s, h, e = compute_trajectory_radial(
                                        rranges,
                                        pt+elevation,
                                        radar_pos,
                                        refraction_method,
                                        N
                                        )
                                        
        list_refraction.append((s, h, e))
    
    for i in range(len(pts_hor)):
        for j in range(len(pts_ver)):

            # GH coordinates
            pt = [pts_hor[i]+azimuth, pts_ver[j]+elevation]
            # Interpolate beam
            lats,lons,list_vars = trilin_interp_radial_WRF(list_variables,
                                                pts_hor[i]+azimuth,
                                                list_refraction[j][0],
                                                list_refraction[j][1])

            weight = weights[i,j]

            # Create dictionary of beams
            dic_beams={}
            # Loop on interpolated variables
            for k, bi in enumerate(list_vars):
                # Do this only for the first variable
                # (same mask for all variables)
                if k == 0:
                    '''
                    mask = 1 : interpolated pt is above COSMO top
                    mask = -1 : intepolated pt is below topography
                    mask = 0 : interpolated pt is ok
                    '''
                    mask_beam = np.zeros((len(bi)))
                    mask_beam[bi == -9999] = 1
                    mask_beam[np.isnan(bi)] = -1
                bi[mask_beam!=0] = np.nan # Assign NaN to all missing data
                dic_beams[keys[k]] = bi # Create dictionary

            subradial = Radial(dic_beams, mask_beam, lats, lons,
                                    list_refraction[j][0],
                                    list_refraction[j][1],
                                    list_refraction[j][2],
                                    pt, weight)

            list_subradials.append(subradial)
    
    return list_subradials


def trilin_interp_radial_WRF(list_vars, azimuth, distances_profile, heights_profile):
    """
    Interpolates a radar radial using a specified quadrature and outputs
    a list of subradials (WRF interface)
    Args:
        list_vars: list of WRF variables to be interpolated
        azimuth: the azimuth angle in degrees (phi) of the subradial
        distances_profile: vector of distances in meters of all gates
            along the subradial (computed with the atmospheric refraction
            scheme)
        heights_profile: vector of heights above ground in meters of all
            gates along the subradial (computed with the atmospheric refraction
            scheme)
    Returns:
        lats_rad: vector of all latitudes along the subradial
        lons_rad: vector of all longitudes along the subradial
        interp_data: dictionary containing all interpolated variables along
            the subradial
    """

    # Get position of virtual radar from user configuration
    from cosmo_pol.config.cfg import CONFIG
    radar_pos = CONFIG['radar']['coords']


    # Initialize WGS84 geoid
    g = pyproj.Geod(ellps='WGS84')

    # Get radar bins coordinates
    lons_rad=[]
    lats_rad=[]
    # Using the distance on ground of every radar gate, we get its latlon coordinates
    for d in distances_profile:
        # Note that pyproj uses lon/lat whereas I used lat/lon
        lon,lat,ang=g.fwd(radar_pos[1], radar_pos[0], azimuth,d)
        lons_rad.append(lon)
        lats_rad.append(lat)

    # Convert to numpy array
    lons_rad = np.array(lons_rad)
    lats_rad = np.array(lats_rad)

    # Initialize interpolated variables
    interp_data = []
    isbot_rad = np.ones(len(distances_profile), dtype=bool)
    istop_rad = np.ones(len(distances_profile), dtype=bool)

    for i,var in enumerate(list_vars):

        # Get model heights and WRF proj from that variable
        model_heights = var.attributes['z-levels']
        model_topo = var.attributes['topograph']

        # WRF index coordinates on lambert projection face
        proj_WRF=var.attributes['proj_info']

        # Get lower left corner of WRF domain in local index coordinates
        llc_WRF=(float(0), float(0))
        llc_WRF=np.asarray(llc_WRF).astype('float32')

        # Get upper right corner of WRF domain in local index coordinates
        urc_WRF=(float(proj_WRF['nI'] - 1), float(proj_WRF['nJ'] - 1))
        urc_WRF=np.asarray(urc_WRF).astype('float32')
        
        # Get resolution
        res_WRF = [1., 1.]

        # Transform radar gate coordinates into local wrf coordinates
        coords_rad_loc = pw.WGS_to_WRF((lats_rad,lons_rad), proj_WRF)

        # Check if all points are within WRF domain
        if np.any(coords_rad_loc[:,1]<llc_WRF[0]) or\
            np.any(coords_rad_loc[:,0]<llc_WRF[1]) or \
                np.any(coords_rad_loc[:,1]>urc_WRF[0]) or \
                    np.any(coords_rad_loc[:,0]>urc_WRF[1]):
                        msg = """
                        ERROR: RADAR DOMAIN IS NOT ENTIRELY CONTAINED IN WRF
                        SIMULATION DOMAIN: ABORTING
                        """
                        raise(IndexError(dedent(msg)))

        # Now we interpolate all variables along beam using C-code file
        ###########################################################################
        rad_interp_values = np.zeros(len(distances_profile),)*float('nan')
        model_data = var.data

        # do some transpose and reverse to fit the cosmo data format
        model_heights = np.transpose(model_heights[::-1,...], axes=(0, 2, 1))
        model_data = np.transpose(model_data[::-1,...], axes=(0, 2, 1))
        model_topo = np.transpose(model_topo, axes=(1, 0))

        arguments_c_code = (len(distances_profile),
                            coords_rad_loc,
                            heights_profile,
                            model_data,
                            model_heights,
                            model_topo,
                            llc_WRF,
                            res_WRF)
        
        rad_interp_values = get_all_radar_pts(*arguments_c_code)

        # if var.name == 'QI_v':
        #     print('QI_v')
        #     print(rad_interp_values[1][:])
        # if var.name == 'QR_v':
        #     print('QR_v')
        #     print(rad_interp_values[1][:])

        interp_data.append(rad_interp_values[1][:])

        isbot_rad &= ~np.isnan(rad_interp_values[1][:])
        istop_rad &= ~(rad_interp_values[1][:]==-9999)

    # remove those partially valid radar gates
    for i,var in enumerate(list_vars):
        interp_data[i][~isbot_rad] = np.nan
        interp_data[i][~istop_rad] = -9999

    return lats_rad, lons_rad, interp_data

# some unit tests
if __name__ == "__main__":
    pass
