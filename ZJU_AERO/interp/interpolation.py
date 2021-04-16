'''
Description: interpolatiom.py: Provides routines for the trilinear interpolation of
model variables to the radar gates
Author: Hejun Xie
Date: 2020-08-15 11:07:01
LastEditors: Hejun Xie
LastEditTime: 2021-04-16 17:39:53
'''

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
from . import Radial, get_all_radar_pts
from ..const import global_constants as constants
from ..utils import nansum_arr, sum_arr
from ..beam import compute_trajectory_radial, compute_trajectory_spaceborne

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
    mask == 1 : all subradials are above model top
    mask > 0 : at least one subradials is above model top
    We will keep only gates where no subradials is above model top and at least
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


def get_interpolated_radial(ds_model, azimuth, elevation, N = None,
                            list_refraction = None):
    '''
    Interpolates a radar radial using a specified quadrature and outputs
    a list of subradials
    Args:
        ds_model: A NWP xarray Dataset to be interpolated
        azimuth: the azimuth angle in degrees (phi) of the radial
        elevation: the elevation angle in degrees (theta) of the radial
        N : if the differential refraction scheme by Zeng and Blahak (2014) is
            used, the refractivity of the atmosphere must be provided as an
            additional model variable
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

    keys = list(ds_model.data_vars.keys())

    # Get options
    from ..config.cfg import CONFIG
    bandwidth_3dB = CONFIG['radar']['3dB_beamwidth']
    integration_scheme = CONFIG['integration']['scheme']
    refraction_method = CONFIG['refraction']['scheme']

    # print('refraction: {}'.format(refraction_method))
    # print('integration: {}'.format(integration_scheme))
    
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

    # Get rrange and coordinates of the radar 
    if CONFIG['radar']['type'] in ['ground']:
        rranges = constants.RANGE_RADAR
        radar_pos = CONFIG['radar']['coords']

    # Get the subradial trajaxtory of the radar
    list_refraction = []

    if CONFIG['radar']['type'] in ['spaceborne']:
        trajectoryListType = 'zen'
        for pt in pts_ver:
            s, h, e = compute_trajectory_spaceborne(pt+elevation)
            list_refraction.append((s, h, e))

    elif CONFIG['radar']['type'] in ['ground']:
        if refraction_method in [1,2]:
            trajectoryListType = 'zen'
            for pt in pts_ver:
                s, h, e = compute_trajectory_radial(
                                                rranges,
                                                pt+elevation,
                                                radar_pos,
                                                refraction_method,
                                                N
                                                )                                    
                list_refraction.append((s, h, e))     
        elif refraction_method in [3]:
            trajectoryListType = 'az,zen'
            for pt_hor in pts_hor:
                s, h, e = compute_trajectory_radial(
                                                    rranges,
                                                    pts_ver + elevation,
                                                    radar_pos,
                                                    refraction_method,
                                                    N,
                                                    pt_hor + azimuth
                                                    )
                list_refraction.append((s, h, e))
                # exit()
    else:
        raise KeyError('No such radar type:{}'.format(CONFIG['radar']['type']))
    

    for i in range(len(pts_hor)):
        for j in range(len(pts_ver)):

            # GH coordinates
            pt = [pts_hor[i]+azimuth, pts_ver[j]+elevation]
            # Interpolate beam
            if trajectoryListType == 'zen':
                lats,lons,list_vars = trilin_interp_radial(ds_model,
                                                    pts_hor[i]+azimuth,
                                                    list_refraction[j][0],
                                                    list_refraction[j][1])
            elif trajectoryListType == 'az,zen':
                lats,lons,list_vars = trilin_interp_radial(ds_model,
                                                    pts_hor[i]+azimuth,
                                                    list_refraction[i][0][j],
                                                    list_refraction[i][1][j])

            weight = weights[i,j]

            # Create dictionary of beams
            dic_beams={}
            # Loop on interpolated variables
            for k, bi in enumerate(list_vars):
                # Do this only for the first variable
                # (same mask for all variables)
                if k == 0:
                    '''
                    mask = 1 : interpolated pt is above model top
                    mask = -1 : intepolated pt is below topography
                    mask = 0 : interpolated pt is ok
                    '''
                    mask_beam = np.zeros((len(bi)))
                    mask_beam[bi == -9999] = 1
                    mask_beam[np.isnan(bi)] = -1
                bi[mask_beam!=0] = np.nan # Assign NaN to all missing data
                dic_beams[keys[k]] = bi # Create dictionary

            if trajectoryListType == 'zen':
                subradial = Radial(dic_beams, mask_beam, lats, lons,
                                        list_refraction[j][0],
                                        list_refraction[j][1],
                                        list_refraction[j][2],
                                        pt, weight)
                list_subradials.append(subradial)
                
            elif trajectoryListType == 'az,zen':
                subradial = Radial(dic_beams, mask_beam, lats, lons,
                                        list_refraction[i][0][j],
                                        list_refraction[i][1][j],
                                        list_refraction[i][2][j],
                                        pt, weight)
                list_subradials.append(subradial)
    
    return list_subradials

def trilin_interp_radial(ds_model, azimuth, distances_profile, heights_profile):
    """
    Interpolates a radar radial using a specified quadrature and outputs
    a list of subradials (all NWP model interface)
    Args:
        ds_model: A NWP xarray Dataset to be interpolated
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
    from ..config.cfg import CONFIG
    radar_pos = CONFIG['radar']['coords']

    if CONFIG['nwp']['name'] == 'grapes':
        from ..nwp.grapes import WGS_to_GRAPES as WGS_to_MODEL
    elif CONFIG['nwp']['name'] == 'wrf':
        from ..nwp.wrf import WGS_to_WRF as WGS_to_MODEL

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

    # Get MODEL heights and topograph from the MODEL dataset
    model_heights = ds_model.coords['z-levels']
    model_topo = ds_model.coords['topograph']

    # Do some transpose jobs
    model_heights = np.transpose(model_heights.data[::-1,...], axes=(0, 2, 1))
    model_topo = np.transpose(model_topo.data, axes=(1, 0))

    # get MODEL projection paramters
    proj_MODEL = ds_model.attrs
    # print(proj_MODEL)
    # print(ds_model['T'].data.shape)
    # exit()

    # Get lower left corner of MODEL domain in local index coordinates
    llc_MODEL=(float(0), float(0))
    llc_MODEL=np.asarray(llc_MODEL).astype('float32')

    # Get upper right corner of MODEL domain in local index coordinates
    urc_MODEL=(float(proj_MODEL['nI'] - 1), float(proj_MODEL['nJ'] - 1))
    urc_MODEL=np.asarray(urc_MODEL).astype('float32')
    
    # Get resolution
    res_MODEL = [1., 1.]

    # Transform radar gate WGS coordinates into MODEL grid coordinates 
    coords_rad_loc = WGS_to_MODEL((lats_rad,lons_rad), proj_MODEL)

    # Check if all points are within MODEL domain
    # TODO: Check earlier results
    if np.any(coords_rad_loc[:,0]<llc_MODEL[0]) or\
        np.any(coords_rad_loc[:,1]<llc_MODEL[1]) or \
            np.any(coords_rad_loc[:,0]>urc_MODEL[0]) or \
                np.any(coords_rad_loc[:,1]>urc_MODEL[1]):
                    msg = """
                    ERROR: RADAR DOMAIN IS NOT ENTIRELY CONTAINED IN NWP MODEL
                    SIMULATION DOMAIN: ABORTING
                    """
                    raise(IndexError(dedent(msg)))

    # Initialize interpolated variables
    interp_data = []

    for i,var in enumerate(list(ds_model.data_vars.keys())):

        # Now we interpolate all variables along beam using C-code file
        ###########################################################################
        rad_interp_values = np.zeros(len(distances_profile),)*float('nan')
        model_data = ds_model[var]

        # do some transpose and reverse to fit the model data format
        model_data = np.transpose(model_data.data[::-1,...], axes=(0, 2, 1))

        arguments_c_code = (len(distances_profile),
                            coords_rad_loc,
                            heights_profile,
                            model_data,
                            model_heights,
                            model_topo,
                            llc_MODEL,
                            res_MODEL)
        
        rad_interp_values = get_all_radar_pts(*arguments_c_code)

        # np.set_printoptions(threshold=30)
        # if var in ['QI_v', 'QR_v', 'QS_v', 'T', 'U']:
        #     print(var)
        #     print(rad_interp_values[1][:])
        #     print(model_data.max())
        #     print(model_data.min())

        interp_data.append(rad_interp_values[1][:])
    

    return lats_rad, lons_rad, interp_data
