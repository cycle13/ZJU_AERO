'''
Description: Provides routines to convert radar operator outputs
to the format used in spaceborne files, as well as routines to retrieve
all azimuth and elevation angles from a spaceborne file, which needed before
simulating the scan.
(only for spaceborne radar pre-research purpose now...)
Author: Hejun Xie
Date: 2020-10-09 15:56:38
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 12:25:53
'''

# Global imports
import sys
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error
import pyproj as pp
import copy

# Local imports
from ..const import global_constants as constants
from ..utils import vlinspace


def get_spaceborne_angles(swath_file, swath_slice=None):
    '''
    From as specified spaceborne radar product file, gets the azimuth (phi), elevation
    (theta) angles and range bins for every radar radial
    (there is one radial for every measurement at the ground),
    based on the position of the satellite given in the file
    Args:
        swath_filename: path of the corresponding spaceborne radar file
        swath_slice: slice of the spaceborne radar pixels in model domain (should be a 2D slice)
    Returns:
        azimuths: all azimuths in the form of a 2D array [degrees]
        elevations: all elevation angles in the form of a 2D array [degrees]
        ranges: all range in the form of a 2D array [m]
    
    !!Attention: Now it is only used for pre-research tests.
    '''
    # Initialize geoid for inverse distance computations
    geoid = pp.Geod(ellps = 'WGS84')

    # Projection from lat/long/alt to eced
    ecef = pp.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pp.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    if sys.version_info[0] >= 3:
        transformer = pp.Transformer.from_proj(lla, ecef)

    # get lat_2D, lon_2D, center_lon_sc, center_lat_sc, pos_sc by makeup or 
    # from real spaceborne radar products (not implemented yet...)
    if swath_file == 'test':
        (nscan, npixel) = (100, 100)
        lat_2D = np.empty((nscan, npixel), dtype=float)
        lon_2D = np.empty((nscan, npixel), dtype=float)
        center_lon_sc = np.empty((nscan), dtype=float)
        center_lat_sc = np.empty((nscan), dtype=float)
        pos_sc = np.empty((nscan, 3), dtype=float)
        
        # now make up the test swath info
        sc_alt_makeup = 407000 # [m] = 407km PolarObiting height
        lon_frame_makeup = [114.0, 117.5, 115.0, 118.5] # [NW, NE, SW, SE]
        lat_frame_makeup = [42.,   42.,  38.0,  38.0] # [NW, NE, SW, SE]

        sc_alt = np.empty((nscan), dtype=float)
        sc_alt[:] = sc_alt_makeup #[m]

        center_lon_sc = np.linspace(0.5*(lon_frame_makeup[0]+lon_frame_makeup[1]),
        0.5*(lon_frame_makeup[2]+lon_frame_makeup[3]), nscan)
        center_lat_sc = np.linspace(0.5*(lat_frame_makeup[0]+lat_frame_makeup[1]),
        0.5*(lat_frame_makeup[2]+lat_frame_makeup[3]), nscan)

        # (nscan, npixels)
        lon_2D = vlinspace(
            np.linspace(lon_frame_makeup[0], lon_frame_makeup[1], npixel),
            np.linspace(lon_frame_makeup[2], lon_frame_makeup[3], npixel),
            nscan
        ).transpose()
        lat_2D = vlinspace(
            np.linspace(lat_frame_makeup[0], lat_frame_makeup[1], npixel),
            np.linspace(lat_frame_makeup[2], lat_frame_makeup[3], npixel),
            nscan
        ).transpose()

        # get pos_sc by pyproj transformation
        for iscan in range(nscan):
            if sys.version_info[0] >= 3:
                sc_x,sc_y,sc_z = transformer.transform(center_lon_sc[iscan], center_lat_sc[iscan], sc_alt[iscan])
            else:
                sc_x,sc_y,sc_z = pp.transform(lla, ecef, center_lon_sc[iscan], center_lat_sc[iscan], sc_alt[iscan])
            pos_sc[iscan, 0], pos_sc[iscan, 1], pos_sc[iscan, 2] = sc_x, sc_y, sc_z

    else:
        raise NotImplementedError('Swath file not implemented yet')

    azimuths = np.zeros(lat_2D.shape)
    ranges = np.zeros(lat_2D.shape)
    elevations = np.zeros(lon_2D.shape)

    # make full slice
    if swath_slice is not None:
        lon_2D = lon_2D[swath_slice]
        lat_2D = lat_2D[swath_slice]
        center_lon_sc = center_lon_sc[swath_slice[0]]
        center_lat_sc = center_lat_sc[swath_slice[0]]
        pos_sc = pos_sc[swath_slice[0]]
    
    [nscan,npixel]=lat_2D.shape

    for i in range(nscan):
        for j in range(npixel):

            a,b,d = geoid.inv(center_lon_sc[i], center_lat_sc[i],
                              lon_2D[i,j], lat_2D[i,j])
            azimuths[i,j]=a

            if sys.version_info[0] >= 3:
                surf_x,surf_y,surf_z = transformer.transform(lon_2D[i,j], lat_2D[i,j], 0)
            else:
                surf_x,surf_y,surf_z = pp.transform(lla,ecef,lon_2D[i,j] ,lat_2D[i,j], 0)

            range_targ=np.sqrt((surf_x-pos_sc[i,0])**2+(surf_y-pos_sc[i,1])**2+(surf_z-pos_sc[i,2])**2)
            ranges[i,j]=range_targ

            H=np.sqrt((pos_sc[i,0])**2+(pos_sc[i,1])**2+(pos_sc[i,2])**2)
            RE=H-sc_alt[i]

            theta=-np.arcsin((H**2+range_targ**2-RE**2)/(2*H*range_targ))/np.pi*180.

            if np.isnan(theta): # Can happen for angles very close to pi
                theta=-90
            elevations[i,j]=-theta # Flip sign since elevations are defined positively in lut

    coords_spaceborne=np.vstack((center_lat_sc, center_lon_sc, sc_alt)).T

    return azimuths, elevations, ranges, coords_spaceborne


class SimulatedSpaceborne():
    '''
    The output class of simulated spaceborne radar swaths in the radar operator
    '''
    def __init__(self, list_radials, swath_dim):
        '''
        Returns a SimulatedSpaceborne Class instance
        Args:
            list_beams: the list of simulated radials as returned in the
                main RadarOperator class
            swath_dim: the horizontal dimension of the simulated spaceborne radar swath, this
                is needed because it can not be guessed from the data
                ex: [88, 49]

        Returns:
            a SimulatedSpaceborne with five fields:
                bin_surface: a 2D array giving the indexes of  radar bins
                    that correspond to the ground level
                lats: a 3D array containing the latitudes at the all spaceborne radar
                    gates (also in vertical)
                lons: a 3D array containing the longitudes at the all spaceborne radar
                    gates (also in vertical)
                heights: a 3D array containing the heights at the all spaceborne radar
                    gates (also in vertical)
                data: a dictionary of 3D arrays with all simulated variables
        '''
        
        # Reorganize output to make it easier to compare with spaceborne radar product format (like GPM-DPR...)
        # The idea is to keep only the bins that are above ground and to put
        # every beam into a common matrix

        [N,M] = swath_dim

        pol_vars = list_radials[0].values.keys() # List of simulated variables

        # make a deepcopy to avoid overwriting the original one
        list_beams_cp = copy.deepcopy(list_radials)

        bin_surface = np.zeros((N,M))
        # Here we remove all points that are below the topography COSMO
        for idx in range(len(list_beams_cp)):
            i = int(np.floor(idx / M))
            j = idx - i * M
            try:
                bin_surface[i,j] = len(list_beams_cp[idx].mask) - \
                    np.where(list_beams_cp[idx].mask<0)[0]
            except:
                bin_surface[i,j]=0
            for k in pol_vars:
                # Remove points that are below topo
                list_beams_cp[idx].values[k] = \
                    list_beams_cp[idx].values[k][list_beams_cp[idx].mask > -1]
            # Take only lat and lon profile where data is valid (above topo)
            list_beams_cp[idx].lats_profile = \
                list_beams_cp[idx].lats_profile[list_beams_cp[idx].mask > -1]

            list_beams_cp[idx].lons_profile =  \
                list_beams_cp[idx].lons_profile[list_beams_cp[idx].mask > -1]
            
            list_beams_cp[idx].heights_profile =  \
                list_beams_cp[idx].heights_profile[list_beams_cp[idx].mask > -1]

        # Length of longest beam
        max_len=np.max([len(r.heights_profile) for r in list_beams_cp])

        # Initalize output dictionary
        list_beams_formatted = {}
        for k in pol_vars:
            list_beams_formatted[k]=np.zeros((N,M,max_len))*float('nan')
        # Initialize lats and lons 3D array

        lats = np.zeros((N,M,max_len))*float('nan')
        lons = np.zeros((N,M,max_len))*float('nan')
        heights = np.zeros((N,M,max_len))*float('nan')

        # Fill the output dictionary starting from the ground
        for idx in range(len(list_beams_cp)):
            i = int(np.floor(idx/M))
            j = idx - i * M
            len_beam=len(list_beams_cp[idx].heights_profile)
            # Flip because we want to start from the ground
            l = list_beams_cp[idx].lats_profile
            ll = list_beams_cp[idx].lons_profile
            h = list_beams_cp[idx].heights_profile

            lats[i, j, 0:len_beam] = l[::-1]
            lons[i, j, 0:len_beam] = ll[::-1]
            heights[i, j, 0:len_beam] = h[::-1]

            # Flip [::-1] because we want to start from the ground
            for k in pol_vars:
                list_beams_formatted[k][i,j,0:len_beam] = \
                    list_beams_cp[idx].values[k][::-1]

        self.bin_surface = bin_surface
        self.lats = lats
        self.lons = lons
        self.heights = heights
        self.data = list_beams_formatted
