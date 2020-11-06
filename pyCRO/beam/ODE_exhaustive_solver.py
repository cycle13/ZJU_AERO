'''
@Description: solving the ODE in an exhaustive way, 
in the method of Zeng and Blahak (2014)
@Author: Hejun Xie
@Date: 2020-08-02 08:31:04
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 22:47:57
'''

# global import
import xarray as xr
import numpy as np
import pyWRF as pw
import pyproj
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp


# Local imports
from ..config.cfg import CONFIG
from ..utilities import get_earth_radius
from ..constants import global_constants as constants


def _beam_solver(range_vec, coords_radar, n, b, topo, elevation_angle, RE):
    '''
    beam propogate ODE:
    dy/dr = beam_ODE(r, y)
    y[0]: s
    y[1]: a = dh/dr
    y[2]: h

    Args:
        range_vec: vector of all ranges along the radial [m]
        coord_radar: radar coordinates as 3D tuple [lat, lon, alt], ex
            [47.35, 12.3, 682]
        n: refraction index of the azimuthal plane (type: xarray DataArray) [-]
        b: the slope of refraction index azimuthal plane (type: xarray DataArray) = dn/dh [-*m^-1]
        topo: interpolated function of topology [m]
        elevation_angle: elevation angle [degrees]
        RE: earth radius at radar latitude [m]
    Returns:
        s: vector of distance at the ground along the radial [m]
        h: vector of heights above ground along the radial [m]
        e: vector of incident elevation angles along the radial [degrees]
    '''

    # define events:
    def hit_topo(r, y):
        h, s = y[2], y[0]
        return h - topo(s)
    
    def hit_model_top(r, y):
        h = y[2]
        return h - constants.MAX_MODEL_HEIGHT
    
    hit_topo.terminal = True
    hit_model_top.terminal = True

    # define the ODE
    def beam_ODE(r, y):
        s, a, h = y[0], y[1], y[2]
        ds2dr = np.sqrt(1-a*a)
        da2dr = (1 - a * a) * (
        b.interp(distance=s, height=h, kwargs={"fill_value": None}) /
        n.interp(distance=s, height=h, kwargs={"fill_value": None}) +
        1. / (RE + h) )
        dh2dr = a
        dy2dr = [ds2dr, da2dr, dh2dr]
        return dy2dr

    # solve the ODE
    y0 = [0., np.sin(np.deg2rad(elevation_angle)), coords_radar[2]] # s, a, h 
    sol = solve_ivp(beam_ODE, [0, range_vec[-1]], y0,
                t_eval=range_vec, events=(hit_topo, hit_model_top))

    s = sol.y[0]
    h = sol.y[2]
    e = np.rad2deg(np.arcsin(sol.y[1]))

    # deal with terminated radials
    pad_num = len(range_vec) - len(h)
    
    if pad_num != 0:
        s = np.pad(s, (0,pad_num), mode='constant', constant_values=0.)
        e = np.pad(e, (0,pad_num), mode='constant', constant_values=1.)

    if len(sol.t_events[0]) != 0: # hit the topo
        h = np.pad(h, (0,pad_num), mode='constant', constant_values=-9999.)
    elif len(sol.t_events[1]) != 0: # hit the model top
        h = np.pad(h, (0,pad_num), mode='constant', constant_values=constants.MAX_MODEL_HEIGHT+9999.)

    # convert from float64 to float32 for the sake of efficiency
    return np.array(s, dtype='float32'), np.array(h, dtype='float32'), np.array(e, dtype='float32')


def ODEZeng2014_exhaustive(range_vec, elevation_angles, azimuth_angle, coords_radar, N):
    '''
    Computes the trajectory of a radar beam along a specified radial with the
    ODE method of Zeng and Blahak (2014) in an exhaustive way
    TODO: adapt for other NWP models
    Args:
        range_vec: vector of all ranges along the radial [m]
        elevation_angles: a list of elevation angles of the beams in the same azimuth_angle [degrees]
        azimuth_angle: azimuth angle of that bunch of beams [degrees]
        coord_radar: radar coordinates as 3D tuple [lat, lon, alt], ex
            [47.35, 12.3, 682]
        N: atmospheric refractivity as a MODEL variable

    Returns:
        s: list of vector of distance at the ground along the radial [m]
        h: list of vector of heights above ground along the radial [m]
        e: list of vector of incident elevation angles along the radial [degrees]
    '''

    from ..config.cfg import CONFIG
    if CONFIG['nwp']['name'] == 'grapes':
        from ..nwp.grapes import WGS_to_GRAPES as WGS_to_MODEL
    elif CONFIG['nwp']['name'] == 'wrf':
        from ..nwp.wrf import WGS_to_WRF as WGS_to_MODEL

    # Get info about NWP coordinate system
    proj_MODEL = N.attrs
    
    RE = get_earth_radius(coords_radar[0])

    # get the mesh of atmosphere refraction index in the azimuth plane
    g = pyproj.Geod(ellps='WGS84')
    
    height_step     = 500
    distance_mesh  = range_vec # sample distances
    height_mesh     = np.arange(0, constants.MAX_MODEL_HEIGHT + height_step, height_step) # sample height
    
    n_data = np.empty((len(distance_mesh), len(height_mesh)), dtype='float32')
    t_data = np.empty((len(distance_mesh)), dtype='float32')

    for idistance, distance in enumerate(distance_mesh):
        lon, lat, ang = g.fwd(coords_radar[1], coords_radar[0], azimuth_angle, distance)

        # Convert WGS84 coordinates to MODEL coordinates model coordinates [x, y]
        coords_aziplane_in_MODEL = WGS_to_MODEL([lat, lon], proj_MODEL)

        # TODO: For model grid only, model grid (I, J)
        llc_MODEL = (0., 0.)
        res_MODEL = [1., 1.]

        # Get index of radar in MODEL coordinates model grids [I, J]
        pos_aziplane_bin = [(coords_aziplane_in_MODEL[0]-llc_MODEL[0]) / res_MODEL[0],
                        (coords_aziplane_in_MODEL[1]-llc_MODEL[1]) / res_MODEL[1]]
        
        # Get refractive index profile from refractivity estimated from NWP variables
        # data (bottom_top, south_north, west_east)
        n_vert_profile = 1 + (N.data[:,int(np.round(pos_aziplane_bin[1])),
                                int(np.round(pos_aziplane_bin[0]))]) * 1E-6
        # Get corresponding altitudes
        h = N.coords['z-levels'][:,int(np.round(pos_aziplane_bin[1])),
                                    int(np.round(pos_aziplane_bin[0]))]

        fn = interp1d(h, n_vert_profile, fill_value='extrapolate')
        n_data[idistance, :] = fn(height_mesh)
        t_data[idistance] = h[0]

    # get DataArray of n, b
    # n: refraction index [-]
    n_data[n_data < 1] = 1
    n = xr.DataArray(n_data, coords=[distance_mesh, height_mesh], dims=['distance', 'height'])
    # b: slope of refraction index [-/m^-1]
    b_data = np.diff(n_data, axis=1) / np.diff(height_mesh)
    height_mesh_b = (height_mesh[:-1] + height_mesh[1:]) / 2.
    b = xr.DataArray(b_data, coords=[distance_mesh, height_mesh_b], dims=['distance', 'height'])
    # get event function: topo [m]
    topo = interp1d(distance_mesh, t_data, fill_value='extrapolate')
    
    s, h, e = [], [], []
    for elevation_angle in elevation_angles:
        ks, kh, ke = _beam_solver(range_vec, coords_radar, n, b, topo, elevation_angle, RE)
        s.append(ks); h.append(kh); e.append(ke)
    
    return s, h, e


if __name__ == "__main__":
    
    FILENAME = '../../../cosmo_pol/pathos/WRF/wsm6/wrfout_d03_2013-10-06_00_00_00'
    file_h = pw.open_file(FILENAME)
    d = file_h.get_variable(['N'], itime=10, assign_heights=True)
    
    range_vec = np.arange(0, 50*1000, 500)
    elevation_angles = [1.0]
    azimuth_angle = 120
    coords_radar = [27.9, 120.8, 200]
    
    s, h, e =  ODEZeng2014_exhaustive(range_vec, elevation_angles, azimuth_angle, coords_radar, d['N'])
    print(s[0])
    print(h[0])
    print(e[0])
    