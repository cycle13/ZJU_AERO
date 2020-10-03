'''
Description: wrapper for pycwr class
Author: Hejun Xie
Date: 2020-10-02 16:32:15
LastEditors: Hejun Xie
LastEditTime: 2020-10-03 19:26:09
'''

# Global import
from pycwr import core
import numpy as np
np.set_printoptions(threshold=np.inf)

# Local import
from ..config import cfg
from ..constants import global_constants as constants

RDOP_to_CINRAD_field_mapping = {
    'ZH'    : 'dBZ',
    'KDP'   : 'KDP',
    'PHIDP' : 'PhiDP',
    'RHOHV' : 'CC',
    'ZDR'   : 'ZDR',
    'RVEL'  : 'V',
}

class PycwrRadop(core.NRadar.PRD):
    def __init__(self, scan_type, scan):
        '''
        Creates a PycwrRadop Class instance
        Args:
            scan_type: the type of scan can be either 'ppi' or 'rhi', for
                vertical scans it doesn't matter, they can be 'ppi' or 'rhi'
            scan: the output of the get_PPI or get_RHI functions of the main
                RadarOperator class
        Returns:
            a PycwrRadop instance which can be used as specified in the
            pycwr doc: https://github.com/YvZheng/pycwr
        '''

        sitename = 'RDOP'
        nsweeps = len(scan['data'])
        rays_per_sweep = np.array([len(scan['data'][isweep]) for isweep in range(nsweeps)])
        nrays = np.sum(rays_per_sweep)
        
        latitude = scan['pos_time']['latitude']
        longitude = scan['pos_time']['longitude']
        altitude = scan['pos_time']['altitude']
        time = [scan['pos_time']['time'],] * nrays

        rrange = constants.RANGE_RADAR
        frequency = cfg.CONFIG['radar']['frequency']

        nyquist_velocity = np.array((nsweeps,), dtype=float)
        unambiguous_range = np.array((nsweeps,), dtype=float)
        bins_per_sweep = np.array((nsweeps,), dtype=int)
        nyquist_velocity[:] = cfg.CONFIG['radar']['nyquist_velocity']
        unambiguous_range[:] = cfg.CONFIG['radar']['range']
        bins_per_sweep[:] = len(rrange)

        if scan_type == 'ppi':
            fixed_angle = np.array(scan['elevations'], dtype=float)
        elif scan_type == 'rhi':
            fixed_angle = np.array(scan['azimuths'], dtype=float)

        # Initialize
        elevation   = []
        azimuth     = []
        sweep_start_ray_index   = []
        sweep_end_ray_index     = []
        idx_start   = 0
        idx_end     = -1

        # Initialize the field
        varnames = scan['data'][0][0].values.keys()
        fields = {}
        for varname in varnames:
            if varname in RDOP_to_CINRAD_field_mapping.keys():
                fields[RDOP_to_CINRAD_field_mapping[varname]] = np.empty((nrays, len(rrange)), dtype=float)

        for isweep in range(nsweeps):
            
            # sweep_start_ray_index and sweep_end_ray_index
            idx_end += rays_per_sweep[isweep]
            sweep_start_ray_index.append(idx_start)
            sweep_end_ray_index.append(idx_end)

            # elevation and azimuth
            if scan_type == 'ppi':
                elevation.extend(list([scan['elevations'][isweep]] * rays_per_sweep[isweep]))
                azimuth.extend(list(scan['azimuths']))
            elif scan_type == 'rhi':
                elevation.extend(list(scan['elevations']))
                azimuth.extend(list([scan['azimuths'][isweep]] * rays_per_sweep[isweep]))
            
            # fields
            for iradial, radial in enumerate(scan['data'][isweep]):

                for varname in varnames:
                    if varname in RDOP_to_CINRAD_field_mapping.keys():

                        data_radial = radial.values[varname]
                        # convert to dB
                        if varname in ['ZH', 'ZV', 'ZDR']:
                            data_radial[data_radial==0]=float('nan')
                            data_radial=10*np.log10(data_radial)

                        fields[RDOP_to_CINRAD_field_mapping[varname]][idx_start+iradial, :] = \
                            data_radial
                        
            # move forward 
            idx_start += rays_per_sweep[isweep]


        elevation = np.array(elevation, dtype=float)
        azimuth = np.array(azimuth, dtype=float)
        sweep_start_ray_index = np.array(sweep_start_ray_index, dtype=int)
        sweep_end_ray_index = np.array(sweep_end_ray_index, dtype=int)
        
        super(PycwrRadop,self).__init__(fields=fields, scan_type=scan_type, time=time,
                range=rrange, azimuth=azimuth, elevation=elevation, 
                latitude=latitude, longitude=longitude, altitude=altitude,
                sweep_start_ray_index=sweep_start_ray_index, sweep_end_ray_index=sweep_end_ray_index, 
                fixed_angle=fixed_angle,
                bins_per_sweep=bins_per_sweep, nyquist_velocity=nyquist_velocity,
                frequency=frequency, unambiguous_range=unambiguous_range,
                nrays=nrays, nsweeps=nsweeps, sitename=sitename)
