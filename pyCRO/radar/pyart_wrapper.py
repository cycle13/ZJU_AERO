'''
Description: wrapper for pyart class
Author: Hejun Xie
Date: 2020-08-24 21:59:44
LastEditors: Hejun Xie
LastEditTime: 2020-10-03 19:53:04
'''

# Global imports
from pyart import graph, filters, config, core, correct
import os, glob, datetime
import numpy as np
import numpy.ma as ma

# Local imports
from ..constants import global_constants as constants
from ..utilities import row_stack

# MODEL

UNITS_SIMUL  = {'ZH':'dBZ',
               'KDP':'deg/km',
               'PHIDP':'deg',
               'RHOHV':'-',
               'ZDR':'dB',
               'RVEL':'m/s',
               'DSPECTRUM':'dBZ',
               'ZV':'dBZ',
               'U':'m/s',
               'V':'m/s',
               'W':'m/s',
               'T':'K',
               'RHO':'kg/m3',
               'QR_v':'kg/m3',
               'QS_v':'kg/m3',
               'QG_v':'kg/m3',
               'QH_v':'kg/m3',
               'ZV':'dBZ',
               'ATT_H':'dBZ',
               'ATT_V':'dBZ'}

VAR_LABELS_SIMUL = {'ZH':'Reflectivity',
                    'KDP':'Specific diff. phase',
                    'RHOHV':'Copolar corr. coeff.',
                    'ZDR':'Diff. reflectivity',
                    'RVEL':'Mean doppler velocity',
                    'DSPECTRUM':'Doppler spectrum',
                    'ZV':'Vert. reflectivity',
                    'U':'U-wind component',
                    'V':'V-wind component',
                    'W':'Vertical wind component',
                    'T':'Temperature',
                    'RHO':'Air density',
                    'QR_v':'Mass density of rain',
                    'QS_v':'Mass density of snow',
                    'QG_v':'Mass density of graupel',
                    'QH_v':'Mass density of hail',
                    'PHIDP':'Diff. phase shift',
                    'ATT_H':'Attenuation at hor. pol.',
                    'ATT_V':'Attenuation at vert. pol.'}

VMIN_SIMUL = {'ZH': 0.,
              'KDP': 0.,
              'PHIDP': 0.,
              'RHOHV': 0.6,
              'ZDR': 0.,
              'RVEL': -25,
              'DSPECTRUM': -50,
              'ZV': 0.,
              'U': -30,
              'V': -30,
              'W': -10,
              'T': 200,
              'RHO': 0.5,
              'QR_v': 0.,
              'QS_v': 0.,
              'QG_v': 0.,
              'QH_v': 0.,
              'ATT_H': 0,
              'ATT_V': 0}

VMAX_SIMUL = {'ZH': 55,
              'KDP': 1,
              'PHIDP': 20,
              'RHOHV': 1,
              'ZDR': 2,
              'RVEL': 25,
              'DSPECTRUM': 30,
              'ZV': 45,
              'U': 30,
              'V': 30,
              'W': 10,
              'T': 300,
              'RHO': 1.4,
              'QR_v': 1E-3,
              'QS_v': 1E-3,
              'QG_v': 1E-3,
              'QH_v': 1E-2,
              'ATT_H': 5,
              'ATT_V': 5}

class PyartRadop(core.Radar):
    """
    This is a class for radar operator outputs that inheritates the
    pyart core Class and thus has the same functionalities
    """
    def __init__(self, scan_type, scan):
        '''
        Creates a PyartRadop Class instance
        Args:
            scan_type: the type of scan can be either 'PPI' or 'RHI', for
                vertical scans it doesn't matter, they can be PPI or RHI
            scan: the output of the get_PPI or get_RHI functions of the main
                RadarOperator class
        Returns:
            a PyartRadop instance which can be used as specified in the
            pyART doc: http://arm-doe.github.io/pyart/
        '''

        N_sweeps=len(scan['data'])

        fields={}
        fixed_angle={}
        fixed_angle['data']=np.zeros(N_sweeps,)

        sweep_start_ray_index={}
        sweep_start_ray_index['data']=[]
        sweep_stop_ray_index={}
        sweep_stop_ray_index['data']=[]

        # Get all variables names
        varnames = list(scan['data'][0][0].values.keys())

        if 'range' in varnames:
            varnames.remove('range')
        for i,k in enumerate(varnames):

            fields[k]={}
            fields[k]['data']=[]
            try: # No info found in header
                fields[k]['long_name'] = VAR_LABELS_SIMUL[k]
            except:
                pass
            try:
                fields[k]['units'] = UNITS_SIMUL[k]
            except:
                pass
            try:
                fields[k]['valid_min'] = VMIN_SIMUL[k]
            except:
                pass
            try:
                fields[k]['valid_max'] = VMAX_SIMUL[k]
            except:
                pass
        # Add latitude and longitude
        fields['Latitude'] = {'data':[],'units':['degrees']}
        fields['Longitude'] = {'data':[],'units':['degrees']}

        # Initialize
        idx_start=0
        idx_stop=0
        elevations=[]
        azimuths=[]

        for i in range(N_sweeps):
            # Convert list of beams to array of data in polar coordinates
            polar_data_sweep={}
            for k in varnames:
                # print(np.array([len(it.values[k]) for it in scan['data'][i]]))
                polar_data_sweep[k] = np.array([it.values[k] for
                                               it in scan['data'][i]])
                if k in ['ZDR','ZV','ZH','DSpectrum']: # Convert to dB
                    polar_data_sweep[k][polar_data_sweep[k]==0]=float('nan')
                    polar_data_sweep[k]=10*np.log10(polar_data_sweep[k])

            # Add also latitude and longitude to variables
            polar_data_sweep['Latitude'] = \
                np.array([it.lats_profile for it in scan['data'][i]])
            polar_data_sweep['Longitude'] = \
                np.array([it.lons_profile for it in scan['data'][i]])


            [N_angles,N_ranges] = polar_data_sweep[varnames[0]].shape

            idx_stop=idx_start+N_angles-1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start = idx_stop + 1

            if scan_type == 'ppi':
                fixed_angle['data'] = scan['elevations'][i]
                elevations.extend(list([scan['elevations'][i]] * N_angles))
                azimuths.extend(list(scan['azimuths']))
            elif scan_type == 'rhi':
                fixed_angle['data'] = scan['azimuths'][i]
                elevations.extend(list(scan['elevations']))
                azimuths.extend(list([scan['azimuths'][i]] * N_angles))

            for k in polar_data_sweep.keys():
                if not len(fields[k]['data']):
                    fields[k]['data'] = polar_data_sweep[k]
                else:
                    fields[k]['data'] = row_stack(fields[k]['data'],
                                                    polar_data_sweep[k])
        

        for k in polar_data_sweep.keys():
            fields[k]['data'] = np.ma.array(fields[k]['data'],
                                mask=np.isnan(fields[k]['data']))
        # Add velocities (for Doppler spectrum) in instrument param
        try:
            instrument_parameters = {'varray': {'data':constants.VARRAY}}
        except:
            instrument_parameters = {}

        '''
        Position and time are obtained from the pos_time field of the
        scan dictionary. Note that these latitude and longitude are the
        coordinates of the radar, whereas the latitude and longitude fields
        are the coords of every gate
        '''

        latitude={'data' : np.array(scan['pos_time']['latitude'],dtype=float)}
        longitude={'data' :  np.array(scan['pos_time']['longitude'],dtype=float)}
        altitude={'data' :  np.array(scan['pos_time']['altitude'],dtype=float)}

        time_units='seconds since '+scan['pos_time']['time']
        time={'data' : np.zeros((N_angles,)),'units': time_units}


        sweep_number={'data' : np.arange(0,N_sweeps, dtype = float)}
        sweep_mode={'data' : [scan_type]*N_sweeps}

        azimuth={'data' : np.array(azimuths)}
        rrange={'data': scan['ranges']}
        elevation={'data' :np.array(elevations, dtype = float)}
        fixed_angle['data'] = np.array([fixed_angle['data']], dtype = float)
        sweep_start_ray_index['data'] = np.array(sweep_start_ray_index['data'], 
                             dtype = int)
        sweep_stop_ray_index['data'] = np.array(sweep_stop_ray_index['data'],
                           dtype = int)
        
        '''
        Finally add ranges as an additional variable, for convenience in order
        to filter gates based on their range
        '''

        fields['rangearray'] = {}
        fields['rangearray']['data'] = np.tile(rrange['data'],(len(elevation['data']),1))

        metadata = {}


        # Create PyART instance
        super(PyartRadop,self).__init__(time, rrange, fields, metadata,
        scan_type, latitude, longitude, altitude, sweep_number, sweep_mode,
        fixed_angle, sweep_start_ray_index, sweep_stop_ray_index, azimuth,
        elevation, instrument_parameters = instrument_parameters)

    def get_field(self,sweep_idx, variable):
        # see the get_field function in the PyART doc, this just flattens
        # the output (squeeze)
        out = super(PyartRadop,self).get_field(sweep_idx, variable)
        return np.squeeze(out)
