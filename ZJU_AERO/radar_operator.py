'''
Description: Defines the main radar operator class RadarOperator
which is the only one you should access as a user. This class allows to
compute PPI scans
Author: Hejun Xie
Date: 2020-08-22 12:45:35
LastEditors: Hejun Xie
LastEditTime: 2021-02-27 22:10:31
'''


# Global imports
from functools import partial
import multiprocess as mp
import numpy as np
import xarray as xr
import copy
import gc
import pickle
from textwrap import dedent

# Local imports
from .radar import PyartRadop, PycwrRadop, get_spaceborne_angles, SimulatedSpaceborne
from .config import createConfig
from .interp import get_interpolated_radial, integrate_radials

from .const import global_constants as constants
from .db import load_all_lut
from .utils import combine_subradials
from .core import get_radar_observables, cut_at_sensitivity

BASE_VARIABLES=['U','V','W','QR_v','QS_v','QG_v','QI_v','RHO','T']

class RadarOperator(object):
    def __init__(self, options_file = None, output_variables = 'all'):
        '''
        Creates a RadarOperator class instance that can be used to compute
        radar profiles (PPI, RHI, spaceborne_swath)
        Args:
            options_file: a .yml file containing the user configuration
                (see examples in the options_files folder)
            output_variables: can be either 'all', 'only_model', or 'only_radar'
                if 'only_model', only model native variables used for the
                radar operator (e.g. temp, concentrations, etc.) will be
                returned at the radar gates, no radar observables will be
                computed. This is fast and can be of use in particular
                circonstances
                if 'only_radar', only radar simulated radar observables will
                be returned at the radar gates (i.e. ZH, ZV, ZDR, KDP,...)
                if 'all', both radar observables and model variables will
                be returned at the radar gates

        Returns:
            A RadarOperator class instance
        '''
        
        createConfig(options_file) # Initialize options with specified file
        from .config.cfg import CONFIG

        self.current_microphys_scheme = '1mom'
        self.dic_vars = None
        self.N = 0 # atmospheric refractivity

        self.lut_sz = None
        self.config = CONFIG

        if output_variables in ['all','only_model','only_radar']:
            self.output_variables = output_variables
        else:
            self.output_variables = 'all'
            msg = """Invalid output_variables input, must be either
            'all', 'only_model' or 'only_radar'
            """
            print(msg)
        
    def close(self):
        '''
        Closes the RadarOperator class instance and deletes its content
        '''
        for hydrom, lut in self.lut_sz.items():
            if lut.type == 'xarray':
                lut.close()
        try:
            del dic_vars, N, lut_sz, output_variables
        except:
            pass
    
    @property
    def config(self):
        return copy.deepcopy(self.__config)

    @config.setter
    def config(self, config_dic):
        '''
        Update the content of the user configuration, applies the necessary
        changes to the constants and if needed reloads the lookup tables
        The best way to use this function is to retrieve the current
        configuration using  deepcopy: copy = copy.deepcopy(radop.config)
        then modify this copy: ex. copy['radar']['frequency'] = 9.41
        and then modify config.
		config = deepcopy.copy(radop.config)
		# change config
		radop.config = config

        Args:
            config_dic: a dictionnary specifying the new configuration to
                use.
        '''
        from .config.cfg import CONFIG

        # 1. get flag of loading scattering property database
        first_config = not hasattr(self, 'config')
        if not first_config:
            lut_needs_reload = CONFIG['radar']['frequency'] != self.config['radar']['frequency'] or \
                        CONFIG['microphysics']['with_melting'] != self.config['microphysics']['with_melting'] or \
                        not self.lut_sz
        flag_load_lut = first_config or lut_needs_reload

        # 2. update radar operator config and CONFIG in config module
        self.__config = copy.deepcopy(CONFIG)

        # 3. update derived constants
        constants.update() # Update constants now that we know user config

        # 4. load scattering property database
        if flag_load_lut:
            self.set_lut()
    
    def set_lut(self):
        '''
        Load a new set of lookup tables for the current radar operator
        based on the user configuration
        '''
        micro_scheme = self.current_microphys_scheme
        has_ice = self.config['microphysics']['with_ice_crystals']
        scattering_method_all = self.config['microphysics']['scattering']
        freq = self.config['radar']['frequency']
        folder_lut = self.config['microphysics']['folder_lut']

        list_hydrom = ['R','S','G']
        if has_ice:
            list_hydrom.extend(['I'])

        def get_method(hydrom_method):
            method = scattering_method_all if hydrom_method == 'default' else hydrom_method 
            return method
        scattering_method = {}
        for hydrom in list_hydrom:
            scattering_method[hydrom] = get_method(self.config['microphysics']['scattering_{}'.format(hydrom)])
            
        lut = load_all_lut(micro_scheme, list_hydrom, freq, scattering_method, folder_lut=folder_lut)

        self.lut_sz = lut

        # exit() # test database
    
    def get_pos_and_time(self):
        '''
        Get the position of the radar and time
        '''
        latitude = self.config['radar']['coords'][0]
        longitude = self.config['radar']['coords'][1]
        altitude = self.config['radar']['coords'][2]
        time=self.dic_vars.attrs['time'] # We could read any variable, T or others

        out={'latitude':latitude,'longitude':longitude,'altitude':altitude,\
        'time':time}

        return out
    
    def define_globals(self):
        '''
        This is used in the parallelization to get all global variables
        '''
        global output_variables
        output_variables = self.output_variables

        global dic_vars
        dic_vars=self.dic_vars

        global N
        N=self.N

        global lut_sz
        lut_sz = self.lut_sz

        return dic_vars, N, lut_sz, output_variables
    
    def load_model_file(self, filename, load_datetime, load_from_file=False, load_file=None):
        '''
        Loads data from a MODEL file, which is a pre-requirement for
        simulating radar observables.

        Args:
            filename: The name of the model file, or a list of model file name in netCDF4 format.
            load_datetime: The datetime we want to simulate, datetime class.
            load_from_file: A bool flag indicating whether to load from files dumped previously.
                Turn it on when execution acceleration is needed, turn it off when model file
                or model datetime is changed.
            load_file: The name of the file we want to load from. This parameter is forced if 
                load_from_file = True.
        Returns:
            No returns, but save the loaded_vars and atmosphere refraction index N 
            in the attributes of this RadarOperator class. 
        '''

        # if single model file is entered, make it a filelist
        if not isinstance(filename, list):
            filename = [filename]

        vars_to_load = copy.deepcopy(BASE_VARIABLES)

        if self.config['nwp']['name'] == 'grapes':
            from .nwp.grapes import check_if_variables_in_file
            from .nwp.grapes import get_grapes_variables as get_model_variables
        elif self.config['nwp']['name'] == 'wrf':
            from .nwp.wrf import check_if_variables_in_file
            from .nwp.wrf import get_wrf_variables as get_model_variables

        if self.config['refraction']['scheme'] in [2,3]:
            if check_if_variables_in_file(['T','P','QV']):
                vars_to_load.extend('N')
            else:
                msg = '''
                Necessary variables for computation of atm. refractivity:
                Pressure, Water vapour mass density and temperature
                were not found in file. Using 4/3 method instead.
                '''
                print(dedent(msg))
                self.config['refraction_method']='4/3'
        
        # check if all vars in file
        vars_ok = check_if_variables_in_file(['P','T','QV','QR','QC','QI','QS','QG','U','V','W'])

        if not vars_ok:
            msg = '''
            Not all necessary variables could be found in file
            For 1-moment scheme, the MODEL file must contain
            Temperature, Pressure, U-wind component, V-wind component,
            W-wind component, and mass-densities (Q) for Vapour, Rain, Snow,
            Graupel, Ice cloud
            '''
            raise ValueError(dedent(msg))
        else:
            print('Using 1-moment scheme')
            self.current_microphys_scheme = '1mom'
        
        # Read variables from netCDF4 file
        print('Reading variables ' + str(vars_to_load) + ' from file')

        if not load_from_file:
            loaded_vars = get_model_variables(filename, vars_to_load, load_datetime)
            loaded_vars.to_netcdf(load_file)
        else:
            loaded_vars = xr.open_dataset(load_file)

        self.dic_vars = loaded_vars #  Assign to class
        
        if 'N' in list(loaded_vars.data_vars.keys()):
            self.N=loaded_vars['N']
            self.dic_vars = self.dic_vars.drop_vars('N') # Remove N from the xarray dataset (we won't need it there)
    
        del loaded_vars
        gc.collect()
        print('-------done------')

        # Check if lookup tables are deprecated
        if not self.output_variables == 'model_only':
            if not self.lut_sz or \
                self.current_microphys_scheme != \
                    self.config['microphysics']['scheme']:
                print('Loading lookup-tables for current specification')
                self.set_lut()
                self.config['microphysics']['scheme'] = self.current_microphys_scheme

    def get_PPI_test(self, elevations, azimuths = None, az_step = None, az_start = 0,
                az_stop = 359):
        '''
        Simulates a PPI scan based on the user configuration (For single thread test)
        Args:
            elevations: a single scalar or a list of elevation angles in
                degrees. If a list is provided, the output will consist
                of several PPI scans (sweeps in the PyART class)
            azimuths: (Optional) a list of azimuth angles in degrees
            az_start az_step, az_stop: (Optional) If 'azimuths' is undefined
                these three arguments will be used to create a list of azimuths
                angles. Defaults are (0, 3dB_beamwidth, 359)
        Returns:
            A PPI profile at the specified elevation(s), in the form of a PyART
            class. To every elevation angle corresponds at sweep
        '''
        # Check if model file has been loaded
        if self.dic_vars=={}:
            print('No model file has been loaded! Aborting...')
            return

        # Check if list of elevations is scalar
        if np.isscalar(elevations):
            elevations=[elevations]

        # Needs to be done in order to deal with Multiprocessing's annoying limitations
        global dic_vars, N, lut_sz, output_variables
        dic_vars, N, lut_sz, output_variables = self.define_globals()
        # Define list of angles that need to be resolved
        if az_step == None:
            az_step=self.config['radar']['3dB_beamwidth']

        # Define list of angles that need to be resolved
        if np.any(azimuths == None):
            # Define azimuths and ranges
            if az_start>az_stop:
                azimuths=np.hstack((np.arange(az_start, 360., az_step),
                                    np.arange(0, az_stop + az_step, az_step)))
            else:
                azimuths=np.arange(az_start, az_stop + az_step, az_step)

        # Define  ranges
        rranges = constants.RANGE_RADAR

        # Initialize computing pool
        def worker(elev, azimuth):#
            print(azimuth)
            list_subradials = get_interpolated_radial(dic_vars,
                                                    azimuth,
                                                    elev,
                                                    N)
            if output_variables in ['all','only_radar']:
                output = get_radar_observables(list_subradials, lut_sz)
            if output_variables == 'only_model':
                output =  integrate_radials(list_subradials)
            elif output_variables == 'all':
                output = combine_subradials((output,
                            integrate_radials(list_subradials)))

            return output

        list_beams = worker(elevations[0], azimuths[0])

        np.set_printoptions(threshold=np.inf)

        list_beams = cut_at_sensitivity([list_beams])

        print(list_beams[0].values['ZH'])
        print(list_beams[0].values['RVEL'])

        self.close()
        
        del dic_vars
        del N
        del lut_sz
        gc.collect()

        exit()    

    def get_PPI(self, elevations, azimuths = None, az_step = None, az_start = 0,
                az_stop = 359, plot_engine='pyart'):
        '''
        Simulates a PPI scan based on the user configuration
        Args:
            elevations: a single scalar or a list of elevation angles in
                degrees. If a list is provided, the output will consist
                of several PPI scans (sweeps in the PyART class)
            azimuths: (Optional) a list of azimuth angles in degrees
            az_start az_step, az_stop: (Optional) If 'azimuths' is undefined
                these three arguments will be used to create a list of azimuths
                angles. Defaults are (0, 3dB_beamwidth, 359)
        Returns:
            A PPI profile at the specified elevation(s), in the form of a PyART
            class. To every elevation angle corresponds at sweep
        '''
        # Check if model file has been loaded
        if self.dic_vars=={}:
            print('No model file has been loaded! Aborting...')
            return

        # Check if list of elevations is scalar
        if np.isscalar(elevations):
            elevations=[elevations]

        # Needs to be done in order to deal with Multiprocessing's annoying limitations
        global dic_vars, N, lut_sz, output_variables
        dic_vars, N, lut_sz, output_variables = self.define_globals()
        # Define list of angles that need to be resolved
        if az_step == None:
            az_step=self.config['radar']['3dB_beamwidth']

        # Define list of angles that need to be resolved
        if np.any(azimuths == None):
            # Define azimuths and ranges
            if az_start>az_stop:
                azimuths=np.hstack((np.arange(az_start, 360., az_step),
                                    np.arange(0, az_stop + az_step, az_step)))
            else:
                azimuths=np.arange(az_start, az_stop + az_step, az_step)

        # Define  ranges
        rranges = constants.RANGE_RADAR

        # Initialize computing pool
        pool = mp.Pool(processes = mp.cpu_count(), maxtasksperchild=1)
        m = mp.Manager()
        event = m.Event()

        list_sweeps=[]
        def worker(event,elev, azimuth):#
            print('Azimuth: {:7.2f}'.format(azimuth))
            try:
                if not event.is_set():
                    list_subradials = get_interpolated_radial(dic_vars,
                                                          azimuth,
                                                          elev,
                                                          N)
                    if output_variables in ['all','only_radar']:
                        output = get_radar_observables(list_subradials, lut_sz)
                    if output_variables == 'only_model':
                        output =  integrate_radials(list_subradials)
                    elif output_variables == 'all':
                        output = combine_subradials((output,
                                 integrate_radials(list_subradials)))

                    return output
            except:
                # Throw signal back
                raise
                event.set()

        for e in elevations: # Loop on the elevations
            func = partial(worker,event,e)
            list_beams = pool.map(func,azimuths)
            list_sweeps.append(list(list_beams))

        pool.close()
        pool.join()

        del dic_vars
        del N
        del lut_sz
        gc.collect()

        # exit()
        if not event.is_set():
            # Threshold at given sensitivity
            if output_variables in ['all','only_radar']:
                list_sweeps = cut_at_sensitivity(list_sweeps)

            simulated_sweep={'elevations':elevations,'azimuths':azimuths,
            'ranges':rranges,'pos_time':self.get_pos_and_time(),
            'data':list_sweeps}

            if plot_engine == 'pyart':
                plot_instance = PyartRadop('ppi',simulated_sweep)
            elif plot_engine == 'pycwr':
                plot_instance = PycwrRadop('ppi',simulated_sweep)

            return plot_instance

            # return simulated_sweep

    def get_RHI(self, azimuths, elevations = None, elev_step = None,
                                            elev_start = 0, elev_stop = 90,
                                            plot_engine='pyart'):
        '''
        Simulates a RHI scan based on the user configuration
        Args:
            azimuths: a single scalar or a list of azimuth angles in
                degrees. If a list is provided, the output will consist
                of several RHI scans (sweeps in the PyART class)
            elevations: (Optional) a list of elevation angles in degrees
            elev_start elev_step, elev_stop: (Optional) If 'elevations' is
                undefinedthese three arguments will be used to create a list
                of elevations angles. Defaults are (0, 3dB_beamwidth, 359)
        Returns:
            A RHI profile at the specified elevation(s), in the form of a PyART
            class. To every azimuth angle corresponds at sweep
        '''
        # Check if model file has been loaded
        if self.dic_vars=={}:
            print('No model file has been loaded! Aborting...')
            return

        # Check if list of azimuths is scalar
        if np.isscalar(azimuths):
            azimuths=[azimuths]


        # Needs to be done in order to deal with Multiprocessing's annoying limitations
        global dic_vars, N, lut_sz, output_variables
        dic_vars, N, lut_sz, output_variables=self.define_globals()

        # Define list of angles that need to be resolved
        if np.any(elevations == None):
            if elev_step == None:
                elev_step = self.config['radar']['3dB_beamwidth']
            # Define elevation and ranges
            elevations = np.arange(elev_start, elev_stop + elev_step,
                                   elev_step)


        # Define  ranges
        rranges = constants.RANGE_RADAR

        # Initialize computing pool
        pool = mp.Pool(processes = mp.cpu_count(),maxtasksperchild=1)
        m = mp.Manager()
        event = m.Event()

        list_sweeps=[]

        def worker(event, azimuth, elev):
            print('Elevation: {:7.2f}'.format(elev))
            try:
                if not event.is_set():
                    list_subradials = get_interpolated_radial(dic_vars,
                                                          azimuth,
                                                          elev,
                                                          N)

                    if output_variables in ['all','only_radar']:
                        output = get_radar_observables(list_subradials, lut_sz)
                    if output_variables == 'only_model':
                        output =  integrate_radials(list_subradials)
                    elif output_variables == 'all':
                        output = combine_subradials((output,
                                 integrate_radials(list_subradials)))

                    return output
            except:
                # Throw signal back
                raise
                event.set()

        for a in azimuths: # Loop on the o
            func = partial(worker, event, a) # Partial function
            list_beams = pool.map(func,elevations)
            list_sweeps.append(list(list_beams))

        pool.close()
        pool.join()

        del dic_vars
        del N
        del lut_sz

        if not event.is_set():
            # Threshold at given sensitivity
            if output_variables in ['all','only_radar']:
                list_sweeps = cut_at_sensitivity(list_sweeps)

            simulated_sweep={'elevations':elevations,'azimuths':azimuths,
            'ranges':rranges,'pos_time':self.get_pos_and_time(),
            'data':list_sweeps}

            if plot_engine == 'pyart':
                plot_instance = PyartRadop('rhi',simulated_sweep)
            elif plot_engine == 'pycwr':
                plot_instance = PycwrRadop('rhi',simulated_sweep)

            return plot_instance

    def get_spaceborne_swath_test(self, swath_file, slice=None):
        '''
        Simulates a spaceborne radar swath, for single thread test
        Args:
            swath_file: a spaceborne file, not implemented yet, only for test.
            slice: slice of the spacebrone radar pixels in model domain. If slice is None,
            then the total swath is simulated.

        Returns:
            An instance of the SimulatedSpaceborne class (see spaceborne_wrapper.py) which
            contains the simulated radar observables.
        '''
        from .config import cfg

        # Check if model file has been loaded
        if self.dic_vars=={}:
            print('No model file has been loaded! Aborting...')
            return

        # Needs to be done in order to deal with Multiprocessing's annoying limitations
        global dic_vars, N, lut_sz, output_variables
        dic_vars, N, lut_sz, output_variables = self.define_globals()

        # get spaceborne radar observation angles
        az,elev,rang,coords_spaceborne = get_spaceborne_angles(swath_file, slice)

        def worker(params):
            # print(params)
            azimuth=params[0]
            elev=params[1]

            # Update spaceborne position and range vector for each beam
            cfg.CONFIG['radar']['range'] = params[2]
            cfg.CONFIG['radar']['coords'] = [params[3],
                                            params[4],
                                            params[5]]

            list_subradials = get_interpolated_radial(dic_vars,
                                                        azimuth,
                                                        elev,N = N)

            # print(integrate_radials(list_subradials).values['T'])

            if output_variables in ['all','only_radar']:
                output = get_radar_observables(list_subradials,lut_sz)
            if output_variables == 'only_model':
                output =  integrate_radials(list_subradials)
            elif output_variables == 'all':
                output = combine_subradials((output,
                        integrate_radials(list_subradials)))

            return output
        
        output = worker((az[0,0], elev[0,0], rang[0,0], 
        coords_spaceborne[0,0], coords_spaceborne[0,1], coords_spaceborne[0,2]))

        (nscan, npixel) = az.shape

        np.set_printoptions(threshold=np.inf)
        print(output.values['ZH'])
        
        self.close()

        del dic_vars
        del N
        del lut_sz
        gc.collect()
        
        # return output
        exit()

    def get_spaceborne_swath(self, swath_file, slice=None):
        '''
        Simulates a spaceborne radar swath
        Args:
            swath_file: a spaceborne file, not implemented yet, only for test.
            slice: slice of the spacebrone radar pixels in model domain. If slice is None,
            then the total swath is simulated.

        Returns:
            An instance of the SimulatedSpaceborne class (see spaceborne_wrapper.py) which
            contains the simulated radar observables.
        '''

        '''
        Attention for this import:
        from .config.cfg import CONFIG 
        is ok for single thread test,
        but in multiprocessing, CONFIG cannot be updated. 
        '''
        from .config import cfg

        # Check if model file has been loaded
        if self.dic_vars=={}:
            print('No model file has been loaded! Aborting...')
            return

        # Needs to be done in order to deal with Multiprocessing's annoying limitations
        global dic_vars, N, lut_sz, output_variables
        dic_vars, N, lut_sz, output_variables = self.define_globals()

        # get spaceborne radar observation angles
        az,elev,rang,coords_spaceborne = get_spaceborne_angles(swath_file, slice)

        # Initialize computing pool
        pool = mp.Pool(processes = mp.cpu_count())
        m = mp.Manager()
        event = m.Event()

        def worker(event,params):
            try:
                if not event.is_set():
                    # print(params)
                    azimuth=params[0]
                    elev=params[1]

                    # Update spaceborne position and range vector for each beam
                    cfg.CONFIG['radar']['range'] = params[2]
                    cfg.CONFIG['radar']['coords'] = [params[3],
                                                    params[4],
                                                    params[5]]

                    list_subradials = get_interpolated_radial(dic_vars,
                                                                azimuth,
                                                                elev,N = N)

                    if output_variables in ['all','only_radar']:
                        output = get_radar_observables(list_subradials,lut_sz)
                    if output_variables == 'only_model':
                        output =  integrate_radials(list_subradials)
                    elif output_variables == 'all':
                        output = combine_subradials((output,
                                integrate_radials(list_subradials)))

                    return output
            except:
                # Throw signal back
                raise
                event.set()

        (nscan, npixel) = az.shape
        
        list_beams=[]
        for iscan in range(nscan):
            print('running scan '+str(iscan))
            # Update radar position
            c0 = np.repeat(coords_spaceborne[iscan,0],len(az[iscan]))
            c1 = np.repeat(coords_spaceborne[iscan,1],len(az[iscan]))
            c2 = np.repeat(coords_spaceborne[iscan,2],len(az[iscan]))
            worker_partial = partial(worker,event)
            list_beams.extend(pool.map(worker_partial,zip(az[iscan],
                                                     elev[iscan],
                                                     rang[iscan],
                                                     c0,
                                                     c1,
                                                     c2)))

        pool.close()
        pool.join()

        del dic_vars
        del N
        del lut_sz
        gc.collect()

        if not event.is_set():
            # Threshold at given sensitivity
            if output_variables in ['all','only_radar']:
                list_beams = cut_at_sensitivity(list_beams)

            list_beams_formatted = SimulatedSpaceborne(list_beams,
                                                    (nscan,npixel))

            return list_beams_formatted

    def get_VPROF(self):
        '''
        Simulates a 90° elevation profile based on the user configuration

        Returns:
            A vertical profile at 90° elevation
        '''

        # Check if model file has been loaded
        if self.dic_vars=={}:
            print('No model file has been loaded! Aborting...')
            return

        # Needs to be done in order to deal with Multiprocessing's annoying limitations
        global dic_vars, N, lut_sz, output_variables
        dic_vars, N, lut_sz, output_variables = self.define_globals()

        # Radar ranges
        rranges = constants.RANGE_RADAR

        list_subradials = get_interpolated_radial(dic_vars, 0.0, 90.0, N)

        output = get_radar_observables(list_subradials, lut_sz)

        # Cut at given sensitivity
        output = cut_at_sensitivity([output])[0]

        if output_variables == 'all':
            output = combine_subradials((output, integrate_radials(list_subradials)))

        del dic_vars
        del N
        del lut_sz
        gc.collect()

        return output
