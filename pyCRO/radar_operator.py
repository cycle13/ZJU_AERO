'''
Description: Defines the main radar operator class RadarOperator
which is the only one you should access as a user. This class allows to
compute PPI scans
Author: Hejun Xie
Date: 2020-08-22 12:45:35
LastEditors: Hejun Xie
LastEditTime: 2020-10-01 16:55:06
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')


# Global imports
from functools import partial
import multiprocess as mp
import numpy as np
import copy
import pyWRF as pw
import gc
import pickle
from textwrap import dedent

# Local imports
from pyCRO.radar import PyartRadop
from pyCRO.config import cfg
from pyCRO.interpolation import get_interpolated_radial, integrate_radials

from pyCRO.constants import global_constants as constants
from pyCRO.lookup import load_all_lut
from pyCRO.utilities import combine_subradials
from pyCRO.scatter import get_radar_observables, cut_at_sensitivity

BASE_VARIABLES=['U','V','W','QR_v','QS_v','QG_v','QI_v','RHO','T']

class RadarOperator(object):
    def __init__(self, options_file = None, output_variables = 'all'):
        '''
        Creates a RadarOperator class instance that can be used to compute
        radar profiles (PPI, RHI, GPM)
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


        # delete the module's globals
        print('Reading options defined in options file')
        cfg.init(options_file) # Initialize options with specified file

        self.current_microphys_scheme = '1mom'
        self.dic_vars = None
        self.N = 0 # atmospheric refractivity

        self.lut_sz = None
        self.config = cfg.CONFIG

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
        pass
        try:
            del dic_vars, N, lut_sz, output_variables
        except:
            pass
        # self.config = None
        # cfg.CONFIG = None
    
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
            check: boolean to state if the specified new configuration
                dictionary should be checked and missing/invalid values should
                be replaced (as when a new .yml file is loaded). True by
                default
        '''
        # if check == False, no sanity check will be done, do this only
        # if you are sure of what you are doing

        print('Loading new configuration...')
        checked_config = cfg.sanity_check(config_dic)

        if hasattr(self, 'config'):
            # If frequency was changed or if melting is considered and not before
            #  reload appropriate lookup tables
            if checked_config['radar']['frequency'] != \
                 self.config['radar']['frequency'] or \
                 checked_config['microphysics']['with_melting'] != \
                 self.config['microphysics']['with_melting'] \
                 or not self.lut_sz:

                 print('Reloading lookup tables...')
                 self.__config = checked_config
                 cfg.CONFIG = checked_config
                 # Recompute constants
                 constants.update() # Update constants now that we know user config
                 self.set_lut()

            else:
                self.__config = checked_config
                cfg.CONFIG = checked_config
                # Recompute constants
                constants.update() # Update constants now that we know user config
        else:
            self.__config = checked_config
            cfg.CONFIG = checked_config
            constants.update()
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
        # freq = 5.6
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
        time=self.dic_vars['T'].attributes['time'] # We could read any variable, T or others

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
    
    def load_model_file(self, filename, itime=0, load_pickle=False, pickle_file=None):
        '''
        Loads data from a MODEL (currently only implemented for WRF model) file, which is a pre-requirement for
        simulating radar observables

        Args:
            filename: the name of the WRF netCDF file
        '''

        file_h = pw.open_file(filename)
        vars_to_load = copy.deepcopy(BASE_VARIABLES)

        vars_ok = file_h.check_if_variables_in_file(['P','T','QV','QR','QC','QI','QS','QG','U','V','W'])

        if self.config['refraction']['scheme'] in [2,3]:
            if file_h.check_if_variables_in_file(['T','P','QV']):
                vars_to_load.extend('N')
            else:
                msg = '''
                Necessary variables for computation of atm. refractivity:
                Pressure, Water vapour mass density and temperature
                were not found in file. Using 4/3 method instead.
                '''
                print(dedent(msg))
                self.config['refraction_method']='4/3'

        if not vars_ok:
            msg = '''
            Not all necessary variables could be found in file
            For 1-moment scheme, the WRF file must contain
            Temperature, Pressure, U-wind component, V-wind component,
            W-wind component, and mass-densities (Q) for Vapour, Rain, Snow,
            Graupel, Ice cloud
            '''
            raise ValueError(dedent(msg))
        else:
            print('Using 1-moment scheme')
            self.current_microphys_scheme = '1mom'
        
        # Read variables from GRIB file
        print('Reading variables ' + str(vars_to_load) + ' from file')

        if not load_pickle:
            loaded_vars = file_h.get_variable(vars_to_load, itime=itime, get_proj_info=True,
                                            shared_heights=False,assign_heights=True)
            
            # To deal with annoying issues with pickle
            for var in loaded_vars.values():
                var.file = None
            with open(pickle_file, "wb") as f:
                pickle.dump(loaded_vars, f)
        else:
            with open(pickle_file, "rb") as f:
                loaded_vars = pickle.load(f)

        self.dic_vars = loaded_vars #  Assign to class
        if 'N' in loaded_vars.keys():
            self.N=loaded_vars['N']
            self.dic_vars.pop('N',None) # Remove N from the variable dictionnary (we won't need it there)
        file_h.close()
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
        del loaded_vars

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

        # print(list_beams.values)

        del dic_vars
        del N
        del lut_sz
        gc.collect()

    def get_PPI(self, elevations, azimuths = None, az_step = None, az_start = 0,
                az_stop = 359):
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

            pyrad_instance = PyartRadop('ppi',simulated_sweep)

            return pyrad_instance

    def get_RHI(self, azimuths, elevations = None, elev_step = None,
                                            elev_start = 0, elev_stop = 90):
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

            pyrad_instance = PyartRadop('rhi',simulated_sweep)
            return  pyrad_instance
    