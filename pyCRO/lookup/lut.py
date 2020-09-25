'''
Description: defines the Lookup Class as well as a set of functions
used to load and save scattering lookup tables
Author: Hejun Xie
Date: 2020-08-19 22:09:15
LastEditors: Hejun Xie
LastEditTime: 2020-09-25 23:58:09
'''


# Global imports
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error

from  scipy import ndimage
import tempfile
import os
import tarfile
import glob
import shutil
import pickle
import xarray as xr
from io import BytesIO
from textwrap import dedent


def load_all_lut(scheme, list_hydrom, frequency, scattering_method, folder_lut=None):
    '''
    Loads all scattering lookup tables for the user specified parameters
    Args:
        scheme: the microphysical scheme, either '1mom' or '2mom'(not implemented yet)
        list_hydrom: the list of hydrometeors for which scattering properties
            should be obtained: 'R': rain, 'S': snow, 'G': graupel. Ex: ['R','S','G']
        frequency: the frequency in GHz, make sure the lookup tables have
            been previously computed for the corresponding frequency!
        scattering_method: A dictionary describing the scattering method that is used by every
            hydrometeors in list_hydrom, can be either 'tmatrix_masc' or 'iitm_masc',
            which correspond to subfolders in the lookup folder. You could add more...
            Ex: {'S':'iitm_masc', 'R':'tmatrix_masc', ...}
        folder_lut: The folder that contains lookup table

    Returns:
        lut_sz: dictionary containing the lookup table for every hydrometeor
            type given in 'list_hydrom', the lookup tables are instances of
            the class Lookup_table (see below)
    '''

    # Get current directory
    if folder_lut is None:
        folder_lut = os.path.dirname(os.path.realpath(__file__))+'/'
    lut_sz = {}

    for h in list_hydrom:
        if scattering_method[h] == 'tmatrix_masc':
            folder_lut_method = folder_lut + 'tmatrix_masc/'
        elif scattering_method[h] == 'iitm_masc':
            folder_lut_method = folder_lut + 'iitm_masc/'

        freq_str = str(frequency).replace('.','_')
        if scattering_method[h] == 'iitm_masc':
            name = 'lut_SZ_' + h + '_' + freq_str + '_' + scheme + '_' + 'LevelB' + '.nc'
        else:
            name = 'lut_SZ_' + h + '_' + freq_str + '_' + scheme + '.lut'
        print(folder_lut_method + name)
        try:
            engine = 'xarray' if scattering_method[h]=='iitm_masc' else 'numpy'
            lut_sz[h] = load_lut(folder_lut_method + name, engine=engine)
        except:
            raise
            msg = """
            Could not find lookup table for scheme = {:s}, hydrometeor =
            {:s}, frequency = {:f} and scattering method = {:s}
            """.format(scheme, h, frequency, scattering_method)
            raise IOError(dedent(msg))

    return lut_sz

def load_lut(filename, engine):
    if engine == 'numpy':
        return load_lut_numpy(filename)
    elif engine == 'xarray':
        return load_lut_xarray(filename)

def load_lut_numpy(filename):
    '''
    Loads an instance of the Loookup_table class, previously saved to the
    drive with the save_lut function in cosmo_pol, (copy from cosmo_pol...) 
    Args:
        filename: the complete filename (with path), indicating where the
            lookup table is stored

    Returns:
        lut: the lookup table as an instance of the class Lookup_table Class
        (see below)
    '''
    tar = tarfile.open(filename, "r")
    
    lut = Lookup_table()
    for member in tar.getmembers():
        array_file = BytesIO()
        array_file.write(tar.extractfile(member).read())
        name = member.name.replace('.npy','')
        array_file.seek(0)
        data = np.load(array_file, allow_pickle = True, encoding = 'latin1')

        if name == 'axes_names':
            data = data.all()
        setattr(lut, name, data)

    tar.close()
    return lut

def load_lut_xarray(filename):
    '''
    Loads an instance of the Loookup_table class by loading the LevelB xarray database,
    previously saved to the drive with the sz_lut function in compute_lut_sz.py
    Args:
        filename: the complete filename (with path), indicating where the
            lookup table is stored

    Returns:
        lut: the lookup table as an instance of the class Lookup_table Class
        (see below)
    '''

    lut = Lookup_table_xarray(filename)
    return lut

class Lookup_table_xarray:
    '''
    The Lookup_table class used to store scattering properties of hydrometeors
    and perform queries of the scattering properties, this class assumes all
    stored data to be defined on regular grids, actually a wrapper for xarray
    '''
    
    def __init__(self, filename):
        self.type = 'xarray'
        self._ds = xr.open_dataset(filename, engine='h5netcdf')
    
    def close(self):
        self._ds.close()
    
    def get_axis_values(self, axis_name):
        return self._ds.coords[axis_name].values
    
    def lookup_line(self, **kwargs):
        '''
        Query the xarray levelB database with interp method.
        get the interpolated and flattened data table
        '''
        interped_data = {}
        # Adapt the kwargs for iitm xarray database
        iitm_kwargs = {}
        for key, value in kwargs.items():
            if key == 'e':
                iitm_kwargs['elevation'] = xr.DataArray(value, dims='line')
            elif key == 't':
                iitm_kwargs['temperature'] = xr.DataArray(value, dims='line')
            else:
                iitm_kwargs[key] = xr.DataArray(value, dims='line')
        # interp the xarray data into queried hyperplane
        for var in self._ds.data_vars.keys():
            interped_data[var] = self._ds[var].interp(**iitm_kwargs).values

        flattened_data = self._flatten_SZ_matrix(interped_data)
        return flattened_data
    
    def _flatten_SZ_matrix(self, interped_data):
        '''
        Flatten the SZ matrix interpolated from xarray levelB matrix
        return a table with dimension (n_valid_gates, nbins_D, 12)
        '''
        shape_single = list(interped_data['p11_bw'].shape)
        shape_single.append(12)
        flattened_shape = tuple(shape_single)
        flattened_data = np.empty(flattened_shape, dtype='float32')
        
        flattened_map_real = {  'p11_bw':0,     'p12_bw':1,     'p13_bw':-1,    'p14_bw':-1,
                                'p21_bw':2,     'p22_bw':3,     'p23_bw':-1,    'p24_bw':-1,
                                'p31_bw':-1,    'p32_bw':-1,    'p33_bw':4,     'p34_bw':5,
                                'p41_bw':-1,    'p42_bw':-1,    'p43_bw':6,     'p44_bw':7}
        flattened_map_cplx = {  's11_fw':8,     's12_fw':-1,    's21_fw':-1,    's22_fw':10}

        # flatten the data
        for key,value in flattened_map_real.items():
            if value != -1:
                flattened_data[...,value] = interped_data[key]

        for key,value in flattened_map_cplx.items():
            if value != -1:
                flattened_data[...,value] = interped_data[key].real
                flattened_data[...,value+1] = interped_data[key].imag
        
        return flattened_data

class Lookup_table:
    '''
    The Lookup_table class used to store scattering properties of hydrometeors
    and perform queries of the scattering properties, this class assumes all
    stored data to be defined on regular grids
    '''
    def __init__(self):

        # Contains the independent variable values known along
        # each axis of the data [[x0, x1, ..., xn], [y0, y1, ..., yn], ...]
        self.type = 'numpy'
        self.axes = []

        # Dictionary whose keys are the name of all axes and its values are
        # the index of the axis, ex: {'x': 0, 'y': 1, 'z': 2}, for 3D
        # Cartesian coordinates
        self.axes_names = {}

        # The limits for every axis (min and max values)
        self.axes_limits= []
         # the step for every axis: x1 - x0
        self.axes_step = []

        # to avoid the array shape check error in set_value_table for axis 'd',
        # which has different axes value for different 'wc'
        self.axes_len = []

        # Stores the dependent data for each point of the regular grid defined
        # the axes
        self.value_table = []


    def add_axis(self, name, axis_values=None):
        '''
        Add an axis to the lookup table. Axis correspond to the dimension of
        the data, for example for 3D coordinates, the axis would be 'x', 'y'
        and 'z'
        Args:
            name: name of the axis, for example 'd' for diameters
            axis_values: 'the values corresponding to the specified new axis
                ex for diameters: np.linspace(0.2, 8, 1024)
        '''
        if self.axes_names.has_key(name):
            raise Error("Axis already exists with name: '%s'" % name)
        axis_i = len(self.axes)

        self.axes_names[name] = axis_i
        axis_values = np.asarray(axis_values).astype('float32')

        self.axes_limits.append([np.min(axis_values), np.max(axis_values)])
        self.axes_step.append(axis_values[1]-axis_values[0])
        self.axes_len.append(len(axis_values[0]) if isinstance(axis_values[0], np.ndarray) else len(axis_values))
        self.axes.append(axis_values)

    def set_axis_values(self, axis_name, axis_values):
        '''
        Set the axis values for the specified axis.
        Axis values define points along the axis at which measurements
        were taken.

        Args:
            axis_name: name of the axis, for example 'd' for diameters
            axis_values: 'the values corresponding to the specified new axis
                ex for diameters: np.linspace(0.2, 8, 1024)
        '''

        axis_i = self.axes_names[axis_name]
        axis_values=np.asarray(axis_values).astype('float32')

        self.axes_limits[axis_i]=[np.min(axis_values),np.max(axis_values)]

        self.axes[axis_i] = axis_values


    def set_value_table(self, value_table):
        """Set the value table to the specified axes
        The shape of the data should correspond to the length of all axis
        i.e. value_table.shape = (len(axis[0]), len(axis[1]),...,len(axis[n]))

        Args:
            value_table: the independent data as a list of lists of numpy
            array
        """
        if not isinstance(value_table,np.ndarray):
            value_table = np.array(value_table)

        # Check dimensions
        tuple_axes_len = tuple(self.axes_len)
        if value_table.shape != tuple_axes_len:
            msg = '''
            The shape of the specified data does not match with the length
            of all axes, please ensure that
            data.shape = (len(axis[0]), len(axis[1]),...,len(axis[n]))
            '''
            return ValueError(dedent(msg))
        else:
            self.value_table = value_table


    def get_axis_name(self, axis_i):
        """Returns the name of an axis given its index

        Args:
            axis_i: the index of the axis (zero-based)
        Returns:
            The name of the axis for example 'd' or 'x' or 'y'
        """
        result = None
        for name, i in self.axes_names.items():
            if i == axis_i:
                result = name
                break
        return result


    def lookup_pts(self, coords, order = 1):
        """ Lookup the interpolated value for given axis values

        Args:
            coords: the coordinates where you want to interpolate the data
                must be contained within the bounds for every axis
                (no extrapolation). The coordinates are in the form of
                 a N x M array, where N is the number of points
                you want to interpolate and M is the number of axes
            order: interpolation order, must be in the range 0 to 5
                0 = nearest neighbour, 1 = linear interpolation,...
        Returns:
            out: the interpolated values at the specified coordinates
        """
        # Check that a value table exists.

        if not len(self.value_table):
            raise Error("No values set for lookup table")

        if coords.shape[0] != len(self.axes_names):
            raise Error("Invalid coordinates dimension")

        coords = [(c - lo) * (n - 1) / (hi - lo) for (lo, hi), c, n in zip(self.axes_limits,
                  coords, self.value_table.shape)]

        out = ndimage.map_coordinates(self.value_table, coords, order=order,
                                      mode = 'constant', cval = np.nan )

        if self.int_transform:
            return self.get_float_data(out)
        else:
            return out

    def lookup_line(self,**kwargs):
        """ Lookup a N-D plane from the data, by fixing one axis or more
            Currently only performs nearest neighbour interpolation
        Args:
            The number of arguments is variable. The function is used like
            this lookup_line(name_0 = val_0,name_1 = al_1,...)
            where 'name_0' and 'name_1' are the names of two axes in the
            lookup table, and val_0 and val_1, are the values used for these
            axes. Basically you extract a slice from the value_table by
            setting a certain number of axes to a fixed value
            For example if you have three axes: 'x', 'y', 'z',
            lookup_line(x = 1) will return the plane corresponding to x = 1
            lookup_line(x = 1, y = 2) will return the line corresponding
            to x = 1 and y = 2
            lookup_line(x = 1, y = 2, z = 4) will return a single line
            THE ORDER OF THE ARGUMENTS IS NOT IMPORTANT
        Returns:
            out: the values at the specified hyperplane
        """
        v = self.value_table
        dim = self.value_table.shape

        I = [slice(None)]*v.ndim
        for i,k in enumerate(kwargs.keys()):
            if k in self.axes_names.keys():
                ax_idx = self.axes_names[k]

                closest = np.floor((kwargs[k]-
                                    self.axes_limits[ax_idx][0])
                                    /self.axes_step[ax_idx])
                closest = np.array(closest, dtype= int)
                closest[closest < 0] = 0
                closest[closest >= dim[ax_idx]] = dim[ax_idx] - 1
                I[ax_idx] = closest

        return v[tuple(I)]


if __name__ == "__main__":
    lut = load_lut_xarray('/home/xhj/wkspcs/Radar-Operator/cosmo_pol/pathos/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelB.nc')
    
    Dmax = lut.get_axis_values('Dmax')
    print(Dmax)

    interped = lut.lookup_line(e=[1.5,2.5], t=[213, 223])
    print(interped.shape)

    lut.close()
