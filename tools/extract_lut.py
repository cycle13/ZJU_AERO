'''
Description: Extract Lookup Table from T-matrix output root directory and save it
as a netCDF format file (LevelA-database)
Author: Hejun Xie
Date: 2020-09-17 10:04:49
LastEditors: Hejun Xie
LastEditTime: 2020-11-20 17:04:30
'''

import os
import sys
import glob
import pandas as pd
import numpy as np
import xarray as xr
from scipy.integrate import quad  

RAW_DATA_ROOT = '/mnt/d/ZJU_AERO/radardb_hexgon'
LEVELA_DATA_TARGET = '/mnt/d/COSMO_POL/lut/iitm_masc/lut_SZ_S_9_41_1mom_LevelA.nc'
NBETA = 91
NELEVATION = 91
ANAME = 'amplmatrix_fw_test.dat'
PNAME = 'phasematrix_bw_test.dat'

def get_amplmatrix(filename):

    dataframe = pd.read_table(filename, delim_whitespace=True, skiprows=0, skipfooter=0, 
    names=[ 'beta', 'elevation',
            's11_r', 's11_i', 's12_r', 's12_i',
            's21_r', 's21_i', 's22_r', 's22_i'])

    dataframe['s11'] = dataframe['s11_r'] + dataframe['s11_i']*1j
    dataframe['s12'] = dataframe['s12_r'] + dataframe['s12_i']*1j
    dataframe['s21'] = dataframe['s21_r'] + dataframe['s21_i']*1j
    dataframe['s22'] = dataframe['s22_r'] + dataframe['s22_i']*1j

    return dataframe

def get_phasematrix(filename):

    dataframe = pd.read_table(filename, delim_whitespace=True, skiprows=0, skipfooter=0, 
    names=['beta', 'elevation',   
            'p11', 'p12', 'p13', 'p14',
            'p21', 'p22', 'p23', 'p24',
            'p31', 'p32', 'p33', 'p34',
            'p41', 'p42', 'p43', 'p44'])

    return dataframe

def get_coords(directory):
    subdirs = os.listdir(directory)
    coords = []
    for subdir in subdirs:
        coords.append(float(subdir.split('_')[-1]))
    coords.sort()
    return coords

def gaussian_pdf(std=10.0, mean=0.0):
    """Gaussian PDF for orientation averaging.

    Args:
        std: The standard deviation in degrees of the Gaussian PDF
        mean: The mean in degrees of the Gaussian PDF.  This should be a number
          in the interval [0, 180)

    Returns:
        pdf(x), a function that returns the value of the spherical Jacobian- 
        normalized Gaussian PDF with the given STD at x (degrees). It is 
        normalized for the interval [0, 180].
    """
    norm_const = 1.0
    def pdf(x):
        return norm_const*np.exp(-0.5 * ((x-mean)/std)**2) * \
            np.sin(np.pi/180.0 * x)
    norm_dev = quad(pdf, 0.0, 180.0)[0]
    # ensure that the integral over the distribution equals 1
    norm_const /= norm_dev 
    return pdf


if __name__ == '__main__':

    dims = ['aspect_ratio', 'temperature', 'Dmax', 'beta', 'elevation']
    dir_asp = RAW_DATA_ROOT
    dir_temp = os.path.join(dir_asp,os.listdir(dir_asp)[0])
    dir_size = os.path.join(dir_temp,os.listdir(dir_temp)[0])

    # add geometry size coordinate
    C = 299792458.
    # some band frequency for radar [GHZ]
    Frequencies = {"S":2.7 ,"C":5.6 ,"X":9.41 ,"Ku":13.6 ,"Ka":35.6}
    Lambda = C / Frequencies['X'] * 1e-9 * 1e+3 # [mm]
    Dmax = np.array(get_coords(dir_size)) * Lambda / np.pi

    coords = {'aspect_ratio': get_coords(dir_asp), 
    'temperature': get_coords(dir_temp),
    'Dmax': Dmax,
    'beta': np.linspace(0, 90, NBETA),
    'elevation': np.linspace(0, 90, NELEVATION)}

    size = tuple([len(coords[dim]) for dim in dims])

    real_variables      = [ 'p11_bw', 'p12_bw', 'p13_bw', 'p14_bw',
                            'p21_bw', 'p22_bw', 'p23_bw', 'p24_bw',
                            'p31_bw', 'p32_bw', 'p33_bw', 'p34_bw',
                            'p41_bw', 'p42_bw', 'p43_bw', 'p44_bw']
    complex_variables   = [ 's11_fw', 's12_fw', 's21_fw', 's22_fw']

    datadic = {}
    for real_variable in real_variables:
        datadic[real_variable] = (dims, np.empty(size, dtype='float32'))
    for complex_variable in complex_variables:
        datadic[complex_variable] = (dims, np.empty(size, dtype='complex64'))

    ds = xr.Dataset(datadic, coords)
    

    nodes = glob.glob('{}/asp_*/temp_*/size_*/'.format(RAW_DATA_ROOT))
    for node in nodes:
        asp = float(node.split('/')[-4].split('_')[-1])
        temp = float(node.split('/')[-3].split('_')[-1])
        size = float(node.split('/')[-2].split('_')[-1])
        iasp = coords['aspect_ratio'].index(asp)
        itemp = coords['temperature'].index(temp)
        isize = get_coords(dir_size).index(size)
        
        file_phasematrix_bw = os.path.join(node, PNAME)
        file_amplmatrix_fw = os.path.join(node, ANAME)
        
        phasematrix_bw = get_phasematrix(file_phasematrix_bw) # * (Lambda / (2*np.pi))**2 # [mm2]
        amplmatrix_fw = get_amplmatrix(file_amplmatrix_fw) # * (Lambda / (2*np.pi)) # [mm]
        
        for real_variable in real_variables:
            ds[real_variable][iasp, itemp, isize, ...] = \
            phasematrix_bw[real_variable.split('_')[0]].to_numpy().reshape(NBETA, NELEVATION) * \
                (Lambda / (2*np.pi))**2 # [mm2]

        for complex_variable in complex_variables:
            ds[complex_variable][iasp, itemp, isize, ...] = \
            amplmatrix_fw[complex_variable.split('_')[0]].to_numpy().reshape(NBETA, NELEVATION) * \
                (Lambda / (2*np.pi)) # [mm]    
    
    ds.coords["size_parameter"] = (("Dmax"), get_coords(dir_size))
    
    print(ds)

    ds.to_netcdf(LEVELA_DATA_TARGET, engine="h5netcdf")

    # with xr.open_dataset("./test.nc", engine="h5netcdf") as ds:
    #     print(ds)
