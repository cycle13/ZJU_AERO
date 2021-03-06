'''
Description: compute scattering lookup tables, with the IITM/TM levelA database
FOR OPERATIONAL USE
Author: Hejun Xie
Date: 2020-09-18 10:16:55
LastEditors: Hejun Xie
LastEditTime: 2021-09-29 15:45:06
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import os
import numpy as np
import xarray as xr
import multiprocessing as mp
import yaml
from scipy import stats
from scipy.integrate import quad

# Local imports
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.hydro import create_hydrometeor_db

# Whether to regenerate  all lookup tables, even if already present
FORCE_REGENERATION_SCATTER_TABLES = True

with open('db.yaml', 'r') as yamlfile:
    db_cfg = yaml.load(yamlfile, Loader=yaml.SafeLoader)

# Load DBS settings
for key in db_cfg['DBS']:
    locals()[key] = db_cfg['DBS'][key]

# Load DIMENSIONS settings:
for key in db_cfg['DIMENSIONS']:
    if key in ['ELEVATIONS', 'BETAS']:
        params = db_cfg['DIMENSIONS'][key]
        locals()[key] = np.arange(params[0], params[1], params[2])
    else:
        locals()[key] = db_cfg['DIMENSIONS'][key]
    
    if isinstance(locals()[key], list):
        locals()[key] = np.array(locals()[key])


# the key in levelA and levelB dictionaries
real_variables      = [ 'p11_bw', 'p12_bw', 'p13_bw', 'p14_bw',
                        'p21_bw', 'p22_bw', 'p23_bw', 'p24_bw',
                        'p31_bw', 'p32_bw', 'p33_bw', 'p34_bw',
                        'p41_bw', 'p42_bw', 'p43_bw', 'p44_bw']
complex_variables   = [ 's11_fw', 's12_fw', 's21_fw', 's22_fw']


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


def _get_weights_asp(hydrom, list_D, list_asp):
    '''
        Computes all probability distribution function (weights) against Dmax and aspect ratio
        Args:
            hydrom: a Hydrometeor Class instance (see hydrometeors.py)
            list_D: list of diameters in mm for which to compute the scattering properties.
            list_asp: list of aspect ratio (>1, Horizontal vs Rotational axis), 
            as level A database provided. But PDF are computed by inverse of aspect_ratio (>1, Horizontal vs Rotational axis), 
            assuming that the list_asp is equally spaced in 'inverse of aspect_ratio' domain.
            
        Returns:
            weights_asp: a xarray of weights against list_D and list_asp (>1, 
            Horizontal vs Rotational axis), for integration in _integration_over_pdf.
    '''
    ar_lambda, ar_loc, ar_mu = hydrom.get_aspect_ratio_pdf_masc(list_D)
    asp_wgt = np.empty((len(list_D), len(list_asp)), dtype='float32')
    
    for l in zip(range(len(list_D)), ar_lambda, ar_loc, ar_mu):
        gamm = stats.gamma(l[1],l[2],l[3])
        wei = gamm.pdf(list_asp)
        wei /= np.sum(wei) # renormalization
        asp_wgt[l[0], :] = wei

    coords = [("Dmax", list_D),("aspect_ratio", list_asp)]
    weights_asp = xr.DataArray(asp_wgt, coords=coords) 

    return weights_asp

def _get_weights_beta(hydrom, list_D, list_beta):
    '''
        Computes all probability distribution function (weights) against Dmax and Euler angle beta
        Args:
            hydrom: a Hydrometeor Class instance (see hydrometeors.py)
            list_D: list of diameters in mm for which to compute the scattering
                properties
            list_beta: list of beta, in range of [0, pi/2.], as level A database provided 
            
        Returns:
            weights_beta: a xarray of weights against list_D and list_beta [0, pi]
            for integration in _integration_over_pdf
    '''
    list_beta_padded = np.pad(list_beta, (0,len(list_beta)-1), 'reflect', reflect_type='odd')
    if hasattr(hydrom, 'get_canting_angle_std_masc'):
        beta_std = hydrom.get_canting_angle_std_masc(list_D)
    else:
        beta_std = hydrom.canting_angle_std * np.ones((len(list_D,)))
    beta_wgt = np.empty((len(list_D), len(list_beta_padded)), dtype='float32')
    
    for iD, D in enumerate(list_D):
        guassian = gaussian_pdf(std=beta_std[iD])
        wei = guassian(list_beta_padded)
        wei /= np.sum(wei) # renormalization
        beta_wgt[iD,:] = wei
    
    coords = [("Dmax", list_D),("beta", list_beta_padded)]
    weights_beta = xr.DataArray(beta_wgt, coords=coords) 

    return weights_beta

def _integrate_over_pdf(hydrom, hydrom_type, frequency, elevation, temperature, list_D, list_beta, list_asp):
    '''
        Computes all scattering properties for a given set of parameters and
        for a given hydrometeor, according its probability distribution function of
        orientation (orientation Euler angle beta and aspect ratio).
        This is to be used for all hydrometeors except melting snow and
        melting graupel (currently only snow implemented)
        Args:
            hydrom: a Hydrometeor Class instance (see hydrometeors.py)
            freqency: the frequency in GHz
            elevation: incident elevation angle in degrees
            temperature: the temperature in K
            list_D: list of diameters in mm for which to compute the scattering
                properties
            list_beta: list of Euler angle beta ([0, pi/2.]), as level A database providing
            list_asp: list of aspect ratio (<1.), as level A database providing
            
        Returns:
            ds_out: the dataset piece to be filled into levelB database 
            temperature: to be used in register callback
            frequency: to be used in register callback
    '''
    
    ds_in = levela_lut.sel(temperature=temperature, elevation=elevation)

    # get weights xarray
    if ARS[hydrom_type] != 'SINGLE':
        weights_asp = _get_weights_asp(hydrom, list_D, list_asp)
    weights_beta = _get_weights_beta(hydrom, list_D, list_beta)

    # start formulating the ds_out
    coords = {'Dmax': list_D}
    datadic = {}
    for var in ds_in.data_vars.keys():
        # Padding the DataArray in beta axis
        padded_da = ds_in[var].pad(beta=(0,len(list_beta)-1), mode='reflect')
        padded_da.coords['beta'] = np.pad(list_beta, (0, len(list_beta)-1), 'reflect', reflect_type='odd')
        # perform intergration over canting angle and aspect ratio probability distribution
        if ARS[hydrom_type] != 'SINGLE':
            datadic[var] = (padded_da*weights_asp*weights_beta).sum(dim='aspect_ratio').sum(dim='beta')
        else:
            datadic[var] = (padded_da*weights_beta).sum(dim='beta')
    
    ds_out = xr.Dataset(datadic, coords=coords)

    # print(ds_out.data_vars['p12_bw'])
    # print(ds_out)
    
    return (ds_out, temperature, elevation)

def assign_dataset_pieces(pack):
    ds, t, e = pack[0], pack[1], pack[2]
    print('temperature={}, elevation={}'.format(t, e))
    for var in ds.data_vars.keys():
        levelb_lut[var].loc[dict(temperature=t, elevation=e)] = ds[var]
    return

def sz_lut(hydrom_type, frequency, levela_name_lut, levelb_name_lut):
    """
        Computes and saves a scattering lookup table (level B database) for a given
        hydrometeor type (non melting) and various frequencies from level A database.
        Args:
            hydrom_type: the hydrometeor type, currently only snow is implemented.
            frequency: list of frequencies for which to obtain the
                lookup tables, in GHz.
            levela_name_lut: the filename of level A database (Input).
            levelb_name_lut: the filename of level B database (Output).
        Returns:
            No output but saves a lookup table (levelb_name_lut).
    """

    scheme = '1mom'
    
    hydrom = create_hydrometeor_db(hydrom_type,'1mom')	
    
    # list_D = np.linspace(hydrom.d_min,hydrom.d_max,NUM_DIAMETERS).astype('float32')
    
    global levela_lut, levelb_lut

    with xr.open_dataset(levela_name_lut, engine="h5netcdf") as levela_lut:
        # Dimensions that formulate the SZ_matices
        list_elevation = levela_lut.coords["elevation"].values
        list_temperature = levela_lut.coords['temperature'].values
        list_D = levela_lut.coords["Dmax"].values
        # Note: the list_D sho\uld meet the hydrometeor list_D for fitted IITM database

        # Dimensions that should be integrated over probability distribution 
        # in levelA to levelB transformation
        if ARS[hydrom_type] != 'SINGLE':
            list_asp = levela_lut.coords['aspect_ratio'].values
        else:
            list_asp = None
        list_beta = levela_lut.coords['beta'].values

        # start formulating the levelB database
        dims = ['temperature', 'Dmax', 'elevation']
        coords = {'temperature': list_temperature,
        'Dmax': list_D, 'elevation': list_elevation}
        size = tuple([len(coords[dim]) for dim in dims])
        datadic = {}
        for real_variable in real_variables:
            datadic[real_variable] = (dims, np.empty(size, dtype='float32'))
        for complex_variable in complex_variables:
            datadic[complex_variable] = (dims, np.empty(size, dtype='complex64'))
        
        levelb_lut = xr.Dataset(datadic, coords)

        '''
        Place to hold the integration over probability distribution
        '''

        # pool = mp.Pool(processes=mp.cpu_count(), maxtasksperchild=1)
        pool = mp.Pool(processes=mp.cpu_count())
        for e in list_elevation:
            for t in list_temperature:
                # pack = _integrate_over_pdf(hydrom, hydrom_type, frequency, e, t, list_D, list_beta, list_asp)
                # assign_dataset_pieces(pack)
                # exit()
                args = (hydrom, hydrom_type, frequency, e, t, list_D, list_beta, list_asp)
                pool.apply_async(_integrate_over_pdf, args=args, callback=assign_dataset_pieces)

        # Gather processes of multiprocess
        pool.close()
        pool.join()

        # pack = _integrate_over_pdf(hydrom, hydrom_type, frequency, 1., 253., list_D, list_beta, list_asp)

        levelb_lut.to_netcdf(levelb_name_lut, engine="h5netcdf")
        
        del levelb_lut

if __name__ == "__main__":
    '''
       Create all lookup tables for the specified hydrometeor types and
        microphysical schemes (currently only 1mom scheme implemented...)
    '''

    scheme = '1mom'
    
    for frequency in FREQUENCIES:
        for hydrom_type in HYDROM_TYPES:

            levela_name_lut = (FOLDER_LUT+"lut_SZ_"+hydrom_type+'_'+
                    str(frequency).replace('.','_')+'_'+scheme+"_LevelA"+".nc")

            levelb_name_lut = (FOLDER_LUT+"lut_SZ_"+hydrom_type+'_'+
                    str(frequency).replace('.','_')+'_'+scheme+"_LevelB"+".nc")
            
            if (FORCE_REGENERATION_SCATTER_TABLES
                or not os.path.exists(levelb_name_lut)):
                msg = '''
                Generating scatter table for 1 moment scheme,
                hydrometeor = {:s}
                freq = {:s}
                '''.format(hydrom_type, str(frequency))
                print(msg)

                sz_lut(hydrom_type, frequency, levela_name_lut, levelb_name_lut)
