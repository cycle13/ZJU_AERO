'''
Description: Define some wrf constant
Author: Hejun Xie
Date: 2020-11-06 10:27:47
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 21:33:02
'''

wrf_constant = {
    'g': 9.81,      # [g*m-2]
    'kappa': 2./7., # [-]
    'Rd': 287.05,   # [J*kg-1*K-1]
    'Rv': 451.51,   # [J*kg-1*k-1]
    'Rdv': 0.6357,  # [-] Rdv = Rd / Rv
    'Rvd': 1.573,   # [-] Rvd = Rv / Rd
}


derived_var_unit = {
    'U':'m s-1',
    'V':'m s-1',
    'W':'m s-1',
    'Z':'m',
    'T':'K',
    'P':'Pa',
    'Pw':'Pa',
    'Qhydro':'kg kg-1',
    'RHO':'kg m-3',
    'QV_v':'kg m-3',
    'QR_v':'kg m-3',
    'QS_v':'kg m-3',
    'QG_v':'kg m-3',
    'QC_v':'kg m-3',
    'QI_v':'kg m-3',
    'N':'-; 1E-6',
}

derived_var_long_name = {
    'U':'U component wind',
    'V':'V component wind',
    'W':'W component wind',
    'Z':'Height',
    'T':'Temperature',
    'P':'Pressure',
    'Pw':'Water vapour pressure',
    'Qhydro':'All hydormeteors mixing ratio',
    'RHO':'Air Density',
    'QV_v':'Vapour mass density',
    'QR_v':'Rain mass density',
    'QS_v':'Snow mass density',
    'QG_v':'Grauple mass density',
    'QC_v':'Cloud mass density',
    'QI_v':'Ice crystals mass density',
    'N':'Refraction index',
}

raw_var_map = {
    'QV':'QVAPOR',
    'QC':'QCLOUD',
    'QR':'QRAIN',
    'QI':'QICE',
    'QS':'QSNOW',
    'QG':'QGRAUP',
    'TP': 'HGT', # topograph [m]
    'U': 'U',
    'V': 'V',
    'W': 'W',
    'PTP': 'T', # perturbation potential temperature [K]
    'PTB': 'T00', # base potential temperature [K]
    'GPP': 'PH', # perturbation geopotential [m2*s-2]
    'GPB': 'PHB', # base geopotential [m2*s-2]
    'PP': 'P', # pressure [Pa]
    'PB': 'PB', # base pressure [Pa]
    'P0': 'P00', # P0 for potential temperature definition [Pa]
}
