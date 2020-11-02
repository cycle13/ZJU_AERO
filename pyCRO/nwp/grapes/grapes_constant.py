'''
Description: Define some grapes constant
Author: Hejun Xie
Date: 2020-11-02 10:02:58
LastEditors: Hejun Xie
LastEditTime: 2020-11-02 12:42:51
'''

grapes_constant = {
    'g': 9.81,      # [g*m-2]
    'P0': 1e5,      # [Pa]
    'kappa': 2./7., # [-]
    'Rd': 287.05,   # [J*kg-1*K-1]
    'Rv': 451.51,   # [J*kg-1*k-1]
    'Rdv': 0.6357,  # [-] Rdv = Rd / Rv
    'Rvd': 1.573,   # [-] Rvd = Rv / Rd
}


derived_var_unit = {
    'U':'m*s-1',
    'V':'m*s-1',
    'W':'m*s-1',
    'T':'K',
    'P':'Pa',
    'Pw':'Pa',
    'Qhydro':'kg*kg-1',
    'RHO':'kg*m-3',
    'QV_v':'kg*m-3',
    'QR_v':'kg*m-3',
    'QS_v':'kg*m-3',
    'QG_v':'kg*m-3',
    'QC_v':'kg*m-3',
    'QI_v':'kg*m-3',
    'N':'-;*1E-6'
}

derived_var_long_name = {
    'U':'U component wind',
    'V':'V component wind',
    'W':'W component wind',
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
    'N':'Refraction index'
}

raw_var_unit = {
    'zz':'m',
    'pi':'-',
    'th':'K',
    'u':'m/s',
    'v':'m/s',
    'w':'m/s',
    'Qv':'kg/kg',
    'Qc':'kg/kg',
    'Qr':'kg/kg',
    'Qi':'kg/kg',
    'Qs':'kg/kg',
    'Qg':'kg/kg',
    'dbz':'dBZ',
}
