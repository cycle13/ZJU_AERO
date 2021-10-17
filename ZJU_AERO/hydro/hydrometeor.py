# -*- coding: utf-8 -*-

'''
Description: hydrometeors.py: Provides classes with relevant functions for all
hydrometeor types considered in the radar operator.
Computes all diameter dependent properties (orientation, aspect-ratio,
dielectric constants, velocity, mass... 
Author: Hejun Xie
Date: 2020-08-18 09:37:31
LastEditors: Hejun Xie
LastEditTime: 2021-10-17 16:19:09
'''

# Global import 
from textwrap import dedent

# Local import
from ._graupel import Graupel, NonsphericalGraupel
from ._ice import IceParticle
from ._snow import Snow, NonsphericalSnow
from ._rain import Rain

def create_hydrometeor(hydrom_type, scheme = '1mom'):
    """
    Creates a hydrometeor class instance, for a specified microphysical
    scheme
    Args:
        hydrom_type: the hydrometeor types, can be either
            'R': rain, 'S': snow aggregates, 'G': graupel
        scheme: microphysical scheme to use, can be either '1mom' (operational
           one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
    Returns:
        A hydrometeor class instance (see below)
    """

    from ..config.cfg import CONFIG

    mp_scheme_name = CONFIG['microphysics']['scheme_name']

    if  hydrom_type == 'R':
       return Rain(scheme, scheme_name=mp_scheme_name)
    elif hydrom_type == 'S':
        if CONFIG['microphysics']['psd_new_solver_S']:
            shape = CONFIG['microphysics']['shape_S']
            param = get_param(CONFIG['microphysics']['scattering_S'])
            return NonsphericalSnow(scheme, shape, scheme_name=mp_scheme_name, param=param)
        else:
            return Snow(scheme, scheme_name=mp_scheme_name)
    elif hydrom_type == 'G':
        if CONFIG['microphysics']['psd_new_solver_G']:
            shape = CONFIG['microphysics']['shape_G']
            return NonsphericalGraupel(scheme, shape, scheme_name=mp_scheme_name)
        else:
            return Graupel(scheme, scheme_name=mp_scheme_name)
    elif hydrom_type == 'I':
        return IceParticle(scheme, scheme_name=mp_scheme_name)
    else:
        msg = """
        Invalid hydrometeor type, must be R, S, G, I
        """
        return ValueError(dedent(msg))

def create_hydrometeor_db(hydrom_type, scheme = '1mom'):
    """
    Creates a hydrometeor class instance, for a specified microphysical
    scheme, only used for LUT generator
    Args:
        hydrom_type: the hydrometeor types, can be either
            'R': rain, 'S': snow aggregates, 'G': graupel
        scheme: microphysical scheme to use, can be either '1mom' (operational
           one-moment scheme) or '2mom' (non-operational two-moment scheme, not implemented yet)
    Returns:
        A hydrometeor class instance (see below)
    """

    if  hydrom_type == 'R':
        return Rain(scheme)
    elif hydrom_type == 'S':
        return Snow(scheme)
    elif hydrom_type == 'G':
        return Graupel(scheme)
    elif hydrom_type == 'I':
        return IceParticle(scheme)
    else:
        msg = """
        Invalid hydrometeor type, must be R, S, G, I
        """
        return ValueError(dedent(msg))

def get_param(scattering_method):
    try:
        param_str = scattering_method.split('_')[-1]
        param = float(param_str)
        return param
    except:
        return None
