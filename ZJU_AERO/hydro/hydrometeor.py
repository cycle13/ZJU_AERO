# -*- coding: utf-8 -*-

'''
Description: hydrometeors.py: Provides classes with relevant functions for all
hydrometeor types considered in the radar operator.
Computes all diameter dependent properties (orientation, aspect-ratio,
dielectric constants, velocity, mass... 
Author: Hejun Xie
Date: 2020-08-18 09:37:31
LastEditors: Hejun Xie
LastEditTime: 2021-06-13 18:26:33
'''

# Global import 
from textwrap import dedent

# Local import
from ._grauple import Graupel
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

