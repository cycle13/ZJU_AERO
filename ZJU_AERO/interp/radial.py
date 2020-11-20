# -*- coding: utf-8 -*-

'''
Description: radial.py: Defines the Radial class, which is used to interpolated data
along a radar radial or subradial 
Author: Hejun Xie
Date: 2020-08-15 10:44:13
LastEditors: Hejun Xie
LastEditTime: 2020-08-15 16:39:33
'''


class Radial():
    def __init__(self, dic_values, mask, lats_profile,lons_profile,
                 dist_ground_profile,heights_profile,
                 elev_profile=[], quad_pt = [], quad_weight=1):
        """
        Creates a Radial instance, which contains all info about the
        interpolated radar radial
        Args:
            dic_values: dictionary containing all interpolated variables
            mask: mask for all radar gates, 1 = above model top, -1 = below
                  topo, 0 = OK
            lats_profile: vector containing the latitude at every gate (WGS84)
            lons_profile: vector containing the longitude at every gate (WGS84)
            dist_ground_profile: vector containing the dist. at ground to
                                 every gate (m)
            elev_profile: vector containing the incident elevation angle at
                                 every gate (degree)
            quad_pt: tuple (phi,theta) giving the azimuth and elevation
                     coordinates of the quadrature point that was used
            quad_weight : float giving the weight of the corresponding
                          quadrature weight

        Returns:
            A Radial class instance

        """
        self.mask = mask
        self.quad_pt = quad_pt
        self.quad_weight = quad_weight
        self.lats_profile = lats_profile
        self.lons_profile = lons_profile
        self.dist_profile = dist_ground_profile
        self.heights_profile = heights_profile
        self.elev_profile = elev_profile
        self.values = dic_values

        # This two are only for the melting scheme and are updated directly
        # in melting.py
        self.has_melting = False
        self.mask_ml = None
