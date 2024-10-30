# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 13:44:33 2021

@author: bjoerk_c


This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.


"""

import numpy as np
import math


class Ray:
    def __init__(self, S, T):
        self.S = np.asarray(S)
        self.T = np.asarray(T)
        self.D12 = np.sqrt(np.sum((self.S - self.T)**2))
        self.vec = self.T - self.S
        self.vec_norm = self.vec/self.D12

        self.zes = np.asarray([self.S[2], self.T[2]])
        self.yes = np.asarray([self.S[1], self.T[1]])
        self.xes = np.asarray([self.S[0], self.T[0]])

        self.acc_depth = 0  # [mm]
        self.acc_depth_before_skin_plane = 0  # [mm]
        self.surface_behind_mesh_by = 0 # [mm]
        self.passed_collimation = False
        self.retracted_by_cloud = False
        self.surface_behind_skin_plane = False
        self.surface_intersection_identified = False
        self.surface_behind_mesh = False
        self.surface_case = ""
        
        self.idx = []
        self._surface_interaction_point = None
        self._skin_plane_point = None
        self.aperture_intersect = None
        self.surface_behind_skin_plane_by = 0
        
        self.surface_interaction_point_candidates = []

        self.title = "From {} to {} ".format(self.S, self.T)

    def add_depth(self, additional_depth):
        self.acc_depth += additional_depth
        
    def add_depth_before_skin_plane(self, additional_depth):
        self.acc_depth_before_skin_plane += additional_depth

    def __repr__(self):
        return self.title

    @property
    def surface_interaction_point(self):
        return self._surface_interaction_point

    @surface_interaction_point.setter
    def surface_interaction_point(self, point):
        self._surface_interaction_point = point
        
    @property
    def skin_plane_point(self):
        return self._skin_plane_point

    @skin_plane_point.setter
    def skin_plane_point(self, point):
        self._skin_plane_point = point        

    def find_intersect_with_plane(self, planeNormal, planePoints):

        if isinstance(planePoints, list):
            planePoint = planePoints[0]
        else:
            planePoint = planePoints

        ndotu = planeNormal.dot(self.vec_norm)

        if ndotu != 0:
            w = self.S - planePoint
            si = -planeNormal.dot(w) / ndotu
            intersect = w + si * self.vec_norm + planePoint

        return intersect
    
    def interaction_with_shell(self, shell):
        points, ind = shell.ray_trace(self.S, self.T)
        
        shortest_distance = 1000000
        first_intersect = []
        for p in points:
            dist = math.dist(p, self.S)
            if dist < shortest_distance:
                shortest_distance = dist
                first_intersect = p
                
        cloud_surface_interaction_point = np.asarray(first_intersect)
        return cloud_surface_interaction_point
    
    
    def find_surface_intersect(self, shell, planeNormal, planePoints ):
        
        skin_plane_point = self.find_intersect_with_plane( planeNormal , planePoints)

        cloud_point = self.interaction_with_shell(shell)

        
        self.skin_plane_point = skin_plane_point
        
        if len(cloud_point) == 0:
            self.surface_interaction_point = skin_plane_point
            return

        if math.dist(cloud_point, self.S) < math.dist(skin_plane_point, self.S):
            self.surface_interaction_point = cloud_point
        else:
            self.surface_interaction_point = skin_plane_point
        