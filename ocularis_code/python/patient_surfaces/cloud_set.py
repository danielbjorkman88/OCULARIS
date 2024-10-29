# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 11:52:50 2023

@author: bjoerk_c
"""


import numpy as np
import patient_surfaces

from patient_surfaces.cloud import shell_ray_first_intersect



class PointCloudSet(patient_surfaces.PatientSurface):
    """
    Finds the average ray intersection point between two point clouds describing the same surface.
    
    """
    def __init__(self, points1, points2, alpha=12., irregular_cloud = False):
        self.type = "cloud set"
        self.points1 = points1
        self.points2 = points2
        self.clouds = [patient_surfaces.PointCloud(points1, alpha, irregular_cloud), patient_surfaces.PointCloud(points2, alpha, irregular_cloud)]
        self.new_points = []
        self.cloud1_points = []
        self.cloud2_points = []
        
    def __repr__(self):
        return "Point cloud set"

    def find_ray_interaction(self, ray):
        
        # intersection_points = []
        
        # for cloud in self.clouds:
        #     p = cloud.find_ray_intersection_with_shell(ray)
        #     intersection_points.append(p)
        # p1 = self.clouds[0].find_ray_intersection_with_shell(ray)
        # p2 = self.clouds[1].find_ray_intersection_with_shell(ray)
        
        p1 = shell_ray_first_intersect(self.clouds[0].shell, ray) # find_ray_intersection_with_shell(ray)
        p2 = shell_ray_first_intersect(self.clouds[1].shell, ray) 
        

        # p_mid = (intersection_points[0] + intersection_points[1]) / 2
        
        if type(p1) != np.ndarray:
            if p1 == None:
                self.clouds[1].find_ray_interaction(ray)
                return
            
        if type(p2) != np.ndarray:
            if p2 == None:
                self.clouds[0].find_ray_interaction(ray)
                return
        
        p_mid = (p1 + p2)/2
        assert type(p_mid) == np.ndarray
        ray.surface_interaction_point_candidates.append(p_mid)
        
        self.cloud1_points.append(p1)
        self.cloud2_points.append(p2)
        self.new_points.append(p_mid)
        