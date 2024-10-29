# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:45:41 2022

@author: Daniel Björkman
"""




import patient_surfaces

class SkinPlane(patient_surfaces.PatientSurface):
    def __init__(self, skin_plane_point, skin_plane_normal):
        self.type = "plane"
        self.skin_plane_point = skin_plane_point
        self.skin_plane_normal = skin_plane_normal
        
    def find_ray_interaction(self, ray):
        skin_plane_intersect = ray.find_intersect_with_plane( self.skin_plane_normal , self.skin_plane_point )
        ray.surface_interaction_point_candidates.append(skin_plane_intersect)

    def __repr__(self):
        return f"Plane, p = {self.skin_plane_point}, normal = {self.skin_plane_normal}"