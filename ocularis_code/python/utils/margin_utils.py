# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 10:29:57 2022

@author: bjoerk_c
"""


import math
import numpy as np
import warnings





class Displacement:
    def __init__(self, vec, magnitude, gaze_vector = [0,0,0]):
        
        # vec is the translation direction
        vec = np.asarray(vec)
        
        norm = np.linalg.norm(vec)
        
        if norm!= 0:
            self.vec = vec/np.linalg.norm(vec)
        else:
            self.vec = vec
        
        self.magnitude = magnitude
        
        self.gaze_vector = np.asarray(gaze_vector)
        
        self.D95 = 0
        self.V90 = 0
        
        self.dose_target_ref = None
        self.dose_target_displaced= None
        self.distal_spare_distance = 0
        self.proximal_spare_distance = 0
        self.distal_R90_distance = 0
        self.distal_point_contained = True
        self.proximal_point_contained = True
        self.target_prior_to_skin_plane = False
        self.min_dose_ref = 0
        self.min_dose_displaced = 0
        self.dose_proximal_end = 0
        self.delta_pupil_pos = None
        self.new_target_com = None
        self.n_elements_mask = 0
        
    def __repr__(self):
        return f"Disp translation = {self.vec}, gaze vector = {self.gaze_vector}"
        
    @property
    def as_dict(self):
        
        out = {}
        
        out["vec"] = self.vec
        out["magnitude"] = self.magnitude
        out["gaze_vector"] = self.gaze_vector
        out["dose_target_ref"] = self.dose_target_ref
        out["dose_target_displaced"] = self.dose_target_displaced
        out["distal_spare_distance"] = self.distal_spare_distance
        out["proximal_spare_distance"] = self.proximal_spare_distance
        out["distal_point_contained"] = self.distal_point_contained
        out["proximal_point_contained"] = self.proximal_point_contained
        out["min_dose_ref"] = self.min_dose_ref
        out["min_dose_displaced"] = self.min_dose_displaced
        out["distal_R90_distance"] = self.distal_R90_distance
        out["target_prior_to_skin_plane"] = self.target_prior_to_skin_plane
        out["dose_proximal_end"] = self.dose_proximal_end
        out["delta_pupil_pos"] = self.delta_pupil_pos
        out["new_target_com"] = self.new_target_com
        out["V90"] = self.V90
        out["n_elements_mask"] = self.n_elements_mask
        return out
        
        