# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 13:44:33 2021

@author: bjoerk_c
"""

import reference_frames
import numpy as np
import math
import os
import sys
from pathlib import Path
import inspect
from mpmath import radians
from scipy.spatial.transform import Rotation

currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)


class OptisGeometry:

    """
    Class containing points, planes and patches needed for raytracing and plotting
    
    """

    def __init__(self, config):
        
        if not type(config) == dict:
            raise ValueError

        self.config = config

        self.aperture_plane = []
        self.aperture_point = []

        self.isocentre = np.asarray([0, 0, 0])

        self.wedge_insert_point = []
        self.wedge_d = []
        #self.wedge_point = []
        self.wedge_vector = []
        self.wedge_top_point = []
        self.wedge_plane = []

        self.target_plane_p1x = []
        self.target_plane_p2x = []
        self.aperture_plane_p1x = []
        self.aperture_plane_p2x = []

        self.wedge_plane_p1x = []
        self.wedge_plane_p2x = []
        self.wedge_plane_p3x = []

        # patches
        self.collimator_top = []
        self.collimator_bot = []
        self.wedge_patch = []

        self.config_geo()

    def config_geo(self):

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')

        self.central_axis = ref.central_axis

        self.aperture_plane = self.central_axis.vec_norm
        self.aperture_point = ref.origin_ref - \
            self.config["AID"] * self.central_axis.vec_norm

        self.fixation_light_plane = self.central_axis.vec_norm
        self.fixation_light_point = ref.origin_ref - \
            self.config["FID"] * self.central_axis.vec_norm

        self.target_plane = self.central_axis.vec_norm
        self.target_point = self.central_axis.T

        half_aperture = self.config["Aperture"]/2

        wedge_extension = 10  # Wedge extends past aperture by 10 mm for plotting

        self.wedge_insert_point = self.aperture_point + ref.xvec_in_ref * \
            math.cos(radians(self.config["wedge_insert_angle"]))*half_aperture + ref.yvec_in_ref*math.sin(
                radians(self.config["wedge_insert_angle"])) * half_aperture
        self.wedge_insert_direction = self.aperture_point - self.wedge_insert_point
        self.wedge_insert_direction = self.wedge_insert_direction / \
            np.linalg.norm(self.wedge_insert_direction)

        self.wedge_apex_point = self.wedge_insert_point + \
            self.wedge_insert_direction*(self.config["wedge_cover"])

        self.wedge_top_point = self.wedge_insert_point + self.central_axis.vec_norm * math.tan(math.radians(self.config["wedge_angle"]))*self.config["wedge_cover"]
            
        
        # self.wedge_insert_point + self.central_axis.vec_norm * \
        #     math.sin(self.config["wedge_angle"]*math.pi/180) * \
        #     (self.config["wedge_cover"] + wedge_extension)
        self.wedge_base_point = self.wedge_insert_point - \
            self.wedge_insert_direction*wedge_extension

        self.skin_plane = self.central_axis.vec_norm
        self.skin_plane_point = ref.origin_ref - \
            self.config["skin_plane"] * self.central_axis.vec_norm
            
            
        self.p_axis_light = ref.origin_ref - \
            self.config["FID"] * self.central_axis.vec_norm

        # self.p_polar = self.p_axis_light + np.asarray([-1, 0, 0])*math.tan(
        #     math.radians(self.polar_angle))*self.tps_config["config"]["FID"]

        # self.light_vector = self.p_polar - self.p_axis_light

        # r = Rotation.from_euler(
        #     'xyz', (0, 0, self.azimuth_angle), degrees=True)

        # self.light_vector = r.apply(self.light_vector)

        # self.light_point = self.p_axis_light + self.light_vector

        p0, p1 = self.wedge_apex_point, self.wedge_top_point

        if (p0 == p1).all():
            """
            In case no wedge
            """
            self.wedge_plane = self.aperture_plane
        else:

            p2 = self.wedge_apex_point + \
                np.cross(self.central_axis.vec_norm,
                         self.wedge_insert_direction)

            x0, y0, z0 = p0
            x1, y1, z1 = p1
            x2, y2, z2 = p2

            ux, uy, uz = u = [x1-x0, y1-y0, z1-z0]
            vx, vy, vz = v = [x2-x0, y2-y0, z2-z0]

            u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]

            normal = np.array(u_cross_v)
            normal = normal/np.linalg.norm(normal)

            self.wedge_plane = normal

        # self.wedge_d = []
        # self.wedge_point = []
        # self.wedge_vector = []

        # self.target_plane_vis_point1 = []
        # self.target_plane_p2x = []
        # self.aperture_plane_p1x = []
        # self.aperture_plane_p2x = []

        # self.wedge_plane_p1x = []
        # self.wedge_plane_p2x = []
        # self.wedge_plane_p3x = []

    def __repr__(self):
        pass
