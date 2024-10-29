#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 11:00:36 2021

@author: bjoerk_c
"""

import os
import logging
import numpy as np
from scipy.interpolate import interp1d
from pathlib import Path
import matplotlib.patches as patches
from scipy import ndimage
import math
from decorators import timer
from shapely.geometry.polygon import Point
import pyvista as pv
from abc import ABC

from utils.normalization_utils import find_R90


import reference_frames
import dose_engine
import geometry

# import jax.numpy as jnp

try:
    # Code disables unimportant error message originating from vtk depricated triangles
    import vtk
    vtk.vtkLogger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_OFF)
except:
    pass


class Algo(ABC):
    def __init__(self, basedata, config):
        """

        Parameters
        ----------
        basedata : dict
            Base data from configuration.py
        config : dict
            Dose engine configuration file.

        Returns
        -------
        None.

        """

        self.logger = logging.getLogger("Algo.log")
        if not self.logger.handlers:
            self.logger.setLevel(logging.INFO)
            self.logger.addHandler(logging.StreamHandler())

        self.curr_path = Path(os.getcwd())
        absolute_path = Path(os.path.dirname(__file__))
        self.data_path = absolute_path.parent.parent / "Data"
        self.test_data_path = absolute_path.parent.parent / "Data" / "TestData"
        self.support_data_path = absolute_path.parent.parent / "Data" / "Support"
        self.plot_path = absolute_path.parent.parent / "plots"
        if not os.path.isdir(self.plot_path):
            os.makedirs(self.plot_path)

        self.data = basedata
        self.config = config

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        
        self.central_axis = ref.central_axis

        self.Nvoxels = 0

        # Orhogonal distance to central axis, from which rays are defined for calculation of radiological depth
        self.orthogonal_limit = 0

        self.dose = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])
        self.target_matrix = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])

        self.radii = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])
        self.dist_aperture = -10 *np.ones(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])
        self.first_quadrant_bev = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])

        # Binary mask with voxels in relevant slice
        self.binary_mask = np.zeros(self.dose.shape)
        if self.config["Image"] == "xslice":  # ZX plane
            #"XZ slice"
            y = config["Slice"][1]
            self.binary_mask[0:, y, 0:] = 1
        elif self.config["Image"] == "yslice":  # YX plane
            #"YZ slice"
            x = config["Slice"][0]
            self.binary_mask[x, 0:, 0:] = 1
        elif self.config["Image"] == "zslice":  # XY plane
            #"YX slice"
            z = config["Slice"][2]
            self.binary_mask[0:, 0:, z] = 1
        elif self.config["Image"] == "3D":
            #"3D volume"
            self.binary_mask[0:, 0:, 0:] = 1
        else:
            assert 0, "Other formats not supported"

        self.nozzle = dose_engine.Nozzle(basedata, self.data_path) #, self.config
        self.DD = dose_engine.DepthProfile()
        self.collimator = dose_engine.Collimator(self.config)

        
        self.wedge = dose_engine.Wedge(self.config)
        self.compensator = None
        
        if 'skin_plane_point' in self.config:
            self.skin_plane = self.config['skin_plane_point']
        else:
            self.skin_plane = []

        self.rays = []
        self._target_rays = []
        self.indices_target = []
        self.surfaces = self.config["surfaces"]

        
        self.target_prior_to_skin_plane = False
        self.depth_dose_configured = False
        self.collimator_configured = False
        self.medium_configued = False
        self._surface_interaction_configured = False
        self.target_configured = False
        self.nozzle_configured = False
        self.rays_defined = False
        self.radii_configured = False
        self.raytracer_configured = False
        self.distances_to_aperture_configured = False
        self.dose_calculated = False
        
        self.deepest_ray = None
        self.shallowest_ray = None
        self.depthprofile_filename = None
        self.target_radiological_depth = None
        self.expected_range = 0
        self.expected_modulation = 0


        self.braggpeaks = []
        self.SOBP_xes = []
        self.SOBP_yes = []
        self.SOBP_f = []
        self.slope_vps_aperture = []
        self.phantom = []
        self.aperture_line_xes = []
        self.aperture_line_yes = []
        self.collimator_top = []
        self.collimator_bot = []
        self.inflection = []
        self.inflection_xes = []
        self.inflection_yes = []
        self.wedge_angle = []

        if self.config["Image"] == "xslice":
            slice_label = "XZ slice"
        elif self.config["Image"] == "yslice":
            slice_label = "YZ slice"
        elif self.config["Image"] == "zslice":
            slice_label = "YX slice"
        elif self.config["Image"] == "3D":
            slice_label = "3D volume"
        else:
            assert 0, "Other formats not supported"
        

        self.config_medium()
        self.raytracer = dose_engine.RayTracer(self.config, self.medium)
        
        self.logger.info('Dose Engine initialized for dose calculation of: ')
        self.logger.info(
            f"  {slice_label} ")  # "with resolution of {self.medium.resolution[0]}, {self.medium.resolution[1]}, {self.medium.resolution[2]} mm")
        
        
    def __repr__(self):
        
        try:
            return self.nozzle.foil + " " + self.nozzle.modulator_wheel    
        except:
            return "Broadbeam instance"
        

    def define_geo_relationships(self):

        #ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')
        self.geo = geometry.OptisGeometry(self.config)

        # self.central_axis = ref.central_axis

        aperture_to_iso = self.config["AID"]
        self.wedge_angle = self.config["wedge_angle"]

        self.slope_vps_aperture = 0.5 * \
            self.config["Aperture"] / \
            (abs(self.config["VPS"]) - self.config["AID"])
        self.phantom = patches.Rectangle((self.medium.mesh_origin[2], self.medium.mesh_origin[0]), self.config["Mesh dimensions"][2], self.config["Mesh dimensions"][0], edgecolor='none',
                                         facecolor='b', label="Phantom size x,y,z= {},{},{} mm".format(self.config["Mesh dimensions"][0], self.config["Mesh dimensions"][1], self.config["Mesh dimensions"][2]), alpha=0.05)

        self.inflection_xes = np.asarray(
            [- self.config["SID"], 0,  - self.config["SID"] + self.config["Mesh dimensions"][2]])
        self.inflection_xes = np.asarray(
            [0, 0,   self.config["Mesh dimensions"][2]])       
        
        
        self.inflection_yes = self.nozzle.inflection_shift_ratio[0] * np.asarray([self.slope_vps_aperture*(abs(self.config["VPS"]) + self.inflection_xes[0]),
                                                                                  self.slope_vps_aperture *
                                                                                  (abs(
                                                                                      self.config["VPS"]) + self.inflection_xes[1]),
                                                                                  self.slope_vps_aperture*(abs(self.config["VPS"]) + self.inflection_xes[2])])



        self.aperture_line_xes = np.asarray([self.central_axis.S[2],  - self.config["AID"],  -
                                            self.config["SID"],  - self.config["SID"] + self.config["Mesh dimensions"][2]])
        self.aperture_line_yes = np.asarray([self.central_axis.S[0], self.config["Aperture"]/2, self.slope_vps_aperture*abs(
            self.config["VPS"]), self.slope_vps_aperture*(abs(self.config["VPS"]) - self.config["SID"] + self.config["Mesh dimensions"][2])])

        self.aperture_point1 = self.central_axis.S
        self.aperture_point2a = self.central_axis.T + bev.xvec_in_ref*self.slope_vps_aperture * \
            math.dist(self.central_axis.T, self.central_axis.S) + bev.yvec_in_ref * \
            self.slope_vps_aperture * \
            math.dist(self.central_axis.T, self.central_axis.S)
        self.aperture_point2b = self.central_axis.T - bev.xvec_in_ref*self.slope_vps_aperture * \
            math.dist(self.central_axis.T, self.central_axis.S) - bev.yvec_in_ref * \
            self.slope_vps_aperture * \
            math.dist(self.central_axis.T, self.central_axis.S)

        # self.aperture_line_xes = np.asarray([self.central_axis.S[2],  - self.config["AID"],  - self.config["SID"],  - self.config["SID"] + self.config["Mesh dimensions"][2] ])
        # self.aperture_line_yes = np.asarray([self.central_axis.S[0], self.config["Aperture"]/2, self.slope_vps_aperture*abs(self.config["VPS"]), self.slope_vps_aperture*(abs(self.config["VPS"]) - self.config["SID"] + self.config["Mesh dimensions"][2]) ]) [0]]], facecolor='gold')

        self.beta_xes = self.inflection_xes
        self.beta_yes = np.asarray([self.nozzle.beta_m - self.nozzle.beta_k*(self.inflection_xes[1] - self.inflection_xes[0]),
                                   self.nozzle.beta_m, self.nozzle.beta_m + self.nozzle.beta_k*(self.inflection_xes[2] - self.inflection_xes[1])])

        # TODO:
        # Need getter and setters which should return 0 interpolating outside of range
        frame_shift = - self.config["SID"]
        self.f_half = interp1d(self.aperture_line_xes - frame_shift, self.aperture_line_yes *
                               self.nozzle.inflection_shift_ratio, bounds_error=False, fill_value="extrapolate")
        self.f_beta = interp1d(self.beta_xes, self.beta_yes,
                               bounds_error=False, fill_value="extrapolate"
                               )
        if "penumbra_factor" in self.config.keys():
            self.penumbra_factor = self.config["penumbra_factor"]
            self.logger.info('Lateral penumbra from config')
        else:
            self.penumbra_factor = 1
            self.logger.info('Lateral penumbra default')

        self.target_plane_p1x = self.central_axis.T + \
            bev.xvec_in_ref*self.config["Mesh dimensions"][0]*1.3/2
        self.target_plane_p2x = self.central_axis.T - \
            bev.xvec_in_ref*self.config["Mesh dimensions"][0]*1.3/2
        self.aperture_plane_p1x = self.geo.aperture_point + \
            bev.xvec_in_ref*self.config["Mesh dimensions"][0]*1.3/2
        self.aperture_plane_p2x = self.geo.aperture_point - \
            bev.xvec_in_ref*self.config["Mesh dimensions"][0]*1.3/2
        self.wedge_apex_point = self.geo.wedge_apex_point
        self.wedge_top_point = self.geo.wedge_top_point

    # def define_geo_relationships(self):

    #     ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')

    #     aperture_to_iso = self.config["AID"]
    #     self.wedge_angle = self.config["wedge_angle"]

    #     self.slope_vps_aperture = 0.5*self.config["Aperture"]/ (abs(self.config["VPS"]) - self.config["AID"] )
    #     self.phantom = patches.Rectangle((  self.medium.mesh_origin[2], self.medium.mesh_origin[0]), self.config["Mesh dimensions"][2],self.config["Mesh dimensions"][0], edgecolor='none', facecolor='b', label = "Phantom size x,y,z= {},{},{} mm".format(self.config["Mesh dimensions"][0], self.config["Mesh dimensions"][1], self.config["Mesh dimensions"][2]), alpha = 0.05)

    #     self.inflection_xes = np.asarray([ - self.config["SID"], 0,  - self.config["SID"] + self.config["Mesh dimensions"][2]])
    #     self.inflection_yes = self.nozzle.inflection_shift_ratio[0] *np.asarray([self.slope_vps_aperture*(abs(self.config["VPS"]) + self.inflection_xes[0]),
    #                                                                              self.slope_vps_aperture*(abs(self.config["VPS"]) + self.inflection_xes[1]),
    #                                                                              self.slope_vps_aperture*(abs(self.config["VPS"]) + self.inflection_xes[2])])

    #     #TODO
    #     # Transformation matrix from treamtent room to BEV
    #     self.T_treatment_bev = np.linalg.solve(self.treatmentroom,self.bev).T

    #     # Transformation matrix from BEV to treatment room
    #     T_bev_treatment = np.linalg.solve(self.bev, self.treatmentroom).T

    #     assert (self.T_treatment_bev * self.treatmentroom.T == self.bev.T).all()
    #     #assert (self.T_treatment_bev * self.treatmentroom_x == self.bev_x).all()

    #     # print(self.T_treatment_bev * self.treatmentroom_x , self.bev_x)

    #     self.aperture_plane = self.central_axis.vec_norm
    #     self.aperture_point = [aperture_to_iso*math.sin(self.config["Gantry rot theta"]*math.pi/180), 0, -aperture_to_iso*math.cos(self.config["Gantry rot theta"]*math.pi/180)]

    #     half_aperture = self.config["Aperture"]/2

    #     # x and y only
    #     insert_point = half_aperture *np.asarray([math.cos(self.config["wedge_insert_angle"]*math.pi/180), math.sin(self.config["wedge_insert_angle"]*math.pi/180)])
    #     insert_direction = -insert_point
    #     insert_direction = insert_direction /np.linalg.norm(insert_direction)

    #     #self.wedge_insert_point = self.aperture_point + half_aperture*self.treatmentroom_x*math.cos(self.config["wedge_insert_angle"]*math.pi/180) + half_aperture*self.treatmentroom_y*math.sin(self.config["wedge_insert_angle"]*math.pi/180)

    #     # self.wedge_insert_point = self.aperture_point - insert_direction *half_aperture

    #     wedge_apex_xy = insert_point + insert_direction*self.config["wedge_cover"]
    #     #self.wedge_apex_point = [wedge_apex_xy[0], wedge_apex_xy[1], -aperture_to_iso * math.cos(self.config["wedge_insert_angle"])]

    #     self.wedge_apex_point = self.aperture_point + azimut_x*math.cos(self.config["wedge_insert_angle"]*math.pi/180)*(self.config["wedge_cover"] - half_aperture)

    #     self.insert_direction = self.wedge_apex_point - self.aperture_point
    #     self.insert_direction = self.insert_direction/np.linalg.norm(self.insert_direction)
    #     #self.wedge_insert_point = self.aperture_point - self.insert_direction *half_aperture

    #     self.wedge_insert_point = self.aperture_point - azimut_x*(half_aperture + 10)*math.cos(self.config["wedge_insert_angle"]*math.pi/180)
    #     self.wedge_insert_direction = self.aperture_point - self.wedge_insert_point

    #     # distance between aperture point and wedge point in z
    #     wedge_apex_r = np.hypot(wedge_apex_xy[0], wedge_apex_xy[1]) # not global
    #     self.wedge_d = wedge_apex_r * math.tan(self.config["wedge_angle"]*math.pi/180) #not global
    #     self.wedge_point = self.aperture_point + self.central_axis.vec_norm*self.wedge_d

    #     # distance_wedgepoint_to_top = self.config["wedge_cover"] * math.tan(self.config["wedge_angle"]*math.pi/180) - self.wedge_d
    #     #self.collimator_end_point = self.aperture_point - self.insert_direction *(half_aperture + 10)

    #     wedge_vector = self.wedge_point - self.wedge_apex_point
    #     self.wedge_vector = wedge_vector/np.linalg.norm(wedge_vector)

    #     #self.wedge_top_point = self.wedge_apex_point + self.wedge_vector *(10 + self.config["wedge_cover"])/math.cos(self.config["wedge_angle"]*math.pi/180)

    #     self.wedge_top_point = self.wedge_insert_point + self.central_axis.vec_norm*math.sin(self.config["wedge_angle"]*math.pi/180)*(self.config["wedge_cover"] + 10)

    #     #------------------------------------------

    #     # point_y = self.wedge_point + azimut_y

    #     # points = [[self.wedge_top_point],
    #     #            [self.wedge_apex_point],
    #     #            [point_y]]

    #     p0, p1 = self.wedge_apex_point, self.wedge_top_point
    #     p2 = self.wedge_apex_point + np.cross(self.central_axis.vec_norm, self.wedge_insert_direction)

    #     # p0, p1, p2 = points

    #     x0, y0, z0 = p0
    #     x1, y1, z1 = p1
    #     x2, y2, z2 = p2

    #     ux, uy, uz = u = [x1-x0, y1-y0, z1-z0]
    #     vx, vy, vz = v = [x2-x0, y2-y0, z2-z0]

    #     u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]

    #     point  = np.array(p0)
    #     normal = np.array(u_cross_v)
    #     normal = normal/np.linalg.norm(normal)

    #     self.wedge_plane = normal
    #     #-------------------------------------------

    #     self.aperture_line_xes = np.asarray([self.central_axis.S[2],  - self.config["AID"],  - self.config["SID"],  - self.config["SID"] + self.config["Mesh dimensions"][2] ])
    #     self.aperture_line_yes = np.asarray([self.central_axis.S[0], self.config["Aperture"]/2, self.slope_vps_aperture*abs(self.config["VPS"]), self.slope_vps_aperture*(abs(self.config["VPS"]) - self.config["SID"] + self.config["Mesh dimensions"][2]) ])

    #     self.aperture_line_xes = np.asarray([self.central_axis.S[2],  - self.config["AID"],  - self.config["SID"],  - self.config["SID"] + self.config["Mesh dimensions"][2] ])
    #     self.aperture_line_yes = np.asarray([self.central_axis.S[0], self.config["Aperture"]/2, self.slope_vps_aperture*abs(self.config["VPS"]), self.slope_vps_aperture*(abs(self.config["VPS"]) - self.config["SID"] + self.config["Mesh dimensions"][2]) ])

    #     colli_depth = 50
    #     # self.collimator_top = patches.Rectangle(( - 70 -colli_depth, self.config["Aperture"]/2), colli_depth,10, edgecolor='none', facecolor='k', label = "Collimator Aperture = {} mm".format(self.config["Aperture"]))
    #     # self.collimator_bot = patches.Rectangle(( - 70 -colli_depth, -self.config["Aperture"]/2), colli_depth,-10, edgecolor='none', facecolor='k')

    #     corner1 = self.aperture_point + azimut_x *self.config["Aperture"]/2 - self.central_axis.vec_norm*colli_depth
    #     corner2 = self.aperture_point - azimut_x *self.config["Aperture"]/2 - self.central_axis.vec_norm*colli_depth

    #     self.collimator_top = patches.Rectangle(( corner2[2], corner2[0]), colli_depth, 10, edgecolor='none', facecolor='k', label = "Collimator Aperture = {} mm".format(self.config["Aperture"]), angle= -self.config["Gantry rot theta"])
    #     self.collimator_bot = patches.Rectangle(( corner1[2], corner1[0]), colli_depth,-10, edgecolor='none', facecolor='k', angle= - self.config["Gantry rot theta"])

    #     # wedge_length = 30
    #     #self.wedge_patch = patches.Wedge(( corner_wedge[2], corner_wedge[0]), wedge_length,90 + self.config["Gantry rot theta"] - self.config["wedge_angle"], 90 + self.config["Gantry rot theta"] ,edgecolor='none', facecolor='gold')
    #     self.wedge_patch = patches.Polygon([[self.wedge_apex_point[2], self.wedge_apex_point[0]], [self.wedge_top_point[2], self.wedge_top_point[0]], [self.wedge_insert_point[2],self.wedge_insert_point
    #                                                                                                                                                    [0]]], facecolor='gold')

    #     self.beta_xes = self.inflection_xes
    #     self.beta_yes = np.asarray([self.nozzle.beta_m - self.nozzle.beta_k*(self.inflection_xes[1] - self.inflection_xes[0]), self.nozzle.beta_m, self.nozzle.beta_m + self.nozzle.beta_k*(self.inflection_xes[2] - self.inflection_xes[1])])

    #     # TODO:
    #     # Need getter and setters which should return 0 interpolating outside of range
    #     frame_shift =  - self.config["SID"]
    #     self.f_half = interp1d(self.aperture_line_xes - frame_shift, self.aperture_line_yes* self.nozzle.inflection_shift_ratio , bounds_error=False, fill_value="extrapolate")
    #     self.f_beta = interp1d(self.beta_xes - frame_shift, self.beta_yes , bounds_error=False, fill_value="extrapolate")

    #     #For plotting
    #     self.target_plane_p1x = self.central_axis.T + azimut_x*self.config["Mesh dimensions"][0]*1.3/2
    #     self.target_plane_p2x = self.central_axis.T - azimut_x*self.config["Mesh dimensions"][0]*1.3/2
    #     self.aperture_plane_p1x = self.aperture_point + azimut_x*self.config["Mesh dimensions"][0]*1.3/2
    #     self.aperture_plane_p2x = self.aperture_point - azimut_x*self.config["Mesh dimensions"][0]*1.3/2
    #     # self.wedge_plane_p1x = [wedge_apex_xy[0], wedge_apex_xy[1], self.aperture_point[2]]
    #     # #self.wedge_point - self.wedge_vector * #
    #     # # self.wedge_plane_p2x = self.wedge_point
    #     # self.wedge_plane_p2x = [self.config["Aperture"]/2, 0 , self.aperture_point[2] + self.config["wedge_cover"]*math.tan(self.config["wedge_angle"]*math.pi/180)]
    #     self.wedge_plane_p1x = self.wedge_apex_point
    #     self.wedge_plane_p2x = self.wedge_point
    #     self.wedge_plane_p3x = self.wedge_top_point
    # @timer
    def config_target_depth(self):
        
        if self.raytracer_configured == False:
            print("Configure raytracer first")
            return
        
        min_depth = 100000
        max_depth = 0

        indicis = np.where(self.target_matrix)

        for x, y, z in zip(indicis[0], indicis[1], indicis[2]):
            if self.raytracer.depth[x, y, z] < min_depth and [x, y, z] in self.raytracer.traced_indices:
                min_depth = self.raytracer.depth[x, y, z]

            if self.raytracer.depth[x, y, z] > max_depth:
                max_depth = self.raytracer.depth[x, y, z]

        self.config["Target_range"] = max_depth  # mm
        self.config["Modulation_range"] = max_depth - min_depth  # mm        
    
        self.logger.info(
            f'Target depth configured from a max depth of {max_depth} and a modulation range of {self.config["Modulation_range"]}')
    # @timer
    def config_nozzle(self):

        self.nozzle.config_nozzle(self.config)

        self.logger.info(
            f'Nozzle configured for modulator wheel {self.nozzle.modulator_wheel} and scattering foil {self.nozzle.foil}')

        self.logger.info('Nozzle configured for a target range of {} mm and a modulation range of {} mm'.format(
            round(self.config["Target_range"], 3), round(self.config["Modulation_range"], 3)))

        self.define_geo_relationships()
        self.nozzle_configured = True
    # @timer
    def config_nozzle_by_mw(self, wheelNo):
        self.nozzle.config_nozzle_by_mw(self.config, wheelNo)
        self.logger.info(
            f'Nozzle configured for modulator wheel {wheelNo} and scattering foil {self.nozzle.foil}')

        self.define_geo_relationships()
        self.nozzle_configured = True
    # @timer
    def config_target(self):
        self.target_matrix = self.config["target_points"]

        com_bins = np.asarray(
            ndimage.measurements.center_of_mass(self.target_matrix))
        translation_vector = com_bins*self.medium.resolution + self.medium.mesh_origin
        self.medium.translate(- translation_vector)
        self.translation_vector = - translation_vector

        #self.target_voxel_com = np.asarray(ndimage.measurements.center_of_mass(self.target_matrix))
        #self.target_com = com_bins*self.medium.resolution + self.medium.mesh_origin

        self.logger.info(f'Target centered by {- translation_vector} mm')

        indicis = np.where(self.target_matrix)

        assert len(indicis[0]) == len(indicis[1]) == len(indicis[2])

        for x, y, z in zip(indicis[0], indicis[1], indicis[2]):

            x_pos = self.medium.resolution[0] * \
                (x + 0.5) + self.medium.mesh_origin[0]
            y_pos = self.medium.resolution[1] * \
                (y + 0.5) + self.medium.mesh_origin[1]
            z_pos = self.medium.resolution[2] * \
                (z + 0.5) + self.medium.mesh_origin[2]

            self.indices_target.append(
                [x_pos, y_pos, z_pos])

        # # Depth min and max along central axis
        # min_target_S_distance = 100000
        # min_target_T_distance = 100000
        # self.min_target_point = []
        # self.max_target_point = []

        for indices in self.indices_target:
            x, y, z = indices

            T = [x, y, z]
            self._target_rays.append(dose_engine.Ray(self.central_axis.S, T))

            p = np.asarray(T)

            # normalized tangent vector
            d = self.central_axis.vec_norm

            # signed parallel distance components
            # s = np.dot(self.central_axis.S - p, d)
            t = np.dot(p - self.central_axis.T, d)

            p_on_axis = self.central_axis.T + t*self.central_axis.vec_norm

            # if math.dist(self.central_axis.S, p_on_axis) < min_target_S_distance:
            #     min_target_S_distance = math.dist(self.central_axis.S, p_on_axis)
            #     self.min_target_point = p_on_axis

            # if math.dist(self.central_axis.T, p_on_axis) < min_target_T_distance:
            #     min_target_T_distance = math.dist(self.central_axis.T, p_on_axis)
            #     self.max_target_point = p_on_axis

            dist = math.dist(p, p_on_axis)

            if dist > self.orthogonal_limit:
                self.orthogonal_limit = dist

        # self.target_range = math.dist(self.central_axis.S, self.max_target_point)
        # self.modulation = math.dist(self.central_axis.T, self.max_target_point)

        #expands the orhogonal limitfor taget rays by a factor
        self.orthogonal_limit = self.orthogonal_limit*2

        self.dist_aperture = - self.dist_aperture*self.orthogonal_limit

        self.target_configured = True
        self.logger.info('Target rays defined')

    def define_artifical_target(self, case=1, center_target = 1):

        #mid_point = np.asarray([10, 0, -34])

        mid_point = self.medium.mesh_origin + self.medium.resolution*self.dose.shape/2

        mid_point2 = np.asarray([10, 7, -34])

        radius = 5

        for x in range(self.target_matrix.shape[0]):
            for y in range(self.target_matrix.shape[1]):
                for z in range(self.target_matrix.shape[2]):

                    x_pos = self.medium.resolution[0] * \
                        (x + 0.5) + self.medium.mesh_origin[0]
                    y_pos = self.medium.resolution[1] * \
                        (y + 0.5) + self.medium.mesh_origin[1]
                    z_pos = self.medium.resolution[2] * \
                        (z + 0.5) + self.medium.mesh_origin[2]

                    if case == 1:
                        if math.sqrt((x_pos - mid_point[0])**2 + ((y_pos - mid_point[1]))**2 + (z_pos - mid_point[2])**2) < radius:
                            if math.sqrt((x_pos - mid_point2[0])**2 + (y_pos - mid_point2[1])**2) < radius*1.2:
                                pass
                            else:
                                self.target_matrix[x, y, z] = 1
                                self.indices_target.append(
                                    [x_pos, y_pos, z_pos])
                    elif case == 2:
                        if math.sqrt(((x_pos - mid_point[0]))**2 + ((y_pos - mid_point[1]))**2 + (z_pos - mid_point[2])**2) < radius*1.8:
                            if abs(x_pos - mid_point[0]) < 3:
                                self.target_matrix[x, y, z] = 1
                                self.indices_target.append(
                                    [x_pos, y_pos, z_pos])
                    elif case == 3:
                        if math.sqrt(0.2*((x_pos - mid_point[0]))**2 + 1.8*((y_pos - mid_point[1]))**2 + (z_pos - mid_point[2])**2) < radius:
                            self.target_matrix[x, y, z] = 1
                            self.indices_target.append(
                                [x_pos, y_pos, z_pos])

        assert np.sum(self.target_matrix) != 0
        
        if center_target:
            translation_vector = - mid_point - np.zeros(mid_point.shape)
            self.medium.translate(translation_vector)
            self.translation_vector = - translation_vector
        else:
            self.translation_vector = np.asarray([0,0,0])

        for indices in self.indices_target:
            x, y, z = indices

            T = [x - mid_point[0], y - mid_point[1], z - mid_point[2]]
            self._target_rays.append(dose_engine.Ray(self.central_axis.S, T))

            p = np.asarray(T)

            # normalized tangent vector
            d = self.central_axis.vec_norm

            # signed parallel distance components
            # s = np.dot(self.central_axis.S - p, d)
            t = np.dot(p - self.central_axis.T, d)

            p_on_axis = self.central_axis.T + t*self.central_axis.vec_norm

            dist = math.dist(p, p_on_axis)

            if dist > self.orthogonal_limit:
                self.orthogonal_limit = dist

        #expands the orhogonal limit by a factor
        self.orthogonal_limit = self.orthogonal_limit*2

        self.dist_aperture = - self.dist_aperture*self.orthogonal_limit

        self.target_configured = True
        
        if center_target:
            self.logger.info(f'Target centered by {- self.translation_vector} mm')
        self.logger.info(f'Defined {len(self._target_rays)} target rays')

    def config_medium(self):
        if self.config["Medium"] == "Phantom":
            self.medium = dose_engine.Phantom(self.config)
            self.logger.info('Medium configured as a phantom')
        elif self.config["Medium"] == "Dicom":
            self.medium = dose_engine.Dicom(self.config)
            self.logger.info(
                'Medium configured based on Sclera binary mask from dicom import')
        else:
            assert 0, "Only defined for Phantom or Dicom for now"

        self.medium_configued = True
    
    # @timer
    def config_collimator(self):
        
        if self.nozzle_configured == False:
            raise AssertionError("Config nozzle first")
        

        if len(self._target_rays) == 0:
            self.logger.info('No rays for target. Collimator NOT configured')
            return
        
        expansion = self.config["collimator_expansion"]
        
        resulting_expansion = self.nozzle.inflection_shift_mm[0] + expansion
        
        self.collimator.add_ray(self._target_rays)
        self.collimator.define_aperture_hull()
        self.collimator.ruler_length = self.config["collimator_ruler_length"]
        self.collimator.expand_aperture(resulting_expansion)

        self.collimator_configured = True
        
        self.logger.info(f'Collimator configured and expanded by {round(resulting_expansion, 4)} mm')

    def config_collimator_from_eyeplan_model(self, eyeplan_points_in_plane, expansion=0):

        # if len(self._target_rays) == 0:
        #     self.logger.info('No rays for target. Collimator NOT configured')
        #     return
        
        resulting_expansion = self.nozzle.inflection_shift_mm[0] + expansion
        
        self.collimator.add_ray(self._target_rays)
        self.collimator.define_aperture_hull_from_points(
            eyeplan_points_in_plane)
        self.collimator.ruler_length = self.config["collimator_ruler_length"]
        self.collimator.expand_aperture(resulting_expansion)

        self.collimator_configured = True

        self.logger.info('Collimator configured from eyeplan model')
        self.logger.info(f'Collimator expanded by {resulting_expansion} mm')

    @timer
    def config_distances_to_aperture(self):

        if not self.collimator.is_configured:
            self.logger.info(
                'Distances to aperture NOT calculated. Configure collimator first')
            return

        self.distances_to_aperture_filename = f"distances_to_aperture_{self.config}"

        self.fq_filename = "first_quadrant" + \
            self.distances_to_aperture_filename[21:]

        if self.distances_to_aperture_filename in os.listdir(self.support_data_path):
            self.dist_aperture = np.load(
                self.support_data_path / self.distances_to_aperture_filename)
            self.first_quadrant_bev = np.load(
                self.support_data_path / self.fq_filename)
            self.logger.info('Distances to aperture values loaded from file')
            return

        hull_in_planes = {}

        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        for indices in self.raytracer.traced_indices:
            x, y, z = indices

            x_pos = self.medium.resolution[0] * \
                (x + 0.5) + self.medium.mesh_origin[0]
            y_pos = self.medium.resolution[1] * \
                (y + 0.5) + self.medium.mesh_origin[1]
            z_pos = self.medium.resolution[2] * \
                (z + 0.5) + self.medium.mesh_origin[2]

            p = np.asarray([x_pos, y_pos, z_pos])
            


            # normalized tangent vector
            d = self.central_axis.vec_norm

            # signed parallel distance components
            #s = np.dot(self.central_axis.S - p, d)
            t = np.dot(p - self.central_axis.T, d)

            p_on_axis = self.central_axis.T + t*self.central_axis.vec_norm

            #dist_to_iso = self.config["VPS"] + t

            vector = p - p_on_axis

            x_in_plane = - np.dot(bev.xvec_in_ref, vector)
            y_in_plane = np.dot(bev.yvec_in_ref, vector)

            #hull_filename = f"hull_{self.config}_{p_on_axis}"

            if tuple(p_on_axis) in hull_in_planes:
                hull_in_plane = hull_in_planes[tuple(p_on_axis)]
            else:
                hull_in_plane = self.collimator.project_aperture(p_on_axis)
                hull_in_planes[tuple(p_on_axis)] = hull_in_plane
            
    
            point = Point(x_in_plane, y_in_plane)
            self.dist_aperture[x, y, z] = - \
                hull_in_plane.boundary.distance(point)

            if hull_in_plane.boundary.contains(point):
                self.dist_aperture[x, y, z] = hull_in_plane.boundary.exterior.distance(
                    point)

            if np.sign(vector[0]) == np.sign(bev.xvec_in_ref[0]) and np.sign(vector[1]) == np.sign(bev.yvec_in_ref[1]):
                self.first_quadrant_bev[x, y, z] = 1
                
                

        self.distances_to_aperture_configured = True

        self.logger.info(
            f'Distances to aperture calculated for {self.Nvoxels} voxels')

    #@timer

    def config_radii(self):
        """
        Calculates the shortest distance to the central axis and stores values in radii cube
        
        """

        if self.config['Image'] == '3D':
            self.radii_filename = f"Radii_{self.config['Image']}_Nvoxels_{self.config['Nvoxels']}_Dimensions_{self.config['Mesh dimensions']}_theta_{self.config['Gantry rot theta']}_phi_{self.config['Gantry rot phi']}_{self.target_configured}.npy"
        elif self.config['Image'] == 'xslice':
            self.radii_filename = f"Radii_{self.config['Image']}_{self.config['Slice'][0]}_Nvoxels_{self.config['Nvoxels']}_Dimensions_{self.config['Mesh dimensions']}_theta_{self.config['Gantry rot theta']}_phi_{self.config['Gantry rot phi']}_{self.target_configured}.npy"
        elif self.config['Image'] == 'yslice':
            self.radii_filename = f"Radii_{self.config['Image']}_{self.config['Slice'][1]}_Nvoxels_{self.config['Nvoxels']}_Dimensions_{self.config['Mesh dimensions']}_theta_{self.config['Gantry rot theta']}_phi_{self.config['Gantry rot phi']}_{self.target_configured}.npy"
        elif self.config['Image'] == 'zslice':
            self.radii_filename = f"Radii_{self.config['Image']}_{self.config['Slice'][2]}_Nvoxels_{self.config['Nvoxels']}_Dimensions_{self.config['Mesh dimensions']}_theta_{self.config['Gantry rot theta']}_phi_{self.config['Gantry rot phi']}_{self.target_configured}.npy"

        self.fq_filename = "first_quadrant" + self.radii_filename[5:]

        if self.radii_filename in os.listdir(self.support_data_path):
            self.radii = np.load(self.support_data_path / self.radii_filename)
            self.first_quadrant_bev = np.load(
                self.support_data_path / self.fq_filename)
            self.logger.info('Radii loaded from file')
            return

        origin = self.medium.mesh_origin
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        S = self.central_axis.S
        T = self.central_axis.T

        for indices in self.raytracer.traced_indices:
            x, y, z = indices

            x_pos = self.medium.resolution[0]*(x + 0.5) + origin[0]
            y_pos = self.medium.resolution[1]*(y + 0.5) + origin[1]
            z_pos = self.medium.resolution[2]*(z + 0.5) + origin[2]

            p = np.asarray([x_pos, y_pos, z_pos])

            # normalized tangent vector
            # d = np.divide(T - S, np.linalg.norm(T - S))
            d = self.central_axis.vec_norm

            # signed parallel distance components
            s = np.dot(S - p, d)
            t = np.dot(p - T, d)

            p_axis = self.central_axis.T + t*self.central_axis.vec_norm
            r = math.dist(p, p_axis)

            vector = p - p_axis

            if np.sign(vector[0]) == np.sign(bev.xvec_in_ref[0]) and np.sign(vector[1]) == np.sign(bev.yvec_in_ref[1]):
                self.first_quadrant_bev[x, y, z] = 1

            # # clamped parallel distance
            # h = np.maximum.reduce([s, t, 0])
            # # perpendicular distance component
            # c = np.cross(p - S, d)

            self.radii[x, y, z] = r  # np.hypot(h, np.linalg.norm(c))

        assert np.amax(self.radii) != 0

        np.save(self.support_data_path / self.radii_filename, self.radii)
        np.save(self.support_data_path /
                self.fq_filename, self.first_quadrant_bev)

        self.radii_configured = True

        self.logger.info('Radii calculated and saved to file')



    def define_target_rays_from_points(self):
        
        geo = geometry.OptisGeometry(self.config)

        for p in self.config["target_points"]:
            ray1 = dose_engine.Ray(geo.central_axis.S, p)
            self._target_rays.append(ray1)
        
        return self._target_rays
        


    # @timer
    def define_rays(self,  sparse=0):

        if self.orthogonal_limit == 0:
            self.logger.info('No rays defined')
            return

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        # limit = 0.8*self.config["Mesh dimensions"][0]/2
        # limit = 3*self.orthogonal_limit

        if sparse:
            Nsteps = 6
        else:
            Nsteps = int(2*self.orthogonal_limit/min(self.medium.resolution))

        if self.config["Image"] == "3D" or self.config["Image"] == "zslice":
            for y in np.linspace(-self.orthogonal_limit, self.orthogonal_limit, num=Nsteps):
                for x in np.linspace(-self.orthogonal_limit, self.orthogonal_limit, num=Nsteps):
                    self.rays.append(dose_engine.Ray(
                        self.central_axis.S, self.central_axis.T + x*bev.xvec_in_ref + y*bev.yvec_in_ref))

        elif self.config["Image"] == "xslice":
            y = 0
            for x in np.linspace(-self.orthogonal_limit, self.orthogonal_limit, num=Nsteps):
                self.rays.append(dose_engine.Ray(
                    self.central_axis.S, self.central_axis.T + x*bev.xvec_in_ref + y*bev.yvec_in_ref))
        elif self.config["Image"] == "yslice":
            x = 0
            for y in np.linspace(-self.orthogonal_limit, self.orthogonal_limit, num=Nsteps):
                self.rays.append(dose_engine.Ray(
                    self.central_axis.S, self.central_axis.T + x*bev.xvec_in_ref + y*bev.yvec_in_ref))

        self.logger.info(f'Defined {len(self.rays)} rays')



    # def config_surface_interaction(self):
        
    #     if "eyeglobe_point" not in self.config.keys():
    #         self.logger.info('No point cloud present to be used for range retraction')
    #         return
        
    #     points = self.config["eyeglobe_point"]
        
    #     cloud = pv.PolyData(points)
    #     volume = cloud.delaunay_3d(alpha=12.)
    #     # volume.plot()
    #     shell = volume.extract_geometry()
        
        
    #     surface_shell = None
    #     if "patient_surface_points" in self.config.keys():
    #         surface_points = self.config["patient_surface_points"]
    #         surface_cloud = pv.PolyData(surface_points)
    #         surface_volume = surface_cloud.delaunay_3d(alpha=12.)
    #         # volume.plot()
    #         surface_shell = surface_volume.extract_geometry()
        
    #     self.rays.append(self.central_axis)
        
    #     for ray in self.rays:
    #         points, ind = shell.ray_trace(ray.S, ray.T)
            
    #         if len(points) == 0:
    #             continue
            
    #         shortest_distance = 1000000
    #         first_intersect = []
    #         for p in points:
    #             dist = math.dist(p, ray.S)
    #             if dist < shortest_distance:
    #                 shortest_distance = dist
    #                 first_intersect = p
            
    #         cloud_surface_interaction_point = np.asarray(first_intersect)
            
            
    #         d3 = 10000000
    #         if type(surface_shell) == pv.core.pointset.PolyData:
    #             points, ind = surface_shell.ray_trace(ray.S, ray.T)
                
    #             if len(points) == 0:
    #                 continue
                
    #             shortest_distance = 1000000
    #             first_intersect_surface = []
    #             for p in points:
    #                 dist = math.dist(p, ray.S)
    #                 if dist < shortest_distance:
    #                     shortest_distance = dist
    #                     first_intersect_surface = p    
                    
    #             cloud_surface_interaction_point2 = np.asarray(first_intersect_surface)        
                
    #             d3 = math.dist(cloud_surface_interaction_point2, ray.S)
            
    #         # if "PSI_model" in self.config.keys():
    #         #     p = self.config["skin_plane_point"]
    #         #     p[2] = p[2] - 8
    #         #     skin_plane_intersect = ray.find_intersect_with_plane( self.config["skin_plane_normal"] , self.config["skin_plane_point"] )
    #         # else:
    #         skin_plane_intersect = ray.find_intersect_with_plane( self.config["skin_plane_normal"] , self.config["skin_plane_point"] )
            
    #         d1, d2 = math.dist(cloud_surface_interaction_point, ray.S), math.dist(skin_plane_intersect, ray.S)
            
    #         # if d3 < d2:
    #         #     ray.surface_interaction_point = cloud_surface_interaction_point2
    #         #     ray.add_depth(d3 - d1)
    #         #     ray.retracted_by_cloud = True   
            
    #         if d1 < d2:
    #             ray.surface_interaction_point = cloud_surface_interaction_point
    #             ray.add_depth(d2 - d1)
    #             ray.retracted_by_cloud = True
    #         else:
    #             ray.surface_interaction_point = skin_plane_intersect
        
    #     self._surface_interaction_configured = True
    #     self.logger.info('Range retracted from point cloud')
        
        
        
    def config_surface_interaction(self):
            
            
            # assert self.config["surfaces"][1].type == "cloud"
            
            for ray in self.rays:
                
                
                for surface in self.surfaces:
                    surface.find_ray_interaction(ray)
                
                # if -11 < ray.T[0] < -9 and -2 < ray.T[1] < -1 and self.config["surfaces"][0].type == "cloud":
                #     apa = 1
                
                # self.config["surfaces"][1].find_ray_interaction(ray)
            
                # print(ray.surface_interaction_point_candidates)
                # for surface in self.config["surfaces"]:
                    
                shortest_distance = 1000000
                first_intersect = []
                for p in ray.surface_interaction_point_candidates:
                    
                    if type(p) != np.ndarray:
                        if p == None:
                            continue
                    
                    dist = math.dist(p, ray.S)
                    if dist < shortest_distance:
                        shortest_distance = dist
                        first_intersect = p
                        
                first_surface_intersect = np.asarray(first_intersect)
                
                cloud_surface_intersect = first_surface_intersect
                
                skin_plane_intersect = ray.find_intersect_with_plane( self.config["skin_plane_normal"] , self.config["skin_plane_point"] )
                
                mesh_surface_intersect = ray.find_intersect_with_plane( self.config["skin_plane_normal"] , np.asarray([0,0,self.medium.maxZ]) )
                
                
                if len(first_surface_intersect) == 3:
                    ray.surface_intersection_identified = True
                
    
                
                # In case no interaction found
                elif len(first_surface_intersect) == 0:
                    ray.surface_intersection_identified = False
                    first_surface_intersect = skin_plane_intersect
                else:
                    raise ValueError
                    
                ray.surface_interaction_point = first_surface_intersect
                
                
                
                d1 = math.dist(first_surface_intersect, ray.S) 
                d2 = math.dist(skin_plane_intersect, ray.S)
                
                d3 = math.dist(mesh_surface_intersect, ray.S)
                
                #d4 = math.dist(cloud_surface_intersect, ray.S)
                
                surface_case = ""
                if d1 < d2:
                    surface_case = "cloud"
                    ray.add_depth(d2 - d1)
                    ray.retracted_by_cloud = True
                
                if d1 > d2:
                    # using measured surfaces as the first interaction. Neither plane or cloud first from the perspecitve of the ray
                    surface_case = "measured"
                    ray.surface_behind_skin_plane = True
                    ray.surface_behind_skin_plane_by = abs(d2 - d1) #mm
                
                if d1 == d2:
                    surface_case = "plane"
                
                if d3 < d1:
                    ray.surface_behind_mesh = True
                    ray.surface_behind_mesh_by = abs(d3 - d1) #mm
    
                # if self.config['anterior_aligned']:
                #     surface_case = ""
                ray.surface_case = surface_case
    
            # for ray1 in self._target_rays:
            #     ray1.find_surface_intersect(shell, self.config["skin_plane_normal"] , self.config["skin_plane_point"] )
                
            #     distance = math.dist(ray1.T, ray1.surface_interaction_point)
                
            #     if distance > deepest_range:
            #         deepest_range = distance
            #         deepest_ray = ray1
            #     if distance < shallowest_range:
            #         shallowest_range = distance
                    # shallowest_ray = ray1      
                
            # for surface in self.config["surfaces"]:
            #     surface.find_ray_interaction(ray)
            
            # deepest_range = 0
            # deepest_ray = []
            # shallowest_range = 100000000
            # shallowest_ray = []        
            # for ray1 in self._target_rays:
            #     ray1.find_surface_intersect(shell, self.config["skin_plane_normal"] , self.config["skin_plane_point"] )
                
            #     distance = math.dist(ray1.T, ray1.surface_interaction_point)
                
            #     if distance > deepest_range:
            #         deepest_range = distance
            #         deepest_ray = ray1
            #     if distance < shallowest_range:
            #         shallowest_range = distance
            #         shallowest_ray = ray1                
            
            # self.deepest_ray = deepest_ray
            # self.shallowest_ray = shallowest_ray
            
            
            
                
            self._surface_interaction_configured = True            
            
            self.logger.info(f'Surface interaction configured from {len(self.config["surfaces"])} patient surfaces')



    def determine_target_and_modulation_range_of_target(self):
        
        if "target_points" not in self.config.keys():
            self.logger.info('No point cloud of target present to determine target and modulation range')
            return
        
        points = self.config["eyeglobe_point"]
        
        cloud = pv.PolyData(points)
        volume = cloud.delaunay_3d(alpha=self.config["eyeglobe_mesh_triangle_size"] )
        # volume.plot()
        shell = volume.extract_geometry()
        
        assert len(self._target_rays) > 0
        
        deepest_range = 0
        deepest_ray = []
        shallowest_range = 100000000
        shallowest_ray = []        
        for ray1 in self._target_rays:
            # ray1.find_surface_intersect(shell, self.config["skin_plane_normal"] , self.config["skin_plane_point"] )
            
            
            for surface in self.config["surfaces"]:
                surface.find_ray_interaction(ray1)
                
            shortest_distance = 1000000
            first_intersect = []
            for p in ray1.surface_interaction_point_candidates:
                
                if type(p) != np.ndarray:
                    if p == None:
                        continue
                
                dist = math.dist(p, ray1.S)
                if dist < shortest_distance:
                    shortest_distance = dist
                    first_intersect = p
            
            ray1.surface_interaction_point = first_intersect
            distance = math.dist(ray1.T, ray1.surface_interaction_point)
            
            if distance > deepest_range:
                deepest_range = distance
                deepest_ray = ray1
            if distance < shallowest_range:
                shallowest_range = distance
                shallowest_ray = ray1                
        
        
        self.deepest_ray = deepest_ray
        self.shallowest_ray = shallowest_ray
        
        assert self.deepest_ray != None
        assert self.shallowest_ray != None
        
        if self.shallowest_ray.T[2] > self.config['skin_plane_point'][2]:
            self.target_prior_to_skin_plane = True
        
        self.target_radiological_depth = math.dist(self.deepest_ray.T, self.deepest_ray.surface_interaction_point)
        target_modulation= self.target_radiological_depth - math.dist(shallowest_ray.T, shallowest_ray.surface_interaction_point)        
        
        
        if self.config['Target_range'] != 0 or self.config['Modulation_range'] != 0:
            self.logger.info('Note: Target and modulation range overwritten')
            
        
        self.config['Target_range'] = self.target_radiological_depth + self.config["distal_margin"] #- math.dist(deepest_ray.surface_interaction_point, deepest_ray.skin_plane_point)
        self.config['Modulation_range'] = target_modulation + self.config["proximal_margin"] + self.config["distal_margin"] 
        
        
        self.expected_range = deepest_ray.surface_interaction_point[2] -  math.dist(deepest_ray.T + deepest_ray.vec_norm*self.config["distal_margin"] , deepest_ray.surface_interaction_point)
        self.expected_modulation = shallowest_ray.T[2] + self.config["proximal_margin"] #expected_range + algo.config['Modulation_range']
        
        self.logger.info(
            'Target and modulation range set from point cloud of target')

        return self.target_radiological_depth, target_modulation, deepest_ray, shallowest_ray

    def config_wedge(self):

        if self.config["wedge_angle"] == 0:
            self.logger.info('Wedge not traced')
            return

        self.wedge.add_ray(self.rays)
        self.wedge.trace_rays()

        self.wedge_configured = True

        self.logger.info(f'Wedge traced for {len(self.wedge.rays)} rays')

    def config_compensator(self, sphere_centre, radius):

        self.compensator = dose_engine.Compensator(self.config, sphere_centre, radius, 1)
        self.compensator.add_ray(self.rays)
        self.compensator.trace_rays()

        self.wedge.add_ray(self.rays)
        self.wedge.trace_rays()

        # self.wedge_configured = True

        self.logger.info(f'Compensator traced for {len(self.wedge.rays)} rays')

    def config_skin_plane(self):

        if self.config["skin_plane"] < self.medium.mesh_apex[2] or self.medium_configued == False:
            self.logger.info('Skin plane not traced')
            return

        self.skin_plane = dose_engine.SkinPlane(self.config, self.medium)
        self.skin_plane.add_ray(self.rays)
        self.skin_plane.trace_rays()

        self.skin_plane_configured = True

        self.logger.info(
            f'Skin plane traced for {len(self.skin_plane.rays)} rays')


    @timer
    def config_raytracer(self):
        
        self.raytracer.add_ray(self.rays)
        self.raytracer.trace_rays()
        
        if not self.raytracer.intersects_found:
            raise AssertionError("No intersects for raytracer")

        self.Nvoxels = len(self.raytracer.traced_indices)

        self.raytracer_configured = True

        self.logger.info(
            f'Raytracing completed for {len(self.raytracer.rays)} rays')
        self.logger.info(
            f'{self.Nvoxels} voxels have been indexed for dose calculation')
    # @timer
    def config_depth_dose_profile(self):
        """
        
        Note: Depth profile is currently extrapolated to make a longer profile in depth
              Only a single tune is considered now

        Returns
        -------
        None.

        """
        

        if self.nozzle_configured == False:
            print("Configure nozzle first")
            self.logger.info('Depth profile NOT configured')
            return

        inverse_square_correction = np.loadtxt(
            self.data_path / "BaseData" / f"Inverse_square_correction_{self.nozzle.foil}.out")
        f_isc = interp1d(
                  inverse_square_correction[0], inverse_square_correction[1], fill_value="extrapolate")
        

        if self.nozzle.tune == 'ESS16': #"new_profile" in self.config.keys() and  
                self.depthprofile_filename = self.data["DepthProfile_ESS16"]
                self.DD.load_scan(
                    self.depthprofile_filename, self.data_path / "DepthProfiles")
        elif self.nozzle.tune == 'ESS18':
                self.depthprofile_filename = "ESS18_z" #"km3-5_km5-2_dp03_12_-12_75.dat"
                self.DD.load_scan(
                    self.depthprofile_filename , self.data_path / "DepthProfiles")
        elif self.nozzle.tune == 'ESS20':
                self.depthprofile_filename = "km3-5_km5-2_dp05_16_-16_75.dat"
                self.DD.load_scan(
                    self.depthprofile_filename, self.data_path / "DepthProfiles")                
        else:
            self.depthprofile_filename = self.data["DepthProfile"]
            self.DD.load_scan(
                self.depthprofile_filename, self.data_path / "DepthProfiles")

        targeted_range = self.nozzle.target_range


        if self.DD.R90 - targeted_range < 0:
            self.DD.extrapolate()
            self.DD.calc()
            base_modulation = self.DD.R90 - targeted_range
            self.logger.info(f'Depth profile for tune {self.nozzle.tune} too short. Profile extrapolated')
        else:
            base_modulation = self.DD.R90 - targeted_range          
            
            
        assert base_modulation >= 0, "Target range larger than range of depth dose curve"
        
        
        self.SOBP_xes = self.DD.xes[:-2]
        self.SOBP_yes = np.zeros(len(self.SOBP_xes))
        

        for i in range(len(self.nozzle.MWwatereq)):

            tot_modulation = float(self.nozzle.MWwatereq[i])
            
            w = float(self.nozzle.MWweights[i])

            xes = np.linspace(
                min(self.DD.xes[:-1]), max(self.DD.xes[:-1]), num=1000) + tot_modulation

            yes = self.DD.f(xes)
            yes = w*f_isc(tot_modulation)*yes

            if yes[np.logical_not(np.isnan(yes))].size > 0:

                braggpeak = dose_engine.BraggPeak(
                    xes - tot_modulation, yes, w, tot_modulation)
                braggpeak.calc()

                self.braggpeaks.append(braggpeak)

        for braggpeak in self.braggpeaks:
            f = interp1d(braggpeak.xes, braggpeak.yes)
            self.SOBP_yes = self.SOBP_yes + f(self.SOBP_xes)

        SOBP_f_rev = interp1d(self.SOBP_yes, self.SOBP_xes)
        x_at_max = SOBP_f_rev(max(self.SOBP_yes))

        SOBP_f = interp1d(self.SOBP_xes, self.SOBP_yes)
        normfactor = SOBP_f(x_at_max)
                
        self.SOBP_yes = self.SOBP_yes/normfactor
        for braggpeak in self.braggpeaks:
            braggpeak.yes = braggpeak.yes/normfactor

        SOBP_pre = dose_engine.SOBP(self.SOBP_xes, self.SOBP_yes/normfactor, base_modulation, 0, 0)
        SOBP_pre.calc()
        
        diff_DD_SOBP = self.DD.R90 - SOBP_pre.R90
        
        self.SOBP_xes = self.SOBP_xes - base_modulation + diff_DD_SOBP
        x_deepest_peak = x_at_max - base_modulation + diff_DD_SOBP
        
        
        x_deepest_peak, normfactor_distal_edge = _find_norm_factor_from_the_right(self.SOBP_xes, self.SOBP_yes)
        self.SOBP_yes = self.SOBP_yes/normfactor_distal_edge
        for braggpeak in self.braggpeaks:
            braggpeak.yes = braggpeak.yes/normfactor_distal_edge
                
        # if self.nozzle.tune == 'ESS16' or "mid_normalization" in self.config:
        x_R90 = find_R90(self.SOBP_xes, self.SOBP_yes)
        SOBP_f = interp1d(self.SOBP_xes, self.SOBP_yes)
        x_at_mid = x_R90 - self.nozzle.modulation/2
        
        normfactor = SOBP_f(x_at_mid)
        self.SOBP_yes = self.SOBP_yes/normfactor
        for braggpeak in self.braggpeaks:
            braggpeak.yes = braggpeak.yes/normfactor
        #     braggpeak.xes -= x_shift
        
        x_shift = find_R90(self.SOBP_xes, self.SOBP_yes) - targeted_range
        self.SOBP_xes -= x_shift

        
        self.logger.info('SOBP normalized to centre of modulation range')
            # print("Additional normalization")
            # print("targeted_range", targeted_range)
            # print("Curve R90", find_R90(self.SOBP_xes, self.SOBP_yes))
        # else:
        #     self.logger.info('SOBP normalized to deepest peak')
        
        x_shallowest_peak, _ = _find_first_peak_from_left(self.SOBP_xes, self.SOBP_yes)
        
        
        self.SOBP = dose_engine.SOBP(self.SOBP_xes, self.SOBP_yes, base_modulation, x_shallowest_peak, x_deepest_peak)
        self.SOBP.calc()
        # SOBP_sum = sum(self.SOBP_yes)
        # for braggpeak in self.braggpeaks:
        #     f = interp1d(braggpeak.xes, braggpeak.yes)
            #braggpeak.relative_contribution = sum(f(self.SOBP_xes))/SOBP_sum

        # Adding padding for simpler interpolation
        self.SOBP_xes = np.append(self.SOBP_xes, 1000)
        self.SOBP_yes = np.append(self.SOBP_yes, 0)

        self.SOBP_f = interp1d(self.SOBP_xes, self.SOBP_yes, fill_value="extrapolate")
        

        self.depth_dose_configured = True

        self.logger.info('Depth profile configured')




def _find_norm_factor_from_the_right(xes, yes):
    idx_split = -1000
    
    prev = -1
    for idx in np.flip(range(len(xes))):
        
        val = yes[idx]
        if val < 0.1:
            prev = val
            continue
        
        if val >= prev:
            prev = val
            continue

        idx_split = idx

        break
    
    sub_xes = xes[idx_split:]
    sub_yes = yes[idx_split:]
    normfactor = np.max(sub_yes)
    sub_yes = sub_yes/max(sub_yes)
    

    x_deepest_peak = xes[idx_split + 1]

    return x_deepest_peak, normfactor




def _find_first_peak_from_left(xes, yes):
    idx_split = -1000
    
    prev = -1
    for idx in range(len(xes)):
        
        val = yes[idx]
        if val < 0.95:
            prev = val
            continue
        
        if val >= prev or abs(prev - val) < 0.0001:
            prev = val
            continue
        

        idx_split = idx - 2
        break
    
    sub_xes = xes[:idx_split]
    sub_yes = yes[:idx_split]
    normfactor = np.max(sub_yes)
    sub_yes = sub_yes/max(sub_yes)
    
    x_shallowest_peak = xes[idx_split + 1]
    
    return x_shallowest_peak, normfactor