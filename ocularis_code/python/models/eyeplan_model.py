# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""



import geometry
from pyproton.volume import Grid
from pyproton.structure.organ_data import OrganData, OrganType
from pyproton.structure.array_structure import ArrayStructure
from pyproton.handler import EyeplanStructureHandler
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import math
from scipy.interpolate import griddata
from scipy.spatial import KDTree
from scipy.optimize import curve_fit
import copy

from scipy.interpolate import interp2d

from utils.filters import densify_point_cloud, densify_point_cloud, filter_points_relative_to_plane, project_furthest_point_behind, segment_and_densify_point_cloud, make_spherical_object_hollow_thin, find_voxel_centers_vectorized
from utils.vector_utils import angle
from pyproton.utils.concave_hull import ConcaveHull
from utils.transformations import rotation_matrix_from_vectors
from utils.inflections import find_inflexions

import models
import dose_engine

from scipy.spatial.transform import Rotation
import pyvista as pv
import sklearn.cluster
from utils.transformations import mapping
from configuration import constants


currentdir = Path(os.getcwd())
newdir = currentdir.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))


class EyePlanModel(models.Model):

    def __init__(self, patient_config):
        super().__init__()
        self.patient_config = patient_config
        # self.config = patient_config["config"]
        self._config = None
        self._file_path = self.patient_config["patient_path"]

        self._structure_set = None
        self._structure_set_gaze_centered = None
        self._structure_set_clips_registered = None
        self._structure_set_resampled_registered = None
        self._doseplane_h = None
        self._doseplane_v = None
        # self._dose_h = None
        # self._dose_v = None
        self._dose_h_gaze_centered = None
        self._dose_v_gaze_centered = None
        self._dose_h_clips_registered = None
        self._dose_v_clips_registered = None
        self._polar_angle = None
        self._azimuth_angle = None
        self._target_range_from_ptd = None
        self._modulation_range_from_ptd = None
        self._gaze_vector_light = None
        self._intersection_points = None
        self._clips = None
        self._clips_gaze_centered = None
        self._grid = None
        self._centre_of_model = None
        self._relative_skin_pos = None
        self._proximal_margin = 0
        self._distal_margin = 0
        self._use_of_wedge = None
        self._treatment_eye = None
        self._fixation_eye = None
        self._pupil_to_centre_distance = None
        self._symmetrical_eye = None
        
        self._retina_defined = False
        self._ciliary_body_defined = False
        self.eyelid_case = False
        
        self._eyelid_thickness = None

        self._clip_transformation = None
        self._ptd_filename = self.patient_config["eyeplan ptd filename"]
        
        try:
            self.patient_number = int(self._ptd_filename[0:5])
        except:
            self.patient_number = "Unable to read"
        
        
        self.title = "Eyeplan model " +  self._ptd_filename
        
        self.n_clips = 4

        self._collimator_points_in_plane = []

    @property
    def config(self):
        return self._config
    
    @config.setter
    def config(self, conf):
        self._config = conf
        # self.add_pupil_structure()

    @property
    def grid(self):
        return self._grid
    
    @grid.setter
    def grid(self, grid):
        self._grid = grid
        
        try:
            self.structure_set_clips_registered.grid = self._grid 
        except:
            pass
        
        # self.structure_set_clips_registered.grid = self._grid

    # @property
    # def most_forward_point(self):
    #     "Doesnt work"
    #     raise ValueError

    #     max_z = 0
    #     for name in self.structure_set.name_list:
    #         if name == "baseArea":
    #             continue
    #         maxi = max(self.structure_set[name].contour_coordinates[0:, 2])
    #         if maxi > max_z:
    #             max_z = maxi
    #     return max_z

    @property
    def clip_transformation(self):
        try:
            self._clip_transformation = np.asarray(pd.read_csv(Path(
                self.patient_config["patient_path"]) / self.patient_config["clip transformation"], header=None))
        except:
            print("Clip transformation not loaded")

        return self._clip_transformation

    @property
    def label_vector(self):
        return np.append(self.centre_of_model, (self.polar_angle, self.azimuth_angle, self.light_point[0], self.light_point[1]))

    @property
    def feature_points(self):
        _ = self.feature_vector
        return self._feature_points

    @property
    def feature_distances(self):
        _ = self.feature_vector
        return self._feature_xes, self._feature_yes, self._feature_zes

    @property
    def feature_labels(self):

        # self._feature_labels = _get_feature_labels()

        self._feature_labels = [
            'target com x',
            'target com y',
            'target com z',
            'optdisk com x',
            'optdisk com y',
            'optdisk com z',
            'macula com x',
            'macula com y',
            'macula com z',
            'target max x x',
            'target max x y',
            'target max x z',
            'target min x x',
            'target min x y',
            'target min x z',
            'target max y x',
            'target max y y',
            'target max y z',
            'target min y x',
            'target min y y',
            'target min y z',
            'target max z x',
            'target max z y',
            'target max z z',
            'target min z x',
            'target min z y',
            'target min z z',

            #F10
            'apex closest x',
            'apex closest y',
            'apex closest z',

            #F11
            'apex furthest x',
            'apex furthest y',
            'apex furthest z',

            #f12 target_point_closest_to_optdisk
            'target.closest optdisk x',
            'target.closest optdisk y',
            'target.closest optdisk z',

            #f13 target_point_closest_to_macula
            'target.closest macula x',
            'target.closest macula y',
            'target.closest macula z',


            # Scalars
            'target eye',
            'fixation eye',
            'use of wedge',
            'target basearea',
            'target.com - Optdisk',
            'target.com - Macula',
            'target.closest - Optdisk',
            'target.closest - Macula',
            'apex height']

        return self._feature_labels

    @property
    def label_labels(self):
        self._label_labels = [

            'treatment pos x',
            'treatment pos y',
            'treatment pos z',
            'polar',
            'azimuth',
            'fixation p x',
            'fixation p y']

        return self._label_labels

    @property
    def feature_vector(self):

        # Points ----------------------------------------

        f1 = self.structure_set_gaze_centered["target"].com
        f2 = self.structure_set_gaze_centered["optdisk"].com
        f3 = self.structure_set_gaze_centered["macula"].com

        target_contour_coordinates = self.structure_set_gaze_centered["target"].contour_coordinates
        idx_max = np.where(target_contour_coordinates ==
                           max(target_contour_coordinates[0:, 0]))
        idx_min = np.where(target_contour_coordinates ==
                           min(target_contour_coordinates[0:, 0]))
        idy_max = np.where(target_contour_coordinates ==
                           max(target_contour_coordinates[0:, 1]))
        idy_min = np.where(target_contour_coordinates ==
                           min(target_contour_coordinates[0:, 1]))
        idz_max = np.where(target_contour_coordinates ==
                           max(target_contour_coordinates[0:, 2]))
        idz_min = np.where(target_contour_coordinates ==
                           min(target_contour_coordinates[0:, 2]))

        f4 = target_contour_coordinates[idx_max[0][0], 0:]
        f5 = target_contour_coordinates[idx_min[0][0], 0:]
        f6 = target_contour_coordinates[idy_max[0][0], 0:]
        f7 = target_contour_coordinates[idy_min[0][0], 0:]
        f8 = target_contour_coordinates[idz_max[0][0], 0:]
        f9 = target_contour_coordinates[idz_min[0][0], 0:]

        apex_point_closest = np.asarray([10000, 10000, 10000])
        closest_distance = 100000
        apex_point_furthest = np.asarray([0.1, 0.1, 0.1])
        distances = np.zeros(len(target_contour_coordinates))
        origo = np.asarray([0, 0, 0])
        for i, p in enumerate(target_contour_coordinates):
            dist = math.dist(p, origo)
            if dist < closest_distance:
                closest_distance = dist
                apex_point_closest = p
        

        vec = apex_point_closest - origo

        distant_point = origo + vec*100

        for i, p in enumerate(target_contour_coordinates):
            dist = math.dist(p, distant_point)
            distances[i] = dist

        smallest_distance = 1000
        smallest_distance_idx = 0
        for i, dist in enumerate(distances):
            if dist < smallest_distance:
                smallest_distance = dist
                smallest_distance_idx = i
                
                
        apex_point_furthest = target_contour_coordinates[smallest_distance_idx]

        f10 = apex_point_closest
        f11 = apex_point_furthest
        
        
        closest_distance_to_optdisk = 100000
        target_point_closest_to_optdisk = np.asarray([100, 100, 100])
        opt_disk_com = self.structure_set_gaze_centered['optdisk'].com
        for p in target_contour_coordinates:
            dist = math.dist(p, opt_disk_com)
            if dist < closest_distance_to_optdisk:
                closest_distance_to_optdisk = dist
                target_point_closest_to_optdisk = p
        
        
        f12 = target_point_closest_to_optdisk

        closest_distance_to_macula = 100000
        target_point_closest_to_macula = np.asarray([100, 100, 100])
        macula_com = self.structure_set_gaze_centered['macula'].com
        for p in target_contour_coordinates:
            dist = math.dist(p, macula_com)
            if dist < closest_distance_to_macula:
                closest_distance_to_macula = dist
                target_point_closest_to_macula = p

        f13 = target_point_closest_to_macula
        
        
        # Scalars -----------

        f14 = int(self.treatment_eye == "R")
        f15 = int(self.fixation_eye == "R")
        f16 = self.use_of_wedge
        if 'trgt_baseArea' in self.structure_set_gaze_centered.name_list:
            f17 = self.structure_set_gaze_centered['trgt_baseArea'].contour_coordinates[0][0]
        else:
            f17 = self.structure_set_gaze_centered['baseArea'].contour_coordinates[0][0]
        
        target_com = self.structure_set_gaze_centered["target"].com
        optdisk_com = self.structure_set_gaze_centered["optdisk"].com
        macula_com = self.structure_set_gaze_centered["macula"].com
        f18 = math.dist(target_com, optdisk_com)
        f19 = math.dist(target_com, macula_com)

        f20 = closest_distance_to_optdisk
        f21 = closest_distance_to_macula

        # Height dist
        f22 = math.dist(apex_point_closest, apex_point_furthest)

        # Assignment of feature vector -------------------

        # Points
        feature_vector = np.array(
            [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13])

        # Scalars
        feature_vector = np.append(
            feature_vector, (f14, f15, f16, f17, f18, f19, f20, f21, f22))

        self._feature_points = []
        self._feature_points = self._feature_points + \
            [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13]

        self._feature_xes = []
        self._feature_xes.append((target_com[0], optdisk_com[0]))
        self._feature_xes.append((target_com[0], macula_com[0]))
        self._feature_xes.append(
            (apex_point_closest[0], apex_point_furthest[0]))

        self._feature_yes = []
        self._feature_yes.append((target_com[1], optdisk_com[1]))
        self._feature_yes.append((target_com[1], macula_com[1]))
        self._feature_yes.append(
            (apex_point_closest[1], apex_point_furthest[1]))

        self._feature_zes = []
        self._feature_zes.append((target_com[2], optdisk_com[2]))
        self._feature_zes.append((target_com[2], macula_com[2]))
        self._feature_zes.append(
            (apex_point_closest[2], apex_point_furthest[2]))
        
        

        return feature_vector

    @property
    def feature_vector2(self):

        f1 = self.structure_set_gaze_centered["target"].com

        z_min = min(
            self.structure_set_gaze_centered["target"].contour_coordinates[0:, 2])
        z_max = max(
            self.structure_set_gaze_centered["target"].contour_coordinates[0:, 2])
        points = self.structure_set_gaze_centered["target"].contour_coordinates

        idx_zmin = np.argmin(points[0:, 2])
        idx_zmax = np.argmax(points[0:, 2])

        cloud = pv.PolyData(
            self.structure_set_gaze_centered["target"].contour_coordinates)
        # cloud.plot()

        volume = cloud.delaunay_3d(alpha=10.)
        shell = volume.extract_geometry()

        feature_points = []

        feature_points.append(points[idx_zmin])
        feature_points.append(points[idx_zmax])

        step = (z_max - z_min)/10

        for z in np.linspace(z_min + step, z_max - step, 8):
            resampled_points = _resample_points_in_plane(
                shell, z)

            idx_ymin = np.argmin(resampled_points[0:, 1])
            idx_ymax = np.argmax(resampled_points[0:, 1])

            idx_xmin = np.argmin(resampled_points[0:, 0])
            idx_xmax = np.argmax(resampled_points[0:, 0])

            feature_points.append(resampled_points[idx_xmin])
            feature_points.append(resampled_points[idx_xmax])

            feature_points.append(resampled_points[idx_ymin])
            feature_points.append(resampled_points[idx_ymax])

        feature_vector = np.asarray(feature_points).flatten()
        feature_vector = np.append(f1, feature_vector)
        return feature_vector

    @property
    def skin_plane_most_anterior(self):
        most_forward_point = max(self.structure_set_clips_registered['eyeglobe'].contour_coordinates[0:,2])
        skin_plane_z_pos = most_forward_point
        return skin_plane_z_pos
        
    @property
    def target_and_modulation_range_of_target(self):
        struct = self.structure_set_clips_registered['target']
        target_points = struct.contour_coordinates

        longest = 0
        deepest_ray = []
        shortest_dist = 100000000
        shallowest_ray = []        
        for p in target_points:
            ray1 = dose_engine.Ray([0,0,constants["VPS"]], p)
            
            if ray1.D12 > longest:
                longest = ray1.D12
                deepest_ray = ray1
            if ray1.D12 < shortest_dist:
                shortest_dist = ray1.D12
                shallowest_ray = ray1                

        
        deepest_ray_skin_plane_point = deepest_ray.find_intersect_with_plane( np.asarray([0,0,1]), np.asarray([0,0, self.skin_plane_most_anterior]))
        shallowest_ray_skin_plane_point = shallowest_ray.find_intersect_with_plane( np.asarray([0,0,1]), np.asarray([0,0, self.skin_plane_most_anterior]))
        

        deepest_ray.surface_interaction_point = deepest_ray.find_intersect_with_plane( np.asarray([0,0,1]), np.asarray([0,0, self.skin_plane_most_anterior]))
        shallowest_ray.surface_interaction_point = shallowest_ray.find_intersect_with_plane( np.asarray([0,0,1]), np.asarray([0,0, self.skin_plane_most_anterior]))
        
        target_radiological_depth = math.dist(deepest_ray.T, deepest_ray.surface_interaction_point)
        target_modulation= target_radiological_depth - math.dist(shallowest_ray.T, shallowest_ray.surface_interaction_point)        
        
        return target_radiological_depth, target_modulation, deepest_ray, shallowest_ray

    def __repr__(self):
         return self.title
    # @property
    # def target_depth_of_target(self):
    #     struct = self.structure_set_clips_registered['target']
    #     target_points = struct.contour_coordinates

    #     longest = 0
    #     deepest_ray = []
    #     for p in target_points:
    #         ray1 = dose_engine.Ray([0,0,constants["VPS"]], p)
            
    #         if ray1.D12 > longest:
    #             longest = ray1.D12
    #             deepest_ray = ray1


    #     deepest_ray.surface_interaction_point = deepest_ray.find_intersect( np.asarray([0,0,1]), np.asarray([0,0, self.skin_plane_most_anterior]))

        
    #     target_radiological_depth = math.dist(deepest_ray.T, deepest_ray.surface_interaction_point)
    #     return target_radiological_depth, deepest_ray

    # @property
    # def modulation_range_of_target(self):
    #     struct = self.structure_set_clips_registered['target']
    #     target_points = struct.contour_coordinates

    #     shortest_dist = 100000000
    #     shallowest_ray = []
    #     for p in target_points:
    #         ray1 = dose_engine.Ray([0,0,constants["VPS"]], p)
            
    #         if ray1.D12 < shortest_dist:
    #             shortest_dist = ray1.D12
    #             shallowest_ray = ray1


    #     shallowest_ray.surface_interaction_point = shallowest_ray.find_intersect( np.asarray([0,0,1]), np.asarray([0,0, self.skin_plane_most_anterior]))

        
    #     target_modulation= math.dist(shallowest_ray.T, shallowest_ray.surface_interaction_point)
    #     return target_modulation, shallowest_ray

    # @property
    # def target_rays(self):    
                
    #     geo = geometry.OptisGeometry(self.config)
        
    #     if len(self._target_rays) == 0:
        
    #         struct = self.structure_set_clips_registered['target']
    #         target_points = struct.contour_coordinates
            
    #         for p in target_points:
    #             ray1 = dose_engine.Ray(geo.central_axis.S, p)
    #             self._target_rays.append(ray1)
        
    #     return self._target_rays
    
    
    def generate_structure_set(self, new_points):
        
        """
        Generates a structure set with points from input variable
        dict new_points 
        """
        
        
        print("Generating new structure set")
        
        # It doesnt matter from where the structures are fetched. 
        # The points describing the structures are replaced anyway
        try:
            structure_handler = EyeplanStructureHandler(
                self.patient_config["patient_model_path"])
            structure_set = structure_handler.structure_set
        except:
            structure_handler = EyeplanStructureHandler(
                self.patient_config["patient_model_clips_path"])
            structure_set = structure_handler.structure_set
        
        found_any = False
        name_list = structure_set.name_list
        
        if self._retina_defined and "retina" in new_points.keys():
            struct = ArrayStructure(
                OrganData("retina", OrganType.UNKNOWN), new_points["retina"])
            struct.grid = self.grid
            
            structure_set.add_structure(struct)
            # self._structure_set_clips_registered.add_structure(struct)
            name_list.append("retina")
            
            
        if self._ciliary_body_defined and "ciliary_body" in new_points.keys():
            struct = ArrayStructure(
                OrganData("ciliary_body", OrganType.UNKNOWN), new_points["ciliary_body"])
            struct.grid = self.grid
            
            structure_set.add_structure(struct)
            
            # structure_set.add_structure(self._structure_set_clips_registered["ciliary_body"])
            name_list.append("ciliary_body") 
        
        # print("name_list", name_list)
        for name in name_list:
            
            if name in new_points.keys():
                found_any = True
            
                struct = structure_set[name]
                struct.contour_coordinates = new_points[name]
            else:
                struct = structure_set[name]
                struct.contour_coordinates = np.zeros([0, struct.contour_coordinates.shape[1]])
                structure_set.name_list.remove(name)
                
        if found_any == False:
            raise ValueError
        
        structure_set.grid = self.grid
        
        return structure_set

    @property
    def structure_set_gaze_centered(self):

        first_call = False

        if self._structure_set == None:
            first_call = True
            structure_handler = EyeplanStructureHandler(
                self.patient_config["patient_model_path"])
            self._structure_set_gaze_centered = structure_handler.structure_set
        
        return self._structure_set_gaze_centered

    @property
    def structure_set(self):

        first_call = False

        if self._structure_set == None:
            first_call = True
            structure_handler = EyeplanStructureHandler(
                self.patient_config["patient_model_path"])
            self._structure_set = structure_handler.structure_set

        if self._clips == None and 'clips' in self._structure_set:
            for idx in self.clips:
                struct = ArrayStructure(
                    OrganData("clip" + str(idx), OrganType.UNKNOWN), self.clips[idx].array)
                struct.grid = self.grid

                self.structure_set.add_structure(struct)

        if first_call:
            # self.structure_set.translate_contour_coordinates(
            #     -self.centre_of_model)
            for name in self.structure_set.name_list:
                struct = self._structure_set[name]
                # struct.contour_coordinates = self.rotation_matrix_z_to_gaze_light.dot(
                #     struct.contour_coordinates.T).T
                struct.contour_coordinates = _rotate_eyeplan_model_to_treatment_room_frame(
                    self, struct.contour_coordinates)

        return self._structure_set

    # @property
    # def structure_set(self):
    #     if self._structure_set_gaze_centered == None:
    #         structure_handler = EyeplanStructureHandler(
    #             self.patient_config["patient_model_path"])
    #         self._structure_set_gaze_centered = structure_handler.structure_set

    #     if self._clips_gaze_centered == None and 'clips' in self._structure_set_gaze_centered:
    #         for idx in self.clips:
    #             struct = ArrayStructure(
    #                 OrganData("clip" + str(idx), OrganType.UNKNOWN), self.clips[idx].array)
    #             struct.grid = self.grid

    #             self._structure_set_gaze_centered.add_structure(struct)

    #     return self._structure_set_gaze_centered

    @property
    def structure_set_clips_registered(self):
        """
        Reads structures registered using clips based method

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        if self._structure_set_clips_registered == None:
            structure_handler = EyeplanStructureHandler(
                self.patient_config["patient_model_clips_path"])
            self._structure_set_clips_registered = structure_handler.structure_set
            self._structure_set_clips_registered._loaded_filenames = structure_handler._loaded_filenames
            
            
        if self._clips_gaze_centered == None and 'clips' in self._structure_set_clips_registered:
            for idx in self.clips:
                struct = ArrayStructure(
                    OrganData("clip" + str(idx), OrganType.UNKNOWN), self.clips[idx].array)
                struct.grid = self.grid

                self._structure_set_clips_registered.add_structure(struct)
        

        return self._structure_set_clips_registered


    @property
    def pupil_to_centre_distance(self):
        if self._pupil_to_centre_distance == None:
            self._pupil_to_centre_distance = math.dist(self.structure_set_clips_registered["pupil"].contour_coordinates[0], self.centre_of_model)
        return self._pupil_to_centre_distance
    
    def dvh_filename(self, target_name):
        filenames =  list(filter(lambda x:x.endswith((".csv")), os.listdir(Path(self.patient_config["dvhs"]))) )
        return _filename_from_filenames(target_name, filenames)
        
    def add_retina_structure(self):
        
        
        
        # print("3")
        # if 'retina' not in self._structure_set_clips_registered.name_list and 'eyeglobe' in self._structure_set_clips_registered.name_list:

        # print("struct.contour_coordinates.shape", struct.contour_coordinates.shape)
        # print("struct.structure_points.shape", struct.structure_points.shape)
        
        # points_for_filtering = copy.deepcopy(eyeglobe_points)
        # points_for_filtering[0:,2] = points_for_filtering[0:,2] + 0.4
        
        
        structure_name = 'lens'
        struct = self.structure_set_clips_registered[structure_name]
        # lens_com = struct.com
        lens_points = struct.structure_points
        
        lens_min, lens_max = _find_min_max_projection(lens_points, self.centre_of_model, self.gaze_vector_light, min_lens = 0.01, max_lens = 0.6)
        
        retina_cutoff = lens_min - self.gaze_vector_light * 3
        

        
        # point_threshold = lens_min #- 2* self.gaze_vector_light
        struct = self.structure_set_clips_registered['eyeglobe']
        struct.resample_contour_points_in_grid(1)
        eyeglobe_points = struct.structure_points
        
        retina_mask = make_spherical_object_hollow_thin(struct.binary_mask)
        retina_points = find_voxel_centers_vectorized(retina_mask, self.grid.origin, self.grid.spacing)
        retina_points = _points_behind_axis_point(retina_points, retina_cutoff, self.gaze_vector_light)
        
        
        # retina_points = _points_behind_axis_point(eyeglobe_points, retina_cutoff, self.gaze_vector_light)
        

        
        # p_end = project_furthest_point_behind(eyeglobe_points, self.centre_of_model, self.gaze_vector_light)
        
        # retina_cuttoff_dist_anterior_distal = 0.6*math.dist(p_end, self.centre_of_model)
        
        
      
        
        # retina_points_anterior, retina_points_distal = filter_points_relative_to_plane(retina_points, self.centre_of_model, self.gaze_vector_light, retina_cuttoff_dist_anterior_distal)
        
        # assert len(retina_points_anterior) > 2
        # assert len(retina_points_distal) > 2
        
        # retina_points_distal = densify_point_cloud(retina_points_distal, 0.6)
        # # retina_points_anterior = retina_points[retina_points[:, 2] > z_threshold]
        
        
        
        
        # resolution = 0.1
        # # print("eyeglobe_points.shape", eyeglobe_points.shape)
        # # print("retina_points_anterior.shape", retina_points_anterior.shape)
        # # print("retina_points_distal.shape", retina_points_distal.shape)
        
        
        # retina_points_anterior = segment_and_densify_point_cloud(retina_points_anterior, 24, resolution)
        
        # retina_points = np.vstack([
        #     retina_points_distal,
        #     retina_points_anterior
        # ])



        
        # retina_points = _points_behind_axis_point(eyeglobe_points, point_threshold, self.gaze_vector_light)
        
        # retina_points = densify_point_cloud(retina_points, 0.4)
        # retina_points = densify_point_cloud(retina_points, 0.6)
        # retina_points = densify_point_cloud(retina_points, 0.4)
        
        
        # retina_points = filter_sphere(retina_points, points_for_filtering)
        
        # # retina_points = points_near_each_other(retina_points, 0.5)
        # retina_points = filter_for_neighbors(retina_points, 0.6)
        # retina_points = filter_for_neighbors(retina_points, 0.6)
        
        struct = ArrayStructure(
            OrganData("retina", OrganType.UNKNOWN), retina_points)
        struct.grid = self.grid
        
        self._structure_set_clips_registered.add_structure(struct)   
        
        
        self._retina_defined = True
        print("retina defined")
        
        
    def add_ciliary_structure(self):
        
        struct = self.structure_set_clips_registered['eyeglobe']
        struct.resample_contour_points_in_grid(1)
        # eyeglobe_mask = struct.binary_mask
        
        # eyeglobe_points = struct.structure_points
        # points_for_filtering = copy.deepcopy(eyeglobe_points)
        
        
        
        gaze_vector_light = self.gaze_vector_light
                
        structure_name = 'lens'
        struct = self.structure_set_clips_registered[structure_name]
        # lens_com = struct.com
        lens_points = struct.structure_points
        
        struct = self.structure_set_clips_registered['eyeglobe']
        eyeglobe_points = struct.structure_points
        
        
        # p_min, p_max = _find_min_max_projection(lens_points, self.centre_of_model, gaze_vector_light)
        
        
        
        
        p1,p2 = _find_min_max_projection(lens_points, self.centre_of_model, gaze_vector_light)
        
        dist = math.dist(p1, p2)
        
        p_mid = p1 + gaze_vector_light*dist
        
        lens_axis_points = []
        for d in np.linspace(0,dist, 20):
            new_p = p1 + gaze_vector_light*d
            lens_axis_points.append(new_p)
        
        all_lens_max = []
        for center_p in lens_axis_points:
            radius_min = _calculate_radius_around_point(lens_points, center_p, gaze_vector_light, 1)
            all_lens_max.append(radius_min)
        
        lens_radius = np.max(all_lens_max)
        
        ciliary_body_points = []
        
        # max_rolling = []
        min_rolling = []
        for center_p in lens_axis_points:
            radius_max = _calculate_radius_around_point(eyeglobe_points, center_p, gaze_vector_light, 1)
            radius_min = _calculate_radius_around_point(lens_points, center_p, gaze_vector_light, 1)
            
            
            if len(min_rolling) > 2:
                min_rolling.pop(0)
            
            # if len(max_rolling) > 2:
            #     max_rolling.pop(0)
            
            # max_rolling.append(radius_max)
            
            effective_min_radius = lens_radius + 2 + math.dist(p_mid, center_p) #+ math.dist(p1, center_p) 
            
            if effective_min_radius < lens_radius:
                effective_min_radius = lens_radius
            
            min_rolling.append(effective_min_radius)
            
            for radi in np.linspace(np.mean(min_rolling), radius_max, 20): #np.mean(max_rolling)
                temp_points = _generate_ring_points(center_p, gaze_vector_light, 150, radi)
                #temp_points = _filter_points_by_radius(temp_points, self.centre_of_model, gaze_vector_light, radius_min)
                ciliary_body_points.extend(temp_points) #interpolate_points(points, eyeglobe_points, gaze_vector_light, 100)
        
        
        ciliary_body_points = np.asarray(ciliary_body_points)
        
      
        
        
        struct = ArrayStructure(
            OrganData("ciliary_body", OrganType.UNKNOWN), ciliary_body_points)
        struct.grid = self.grid
        
        self._structure_set_clips_registered.add_structure(struct)   
        
        self._ciliary_body_defined = True
        print("ciliary body defined")        
        
    #     print(struct)
        
    # #     self._structure_set_clips_registered.add_structure(struct)
        
    #     # Code for adding pupil to structure set
    #     p_front = self.centre_of_model + 50*self.gaze_vector_light

    #     struct = self.structure_set_clips_registered["lens"]

    #     points = struct.structure_points
        
    #     smallest_dist = float("inf")

    #     best_points = points[0]
    #     for p in points:
    #         dist = math.dist(p, p_front)
    #         if dist < smallest_dist:
    #             smallest_dist = dist
    #             best_points = p

    #     d = self.gaze_vector_light
    #     t = np.dot(best_points - p_front, d)
    #     p_on_axis = p_front + t*self.gaze_vector_light        
        

    #     new_struct = ArrayStructure(
    #         OrganData("pupil", OrganType.UNKNOWN), p_on_axis.reshape(1, 3))
    #     new_struct.grid = self.grid
    #     self._structure_set_clips_registered.add_structure(new_struct)    
    
        
        
    # def add_pupil_structure(self):
    #     # Code for adding pupil to structure set
    #     p_front = self.centre_of_model + 50*self.gaze_vector_light

    #     struct = self.structure_set_clips_registered["lens"]

    #     points = struct.structure_points
        
    #     smallest_dist = float("inf")

    #     best_points = points[0]
    #     for p in points:
    #         dist = math.dist(p, p_front)
    #         if dist < smallest_dist:
    #             smallest_dist = dist
    #             best_points = p

    #     d = self.gaze_vector_light
    #     t = np.dot(best_points - p_front, d)
    #     p_on_axis = p_front + t*self.gaze_vector_light        
        

    #     new_struct = ArrayStructure(
    #         OrganData("pupil", OrganType.UNKNOWN), p_on_axis.reshape(1, 3))
    #     new_struct.grid = self.grid
    #     self._structure_set_clips_registered.add_structure(new_struct)


    # @property
    # def structure_set_resampled_registered(self):
    #     """
    #     Reads structures registered using clips and resampling

    #     Returns
    #     -------
    #     TYPE
    #         DESCRIPTION.

    #     """

    #     if self._structure_set_resampled_registered == None:
    #         structure_handler = EyeplanStructureHandler(
    #             self.patient_config["patient_model_resampled_path"])
    #         self._structure_set_resampled_registered = structure_handler.structure_set

    #     if self._clips_gaze_centered == None and 'clips' in self._structure_set_resampled_registered:
    #         for idx in self.clips:
    #             struct = ArrayStructure(
    #                 OrganData("clip" + str(idx), OrganType.UNKNOWN), self.clips[idx].array)
    #             struct.grid = self.grid

    #             self._structure_set_resampled_registered.add_structure(struct)

    #     return self._structure_set_resampled_registered

    @property
    def clips(self):

        if self._clips == None:
            self._clips = {}
            for i in range(self.n_clips):
                self._clips[i] = Clip(i)

            self.classify_clips()
        return self._clips

    # @property
    # def com_dose_planes(self):

    #     print("Warning, COM poorly defined")

    #     horizontal_points = self.dose_horizontal
    #     vertical_points = self.dose_vertical

    #     horizontal_points = np.delete(horizontal_points, 3, 1)
    #     vertical_points = np.delete(vertical_points, 3, 1)

    #     com_horizontal = np.asarray([sum(row[i] for row in horizontal_points)/len(
    #         horizontal_points) for i in range(len(horizontal_points[0]))])

    #     com_vertical = np.asarray([sum(row[i] for row in vertical_points)/len(
    #         vertical_points) for i in range(len(vertical_points[0]))])

    #     com = com_horizontal - (com_horizontal - com_vertical)/2

    #     return com

    def read(self):
        # try:

        # Using the transformation scheme
        self._dose_h = np.array(pd.read_csv(
            Path(self._file_path) / self.patient_config["dose horizontal filename"], skiprows=2, header=None))
        self._dose_h[0:, 0:3] = _transform_points_to_treatment_room_frame(
            self, self._dose_h[0:, 0:3])

        # self._dose_v = np.array(pd.read_csv(
        #     Path(self._file_path) / self.patient_config["dose vertical filename"], skiprows=2, header=None))
        # self._dose_v[0:, 0:3] = _transform_points_to_treatment_room_frame(
        #     self, self._dose_v[0:, 0:3])

        # Clip registered
        self._dose_h_clips_registered = np.array(pd.read_csv(
            Path(self._file_path) / self.patient_config["dose horizontal filename"], skiprows=2, header=None))
        self._dose_h_clips_registered[0:, 0:3] = mapping(
            self._dose_h_clips_registered[0:, 0:3], self.clip_transformation)

        # self._dose_v_clips_registered = np.array(pd.read_csv(
        #     Path(self._file_path) / self.patient_config["dose vertical filename"], skiprows=2, header=None))
        # self._dose_v_clips_registered[0:, 0:3] = mapping(
        #     self._dose_v_clips_registered[0:, 0:3], self.clip_transformation)

        # Raw, ie in the Gaze centred reference frame
        self._dose_h_gaze_centered = np.array(pd.read_csv(
            Path(self._file_path) / self.patient_config["dose horizontal filename"], skiprows=2, header=None))
        # self._dose_v_gaze_centered = np.array(pd.read_csv(
        #     Path(self._file_path) / self.patient_config["dose vertical filename"], skiprows=2, header=None))
        # # except:
        #     print("Unable to load dose planes from file")

    @property
    def dose_h_raw(self):
        return np.array(pd.read_csv(
            Path(self._file_path) / self.patient_config["dose horizontal filename"], skiprows=2, header=None))

    @property
    def dose_v_raw(self):
        return np.array(pd.read_csv(
            Path(self._file_path) / self.patient_config["dose vertical filename"], skiprows=2, header=None))
        #
        # try:
        #     self._dose_h = np.array(pd.read_csv(
        #         Path(self._file_path) / self.patient_config["dose horizontal filename"], skiprows=2, header=None))
        #     self._dose_v = np.array(pd.read_csv(
        #         Path(self._file_path) / self.patient_config["dose vertical filename"], skiprows=2, header=None))
        # except:
        #     print("Unable to load dose planes from file")

    @property
    def rotation_matrix_z_to_gaze_light(self):
        mat = rotation_matrix_from_vectors(
            np.asarray([0, 0, 1]), self.gaze_vector_light)
        return mat

    @property
    def rotation_matrix_gaze_light_to_z(self):
        mat = rotation_matrix_from_vectors(
            self.gaze_vector_light, np.asarray([0, 0, 1]))
        return mat

    @property
    def centre_of_model(self):
        self.parse_ptd()
        return self._centre_of_model

    @property
    def relative_skin_pos(self):
        self.parse_ptd()
        return self._relative_skin_pos

    @property
    def doseplane_h(self):
        if self._doseplane_h == None:
            self._doseplane_h = DosePlane(self.patient_config["patient_path"], self.patient_config["dose horizontal filename"], self.clip_transformation, True)
        return self._doseplane_h

    @property
    def doseplane_v(self):
        if self._doseplane_v == None:
            self._doseplane_v = DosePlane(self.patient_config["patient_path"], self.patient_config["dose vertical filename"], self.clip_transformation, False)
        return self._doseplane_v

    # @property
    # def dose_horizontal(self):
    #     self.read()
    #     return self._dose_h

    # @property
    # def dose_vertical(self):
    #     self.read()
    #     return self._dose_v

    # @property
    # def dose_horizontal_gaze_centered(self):
    #     self.read()
    #     return self._dose_h_gaze_centered

    # @property
    # def dose_vertical_gaze_centered(self):
    #     self.read()
    #     return self._dose_v_gaze_centered

    # @property
    # def dose_horizontal_clips_registered(self):
    #     self.read()
    #     return self._dose_h_clips_registered

    # @property
    # def dose_vertical_clips_registered(self):
    #     self.read()
    #     return self._dose_v_clips_registered

    # @property
    # def intersection_points(self):
    #     _ = self.gaze_vector_planes
    #     return self._intersection_points

    # @property
    # def gaze_vector_planes(self):

    #     import copy

    #     translation_vector = self.centre_of_model

    #     horizontal_points = copy.deepcopy(self.dose_horizontal[::10])
    #     vertical_points = copy.deepcopy(self.dose_vertical[::10])

    #     horizontal_points[0:, 0] += translation_vector[0]
    #     horizontal_points[0:, 1] += translation_vector[1]
    #     horizontal_points[0:, 2] += translation_vector[2]

    #     vertical_points[0:, 0] += translation_vector[0]
    #     vertical_points[0:, 1] += translation_vector[1]
    #     vertical_points[0:, 2] += translation_vector[2]

    #     horizontal_points = np.delete(horizontal_points, 3, 1)
    #     vertical_points = np.delete(vertical_points, 3, 1)

    #     horizontal_points = self.rotation_matrix_z_to_gaze_light.dot(
    #         horizontal_points.T).T
    #     vertical_points = self.rotation_matrix_z_to_gaze_light.dot(
    #         vertical_points.T).T

    #     intersection_points = _intersection_points_from_point_planes(
    #         horizontal_points, vertical_points)

    #     self._intersection_points = intersection_points

    #     p_min = intersection_points[0]
    #     p_max = intersection_points[-1]
    #     #p_mid = intersection_points[int(len(intersection_points)/2)]

    #     gaze = p_max - p_min
    #     gaze = gaze / np.linalg.norm(gaze)

    #     return gaze

    @property
    def light_point(self):
        _ = self.gaze_vector_light
        return self._light_point

    @property
    def gaze_vector_light(self):
        
        geo = geometry.OptisGeometry(self.config)

        p_axis_light = geo.p_axis_light

        # p_polar = geo.p_polar

        # light_point = geo.light_point

        # p_axis_light = np.asarray([0, 0, constants["FID"]])

        p_polar = p_axis_light + np.asarray([1, 0, 0])*math.tan(
            math.radians(self.polar_angle))*constants["FID"]

        light_vector = p_polar - p_axis_light

        r = Rotation.from_euler(
            'xyz', (0, 0, self.azimuth_angle), degrees=True)

        light_vector = r.apply(light_vector)

        light_point = p_axis_light + light_vector

        assert (light_point == _find_light_fixation_point(
            self.config, self.polar_angle, self.azimuth_angle)).all()

        r = np.tan(np.radians(self.polar_angle))*132.5

        y = r*np.sin(np.radians(self.azimuth_angle))
        x = r*np.cos(np.radians(self.azimuth_angle))

        light_point2 = np.asarray([x, y, constants["FID"]])
        
        assert math.dist(light_point2 , light_point) < 0.0001

        self._light_point = light_point

        gaze_vector_light = light_point - self.centre_of_model

        gaze_vector_light = gaze_vector_light / \
            np.linalg.norm(gaze_vector_light)

        self._gaze_vector_light = gaze_vector_light
        return self._gaze_vector_light

    @property
    def treatment_eye(self):
        if self._treatment_eye == None:
            self.parse_ptd()
        return self._treatment_eye

    @property
    def fixation_eye(self):
        if self._fixation_eye == None:
            self.parse_ptd()
        return self._fixation_eye

    @property
    def symmetrical_eye(self):
        if self._symmetrical_eye == None:
            self.parse_ptd()

            if self._symmetrical_eye == None:
                raise ValueError

        return self._symmetrical_eye

    @property
    def use_of_wedge(self):
        if self._use_of_wedge == None:
            self.parse_ptd()
        return self._use_of_wedge

    @property
    def polar_angle(self):
        if self._polar_angle == None:
            self.parse_ptd()
        return self._polar_angle

    @property
    def azimuth_angle(self):
        if self._azimuth_angle == None:
            self.parse_ptd()
        return self._azimuth_angle

    @property
    def target_range_from_ptd(self):
        if self._target_range_from_ptd == None:
            self.parse_ptd()
        return self._target_range_from_ptd

    @property
    def modulation_range_from_ptd(self):
        if self._modulation_range_from_ptd == None:
            self.parse_ptd()
        return self._modulation_range_from_ptd

    @property
    def proximal_margin(self):
        if self._proximal_margin == None:
            self.parse_ptd()
        return self._proximal_margin

    @proximal_margin.setter
    def proximal_margin(self, value):
        self._proximal_margin = value

    @property
    def distal_margin(self):
        if self._distal_margin == None:
            self.parse_ptd()        
        return self._distal_margin

    @distal_margin.setter
    def distal_margin(self, value):
        self._distal_margin = value
        
    @property
    def eyelid_thickness(self):
        if self._eyelid_thickness == None:
            self.parse_ptd()  
        return self._eyelid_thickness

    @distal_margin.setter
    def eyelid_thickness(self, value):
        self._eyelid_thickness = value        

    def infer_polar_and_azimuth(self, point):
        """
        Point on eye fixation plane
        """

        if len(point) == 2:
            point = np.asarray([point[0], point[1], 0])

        point = np.asarray(point)

        return _infer_polar_and_azimuth(self.patient_config["config"], point)

   

    def translate_dose_planes(self, translation_vector):

        self.dose_horizontal[0:, 0] += translation_vector[0]
        self.dose_horizontal[0:, 1] += translation_vector[1]
        self.dose_horizontal[0:, 2] += translation_vector[2]

        self.dose_vertical[0:, 0] += translation_vector[0]
        self.dose_vertical[0:, 1] += translation_vector[1]
        self.dose_vertical[0:, 2] += translation_vector[2]

    def rotate_dose_planes_around_axis(self, rot_vector):

        rot_vector = np.asarray(rot_vector)

        r = Rotation.from_euler(
            'xyz', (rot_vector[0], rot_vector[1], rot_vector[2]), degrees=True)

        for idx in range(len(self.dose_horizontal)):
            points = self.dose_horizontal[idx][0:3]
            self.dose_horizontal[idx][0:3] = r.apply(points)

        for idx in range(len(self.dose_vertical)):
            points = self.dose_vertical[idx][0:3]
            self.dose_vertical[idx][0:3] = r.apply(points)

        print("EP model rotated by ", rot_vector)

    def parse_ptd(self):

        self._use_of_wedge = 0

        with open(Path(self._file_path) / self.patient_config["eyeplan ptd filename"]) as f:
            lines = f.readlines()
            
            DIAM, DIAH, DIAV = -1, -1, -1

            for line in lines:
                split_line = line.split()
                if len(split_line) > 1:
                    if split_line[0] == 'POLR':
                        self._polar_angle = float(split_line[1])
                    if split_line[0] == 'AZIM':
                        self._azimuth_angle = float(split_line[1])
                    if split_line[0] == 'RNGE':
                        self._target_range_from_ptd = float(split_line[1])
                    if split_line[0] == 'MODN':
                        self._modulation_range_from_ptd = float(split_line[1])
                    if split_line[0] == 'SKIN':
                        self._relative_skin_pos = float(split_line[1])
                    if split_line[0] == 'ECEN':
                        self._centre_of_model = np.asarray(
                            [float(split_line[1][:-1]), float(split_line[2][:-1]), float(split_line[3])])
                    if split_line[0] == 'TEYE':
                        self._treatment_eye = split_line[1]
                    if split_line[0] == 'FIXN':
                        self._fixation_eye = split_line[1]
                    if split_line[0] == 'WED1':
                        self._use_of_wedge = 1
                    if split_line[0] == 'DTMG':
                        self._distal_margin = float(split_line[1])                   
                    if split_line[0] == 'PXMG':
                        self._proximal_margin = float(split_line[1])      
                    if split_line[0] == 'ULTH':
                        self._eyelid_thickness = float(split_line[1])                           
                    if split_line[0] == 'DIAM':
                        DIAM = float(split_line[1])           
                    if split_line[0] == 'DIAH':
                        DIAH = float(split_line[1])      
                    if split_line[0] == 'DIAV':
                        DIAV = float(split_line[1])      
                        
        if DIAM > 0 and DIAH > 0 and DIAV > 0:
            
            if DIAM == DIAH == DIAV:
                self._symmetrical_eye = True           
            else:
                self._symmetrical_eye = False    
                    
    @property
    def collimator_points_in_plane(self):

        collimator_points_in_plane = []

        key = "eyeplan mil filename"
        if key not in self.patient_config:
            print("mil file not opened")
            return

        try:
            with open(Path(self._file_path) / self.patient_config[key]) as f:
                lines = f.readlines()
                lines = lines[2:-2]

                for line in lines:
                    split_line = line.split()

                    x = float(split_line[2][1:])
                    y = float(split_line[3][1:])
                    collimator_points_in_plane.append([x, y])
        except:
            print("Unable to open mil file")
            return

        # if collimator_points_in_plane[0] == collimator_points_in_plane[-1]:
        #     collimator_points_in_plane = collimator_points_in_plane[0:-1]
        
        
        # collimator_points_in_plane = _unscramble_contour_points(
        #     collimator_points_in_plane)
        
        self._collimator_points_in_plane = np.asarray(collimator_points_in_plane)
        
        return self._collimator_points_in_plane

    def classify_clips(self, n_clips=4):
        """"
        https://stackoverflow.com/questions/66099958/density-clustering-around-a-separate-point-python
        """

        points = self.structure_set["clips"].contour_coordinates

        cluster = sklearn.cluster.AgglomerativeClustering(
            n_clusters=n_clips, affinity='euclidean', linkage='ward')
        cluster.fit_predict(points)

        for point, label in zip(points, cluster.labels_):
            self.clips[label].add_row(point)


class Clip():
    def __init__(self, label):
        self._array = np.empty(shape=[0, 3])
        self.label = label

    def add_row(self, row):
        self._array = np.append(self._array, [np.asarray(row)], 0)

    @property
    def com(self):
        try:
            return np.asarray([sum(row[i] for row in self._array)/len(
                self._array) for i in range(len(self._array[0]))])
        except:
            print("No points")
            return 0

    @property
    def array(self):
        return self._array


# a =  np.empty(shape=[0, 3])


# a = np.append(a, [np.asarray([1,1,1])],0)
# a = np.append(a, [np.asarray([1,1,1])],0)

# print(a)


# clip = Clip(1)

# clip.add_row([1,1,1])
# clip.add_row([1,2,3])
# clip.add_row([5,6,7])

# print(clip.array)


class DosePlane:
    def __init__(self, path, filename, transformation_matrix, isHorizontal = True):
        self.path = path
        self.filename = filename
        self.transformation_matrix = transformation_matrix
        self._isHorizontal = isHorizontal
        
        self.image = None
        self.xs_mean = None
        self.zs_mean = None
        
        self.parse()
        
    def parse(self):
        
        if self.filename == "":
            print("Warning: No dose from eyeplan")
            return
        
        # if self.transformation_matrix == None:
        #     raise ValueError
        
        try:
            file_content =  np.asarray(pd.read_csv(Path(self.path ) / self.filename, header=None, skiprows = 2))
        except:
            file_content =  np.asarray(pd.read_csv(Path(self.path) / self.filename, header=None, skiprows=2, encoding='ISO-8859-1'))
        
        

        
        self.locations_gaze_centered = file_content[0:,0:3]
        
        self.locations = mapping(file_content[0:,0:3], self.transformation_matrix ) 
        
        self.values = file_content[0:,3]
        
        N = int(math.sqrt(len(file_content)))
        
        # if int(math.sqrt(len(file_content))) == math.sqrt(len(file_content)):
        #     N = int(math.sqrt(len(file_content)))
        # else:
        #     return ValueError

        xs_mean = []
        ys_mean = []
        zs_mean = []

        for i in range(N):
            points = self.locations[i*N:59+i*N,2]
            zs_mean.append(np.mean(points))
        zs_mean = np.flip(zs_mean)
        self.zs_mean = zs_mean



        for i in range(N):
            x_pos = []
            for j in range(N):
                x_pos.append(self.locations[j*60 + i,0])
            xs_mean.append(np.mean(x_pos))
        xs_mean = np.flip(xs_mean)
        self.xs_mean = xs_mean


        for i in range(N):
            y_pos = []
            for j in range(N):
                y_pos.append(self.locations[j*60 + i,1])
            ys_mean.append(np.mean(y_pos))
        ys_mean = np.flip(ys_mean)
        self.ys_mean = ys_mean

        
      
        
        dist_between_zs = []
        for i in range(len(zs_mean) - 1):
            dist = abs(zs_mean[i] - zs_mean[i + 1])
            dist_between_zs.append(dist)
            
        res_z = np.mean(dist_between_zs)
        
        resolution = np.asarray([res_z, res_z, res_z])
        self.resolution = resolution
        
        
        # origin = np.asarray([ min(self.locations[0:,0]) - resolution[0]/2,
        #                      min(self.locations[0:,1]) - resolution[1]/2,
        #                      min(self.locations[0:,2]) - resolution[2]/2])
        
        
        
        # origin_plane = np.asarray([ min(self.locations[0:,0]),
        #                      min(self.locations[0:,1]),
        #                      min(self.locations[0:,2])])
        
        # apex_plane = np.asarray([ max(self.locations[0:,0]),
        #                      max(self.locations[0:,1]),
        #                      max(self.locations[0:,2])])
        
        
        origin = np.asarray([ min(xs_mean) - resolution[0]/2,
                             min(ys_mean) - resolution[1]/2,
                             min(zs_mean) - resolution[2]/2])


        
        origin_plane = np.asarray([ min(xs_mean),
                             min(ys_mean),
                             min(zs_mean)])

        apex_plane = np.asarray([ max(xs_mean),
                             max(ys_mean),
                             max(zs_mean)])

        
        
        image_h = np.zeros([60,60])
        image_v = np.zeros([60,60])
        
        
        apex = np.zeros(3)
        
        apex[0] = origin[0] + resolution[0]*N
        apex[1] = origin[1] + resolution[1]*N
        apex[2] = max(zs_mean) + resolution[2]/2  # origin[2] + resolution[2]*60. #
        origin[2] = apex[2] - resolution[2]*N
        
        
        minZ = origin[2]
        maxZ = apex[2]
        widthZ = resolution[2]
        
        minX = origin[0]
        maxX = apex[0]
        widthX = resolution[0]
        
        self.mesh_dimensions = [abs(minX) + maxX, abs(origin[1]) + apex[1], abs(minZ) + maxZ]
        
        # Z and X are used with plt.pcolor for plotting meshes in space
        Z, X = np.meshgrid(np.arange(minZ, maxZ + widthZ,
                                      widthZ), np.arange(minX, maxX + widthX, widthX))
        # Z, X = np.meshgrid(np.arange(minZ + widthZ/2, maxZ + widthZ*1.5,
        #                               widthZ), np.arange(minX, maxX + widthX, widthX))
        
        
        if Z.shape == (N+2, N+2):
            Z = Z[0:-1,0:-1]
        if X.shape == (N+2, N+2):
            X = X[0:-1,0:-1]     
            
        self.X = X
        self.Z = Z
            
        
        # z_voxel_pos = np.linspace(
        #     minZ + widthZ, maxZ, num=image_h.shape[0])
        
        z_voxel_pos = np.linspace(
            origin_plane[2] , apex_plane[2], num=image_h.shape[0])
        
        y_voxel_pos = np.linspace(
            origin_plane[1] , apex_plane[1], num=image_h.shape[1])        
        
        x_voxel_pos = np.linspace(
            origin_plane[0] , apex_plane[0], num=image_h.shape[1])# -
        
        self.z_voxel_pos = z_voxel_pos
        self.y_voxel_pos = y_voxel_pos
        self.x_voxel_pos = x_voxel_pos
        
        x_bin = 0
        z_bin = 59
        for x,y,z, val in zip(self.locations[0:,0], self.locations[0:,1], self.locations[0:,2], self.values):
            image_h[x_bin,z_bin] = val
            
            x_bin += 1
            
            if x_bin >= 60:
                x_bin -= 60
                z_bin -= 1
            
        
        image_h = np.flip(image_h, 0)
        self.image = image_h
        
        x = self.zs_mean
        y = self.xs_mean
        X2, Y2 = np.meshgrid(x, y)
        self.interpolation_2D_f = interp2d(x, y, self.image, kind='cubic')
        
        
        
        
    

    def interpolate_from(self, zes, xes):
        image_intp = self.interpolation_2D_f( zes, xes )
        image_intp[image_intp < 0] = 0
        image_intp[image_intp < 1e-5] = 0
        return image_intp


    @property
    def inflections(self):
        return find_inflexions(self.image, self.xs_mean, self.zs_mean)






def _filename_from_filenames(target_name, filenames):
    for item in filenames:
        if target_name in item:
            return item
    
    return "Not found"

def _unscramble_contour_points(points):

    length_of_ruler = 80

    hull = ConcaveHull()
    hull.loadpoints(points)
    hull.calculatehull(length_of_ruler)

    xes, yes = hull.boundary.boundary.xy

    unscrabled_points = []

    for x, y in zip(xes, yes):
        unscrabled_points.append((x, y))

    return unscrabled_points


def _plane_normal_from_points(points):
    p0 = points[0]
    p1 = points[int(len(points)/2)]
    p2 = points[-1]

    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    ux, uy, uz = u = [x1-x0, y1-y0, z1-z0]
    vx, vy, vz = v = [x2-x0, y2-y0, z2-z0]

    u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]

    normal = np.array(u_cross_v)

    normal = normal / np.linalg.norm(normal)

    return normal


def _intersection_points_from_point_planes(points1, points2):

    normal_points1 = _plane_normal_from_points(points1)
    normal_points2 = _plane_normal_from_points(points2)

    plane_size = 150

    p1 = pv.Plane(center=points1[0], direction=normal_points1,
                  i_size=plane_size, j_size=plane_size).triangulate()
    p2 = pv.Plane(center=points2[0], direction=normal_points2,
                  i_size=plane_size, j_size=plane_size).triangulate()

    intersection = p1.intersection(p2)

    intersection_points = intersection[0].points

    intersection_points = np.array(
        sorted(intersection_points, key=lambda x: x[2]))

    # pl = pv.Plotter()

    # pl.add_mesh(p1, style='wireframe')
    # pl.add_mesh(p2, style='wireframe')
    # pl.add_mesh(intersection[0], color='r', line_width=10)
    # pl.show()

    return intersection_points



def _find_light_fixation_point(config, polar_angle, azimuth_angle):
    geo = geometry.OptisGeometry(config)

    p_axis_light = geo.p_axis_light

    p_polar = p_axis_light + np.asarray([1, 0, 0])*math.tan(
        math.radians(polar_angle))*config["FID"]

    light_vector = p_polar - p_axis_light

    r = Rotation.from_euler(
        'xyz', (0, 0, azimuth_angle), degrees=True)

    light_vector = r.apply(light_vector)

    light_point = p_axis_light + light_vector

    return light_point


def _transform_points_to_treatment_room_frame(my_ep_model, my_xyz):
    t_c = my_ep_model.centre_of_model

    T_t = np.identity(4)
    T_t[0:3, 3] = t_c

    T_ = np.identity(4)
    T_[0:3, 3] = -t_c

    T_r = np.identity(4)
    T_r[0:3, 0:3] = my_ep_model.rotation_matrix_z_to_gaze_light

    T = np.identity(4)
    T[0:3, 3] = t_c

    my_xyzT = mapping(my_xyz, T_t)
    my_xyzT_T = mapping(my_xyzT, T_)
    my_xyzRT_T = mapping(my_xyzT_T, T_r)
    my_xyzf = mapping(my_xyzRT_T, T)
    return my_xyzf


def _rotate_eyeplan_model_to_treatment_room_frame(my_ep_model, points):

    cntEYE = my_ep_model.structure_set["cor"].contour_coordinates

    t_c = my_ep_model.centre_of_model - cntEYE

    T_t = np.eye(4)
    T_t[0:3, 3] = t_c

    T_ = np.eye(4)
    T_[0:3, 3] = -my_ep_model.centre_of_model

    T_r = np.eye(4)
    T_r[0:3, 0:3] = my_ep_model.rotation_matrix_z_to_gaze_light

    T = np.eye(4)
    T[0:3, 3] = my_ep_model.centre_of_model

    new_points = mapping(points, T_t)
    new_points2 = mapping(new_points, T_)
    new_points3 = mapping(new_points2, T_r)
    new_points4 = mapping(new_points3, T)

    return new_points4


def _resample_points_in_plane(shell, z):

    points_in_plane = []

    normal = (0, 0, 1)

    clipped = shell.clip(normal=normal, origin=[0, 0, z])
    # clipped.plot()

    cell_indices = clipped.faces

    cell_node_ids = cell_indices.reshape(-1, 4)[:, 1:4].ravel()

    cell_nodes = clipped.points[cell_node_ids]

    filtered_points = list(filter(lambda p: abs(p[2] - z) < 1e-4, cell_nodes))

    for p in filtered_points:
        points_in_plane.append([p[0], p[1], z])

    if len(filtered_points) == 0:
        return np.empty((0, 3))

    return np.array(points_in_plane)


def _infer_polar_and_azimuth(config, point):

    polar_length = math.sqrt(point[0]**2 + point[1]**2)

    inferred_polar = math.degrees(math.atan(polar_length/config["FID"]))

    #r = np.sqrt(point[0]**2 + point[1]**2)

    # y = r*np.sin(np.radians(50))
    # x = r*np.cos(np.radians(50))

    # Vector in light fixation plane
    # vec = point - np.asarray([0,0,0])

    # inferred_azimuth = angle(vec, np.asarray([1, 0, 0]), 1)
    
    #inferred_azimuth = np.degrees(np.arcsin(point[1]/r))
    inferred_azimuth = np.degrees(math.atan2(point[1], point[0]))

    if inferred_azimuth < 0:
        inferred_azimuth += 360

    return inferred_polar, inferred_azimuth





def _points_behind_axis_point(points, axis_point, axis_direction):
    """
    Find the subset of points in a point cloud that lie behind the axis point.

    Parameters:
    - points (numpy.ndarray): Nx3 array representing the point cloud.
    - axis_point (numpy.ndarray): 3D coordinates of a point.
    - axis_direction (numpy.ndarray): 3D vector representing the direction.

    Returns:
    - numpy.ndarray: Subset of points in the point cloud that lie behind the axis point.
    """
    # Normalize the axis direction vector to ensure consistency
    axis_direction_normalized = axis_direction / np.linalg.norm(axis_direction)
    
    # Vector from the axis point to each point in the point cloud
    vec_to_points = points - axis_point

    # Project each vector onto the axis direction to get their scalar projection
    scalar_projection = np.dot(vec_to_points, axis_direction_normalized)

    # Check if the points are behind the axis point
    points_behind = scalar_projection < 0

    # Return the subset of points behind the axis point
    return points[points_behind]




def _generate_ring_points(center_point, axis_direction, num_points_per_ring, radius):
    """
    Generate points in a ring around a center point, orthogonal to a given vector.

    Parameters:
    - center_point (numpy.ndarray): 3D coordinates of the center point of the ring.
    - axis_direction (numpy.ndarray): 3D vector representing the direction orthogonal to the rings.
    - num_points_per_ring (int): Number of points to be generated per ring.
    - radius (float): Radius of the ring.

    Returns:
    - numpy.ndarray: Nx3 array representing the generated points in the ring.
    """
    # Create an orthonormal basis with the axis direction as one of the basis vectors
    basis_vec_1 = np.array([1, 0, 0])  # Arbitrarily chosen
    basis_vec_2 = np.cross(axis_direction, basis_vec_1)
    basis_vec_1 = np.cross(basis_vec_2, axis_direction)
    basis_vec_1 /= np.linalg.norm(basis_vec_1)
    basis_vec_2 /= np.linalg.norm(basis_vec_2)

    # Generate points in a ring around the center point
    ring_points = []
    for i in range(num_points_per_ring):
        angle = 2 * np.pi * i / num_points_per_ring
        offset = radius * (basis_vec_1 * np.cos(angle) + basis_vec_2 * np.sin(angle))
        ring_points.append(center_point + offset)

    # Convert the list of ring points to a numpy array
    ring_points = np.array(ring_points)

    return ring_points


def _calculate_radius_around_point(point_cloud, center_point, axis_direction, tolerance=0.1):
    """
    Calculate the approximate radius around a point to the point cloud along a vector.

    Parameters:
    - point_cloud (numpy.ndarray): Nx3 array representing the point cloud.
    - center_point (numpy.ndarray): 3D coordinates of the center point.
    - axis_direction (numpy.ndarray): 3D vector representing the direction orthogonal to the rings.
    - tolerance (float): Tolerance for filtering points near the plane. Points within this tolerance are considered part of the plane.

    Returns:
    - float: Approximate radius around the center point to the point cloud.
    """
    # Calculate the signed distances from each point to the plane orthogonal to the axis direction
    axis_direction = axis_direction / np.linalg.norm(axis_direction)
    signed_distances = np.dot(point_cloud - center_point, axis_direction)

    # Filter points near the plane based on the tolerance
    filtered_points = point_cloud[np.abs(signed_distances) <= tolerance]

    # Calculate distances from the center point to the filtered points
    distances = np.linalg.norm(filtered_points - center_point, axis=1)

    # Find the maximum distance (radius)
    radius = np.mean(distances)

    return radius



def _find_min_max_projection(point_cloud, point, axis_direction, min_lens = 0.1, max_lens = 0.6):
    """
    Find the 3D coordinates of the points corresponding to 20% and 80% of the length of the point cloud
    along the axis direction.

    Parameters:
    - point_cloud (numpy.ndarray): Nx3 array representing the point cloud.
    - point (numpy.ndarray): 3D coordinates of the point.
    - axis_direction (numpy.ndarray): 3D vector representing the direction of the axis.

    Returns:
    - tuple: ((x_20, y_20, z_20), (x_80, y_80, z_80))
    """
    projected_points = _project_point_cloud(point_cloud, point, axis_direction)
    
    # Sort projected points along the axis direction
    sorted_points = projected_points[np.argsort(projected_points[:, 0])]

    # Calculate the indices corresponding to 20% and 80% of the length
    index_20 = int(min_lens * len(sorted_points))
    index_80 = int(max_lens * len(sorted_points))

    # Get the points corresponding to the 20% and 80% positions
    point_20 = sorted_points[index_20]
    point_80 = sorted_points[index_80]
    
    
    points = []
    points.append(point_20)
    points.append(point_80)
    
    sorted_points = sorted(points, key=lambda point: point[2])
    

    return sorted_points


def _project_point_cloud(point_cloud, point, axis_direction):
    """
    Project each point in the point cloud onto the given axis vector.

    Parameters:
    - point_cloud (numpy.ndarray): Nx3 array representing the point cloud.
    - point (numpy.ndarray): 3D coordinates of the point.
    - axis_direction (numpy.ndarray): 3D vector representing the direction of the axis.

    Returns:
    - numpy.ndarray: Array of projected points.
    """
    axis_direction = axis_direction / np.linalg.norm(axis_direction)
    return point + np.dot(point_cloud - point, axis_direction)[:, None] * axis_direction




def _filter_points_by_radius(points, axis_point, axis_vector, radius):
    """
    Filter a point cloud to remove points that are less than a given radius from an axis defined by a point and a vector.

    Parameters:
    - points (numpy.ndarray): Nx3 array representing the points in the point cloud.
    - axis_point (numpy.ndarray): 3D coordinates of a point on the axis.
    - axis_vector (numpy.ndarray): 3D vector representing the direction of the axis.
    - radius (float): The radius within which points are to be retained.

    Returns:
    - numpy.ndarray: Filtered points.
    """
    def point_to_line_distance(point, line_point, line_vector):
        """
        Calculate the perpendicular distance from a point to a line defined by a point and a vector.

        Parameters:
        - point (numpy.ndarray): 3D coordinates of the point.
        - line_point (numpy.ndarray): 3D coordinates of a point on the line.
        - line_vector (numpy.ndarray): 3D vector representing the direction of the line.

        Returns:
        - float: Perpendicular distance from the point to the line.
        """
        return np.linalg.norm(np.cross(line_vector, point - line_point)) / np.linalg.norm(line_vector)

    distances = np.array([point_to_line_distance(point, axis_point, axis_vector) for point in points])
    return points[distances >= radius]


# def densify_point_cloud(points, resolution = 0.5):
#     """
#     Densify a point cloud by interpolating new points.

#     Parameters:
#     - points (numpy.ndarray): Nx3 array representing the original point cloud.
#     - resolution (float): The desired resolution for the densified point cloud.

#     Returns:
#     - numpy.ndarray: Mx3 array representing the densified point cloud.
#     """
#     # Define the bounding box of the original point cloud
#     x_min, y_min, z_min = np.min(points, axis=0)
#     x_max, y_max, z_max = np.max(points, axis=0)

#     # Generate a grid of points within the bounding box with the desired resolution
#     x_grid, y_grid, z_grid = np.meshgrid(np.arange(x_min, x_max, resolution),
#                                          np.arange(y_min, y_max, resolution),
#                                          np.arange(z_min, z_max, resolution))

#     # Interpolate the values of the new points based on the original points
#     dense_points = griddata(points[:, :2], points[:, 2], (x_grid, y_grid), method='linear')

#     # Filter out NaN values
#     valid_indices = ~np.isnan(dense_points)
#     dense_points = np.column_stack((x_grid.ravel()[valid_indices.ravel()],
#                                     y_grid.ravel()[valid_indices.ravel()],
#                                     dense_points.ravel()[valid_indices.ravel()]))

#     return dense_points




# def filter_for_neighbors(points_to_be_filtered, R, epsilon=1e-6):
#     """
#     Filter points that have at least 2 neighbors within distance R.

#     Parameters:
#         points_to_be_filtered (numpy.ndarray): Array of shape (n, 3) representing points in 3D space.
#         R (float): Radius within which points are considered near each other.
#         epsilon (float, optional): Small distance to exclude points that are too close. Defaults to 1e-6.

#     Returns:
#         numpy.ndarray: Subset of points that have at least 2 neighbors within distance R.
#     """
#     tree = KDTree(points_to_be_filtered)
#     points_of_interest = np.zeros(len(points_to_be_filtered), dtype=bool)
    
#     for i, query_point in enumerate(points_to_be_filtered):
#         # Find the neighbors within distance R, excluding the point itself
#         neighbors = tree.query_ball_point(query_point, r=R)
        
#         # Exclude neighbors that are too close to the query point
#         filtered_neighbors = [idx for idx in neighbors if np.linalg.norm(query_point - points_to_be_filtered[idx]) > epsilon]
        
#         # Check if the point has at least 2 valid neighbors (including itself) within distance R
#         if len(filtered_neighbors) >= 5:  # 2 because the point itself is included in neighbors
#             points_of_interest[i] = True
    
#     # Return the subset of points that have at least 2 neighbors within distance R
#     return points_to_be_filtered[points_of_interest]




    
# def sphere_func(coords, x0, y0, z0, r):
#     """Defines a sphere function."""
#     x, y, z = coords
#     return (x - x0)**2 + (y - y0)**2 + (z - z0)**2 - r**2

# def fit_sphere(points):
#     """Fits a sphere to a set of points."""
#     x = points[:, 0]
#     y = points[:, 1]
#     z = points[:, 2]
#     # Initial guess for the sphere parameters (center and radius)
#     initial_guess = (np.mean(x), np.mean(y), np.mean(z), np.std(x + y + z))
#     # Fit the sphere parameters using curve_fit
#     popt, pcov = curve_fit(sphere_func, (x, y, z), np.zeros_like(x), p0=initial_guess)
#     # Extract the sphere parameters
#     x0, y0, z0, r = popt
#     return (x0, y0, z0), r

# def filter_sphere(my_points, sphere_points):
#     """
#     Filter points in my_points that are outside the sphere defined by sphere_points.

#     Parameters:
#     - my_points: numpy array of shape (N, 3) representing points to be filtered
#     - sphere_points: numpy array of shape (M, 3) representing points defining the sphere

#     Returns:
#     - outside_points_array: numpy array of shape (P, 3) containing points from my_points
#                             that are outside the sphere
#     """

#     # Fit a sphere to the sphere_points
#     sphere_center, sphere_radius = fit_sphere(sphere_points)

#     # Compute the squared distance of each point from the sphere center
#     squared_distances = np.sum((my_points - sphere_center)**2, axis=1)

#     # Find points that are outside the sphere
#     outside_points = squared_distances > sphere_radius**2

#     # Get the points that are outside the sphere
#     outside_points_array = my_points[outside_points]

#     return outside_points_array


# def _get_feature_labels():

#     return  [
#      'target com x',
#      'target com y',
#      'target com z',
#      'optdisk com x',
#      'optdisk com y',
#      'optdisk com z',
#      'macula com x',
#      'macula com y',
#      'macula com z',
#      'target max x x',
#      'target max x y',
#      'target max x z',
#      'target min x x',
#      'target min x y',
#      'target min x z',
#      'target max y x',
#      'target max y y',
#      'target max y z',
#      'target min y x',
#      'target min y y',
#      'target min y z',
#      'target max z x',
#      'target max z y',
#      'target max z z',
#      'target min z x',
#      'target min z y',
#      'target min z z',

#      #F10
#      'apex closest x',
#      'apex closest y',
#      'apex closest z',

#      #F11
#      'apex furthest x',
#      'apex furthest y',
#      'apex furthest z',

#      #f12 target_point_closest_to_optdisk
#      'target.closest optdisk x',
#      'target.closest optdisk y',
#      'target.closest optdisk z',

#      #f13 target_point_closest_to_macula
#      'target.closest macula x',
#      'target.closest macula y',
#      'target.closest macula z',


#      # Scalars
#      'target eye',
#      'use of wedge',
#      'target basearea',
#      'target.com - Optdisk',
#      'target.com - Macula',
#      'target.closest - Optdisk',
#      'target.closest - Macula',
#      'apex height']
