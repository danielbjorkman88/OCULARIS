# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.
"""

import copy
import numpy as np
import dose_engine

import patient_surfaces
import healthy_tissue_probabilties
from utils.normalization_utils import find_DX
from utils.motion_utils import rotate_points
from scipy.interpolate import interp1d

from pyproton.metrics.dvh import DVH
import pandas as pd
from generators import  dose_engine_generators
from configuration import config, data

from pyproton.volume.grid import Grid

from utils.color_utils import my_colors
from utils.geo_utils import find_maximum_radius, points_outside_cylinder, filter_points_within_radius
from models.eyeplan_model import _find_min_max_projection

import math

def find_intersect_plane(p, vector):
    p1 = p #self.orientation_vec[:3]
    p2 = p1 + vector*50

    ray_eye_axis = dose_engine.Ray(p1,p2)
    return ray_eye_axis.find_intersect_with_plane(np.asarray([0,0,1]), np.asarray([0,0, config['FID']]) )




class Orientation:
    def __init__(self, orientation_config):
        
        self.orientation_config = orientation_config
        # self.orientation_vec = orientation_vec
        self.ref_orientation = None
        self.dose = None
        self.structure_set = None
        self.ep_model = None
        self.patient = None
        self.orientation_vec = None
        self._fixation_point = None
        self._polar_angle = None
        self.patient_alphas = None
        self.algo = None
        
        
        if "ref" in orientation_config.keys():
            self.ref = True
            self.orientation_vec = orientation_config["orientation"]
            self.dose = orientation_config["dose"]
            self.ep_model = orientation_config["ep_model"]
            self.structure_set = self.ep_model.structure_set_clips_registered
            self.patient = orientation_config["patient"]
            self.algo = orientation_config["algo"]
            
        else:
            self.ref = False
            self.ref_orientation = orientation_config["ref_orientation"]
            self.ep_model = copy.deepcopy(orientation_config["ep_model"])
            self.patient = orientation_config["patient"]
            self.gaze_vector = orientation_config["gaze_vector"]
            self.orientation_vec = np.concatenate((self.ep_model.centre_of_model, self.gaze_vector)) 
            
            
        self.skin_plane_ref = self.orientation_config["skin_plane"]
            
        
        self.max_doses_structures = {}
        
            
        self.masks = {}
        
        self.NTCPs = {}
        
        self.important_dose_points = {}
            
        # self.patient_alphas = orientation_config["patient_alphas"]
        
            # orientation_config["gaze_vector"] = new_gaze_vector
            # orientation_config["ref_orientation"] = orientation_start
            # orientation_config["ep_model"] = ep_model
            # orientation_config["patient"] = pat1
            
            
        # orientation_config["orientation"] = orientation_clinical
        # orientation_config["dose"] = algo.dose
        # orientation_config["structure_set"] = ep_model.structure_set_clips_registered
        # orientation_config["patient"] = pat1
        
        
        
        self.instance_configured = False
        self.structure_set_derived = False
        self.algo_defined = False

        self.dvh_curves = {}
        self.metrics =  {}
        self.metric_df = None
        self.structurenames = ["eyeglobe", "lens", "cornea", "target", "macula", "optdisk", "retina", "ciliary_body"]
        self.oar_structure_names = ["eyeglobe", "lens", "cornea", "macula", "optdisk", "retina", "ciliary_body"]
        
        self.new_points = None
        self.target_com_naive = None
        self.target_com_binary = None
        
        self.rotated_ref_target_com = None
        self.intersect_point = None
        
        
        print("orientation_vec initialized", self.orientation_vec)
    
    @property
    def fixation_point(self):
        
        self._fixation_point = find_intersect_plane(self.orientation_vec[0:3] , self.orientation_vec[3:])
        
        return self._fixation_point
    
    @property
    def polar_angle(self):
        
        fixation_point = self.fixation_point
        
        r = math.sqrt( fixation_point[0]**2 + fixation_point[1]**2)
        
        # polar_length = math.sqrt(point[0]**2 + point[1]**2)

        self._polar_angle = math.degrees(math.atan(r/132.5))
        
        
        return self._polar_angle
    
    @property
    def azimuth_angle(self):
        
        fixation_point = self.fixation_point    

        inferred_azimuth = np.degrees(math.atan2(fixation_point[1], fixation_point[0]))

        if inferred_azimuth < 0:
            inferred_azimuth += 360
            
        self._azimuth_angle = inferred_azimuth
        
        return self._azimuth_angle
        
        
        
    def config_orientation(self):
        
        if self.ref == True:
            self.config_ref()
            # self.extract_info()
        else:
            if self.instance_configured == False:
                self.config_offset()
            if self.structure_set_derived == False:
                self.derive_new_structure_set()
            if self.algo_defined == False:
                self.define_algo()
        
        self.extract_info()
        
    
    def config_ref(self):
        new_points = {}
        new_binary_points = {}
        for structure_name in self.structurenames:
            points = self.structure_set[structure_name].structure_points
            new_points[structure_name] = points
            
            # points = self.structure_set[structure_name].points_of_binary_mask #
            # # points = rotate_points(points, ep_model.centre_of_model, ep_model.gaze_vector_light, new_gaze_vector)
            # new_binary_points[structure_name] = points
            
        self.new_points = new_points
        # self.new_binary_points = new_binary_points
        
        target_points = self.new_points["target"]
        self.target_com_naive =  np.asarray([(np.max(target_points[0:, 0]) - np.min(target_points[0:, 0]))/2 + np.min(target_points[0:, 0]),
                                  (np.max(target_points[0:, 1]) - np.min(target_points[0:, 1]))/2 + np.min(target_points[0:, 1]),
                                  (np.max(target_points[0:, 2]) - np.min(target_points[0:, 2]))/2 + np.min(target_points[0:, 2])])
        
        
        # target_points = self.new_binary_points["target"]
        # self.target_com_binary =  np.asarray([(np.max(target_points[0:, 0]) - np.min(target_points[0:, 0]))/2 + np.min(target_points[0:, 0]),
        #                           (np.max(target_points[0:, 1]) - np.min(target_points[0:, 1]))/2 + np.min(target_points[0:, 1]),
        #                           (np.max(target_points[0:, 2]) - np.min(target_points[0:, 2]))/2 + np.min(target_points[0:, 2])])


        
        
        self.rotated_ref_target_com = self.target_com_naive
        
    def config_offset(self):
        
        skin_plane_z = self.orientation_config["skin_plane"]
        anterior_alignment_z = skin_plane_z + 2
        
        new_points = {}
        new_binary_points = {}
        for structure_name in self.structurenames:
            points = self.ref_orientation.structure_set[structure_name].structure_points#
            points = rotate_points(points, self.ep_model.centre_of_model, self.ep_model.gaze_vector_light, self.gaze_vector)
            

            
            new_points[structure_name] = points
            
            points = self.ref_orientation.structure_set[structure_name].points_of_binary_mask#
            points = rotate_points(points, self.ep_model.centre_of_model, self.ep_model.gaze_vector_light, self.gaze_vector)
            new_binary_points[structure_name] = points
            
            
        
        
        # Clean up any retina
        radius = find_maximum_radius(new_points["eyeglobe"], self.ep_model.centre_of_model, self.gaze_vector)
        infront_of = -6
        behind = 6
        new_points["retina"] = points_outside_cylinder(new_points["retina"], self.ep_model.centre_of_model, self.gaze_vector, radius*0.8, infront_of, behind)
        
        
        # Clean up any ciliary_body
        lens_min, lens_max = _find_min_max_projection(new_points["lens"], self.ep_model.centre_of_model, self.gaze_vector, min_lens = 0.01, max_lens = 0.99)
        radius = find_maximum_radius(new_points["lens"], self.ep_model.centre_of_model, self.gaze_vector)
        new_points["ciliary_body"] = filter_points_within_radius(new_points["ciliary_body"], self.ep_model.centre_of_model, self.gaze_vector, radius)
        
            
            
        self.new_points = new_points
        self.new_binary_points = new_binary_points
        
        assert set( list(self.new_points.keys())  ).issubset(set(self.structurenames))
        
        

        target_points = self.new_points["target"]
        self.target_com_naive =  np.asarray([(np.max(target_points[0:, 0]) - np.min(target_points[0:, 0]))/2 + np.min(target_points[0:, 0]),
                                  (np.max(target_points[0:, 1]) - np.min(target_points[0:, 1]))/2 + np.min(target_points[0:, 1]),
                                  (np.max(target_points[0:, 2]) - np.min(target_points[0:, 2]))/2 + np.min(target_points[0:, 2])])

        target_points = self.new_binary_points["target"]
        self.target_com_binary =  np.asarray([(np.max(target_points[0:, 0]) - np.min(target_points[0:, 0]))/2 + np.min(target_points[0:, 0]),
                                  (np.max(target_points[0:, 1]) - np.min(target_points[0:, 1]))/2 + np.min(target_points[0:, 1]),
                                  (np.max(target_points[0:, 2]) - np.min(target_points[0:, 2]))/2 + np.min(target_points[0:, 2])])



        self.rotated_ref_target_com = rotate_points(self.ref_orientation.target_com_naive, self.ep_model.centre_of_model, self.ep_model.gaze_vector_light, self.gaze_vector)
        
        
        self.translation_vector = - self.target_com_naive
        
        # print("Before translation", self.orientation_vec)
        for structure_name in self.structurenames:
            points = new_points[structure_name]
            new_points[structure_name] = points + self.translation_vector
            
        eye_centre = self.orientation_vec[0:3]
        # print(eye_centre , self.translation_vector)
        self.orientation_vec[0:3] = eye_centre + self.translation_vector
            
        # print("After translation", self.orientation_vec)
        
        assert set(list(self.new_points.keys())).issubset(set(self.structurenames))
        
        self.instance_configured = True
            
    
    def define_algo(self):

        
        self.skin_plane_ref = self.orientation_config["skin_plane"]

        assert self.structure_set_derived == True
        
        conf = self.patient.generate_config_copy

        conf["eyeglobe_point"] = self.new_points['eyeglobe']
        conf["target_points"] = self.new_points['target']
      
        # assert len(conf["surfaces"]) == 1
        
        
        cloud_surface = patient_surfaces.PointCloud(conf["eyeglobe_point"])
        
        # conf["surfaces"].append(cloud_surface)
        
        points = self.structure_set["eyeglobe"].structure_points#
        
        index = np.where(points == np.max(points[:, 2]))[0][0]
        most_anterior_z = points[index, 2]
        
        skin_plane_z = most_anterior_z - 2
        
        skin_plane_point = np.asarray([0,0,skin_plane_z])
        skin_plane_normal = np.asarray([0,0,1])
        
        
        diff_skin_planes = self.skin_plane_ref - skin_plane_z
        
        
        origin = self.ep_model.grid.origin
        origin[2] = origin[2] - diff_skin_planes
        
        spacing = self.ep_model.grid.spacing
        
        size = self.ep_model.grid.size
        
        grid = Grid(origin, spacing, size)
        
        self.ep_model.grid = grid
        
        
        self.orientation_config["skin_plane"] = skin_plane_z
        
        
        plane_surface = patient_surfaces.SkinPlane(skin_plane_point, skin_plane_normal)
        #orientation_config["eyelid"] = surface_eyelid
        if "eyelid" in self.orientation_config.keys():
            conf["surfaces"] = [plane_surface, cloud_surface, self.orientation_config["eyelid"]]
        else:
            conf["surfaces"] = [plane_surface, cloud_surface]
        
        
        conf["skin_plane"] = skin_plane_z
        conf["skin_plane_point"][2] = skin_plane_z

        algo2 = dose_engine_generators.generate_reference(data, conf, self.ep_model)
        
        
        self.algo = algo2
        
        self.dose = algo2.dose
        
        self.algo_defined = True
    
    def derive_new_structure_set(self):
        structure_set = self.ep_model.generate_structure_set(self.new_points)
        self.structure_set = structure_set
        
        assert self.structure_set["retina"]._points_resampled == False
        assert self.structure_set["ciliary_body"]._points_resampled == False
        
        assert set( self.structurenames ).issubset(set(self.structure_set.name_list))
        
        
        self.ep_model._structure_set_clips_registered = structure_set
        
        
        self.ep_model.surfaces = self.orientation_config["surfaces"]
        
        self.structure_set_derived = True
        
    def extract_info(self):
        
        # target_com = self.structure_set['target'].com

        # target_com =  np.asarray([(np.max(target_points[0:, 0]) - np.min(target_points[0:, 0]))/2 + np.min(target_points[0:, 0]),
        #                           (np.max(target_points[0:, 1]) - np.min(target_points[0:, 1]))/2 + np.min(target_points[0:, 1]),
        #                           (np.max(target_points[0:, 2]) - np.min(target_points[0:, 2]))/2 + np.min(target_points[0:, 2])])

        
        # print("target_com", target_com)
        # self.target_com_naive = target_com
        
        lens_structure = self.structure_set["lens"]
        lens_structure.resample_contour_points_in_grid(1)
        
        eyeglobe_structure = self.structure_set["eyeglobe"]
        eyeglobe_structure.resample_contour_points_in_grid(1)
        
        dose_fractions = np.linspace(0, 1.1, 1000)
        
        # curves = {}
        
        for structure_name in self.structurenames:
        
            struct = self.structure_set[structure_name]
            
            if structure_name not in ["optdisk",  "retina", "ciliary_body"]:
                struct.resample_contour_points_in_grid(1)


            
            self.metrics[structure_name] = {}
            
            if structure_name in [ "retina", "ciliary_body"]:   
                
                print(f":::::::::::::::::::::: {structure_name} ::::::::::::::::::::::")
                # print("Was", np.sum(self.structure_set[structure_name].binary_mask))
                
                
                # assert self.structure_set["retina"]._points_resampled == False
                
                points = struct.contour_coordinates
                
                
                assert len(points) > 1
                
                # points = filter_for_neighbors(points, 0.6)
                # points = filter_for_neighbors(points, 0.6)
                
                
                # if structure_name == "retina":
                #     points = densify_point_cloud(points, 0.5)
                
                # binary_mask = np.zeros(self.dose.shape)
                # binary_mask[binary_mask == 0] = False
                
                binary_mask = np.zeros(self.dose.shape, dtype=bool)

                # # If you want to explicitly make sure that all elements are False
                # # you can do this:
                # binary_array.fill(False)
                
                
                grid_origin = self.ep_model.grid.origin
                grid_spacing = self.ep_model.grid.spacing
                grid_size = self.ep_model.grid.size
                voxel_indices = ((points - grid_origin) / grid_spacing).astype(int)
                within_grid = np.all((voxel_indices >= 0) & (voxel_indices < grid_size), axis=1)
                
                binary_mask[tuple(voxel_indices[within_grid].T)] = True
                binary_mask[lens_structure.binary_mask] = False
                
                if structure_name in ["ciliary_body"]:   
                    binary_mask = binary_mask * eyeglobe_structure.binary_mask # [lens_structure.binary_mask] = False
                
                # print(binary_mask)
                
                self.structure_set[structure_name].binary_mask = binary_mask
                
                # self.structure_set[structure_name].binary_mask[tuple(voxel_indices[within_grid].T)] = True
                
                # self.structure_set[structure_name].binary_mask[lens_structure.binary_mask] = False
                
                
                
                print("Now", np.sum(self.structure_set[structure_name].binary_mask), "elements")
                
                
                # binary_mask = np.zeros(algo.dose.shape)
                # for point in points:
                #     x_idx = int((point[0] - algo.medium.mesh_origin[0]) / algo.medium.resolution[0]) 
                #     y_idx = int((point[1] - algo.medium.mesh_origin[1]) / algo.medium.resolution[1]) 
                #     z_idx = int((point[2] - algo.medium.mesh_origin[2]) / algo.medium.resolution[2]) 
                #     if 0 <= x_idx < binary_mask.shape[0] and 0 <= y_idx < binary_mask.shape[1] and 0 <= z_idx < binary_mask.shape[2]:
                #         binary_mask[x_idx, y_idx, z_idx] = 1
                # self.structure_set["retina"].binary_mask = binary_mask
            
            # if structure_name == "ciliary_body":
            #     points = struct.contour_coordinates
                # binary_mask = np.zeros(algo.dose.shape)
                # for point in points:
                #     x_idx = int((point[0] - algo.medium.mesh_origin[0]) / algo.medium.resolution[0]) 
                #     y_idx = int((point[1] - algo.medium.mesh_origin[1]) / algo.medium.resolution[1]) 
                #     z_idx = int((point[2] - algo.medium.mesh_origin[2]) / algo.medium.resolution[2]) 
                #     if 0 <= x_idx < binary_mask.shape[0] and 0 <= y_idx < binary_mask.shape[1] and 0 <= z_idx < binary_mask.shape[2]:
                #         binary_mask[x_idx, y_idx, z_idx] = 1
                # self.structure_set["ciliary_body"].binary_mask = binary_mask
                        
            if structure_name == "optdisk":
                struct.resample_contour_points_in_grid(1)
                
                print(":::::::::::::::::::::: optdisk ::::::::::::::::::::::")
                print("Was", np.sum(self.structure_set[structure_name].binary_mask))
                
                optdisk_coordinates = self.structure_set["optdisk"].contour_coordinates
                grid_origin = self.ep_model.grid.origin
                grid_spacing = self.ep_model.grid.spacing
                grid_size = self.ep_model.grid.size
                
                voxel_indices = ((optdisk_coordinates - grid_origin) / grid_spacing).astype(int)
                within_grid = np.all((voxel_indices >= 0) & (voxel_indices < grid_size), axis=1)
                self.structure_set["optdisk"].binary_mask[tuple(voxel_indices[within_grid].T)] = True
                
                print("Now", np.sum(self.structure_set[structure_name].binary_mask))
                
                
                D_mean = np.mean(self.dose[self.structure_set[structure_name].binary_mask])
                
                self.important_dose_points["optdisk_D_mean"] = D_mean
                
                # for x in np.linspace(0.05, 1, 20):
                #     self.metrics[structure_name][x] = 0 
                # self.dvh_curves[structure_name] = (np.linspace(0.05, 1, 20), np.zeros(20))
                # print("Warning. Binary mask is empty")
                # continue
            
            result_mask = struct.binary_mask
            
            if structure_name == "eyeglobe":
                
                struct = self.structure_set["target"]
                struct.resample_contour_points_in_grid(1)
                
                target_mask = struct.binary_mask
                result_mask = result_mask & ~target_mask
            
            if structure_name == "macula":
                self.masks[structure_name] = result_mask
            
            
            if np.sum(result_mask) == 0:
                max_value = 0
            else:
                max_value = np.max(self.dose[result_mask])
            
            self.max_doses_structures[structure_name] = max_value
            
            dvh = DVH(self.dose, result_mask)

            volumes = dvh.V(dose_fractions)
            
            self.dvh_curves[structure_name] = (dose_fractions, volumes)
            
            
            
            
            f = interp1d(dose_fractions, volumes )
            
            for x in np.linspace(0.05, 1, 20):
            
                V_val =  f(x)
                
                self.metrics[structure_name][x] = V_val.item()
                
            
        self.metric_df = pd.DataFrame.from_dict(self.metrics, orient='index')

        # vec = self.metric_df[0.05]
        
        # v5s = np.zeros(4)
        
        # for i, stucture_name in enumerate(["eyeglobe", "lens", "cornea", "macula"]):
        
        #     v5s[i] = float(vec[stucture_name])
            

        # self.v5s = v5s
        
        self.calc_NTCPs()



    @property            
    def weighted_sum(self):
        
        w_sum = 0
        
        
        tissue_weights_df = self.orientation_config["tissue_weights"]
        
        weight_macula = tissue_weights_df['w_macula'].item()
        weight_optic_disc = tissue_weights_df['w_optdisk'].item()
        weight_cornea = tissue_weights_df['w_cornea'].item()
        weight_lens = tissue_weights_df['w_lens'].item()
        weight_retina = tissue_weights_df['w_retina'].item()
        weight_ciliary_body = tissue_weights_df['w_ciliary_body'].item()
        
        
        # weight_macula, weight_optic_disc, weight_cornea, weight_lens = self.patient_alphas
        #weight_lens = 1
        
        for structure_name, color in zip(self.oar_structure_names, my_colors): 
        
            curves = self.dvh_curves[structure_name]
            xes, yes = curves[0], curves[1]
            
            if structure_name == "macula":
                D_02 = find_DX(xes, yes, 0.02)
                w_sum += weight_macula*D_02
                self.important_dose_points["macula_D02"] = D_02
       
            if structure_name == "optdisk":
                D_20 = find_DX(xes, yes, 0.2)
                w_sum += weight_optic_disc*D_20
                self.important_dose_points["optdisk_D20"] = D_20
                
            if structure_name == "cornea":
                D_20 = find_DX(xes, yes, 0.2)
                w_sum += weight_cornea*D_20
                self.important_dose_points["cornea_D20"] = D_20
    
            if structure_name == "lens":
                D_05 = find_DX(xes, yes, 0.05)
                w_sum += weight_cornea*D_05
                self.important_dose_points["lens_D5"] = D_05
                
            if structure_name == "retina":
                f = interp1d(xes, yes )
                retina_V52 = f(0.915) #52*0.96*1.1/60 = 0.915, 0.96 from the different fractionation scheme by Espensen
                w_sum += weight_retina*retina_V52
                self.important_dose_points["retina_V52"] = retina_V52
                
            if structure_name == "ciliary_body":
                f = interp1d(xes, yes )
                ciliary_body_V26 = f(0.457) #26*0.96*1.1/60 = 0.457, 0.96 from the different fractionation scheme by Espensen
                w_sum += weight_ciliary_body*ciliary_body_V26
                
                self.important_dose_points["ciliary_body_V26"] = ciliary_body_V26
                
                

        
        for i, structure_name in enumerate(self.oar_structure_names):
            
            
            if structure_name == "macula" :
                weight = weight_macula  
            elif structure_name == "optdisk":
                weight = weight_optic_disc
            elif structure_name == "cornea":
                weight = weight_cornea
            elif structure_name == "lens":
                weight = weight_lens
            elif structure_name == "retina":
                weight = weight_retina
            elif structure_name == "ciliary_body":
                weight = weight_ciliary_body 
            else:
                weight = 1
                   
    
            bias = 0
                
            row_data = np.asarray(self.metric_df.loc[structure_name])
            
            w_sum += np.sum(weight*(row_data - bias)) / len(row_data)

        return w_sum


  


    @property
    def intersect_with_fixation_plane(self):
        p1 = self.orientation_vec[:3]
        p2 = p1 + self.orientation_vec[3:]*50

        ray_eye_axis = dose_engine.Ray(p1,p2)
        self.intersect_point = ray_eye_axis.find_intersect_with_plane(np.asarray([0,0,1]), np.asarray([0,0,self.algo.config['FID']]) )

        return self.intersect_point
        
        
    def __repr__(self):
        return f"{self.orientation_vec}"
    
    
    def calc_NTCPs(self):
        
        xes, yes = self.dvh_curves["macula"]
        f = interp1d(xes, yes )
        macula_V28 = f(0.28)
        
        
        xes, yes = self.dvh_curves["optdisk"]
        f = interp1d(xes, yes )
        optdisk_V28 = f(0.28)
        
        self.NTCPs["increased_risk_neovascular_glaucoma_Mishra"] = healthy_tissue_probabilties.increased_risk_neovascular_glaucoma_Mishra(macula_V28, optdisk_V28)
        
        structure_name = "ciliary_body"
        curves = self.dvh_curves[structure_name]
        xes, yes = curves[0], curves[1]
        f = interp1d(xes, yes )
        ciliary_body_V26 = f(0.457) #26*0.96*1.1/60 = 0.457, 0.96 from the different fractionation scheme by Espensen
        
        self.important_dose_points["ciliary_body_V26"] = ciliary_body_V26
        self.NTCPs["probability_cataract"] = healthy_tissue_probabilties.probability_cataract(ciliary_body_V26)
        


        structure_name = "retina"
        curves = self.dvh_curves[structure_name]
        xes, yes = curves[0], curves[1]
        f = interp1d(xes, yes )
        retina_V52 = f(0.915) #52*0.96*1.1/60 = 0.915, 0.96 from the different fractionation scheme by Espensen
        
        self.important_dose_points["retina_V52"] = retina_V52
        
        #V5Gy, V10Gy, V20Gy, V30Gy
        retina_V5Gy = f(5/60)
        retina_V10Gy = f(10/60)
        retina_V20Gy = f(20/60)
        retina_V30Gy = f(30/60)
        
        self.important_dose_points["retina_V5Gy"] = retina_V5Gy
        self.important_dose_points["retina_V10Gy"] = retina_V10Gy
        self.important_dose_points["retina_V20Gy"] = retina_V20Gy
        self.important_dose_points["retina_V30Gy"] = retina_V30Gy
        
        retina_D20 = find_DX(xes, yes, 0.20)
        self.important_dose_points["retina_D20"] = retina_D20

        # retina_V52 = self.important_dose_points["retina_V52"]
        self.NTCPs["probability_retinal_detachment"] = healthy_tissue_probabilties.probability_retinal_detachment(retina_V52)
        
                
        
        # xes, yes = curves[0], curves[1]
        
        
        # # plt.plot(100*xes, 100*yes, label = structure_name.capitalize(), color = color, linewidth = 2)
        
        # if structure_name == "macula":
        #     D_02 = find_DX(xes, yes, 0.02)
        #     plt.scatter(100*D_02, 2, color = color, s=100)
            
            
        #     f = interp1d(xes, yes )
        #     V28 = f(0.28)
        #     plt.scatter(100*0.28, 100*V28, s=100, marker='o', facecolors='none', edgecolors=color)
            
        
        
        # max_macula = self.max_doses_structures["macula"]
        # max_optdisk = self.max_doses_structures["optdisk"]
        
        
        
        
        
        
        for structure_name, color in zip(self.oar_structure_names, my_colors): #orientation_list[0].oar_structure_names
        # # for structure_name, color in orientation_list[idx].structurenames:
            
            
            # structure_name = "macula"
            # color = "C0"
            
            # if structure_name == "optdisk":
            #     continue
        
            #orientation_start.weighted_sum
        
            curves = self.dvh_curves[structure_name]
            xes, yes = curves[0], curves[1]
            
            # xes, yes = curves[0], curves[1]
            # plt.plot(100*xes, 100*yes, label = structure_name.capitalize(), color = color)
            
            if structure_name == "macula":
                D_02 = find_DX(xes, yes, 0.02)
                # plt.scatter(100*D_02, 2, color = color)
                
                self.NTCPs["visual_acuity"] = healthy_tissue_probabilties.probability_visual_acuity(D_02)
        
            if structure_name == "cornea":
                D_20 = find_DX(xes, yes, 0.2)
                # plt.scatter(100*D_20, 20, color = color)
                
                self.NTCPs["neovascular_glaucoma_Espensen"] = healthy_tissue_probabilties.probability_neovascular_glaucoma_Espensen(D_20)
                
            if structure_name == "optdisk":
                D_20 = find_DX(xes, yes, 0.2)
                self.NTCPs["optic_neuropathy"] = healthy_tissue_probabilties.probability_optic_neuropathy(D_20)
    
    
            if structure_name == "lens":
                D_05 = find_DX(xes, yes, 0.05)
                # plt.scatter(100*D_20, 20, color = color)            
    
                self.NTCPs["increased_risk_cataract_Thariat"] = healthy_tissue_probabilties.increased_risk_cataract_Thariat(D_05)
        

# def generate_reference(data, config, ep_model, simple = False):

    
#     assert type(simple) == bool
    
#     algo = dose_engine.BroadBeam(data, config)
    
    
#     if config["anterior_aligned"]:
#         # # Aligns mesh with anterior segment of eyeglobe + 0.5 mm
#         # z_alignment = np.max(ep_model.structure_set_clips_registered["eyeglobe"].contour_coordinates[0:,2]) + 0.5
        
   
#         n_voxels_front_of_skinplane = config["n_voxels_front_of_skinplane"] #int(np.ceil(3/ep_model.doseplane_h.resolution[2]))
#         z_alignment = algo.config['skin_plane_point'][2] + n_voxels_front_of_skinplane*ep_model.doseplane_h.resolution[2]
        
#         config["n_voxels_front_of_skinplane"] = n_voxels_front_of_skinplane
        
#         print(f"Dose_engine: Mesh aligned to {n_voxels_front_of_skinplane*ep_model.doseplane_h.resolution[2]} mm in front of skin plane")
#         # print("Dose_engine: Mesh aligned to 2.5 mm in front of skin plane")
#     else:
#         # Aligns mesh with skin plane
#         z_alignment = algo.config["skin_plane_point"][2]
#         print("Dose_engine: Mesh aligned to skin plane")
    
    
#     translation_vector = np.zeros(3)

    
#     translation_vector[2] = z_alignment - algo.medium.maxZ
#     translation_vector[1] = ep_model.doseplane_h.resolution[1]/2
#     translation_vector[0] = ep_model.doseplane_h.resolution[0]/2
#     algo.medium.translate(translation_vector)
#     assert - 0.0001 < (algo.medium.maxZ -
#                        z_alignment)/z_alignment < 0.0001
    
#     if simple:
#         print("Generated simplified dose engine instance")
#         return algo
    
#     if config["anterior_aligned"]:
#         struct = ep_model.structure_set_clips_registered["eyeglobe"]
#         struct.resample_contour_points_in_grid(1)        
#         algo.medium.medium[struct.binary_mask] = 1
        
#         struct = ep_model.structure_set_clips_registered["target"]
#         struct.resample_contour_points_in_grid(1)        
#         algo.medium.medium[struct.binary_mask] = 1
        
#         struct = ep_model.structure_set_clips_registered["cornea"]
#         struct.resample_contour_points_in_grid(1)        
#         algo.medium.medium[struct.binary_mask] = 1
        
#     algo.define_target_rays_from_points()
#     algo.orthogonal_limit = 15
#     algo.define_rays()
#     algo.config_surface_interaction()
#     algo.determine_target_and_modulation_range_of_target()
#     algo.config_wedge()

#     algo.config_raytracer()
#     algo.config_nozzle()
#     algo.config_depth_dose_profile()
#     algo.config_collimator()  
#     algo.config_distances_to_aperture()
#     algo.calc()
    
#     print("Generated full dose engine instance")

#     return algo



# def generate_reference_anterior_aligned(data, config, ep_model, simple = False):

    
#     assert type(simple) == bool
    
    
#     alignment_z = np.max(ep_model.structure_set_clips_registered["eyeglobe"].contour_coordinates[0:,2])
    
#     algo = dose_engine.BroadBeam(data, config)
#     translation_vector = np.zeros(3)

#     # Aligns mesh with skin plane
#     translation_vector[2] = (alignment_z - algo.medium.maxZ)
#     translation_vector[1] = ep_model.doseplane_h.resolution[1]/2
#     translation_vector[0] = ep_model.doseplane_h.resolution[0]/2
#     algo.medium.translate(translation_vector)
#     assert - 0.0001 < (algo.medium.maxZ -
#                        alignment_z)/alignment_z < 0.0001
    
#     if simple:
#         print("Generated simplified dose engine instance")
#         return algo
    
#     algo.define_target_rays_from_points()
#     algo.orthogonal_limit = 15
#     algo.define_rays()
#     algo.config_surface_interaction()
#     algo.determine_target_and_modulation_range_of_target()
#     algo.config_wedge()

#     algo.config_raytracer()
#     algo.config_nozzle()
#     algo.config_depth_dose_profile()
#     algo.config_collimator()  
#     algo.config_distances_to_aperture()
#     algo.calc()
    
#     print("Generated fuill dose engine instance")

#     return algo


# def generate_subfraction(algo, ep_model, new_points):

#     conf = copy.deepcopy(algo.config)

#     assert "target" in new_points.keys()
#     assert "eyeglobe" in new_points.keys()

#     conf["eyeglobe_point"] = new_points["eyeglobe"]
#     conf["target_points"] = new_points["target"]

#     alg = dose_engine.BroadBeam(algo.data, conf)
#     translation_vector = np.zeros(3)
#     # abs(alg.medium.mesh_apex[2]) + 35.6/2 #skin_plane
#     translation_vector[2] = (conf["skin_plane"] - alg.medium.maxZ)
#     translation_vector[1] = ep_model.doseplane_h.resolution[1]/2
#     translation_vector[0] = ep_model.doseplane_h.resolution[0]/2
#     #alg.config_nozzle_by_mw(604)
#     alg.medium.translate(translation_vector)
#     alg.define_target_rays_from_points()
#     alg.orthogonal_limit = 15
#     alg.define_rays()
#     alg.config_surface_interaction()
#     # alg.determine_target_and_modulation_range_of_target()
#     alg.config_wedge()

#     alg.config_raytracer()
#     alg.nozzle = copy.deepcopy(algo.nozzle)
#     alg.define_geo_relationships()
#     alg.nozzle_configured = True
#     # alg.config_nozzle()
#     alg.config_depth_dose_profile()
#     alg.collimator = copy.deepcopy(algo.collimator)
#     # alg.config_collimator() #- alg.nozzle.inflection_shift_mm[0]
#     alg.config_distances_to_aperture()
#     alg.calc()

#     return alg


# def generate_subfractions_from_motion(algo, ep_model, fraction1, slice_idx, path):
    
#     from plot_utils.motion_plot_utils import plot_motion_subfraction_frame
    
#     investigated_indices = range(fraction1.N_resampled) #[5, 6, 9, 12, 14, 18]
    
#     tot_acc_dose = np.zeros(algo.dose.shape)
#     for idx_of_interest in investigated_indices:
        
        
#         print("----- idx_of_interest ", idx_of_interest)
        
#         acc_dose_subfraction = np.zeros(algo.dose.shape)
#         # acc_dose_subfraction2 = np.zeros(algo.dose.shape)
#         # acc_dose_not_eye = np.zeros(algo.dose.shape)
#         # acc_dose_eye = np.zeros(algo.dose.shape)
        
#         algo_moved_model = generate_subfraction(algo, ep_model, fraction1.points_at_idx(idx_of_interest, True))
        
#         structure_set = ep_model.generate_structure_set(fraction1.points_at_idx(idx_of_interest, True))
#         structure_set.grid = ep_model.grid
        
#         struct = structure_set["eyeglobe"]
#         struct.resample_contour_points_in_grid(1)
#         binary_mask = struct.binary_mask
#         negative_mask = binary_mask == False
#         points_of_binary_mask = struct.points_of_binary_mask
#         transformed_points = fraction1.transform_points_by_idx(points_of_binary_mask, idx_of_interest, True)
            
        
#         acc_dose_subfraction[negative_mask] = algo_moved_model.dose[negative_mask]
#         # acc_dose_subfraction2[negative_mask] = algo_moved_model.dose[negative_mask]
        
#         origin_point_x_bins = np.asarray(np.floor((points_of_binary_mask[0:,0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
#         origin_point_y_bins = np.asarray(np.floor((points_of_binary_mask[0:,1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
#         origin_point_z_bins = np.asarray(np.floor((points_of_binary_mask[0:,2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
                    
        
        
#         sampling_x_bins = np.asarray(np.floor((transformed_points[0:,0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
#         sampling_y_bins = np.asarray(np.floor((transformed_points[0:,1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
#         sampling_z_bins = np.asarray(np.floor((transformed_points[0:,2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
        
#         # In case sampling point is outside the grid, then it sampling instead within the grid
#         sampling_x_bins[sampling_x_bins>algo.dose.shape[0]-1] = algo.dose.shape[0] -1
#         sampling_y_bins[sampling_y_bins>algo.dose.shape[1]-1] = algo.dose.shape[1] -1
#         sampling_z_bins[sampling_z_bins>algo.dose.shape[2]-1] = algo.dose.shape[2] -1
#         sampling_x_bins[sampling_x_bins < 0] = 0
#         sampling_y_bins[sampling_y_bins < 0] = 0
#         sampling_z_bins[sampling_z_bins < 0] = 0
        
          
        
#         acc_dose_subfraction[origin_point_x_bins, origin_point_y_bins, origin_point_z_bins] += algo.dose[sampling_x_bins, sampling_y_bins, sampling_z_bins]
        
#         # # Sample each
#         # for idx, (sample_point, origin_point) in enumerate(zip(transformed_points, points_of_binary_mask)):
            
#         #     # Sample points
#         #     x_bin = int((sample_point[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0])
#         #     y_bin = int((sample_point[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1])
#         #     z_bin = int((sample_point[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2])
            
#         #     if x_bin >= algo.dose.shape[0]:
#         #         x_bin = algo.dose.shape[0] - 1 
#         #     if y_bin >= algo.dose.shape[1]:
#         #         y_bin = algo.dose.shape[1] - 1 
#         #     if z_bin >= algo.dose.shape[2]:
#         #         z_bin = algo.dose.shape[2] - 1          
                
#         #     if x_bin < 0:
#         #         x_bin = 0
#         #     if y_bin < 0:
#         #         y_bin = 0
#         #     if z_bin < 0:
#         #         z_bin = 0       
            
#         #     # Origin points
#         #     origin_point_x_bin = int((origin_point[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0])
#         #     origin_point_y_bin = int((origin_point[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1])
#         #     origin_point_z_bin = int((origin_point[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2])        
            
            
#         #     assert [x_bin, y_bin, z_bin] == [sampling_x_bins[idx], sampling_y_bins[idx], sampling_z_bins[idx]]
#         #     assert [origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] == [origin_point_x_bins[idx], origin_point_y_bins[idx], origin_point_z_bins[idx]]
            
        
#         #     assert binary_mask[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] == True
            
#         #     try:
#         #         acc_dose_subfraction[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] += algo.dose[x_bin, y_bin, z_bin]
        
#         #     except:
#         #         #Sampling outside grid. Solved by assuming same value as its origin, ie the boundary value
#         #         acc_dose_subfraction[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] += algo.dose[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin]
#         #         # print("Removed", x_bin, y_bin, z_bin)
#         #         pass
            
#         # assert (acc_dose_subfraction == acc_dose_subfraction2).all()
        
#         np.save(path / f'subfraction_{idx_of_interest}', acc_dose_subfraction)
#         tot_acc_dose += acc_dose_subfraction/len(investigated_indices) #fraction1.N_resampled
        
    
#         plot_motion_subfraction_frame(algo, algo_moved_model, ep_model, idx_of_interest, acc_dose_subfraction, slice_idx,  path)
    