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



def generate_reference(data, config, ep_model, simple = False):

    
    assert type(simple) == bool
    
    algo = dose_engine.BroadBeam(data, config)
    
    
    if config["anterior_aligned"]:
        # # Aligns mesh with anterior segment of eyeglobe + 0.5 mm
        # z_alignment = np.max(ep_model.structure_set_clips_registered["eyeglobe"].contour_coordinates[0:,2]) + 0.5
        
   
        n_voxels_front_of_skinplane = config["n_voxels_front_of_skinplane"] #int(np.ceil(3/ep_model.doseplane_h.resolution[2]))
        z_alignment = algo.config['skin_plane_point'][2] + n_voxels_front_of_skinplane*ep_model.doseplane_h.resolution[2]
        
        config["n_voxels_front_of_skinplane"] = n_voxels_front_of_skinplane
        
        print(f"Dose_engine: Mesh aligned to {n_voxels_front_of_skinplane*ep_model.doseplane_h.resolution[2]} mm in front of skin plane")
        # print("Dose_engine: Mesh aligned to 2.5 mm in front of skin plane")
    else:
        # Aligns mesh with skin plane
        z_alignment = algo.config["skin_plane_point"][2]
        print("Dose_engine: Mesh aligned to skin plane")
    
    
    translation_vector = np.zeros(3)

    
    translation_vector[2] = z_alignment - algo.medium.maxZ
    translation_vector[1] = ep_model.doseplane_h.resolution[1]/2
    translation_vector[0] = ep_model.doseplane_h.resolution[0]/2
    algo.medium.translate(translation_vector)
    assert - 0.0001 < (algo.medium.maxZ -
                       z_alignment)/z_alignment < 0.0001
    
    if simple:
        print("Generated simplified dose engine instance")
        return algo
    
    if config["anterior_aligned"]:
        struct = ep_model.structure_set_clips_registered["eyeglobe"]
        struct.resample_contour_points_in_grid(1)        
        algo.medium.medium[struct.binary_mask] = 1
        
        struct = ep_model.structure_set_clips_registered["target"]
        struct.resample_contour_points_in_grid(1)        
        algo.medium.medium[struct.binary_mask] = 1
        
        struct = ep_model.structure_set_clips_registered["cornea"]
        struct.resample_contour_points_in_grid(1)        
        algo.medium.medium[struct.binary_mask] = 1
        
    algo.define_target_rays_from_points()
    algo.orthogonal_limit = 15
    algo.define_rays()
    algo.config_surface_interaction()
    algo.determine_target_and_modulation_range_of_target()
    algo.config_wedge()

    algo.config_raytracer()
    algo.config_nozzle()
    algo.config_depth_dose_profile()
    algo.config_collimator()  
    algo.config_distances_to_aperture()
    algo.calc()
    
    print("Generated full dose engine instance")

    return algo



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


def generate_subfraction(algo, ep_model, new_points):

    conf = copy.deepcopy(algo.config)

    assert "target" in new_points.keys()
    assert "eyeglobe" in new_points.keys()

    conf["eyeglobe_point"] = new_points["eyeglobe"]
    conf["target_points"] = new_points["target"]

    alg = dose_engine.BroadBeam(algo.data, conf)
    translation_vector = np.zeros(3)
    # abs(alg.medium.mesh_apex[2]) + 35.6/2 #skin_plane
    translation_vector[2] = (conf["skin_plane"] - alg.medium.maxZ)
    translation_vector[1] = ep_model.doseplane_h.resolution[1]/2
    translation_vector[0] = ep_model.doseplane_h.resolution[0]/2
    #alg.config_nozzle_by_mw(604)
    alg.medium.translate(translation_vector)
    alg.define_target_rays_from_points()
    alg.orthogonal_limit = 15
    alg.define_rays()
    alg.config_surface_interaction()
    # alg.determine_target_and_modulation_range_of_target()
    alg.config_wedge()

    alg.config_raytracer()
    alg.nozzle = copy.deepcopy(algo.nozzle)
    alg.define_geo_relationships()
    alg.nozzle_configured = True
    # alg.config_nozzle()
    alg.config_depth_dose_profile()
    alg.collimator = copy.deepcopy(algo.collimator)
    # alg.config_collimator() #- alg.nozzle.inflection_shift_mm[0]
    alg.config_distances_to_aperture()
    alg.calc()

    return alg


def generate_subfractions_from_motion(algo, ep_model, fraction1, slice_idx, path):
    
    from plot_utils.motion_plot_utils import plot_motion_subfraction_frame
    
    investigated_indices = range(fraction1.N_resampled) #[5, 6, 9, 12, 14, 18]
    
    tot_acc_dose = np.zeros(algo.dose.shape)
    for idx_of_interest in investigated_indices:
        
        
        print("----- idx_of_interest ", idx_of_interest)
        
        acc_dose_subfraction = np.zeros(algo.dose.shape)
        # acc_dose_subfraction2 = np.zeros(algo.dose.shape)
        # acc_dose_not_eye = np.zeros(algo.dose.shape)
        # acc_dose_eye = np.zeros(algo.dose.shape)
        
        algo_moved_model = generate_subfraction(algo, ep_model, fraction1.points_at_idx(idx_of_interest, True))
        
        structure_set = ep_model.generate_structure_set(fraction1.points_at_idx(idx_of_interest, True))
        structure_set.grid = ep_model.grid
        
        struct = structure_set["eyeglobe"]
        struct.resample_contour_points_in_grid(1)
        binary_mask = struct.binary_mask
        negative_mask = binary_mask == False
        points_of_binary_mask = struct.points_of_binary_mask
        transformed_points = fraction1.transform_points_by_idx(points_of_binary_mask, idx_of_interest, True)
            
        
        acc_dose_subfraction[negative_mask] = algo_moved_model.dose[negative_mask]
        # acc_dose_subfraction2[negative_mask] = algo_moved_model.dose[negative_mask]
        
        origin_point_x_bins = np.asarray(np.floor((points_of_binary_mask[0:,0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
        origin_point_y_bins = np.asarray(np.floor((points_of_binary_mask[0:,1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
        origin_point_z_bins = np.asarray(np.floor((points_of_binary_mask[0:,2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
                    
        
        
        sampling_x_bins = np.asarray(np.floor((transformed_points[0:,0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
        sampling_y_bins = np.asarray(np.floor((transformed_points[0:,1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
        sampling_z_bins = np.asarray(np.floor((transformed_points[0:,2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
        
        # In case sampling point is outside the grid, then it sampling instead within the grid
        sampling_x_bins[sampling_x_bins>algo.dose.shape[0]-1] = algo.dose.shape[0] -1
        sampling_y_bins[sampling_y_bins>algo.dose.shape[1]-1] = algo.dose.shape[1] -1
        sampling_z_bins[sampling_z_bins>algo.dose.shape[2]-1] = algo.dose.shape[2] -1
        sampling_x_bins[sampling_x_bins < 0] = 0
        sampling_y_bins[sampling_y_bins < 0] = 0
        sampling_z_bins[sampling_z_bins < 0] = 0
        
          
        
        acc_dose_subfraction[origin_point_x_bins, origin_point_y_bins, origin_point_z_bins] += algo.dose[sampling_x_bins, sampling_y_bins, sampling_z_bins]
        
        # # Sample each
        # for idx, (sample_point, origin_point) in enumerate(zip(transformed_points, points_of_binary_mask)):
            
        #     # Sample points
        #     x_bin = int((sample_point[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0])
        #     y_bin = int((sample_point[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1])
        #     z_bin = int((sample_point[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2])
            
        #     if x_bin >= algo.dose.shape[0]:
        #         x_bin = algo.dose.shape[0] - 1 
        #     if y_bin >= algo.dose.shape[1]:
        #         y_bin = algo.dose.shape[1] - 1 
        #     if z_bin >= algo.dose.shape[2]:
        #         z_bin = algo.dose.shape[2] - 1          
                
        #     if x_bin < 0:
        #         x_bin = 0
        #     if y_bin < 0:
        #         y_bin = 0
        #     if z_bin < 0:
        #         z_bin = 0       
            
        #     # Origin points
        #     origin_point_x_bin = int((origin_point[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0])
        #     origin_point_y_bin = int((origin_point[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1])
        #     origin_point_z_bin = int((origin_point[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2])        
            
            
        #     assert [x_bin, y_bin, z_bin] == [sampling_x_bins[idx], sampling_y_bins[idx], sampling_z_bins[idx]]
        #     assert [origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] == [origin_point_x_bins[idx], origin_point_y_bins[idx], origin_point_z_bins[idx]]
            
        
        #     assert binary_mask[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] == True
            
        #     try:
        #         acc_dose_subfraction[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] += algo.dose[x_bin, y_bin, z_bin]
        
        #     except:
        #         #Sampling outside grid. Solved by assuming same value as its origin, ie the boundary value
        #         acc_dose_subfraction[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin] += algo.dose[origin_point_x_bin, origin_point_y_bin, origin_point_z_bin]
        #         # print("Removed", x_bin, y_bin, z_bin)
        #         pass
            
        # assert (acc_dose_subfraction == acc_dose_subfraction2).all()
        
        np.save(path / f'subfraction_{idx_of_interest}', acc_dose_subfraction)
        tot_acc_dose += acc_dose_subfraction/len(investigated_indices) #fraction1.N_resampled
        
    
        plot_motion_subfraction_frame(algo, algo_moved_model, ep_model, idx_of_interest, acc_dose_subfraction, slice_idx,  path)
    