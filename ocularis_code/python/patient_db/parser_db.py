# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 16:30:55 2023

@author: bjoerk_c
"""


import patient_surfaces

import json
from copy import deepcopy
# from pathlib import Path

from configuration import config
from ocularis_tps.patient import Patient
import numpy as np


def parser_database(path , patient_id, anterior_aligned = False, wider = False):
    
    patient_path = path / patient_id
    
    with open(patient_path / f'{patient_id}.json') as file:
        data = json.load(file)
        
    # Access the separate dictionaries
    patient_config2 = data['patient_config']
    dose_engine_config = data['dose_engine_config']
    
    
    patient_config = {}
    patient_config['patient_path'] = patient_path
    patient_config['dose horizontal filename'] = patient_config2['dose horizontal filename']
    patient_config['eyeplan ptd filename'] = patient_config2['eyeplan ptd filename']
    patient_config['eyeplan mil filename'] = patient_config2['eyeplan mil filename']
    patient_config['clip transformation'] = patient_config2['clip transformation']
    patient_config['patient_model_path'] = patient_path / "model"
    patient_config['patient_model_clips_path'] = patient_path / "model_clips"
    patient_config['skinplane_to_anterior_distance'] = patient_config2['skinplane_to_anterior_distance']
    
    patient_config['patient_number'] = patient_config2['patient_number'] #patient_id #
    patient_config['fractions'] = patient_config2['fractions']
    
    
    if "motion" in patient_config2.keys():
        patient_config["motion"] = patient_config2["motion"]    
    
    if "basepath" in patient_config2.keys():
        patient_config["basepath"] = patient_config2["basepath"]
        
    if "plot_path" in patient_config2.keys():
        patient_config["plot_path"] = patient_config2["plot_path"]        

    config3 = deepcopy(config)
    
    config3["patient_id"] = patient_id
    config3["patient_number"] = int(patient_id[1:])
    
    if wider:
        config3["Nvoxels"] =  [dose_engine_config["Nvoxels"][0]*2, dose_engine_config["Nvoxels"][1]*2, dose_engine_config["Nvoxels"][2]] 
        print("Initializing a wider mesh")
    else:
        config3["Nvoxels"] = dose_engine_config["Nvoxels"]
    config3["Image"] = dose_engine_config["Image"]
    config3["Slice"] = dose_engine_config["Slice"]
    config3["proximal_margin"] = dose_engine_config["proximal_margin"]
    config3["distal_margin"] = dose_engine_config["distal_margin"]
    
    if "wedge_insert_angle" in dose_engine_config.keys():
        config3["wedge_insert_angle"] = dose_engine_config["wedge_insert_angle"]
    
    if "wedge_cover" in dose_engine_config.keys():
        config3["wedge_cover"] = dose_engine_config["wedge_cover"]
    
    if "wedge_angle" in dose_engine_config.keys():
        config3["wedge_angle"] = dose_engine_config["wedge_angle"]
    
    if "collimator_ruler_length" in dose_engine_config.keys():
        config3["collimator_ruler_length"] = dose_engine_config["collimator_ruler_length"]
    
    if "lid_point" in dose_engine_config.keys():
        config3["lid_point"] = np.asarray(dose_engine_config["lid_point"])
        
    if "lid_thickness" in dose_engine_config.keys():
        config3["lid_thickness"] = dose_engine_config["lid_thickness"]
    
    pat = Patient(patient_config)
    
    ep_model = pat.patient_model.eyeplan_model
    

    # ep_model.config = config3
    target_radiological_depth, target_modulation, deepest_ray, shallowest_ray = ep_model.target_and_modulation_range_of_target
    
    # skin_plane = ep_model.skin_plane_most_anterior - 2
    
    # target_radiological_depth, target_modulation, deepest_ray, shallowest_ray = pat1.patient_model.eyeplan_model.target_and_modulation_range_of_target
    # ep_model = pat1.patient_model.eyeplan_model
    skin_plane_z = ep_model.skin_plane_most_anterior - patient_config['skinplane_to_anterior_distance']

    config3['Target_range'] = target_radiological_depth + ep_model.distal_margin
    config3['Modulation_range'] = target_modulation + ep_model.proximal_margin
    if wider:
        config3["Mesh dimensions"] = [ep_model.doseplane_h.mesh_dimensions[0]*2, ep_model.doseplane_h.mesh_dimensions[1]*2, ep_model.doseplane_h.mesh_dimensions[2] ]
    else:
        config3["Mesh dimensions"] = ep_model.doseplane_h.mesh_dimensions
    config3["skin_plane_normal"] = np.asarray([0, 0, 1])
    config3["skin_plane_point"] = np.asarray([0, 0, skin_plane_z])
    
    
    eyeglobe_points = pat.patient_model.eyeplan_model.structure_set_clips_registered['eyeglobe'].contour_coordinates
    
    config3["eyeglobe_point"] = eyeglobe_points
    config3["target_points"] = pat.patient_model.eyeplan_model.structure_set_clips_registered['target'].contour_coordinates
    config3["skin_plane"] = skin_plane_z
    
    ep_model.surfaces.append(patient_surfaces.SkinPlane(config3["skin_plane_point"], config3["skin_plane_normal"]))
    ep_model.surfaces.append(patient_surfaces.PointCloud(config3["eyeglobe_point"], config["eyeglobe_mesh_triangle_size"]))
    
    
    
    if "lid_point" in dose_engine_config.keys() and "lid_thickness" in dose_engine_config.keys():
        ep_lid = deepcopy(eyeglobe_points)
        ep_lid[0:,2] += 2.5
        ep_lid = np.asarray(list(filter(lambda p: p[1] > config3["lid_point"][1] , ep_lid)))
        ep_lid = np.asarray(list(filter(lambda p: p[2] > skin_plane_z , ep_lid)))

        ep_model.surfaces.append(patient_surfaces.PointCloud(ep_lid, config["eyelid_mesh_triangle_size"]))
        ep_model.eyelid_case = True
        print("Eyelid defined")

    config3["surfaces"] = ep_model.surfaces
    
    assert len(config3["surfaces"]) > 1
    
    if anterior_aligned:
        config3["anterior_aligned"] = True
    else:
        config3["anterior_aligned"] = False
    
    
    #Grid is defined here
    pat.dose_engine_config = config3
    
    pat.patient_model.eyeplan_model.config = config3
    
    return pat
