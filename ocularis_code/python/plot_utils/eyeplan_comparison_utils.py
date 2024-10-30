# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:44:12 2022

@author: bjoerk_c
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import numpy as np
from scipy.interpolate import interp1d
#import matplotlib.patches as mpatches
import reference_frames
import pandas as pd
import copy
from utils import color_utils
import pymedphys
from sklearn import preprocessing
from pathlib import Path
import pathlib
import geometry
import dose_engine

import matplotlib.patches as patches
    

import matplotlib.cm as cmx
import matplotlib.colors as colors

from utils.normalization_utils import find_R90_increasing, find_R10_from_left

from utils import normalization_utils


from utils.normalization_utils import find_lateral_inflections, find_lateral_sharpness, find_80_20_L
from pyproton.metrics.dvh import DVH

from geometry import OptisGeometry


from plot_utils.PlotProjectedAperture import Arrow3D

def find_R90_from_left(xes, yes):
    
    # xs, ys = xes[0:split], yes[0:split]
    
    tolerance = 0.01  # adjust this as needed
    for i in range(len(yes) - 1):
        
        xs = [xes[i], xes[i + 1] ]
        
        ys = [yes[i], yes[i + 1] ]
        
        f_rev = interp1d( ys, xs)
        
        try:
            return f_rev(0.9)
        except:
            pass


def plot_wedge_use(algo, ep_model, path = ""):
    
    
    if algo.config["wedge_angle"] == 0:
        print("No wedge")
        return
    
    
    geo = OptisGeometry(algo.config)
    
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.axis('equal')
        
    plt.plot([-17.5, 17.5], [-17.5, 17.5], color = "k", alpha = 0.5)
    plt.plot([17.5, -17.5], [-17.5, 17.5], color = "k", alpha = 0.5)
    
    eye_points = ep_model.structure_set_clips_registered['eyeglobe'].structure_points
    plt.scatter(eye_points[0:, 0], eye_points[0:, 1], color = "k", alpha = 0.1)
    
    target_points = ep_model.structure_set_clips_registered['target'].structure_points
    plt.scatter(target_points[0:, 0], target_points[0:, 1], color = "g")
    
    plt.scatter(geo.wedge_insert_point[0], geo.wedge_insert_point[1], color = "C0")
    plt.scatter(geo.wedge_apex_point[0], geo.wedge_apex_point[1], color = "C0")
    plt.plot( [geo.wedge_apex_point[0], geo.wedge_insert_point[0]],[geo.wedge_apex_point[1], geo.wedge_insert_point[1]], color = "C0", label = "Wedge insertion")
    
    plt.quiver(geo.wedge_apex_point[0], geo.wedge_apex_point[1], geo.wedge_insert_direction[0], geo.wedge_insert_direction[1], color = "C0")
    
    circle2 = plt.Circle((0.0, 0.0), 17.5, color='k', fill = None, label = "Aperture limit")
    ax.add_patch(circle2)
    
    plt.legend()
    
    plt.xlim(-19,19)
    plt.ylim(-19,19)
    
    
    plt.xlabel("x [mm]", fontsize = 12)
    plt.ylabel("y [mm]", fontsize = 12)
    
    plt.title(f"Wedge position for P{ep_model.patient_config['patient_number']}", fontsize = 16)
    
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
   # plt.show()
  #  if type(path) == pathlib.WindowsPath:
     #   try:
    plt.savefig(path / f"P{ep_model.patient_config['patient_number']}_wedge_use",
                        bbox_inches='tight', pad_inches=0.1)
      #  except:
       #     pass



def one_dimensional_eyeplan_comparison(algo, ep_model, path = ""):
    
    # slice_idx = algo.config['Slice']
    
    
    x_bin,y_bin, z_bin = int(algo.dose.shape[0]/2), int(algo.dose.shape[1]/2), int(algo.dose.shape[2]/2)
    
    
    skin_plane = algo.config["skin_plane"]
    
    
    x_pos = algo.medium.resolution[0] * \
        (x_bin + 0.5) + algo.medium.mesh_origin[0]
    y_pos = algo.medium.resolution[1] * \
        (y_bin + 0.5) + algo.medium.mesh_origin[1]
    z_pos = algo.medium.resolution[2] * \
        (z_bin + 0.5) + algo.medium.mesh_origin[2]
    
    p = np.asarray([x_pos, y_pos, z_pos])
    

    
     
    # one_dimension(algo, ep_model, algo.dose[algo.config['Slice'][0],algo.config['Slice'][1], 0:], ep_model.doseplane_h.image[29,0:] , skin_plane)
    
    
    
    # x, y, z = 35, 35, algo.config['Slice'][2]
    
    # x, y, z = algo.config['Slice'][0], algo.config['Slice'][1], algo.config['Slice'][2]
    
    
    
    # dose = copy.deepcopy(algo.dose)
    
    
    
    dose_engine_vector = algo.dose[x_bin, y_bin, 0:]
    
    
    index_fraction = ep_model.doseplane_h.image.shape[0]/ algo.dose.shape[0]
    
    ep_x_index = int(index_fraction * x_bin)
    
    ep_model_vector = ep_model.doseplane_h.image[ep_x_index,0:]
    
    
    
    p_deepest = algo.deepest_ray.T
    p_deepest_2d = np.asarray([p_deepest[2], p_deepest[0]])
    
    x_bin = np.asarray(np.floor((p_deepest[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_deepest[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    z_bin = np.asarray(np.floor((p_deepest[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
    
    
    
    # x2, y2, z2 = 35, algo.config['Slice'][1], algo.config['Slice'][2]
    
    x2, y2, z3 = x_bin, y_bin, z_bin
    
    dose_engine_vector2 = algo.dose[x2, y2, 0:]
    
    # ep_model_vector2 = ep_model.doseplane_h.image[x2,0:]
    
    
    
    dose_show = copy.deepcopy(algo.dose)
    dose_show[x_bin,y_bin, 0: ] = 1.2
    dose_show[x2,y2, 0: ] = 1.3
    
    
    
    image_doseengine_normal = dose_show[0:,y_bin,0:]
    
    image_doseengine_distal_point = dose_show[0:,y2,0:]
    
    
    
    

    
    
    split = 17
    
    fig = plt.figure()
    
    ax = fig.add_subplot(221)
    ax.axis('equal')
    
    
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, image_doseengine_normal)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    plt.grid()
    
    
    eye_points = ep_model.structure_set_clips_registered['eyeglobe'].structure_points
    plt.scatter(eye_points[0:, 2], eye_points[0:, 0], color = "k", alpha = 0.1)
    
    target_points = ep_model.structure_set_clips_registered['target'].structure_points
    plt.scatter(target_points[0:, 2], target_points[0:, 0], color = "g", alpha = 0.3)
    
    
    plt.xlim(-16.2, skin_plane)
    plt.ylabel("x [mm]", fontsize = 12)
    
    plt.title(f"y = {y_bin}, sampled in the plane of the beam axis", fontsize = 10)
    
    # plt.title("Dose engine", fontsize = 14)
    # ax = plt.subplot(gs[0:70, 0])
    
    ax = fig.add_subplot(223)
    ax.axis('equal')
    
    
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, image_doseengine_distal_point)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    plt.grid()
    
    
    eye_points = ep_model.structure_set_clips_registered['eyeglobe'].structure_points
    plt.scatter(eye_points[0:, 2], eye_points[0:, 0], color = "k", alpha = 0.1)
    
    target_points = ep_model.structure_set_clips_registered['target'].structure_points
    plt.scatter(target_points[0:, 2], target_points[0:, 0], color = "g", alpha = 0.3)
    
    
    plt.xlim(-16.2, skin_plane)
    plt.ylabel("x [mm]", fontsize = 12)
    plt.xlabel("z [mm]", fontsize = 12)
    
    
    plt.title(f"y = {y2}, sampled at the distal target point", fontsize = 10)
    
    
    
    ax = plt.subplot(222)

    xes, yes = algo.medium.z_voxel_pos, dose_engine_vector
    plt.plot(xes, yes, label = f"Dose engine R90 at z = {round(algo.config['skin_plane_point'][2] - algo.config['Target_range'], 3)} mm")
    
    # f_rev = interp1d( dose_engine_vector[0:split], algo.medium.z_voxel_pos[0:split])
    z = find_R90_from_left(xes, yes)
    
    p_de = np.asarray([z, 0.9])
    
    plt.scatter(p_de[0], p_de[1], label = z, color = "C0")
    
    
    dist = abs(p_deepest[2] - z_bin)
    # print(dist)
    
    xes, yes = ep_model.doseplane_h.zs_mean, ep_model_vector
    plt.plot(xes, yes, label = "Eyeplan", linewidth = 2)
    
    # f_rev = interp1d( ep_model_vector[0:split], ep_model.doseplane_h.zs_mean[0:split])
    z2 = find_R90_from_left(xes, yes)
    plt.scatter(z2, 0.9, label = z2, color = "C1")
    
    
    dist2 = abs(p_deepest[2] - z2)
    # print(dist2)
    
    plt.axvline(x=skin_plane, color = "g", linewidth = 4, label = f"Skin plane at z = {round(algo.config['skin_plane_point'][2], 3)}")
    
    
    plt.axvline(x=algo.expected_range, color = "k", label = f"expected range {round(algo.expected_range,3)} mm")
    plt.axvline(x=algo.expected_modulation, color = "k")
    
    R90_from_target = algo.expected_range
    rect = patches.Rectangle((R90_from_target, -10),algo.config["distal_margin"], 50, linewidth=1, edgecolor=None, facecolor='b', alpha = 0.1, label = f"Distal margin of {algo.config['distal_margin']} mm")
    ax.add_patch(rect)
    
    
    # proximal_target_max_point = max(target_points[0:,2])
    modulation_end = algo.shallowest_ray.T[2]
    rect = patches.Rectangle((modulation_end, -10), algo.config["proximal_margin"], 50, linewidth=1, edgecolor=None, facecolor='c', alpha = 0.3, label = f"Proximal margin of {algo.config['proximal_margin']} mm")
    ax.add_patch(rect)
    
    plt.legend()
    
    plt.grid(linewidth = 0.3)
    
    
    
    plt.ylabel("Dose intensity", fontsize = 12)
    
    
    # plt.title(algo.nozzle)
    
    
    plt.xlabel("Z [mm]", fontsize = 12)
    
    ax2 = ax.twinx()
    
    
    eye_points = ep_model.structure_set_clips_registered['eyeglobe'].structure_points
    plt.scatter(eye_points[0:, 2], eye_points[0:, 0], color = "k", alpha = 0.1)
    
    target_points = ep_model.structure_set_clips_registered['target'].structure_points
    plt.scatter(target_points[0:, 2], target_points[0:, 0], color = "g")
    
    plt.plot(algo.deepest_ray.zes, algo.deepest_ray.xes, color = "C3" )
    
    plt.scatter(algo.deepest_ray.surface_interaction_point[2], algo.deepest_ray.surface_interaction_point[0], color = "C3", s = 40)
    
    
    plt.plot(algo.shallowest_ray.zes, algo.shallowest_ray.xes, color = "C4" )
    
    plt.scatter(algo.shallowest_ray.surface_interaction_point[2], algo.shallowest_ray.surface_interaction_point[0], color = "C4", s = 40)
    
    
    
    
    plt.xlim(-25,20)
    
    
    plt.xlabel("z [cm]", fontsize = 12)
    plt.ylabel("x [cm]", fontsize = 12)
    
    
    
    
    ax = plt.subplot(224)

    
    # fig = plt.figure()
    
    xes, yes = algo.medium.z_voxel_pos, dose_engine_vector2
    plt.plot(xes, yes, label = f"Dose engine R90 at z = {round(algo.config['skin_plane_point'][2] - algo.config['Target_range'], 3)} mm")
    
    z = find_R90_from_left(xes, yes)
    plt.scatter(z, 0.9, label = z, color = "C0")
    
    
    dist = abs(p_deepest[2] - z_bin)
    # print(dist)
    
    # xes, yes = ep_model.doseplane_h.zs_mean, ep_model_vector2
    # plt.plot(xes, yes, label = "Eyeplan", linewidth = 2)
    
    # z2 = find_R90_from_left(xes, yes)
    # plt.scatter(z2, 0.9, label = z2, color = "C1")
    
    # dist2 = abs(p_deepest[2] - z2)
    # print(dist2)
    
    # plt.axvline(x=skin_plane, color = "g", linewidth = 4, label = f"Skin plane at z = {round(algo.config['skin_plane_point'][2], 3)}")
    # plt.axvline(x=algo.expected_range, color = "k", label = f"expected range {round(algo.expected_range,3)} mm")
    # plt.axvline(x=algo.expected_modulation, color = "k")
    
    R90_from_target = algo.expected_range
    rect = patches.Rectangle((R90_from_target, -10),algo.config["distal_margin"], 50, linewidth=1, edgecolor=None, facecolor='b', alpha = 0.1, label = f"Distal margin of {algo.config['distal_margin']} mm")
    ax.add_patch(rect)
    
    
    # proximal_target_max_point = max(target_points[0:,2])
    modulation_end = algo.shallowest_ray.T[2]
    rect = patches.Rectangle((modulation_end, -10), algo.config["proximal_margin"], 50, linewidth=1, edgecolor=None, facecolor='c', alpha = 0.3, label = f"Proximal margin of {algo.config['proximal_margin']} mm")
    ax.add_patch(rect)
    
    plt.legend(loc='center left')
    
    plt.grid(linewidth = 0.3)
    
    
    plt.title(f"R90 diff = {round(z2 - z, 3)}")
    
    
    
    plt.ylabel("Dose intensity", fontsize = 12)
    
    
    
    
    
    plt.xlabel("Z [mm]", fontsize = 12)
    
    plt.ylim(0.05,1.05)
    
    
    ax2 = ax.twinx()
    
    
    eye_points = ep_model.structure_set_clips_registered['eyeglobe'].structure_points
    plt.scatter(eye_points[0:, 2], eye_points[0:, 0], color = "k", alpha = 0.1)
    
    target_points = ep_model.structure_set_clips_registered['target'].structure_points
    plt.scatter(target_points[0:, 2], target_points[0:, 0], color = "g")
    
    plt.plot(algo.deepest_ray.zes, algo.deepest_ray.xes, color = "C3" )
    
    plt.scatter(algo.deepest_ray.surface_interaction_point[2], algo.deepest_ray.surface_interaction_point[0], color = "C3", s = 40)
    
    
    plt.plot(algo.shallowest_ray.zes, algo.shallowest_ray.xes, color = "C4" )
    
    plt.scatter(algo.shallowest_ray.surface_interaction_point[2], algo.shallowest_ray.surface_interaction_point[0], color = "C4", s = 40)
    
    
    
    
    plt.xlim(-25,20)
    
    
    plt.xlabel("z [cm]", fontsize = 12)
    plt.ylabel("x [cm]", fontsize = 12)
    
    
    
    
    plt.suptitle(ep_model.patient_config["patient_number"] +" | " + algo.nozzle.title)
    
    
    xlength = 14
    fig.set_size_inches(xlength, xlength/1.61)
    plt.show()
    
    if type(path) == pathlib.WindowsPath:
    
        try:
            plt.savefig(Path(path) / f"P{ep_model.tps_config['patient_number']}_1d_comp.jpeg",
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
    
    
    
    


def collimator_comparison(algo, ep_model, path = ""):
    # ------ Aperture plane ------------------
    intersects_colli = []
    intersects_vert_colli = []
    for edge in algo.collimator.expanded_hull.boundary_edges:
        xes = edge[0]
        yes = edge[1]
        zes = edge[2]
        
        f = interp1d(xes, yes)
        f_rev = interp1d(yes, xes)
        try:
            intersects_colli.append([f_rev(0), 0])
        except:
            pass
        try:
            intersects_vert_colli.append([0, f(0)])
        except:
            pass    
    
    
    ref = reference_frames.Reference_Frame(algo.config, 'TreatmentRoom')
    geo = geometry.OptisGeometry(algo.config)


    
    p_on_axis = geo.aperture_point + algo.central_axis.vec_norm*70
            
    
    iso_projection = algo.collimator.project_aperture(p_on_axis)
    
    
    Nplanes = 10
    
    colors = cm.rainbow(np.linspace(0, 1, Nplanes))
    
    dist_sifts = np.linspace(-20, 90, Nplanes)
    
    
    my_plane_hulls = []
    
    
    for val in dist_sifts:
        p_on_axis = geo.aperture_point + algo.central_axis.vec_norm*(70 + val)
        projection = algo.collimator.project_aperture(p_on_axis)
        my_plane_hulls.append(projection)
    


    exclude_keys = ['surfaces']
    new_config = {k: algo.config[k] for k in set(list(algo.config.keys())) - set(exclude_keys)}
    config = copy.deepcopy(new_config)
    config['surfaces'] = algo.config['surfaces']
    

    collimator = dose_engine.Collimator(config)
    collimator.add_ray(algo._target_rays)
    collimator.define_aperture_hull()
    collimator.ruler_length = config["collimator_ruler_length"]
    collimator.expand_aperture(algo.config['collimator_expansion'])


    
    fig = plt.figure()
    
    ax = plt.subplot(121)
    
    
    plt.plot(algo.collimator.aperture_hull.boundary_in_plane_xes, algo.collimator.aperture_hull.boundary_in_plane_yes, label = "Projected target shape", linewidth = 2)

    plt.plot(collimator.expanded_hull.boundary_in_plane_xes, collimator.expanded_hull.boundary_in_plane_yes, color = "k", linewidth = 4, label = f"Physical aperture | {round(collimator.expansion,5)} mm expansion") 
    
    
    plt.plot(algo.collimator.expanded_hull.boundary_in_plane_xes, algo.collimator.expanded_hull.boundary_in_plane_yes, label = f"Apparent aperture | {round(algo.collimator.expansion,5)} mm expansion", linewidth = 3)
    

    
    
    try:
        plt.plot(ep_model.collimator_points_in_plane[0:,0], ep_model.collimator_points_in_plane[0:,1], label = "Eyeplan aperture", color = "C2")
        # plt.scatter(ep_model.collimator_points_in_plane[0:,0], ep_model.collimator_points_in_plane[0:,1], color = "C2")
    except:
        print("No mil file?")
    
    plt.scatter(intersects_colli[1][0], intersects_colli[1][1], color = "k")
    plt.scatter(intersects_colli[0][0], intersects_colli[0][1], color = "k")
    
    plt.grid(linewidth = 0.3)
    plt.legend()
    
    plt.xlabel("x [mm]", fontsize = 18)
    plt.ylabel("y [mm]", fontsize = 18)
    

    plt.title("Aperture plane", fontsize = 12)
    

    
    
    
    # fig = plt.figure()
    
    #------------------------------------------
    ax1 = fig.add_subplot(122, projection='3d')
    ax1.view_init(-70,90) 
    plt.title("Treatment Room", fontsize = 12, y=1.15)
    
    ref.add_frame_to_plot(ax1, ref, 10)
    
    
    arrow_prop_dict = dict(mutation_scale=20, arrowstyle='->', shrinkA=0, shrinkB=0)
    
    
    
    # xes = self.algo.collimator.aperture_hull.boundary.boundary.xy[0]
    # yes = self.algo.collimator.aperture_hull.boundary.boundary.xy[1]  
    
    # for i in range(len(colors)):
    #     plt.scatter(xes[i], yes[i], 70, color = colors[i])
    #     ray = self.algo.collimator.aperture_rays[i]
    #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color=colors[i])
    #     ax1.add_artist(a)    
    
    plt.plot([], [], [], color = "k", linestyle = "-", label = "Projected aperture iso plane", linewidth = 1)
    for edge in iso_projection.boundary_edges:
        plt.plot(edge[0], edge[1], edge[2], color = "k", linestyle = "-", linewidth = 1)
    
    
    
    
    plt.plot([],[],[], label = "Apparent aperture", color = "C1", linestyle = "--")
    for edge in algo.collimator.expanded_hull.boundary_edges:
        xes = edge[0]
        yes = edge[1]
        zes = edge[2]
        plt.plot(xes, yes, zes, color = "C1", linestyle = "--")
        
    
    
    
    
    # for i in range(len(colors)):
    #     p = algo.collimator.point_ray_pair_expanded_aperture[i][0]
    #     plt.scatter(p[0], p[1], 70, color = colors[i])
    
        
    # for i in algo.collimator.aperture_hull.boundary_vertices:
        
    #     if i >= len(colors):
    #         break
        
    #     ray = algo.collimator.aperture_rays[i]
    #     # ray = algo.collimator.point_ray_pair_expanded_aperture[i][1]
    
    #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color=colors[i])
    #     ax1.add_artist(a)
    
    
    
    
    
    
        
    plt.plot([],[],[], label = "Projected target shape", color = "C0", linestyle = "--")
    for edge in algo.collimator.aperture_hull.boundary_edges:
        xes = edge[0]
        yes = edge[1]
        zes = edge[2]
        plt.plot(xes, yes, zes, color = "C0", linestyle = "--")    
        
    # plt.plot(algo.central_axis.xes, algo.central_axis.yes, algo.central_axis.zes,  color = "k", label = "Central axis")
    
    
    plt.plot([],[],[], color = "k", label = "Central axis")
    a = Arrow3D(algo.central_axis.xes, algo.central_axis.yes, algo.central_axis.zes , **arrow_prop_dict, color="k")
    ax1.add_artist(a)
    
    
    
    # # plt.plot([],[],[], color = "C0", label = "Aperture rays")
    # for ray in self.algo.collimator.aperture_rays:
    #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color="C0", label = "Aperture rays")
    #     ax1.add_artist(a)
    
    
    # # plt.plot([],[],[], color = "C1", label = "Rays of expanded aperture")
    # for ray in self.algo.collimator.aperture_rays_expanded:
    #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color="C1", label = "Rays of expanded aperture")
    #     ax1.add_artist(a)
    
    
    # for idx in range(len(colors)):
    #     ray = self.algo.collimator.aperture_rays[idx]
    #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color=colors[idx], label = "Aperture rays")
    #     ax1.add_artist(a)
    
    
    #     ray = self.algo.collimator.aperture_rays_expanded[idx]
    #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color=colors[idx], label = "Rays of expanded aperture")
    #     ax1.add_artist(a)        
    
    
    
    
    
    for idx, projection, color in zip(np.linspace(0,len(dist_sifts)-1, len(dist_sifts)), my_plane_hulls, colors):
        # plt.plot([],[],[], label = f"{dist_sifts[int(idx)]} mm")
        for edge in projection.boundary_edges:
            plt.plot(edge[0], edge[1], edge[2], color = color, linestyle = "-", linewidth = 1)


    
    plt.legend()
    
    ax1.text(ref.origin_ref[0], ref.origin_ref[1], ref.origin_ref[2] -0.1, r'$0_{iso}$')
    
    
    # if (self.algo.central_axis.vec_norm == [0,0,-1]).all():
    #     ax1.set_xlim(-10, 10)            
    #     ax1.set_ylim(-10, 10)
    # else:
    ax1.set_xlim(-20, 20)            
    ax1.set_ylim(-20, 20)
        
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_zlim(-10, 80)
    
    

    
    # plt.suptitle(f"Collimator definition comparison | ruler {algo.config['collimator_ruler_length']} | P{ep_model.patient_config['patient_number']}", fontsize = 16)
    plt.suptitle(f"Collimator definition comparison", fontsize = 16)
        
    
    
    xlength = 12
    fig.set_size_inches(xlength, xlength/2)
   # plt.show()
    #plt.show()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f"P{ep_model.patient_config['patient_number']}_collimator_comparison",
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass





def compare_data(algo, ep_model):
   
     if ep_model.symmetrical_eye == True:
         sym_eye = 'True'
     else:
         sym_eye = 'False'
         
     struct = ep_model.structure_set_clips_registered["target"]
     target_dose = algo.dose[struct.binary_mask]
     target_max = np.max(target_dose)
     target_min = np.min(target_dose)
     target_homo = target_max/target_min
     
       
     out = {}
     out["Eyeplan target range"] = ep_model.target_range_from_ptd
     out["Dose engine target range"] = algo.nozzle.target_range
     out["Eyeplan modulation range"] = ep_model.modulation_range_from_ptd
     out["Dose engine modulation range"] = algo.nozzle.modulation
     out["Symmetrical eye"] = sym_eye
     out["Target homogeneity"] = target_homo
     out["Max dose in target"] = target_max
     out["Min dose in target"] = target_min
     out["Dose engine proximal margin"] = algo.config['proximal_margin']
     out["Dose engine distal margin"] = algo.config['distal_margin']
     out["Eyeplan proximal margin"] = ep_model.proximal_margin
     out["Eyeplan distal margin"] = ep_model.distal_margin
    
     return out      
   
    
def compare_lateral_profile_at_iso(algo, ep_model, path = ""):

    p_iso = np.asarray([0,0,0])
    x_bin = np.asarray(np.floor((p_iso[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_iso[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    z_bin = np.asarray(np.floor((p_iso[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
 
    eyeplan_image = ep_model.doseplane_h.interpolate_from(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos)
    de_lateralsharpness = find_lateral_sharpness(algo.medium.x_voxel_pos,algo.dose[0:,y_bin,z_bin])
    ep_lateralsharpness = find_lateral_sharpness(ep_model.doseplane_h.x_voxel_pos,eyeplan_image[0:,z_bin])
    de_deltaxL = de_lateralsharpness[1] - de_lateralsharpness[0]
    de_deltaxR = de_lateralsharpness[3] - de_lateralsharpness[2]
    ep_deltaxL = ep_lateralsharpness[1] - ep_lateralsharpness[0]
    ep_deltaxR = ep_lateralsharpness[3] - ep_lateralsharpness[2]
    compare_20L = de_lateralsharpness[0] - ep_lateralsharpness[0]
    compare_80L = de_lateralsharpness[1] - ep_lateralsharpness[1]
    compare_20R = de_lateralsharpness[2] - ep_lateralsharpness[2]
    compare_80R = de_lateralsharpness[3] - ep_lateralsharpness[3]
    
    de_lateralinflections = find_lateral_inflections(algo.medium.x_voxel_pos,algo.dose[0:,y_bin,z_bin])
    ep_lateralinflections = find_lateral_inflections(ep_model.doseplane_h.x_voxel_pos,eyeplan_image[0:,z_bin])
    compare_50L = de_lateralinflections[0] - ep_lateralinflections[0]
    compare_50R = de_lateralinflections[1] - ep_lateralinflections[1]
    de_fieldsize = de_lateralinflections[1] - de_lateralinflections[0]
    ep_fieldsize = ep_lateralinflections[1] - ep_lateralinflections[0]
    
    comparison_deltaxL = de_deltaxL - ep_deltaxL
    comparison_deltaxR = de_deltaxR - ep_deltaxR
    comparison_fieldsize = de_fieldsize - ep_fieldsize
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ep_model.doseplane_h.x_voxel_pos,eyeplan_image[0:,z_bin],color='C0', label = "Eyeplan model")
    ax.plot(algo.medium.x_voxel_pos,algo.dose[0:,y_bin,z_bin],color='C1', label = "Dose engine")
    
    ax.scatter(ep_lateralsharpness[0], 0.2,color='C0',label= f'LP 80-20% (left) = {round(ep_deltaxL,3)} mm')
    ax.scatter(ep_lateralsharpness[1], 0.8,color='C0')
    ax.scatter(ep_lateralsharpness[2], 0.8,color='C0')
    ax.scatter(ep_lateralsharpness[3], 0.2,color='C0',label= f'LP 80-20% (right) = {round(ep_deltaxR,3)} mm')
    ax.scatter(ep_lateralinflections[0], 0.5,color='C0')
    ax.scatter(ep_lateralinflections[1], 0.5,color='C0',label= f'Field Size = {round(ep_fieldsize,3)} mm')
    ax.scatter(de_lateralsharpness[0], 0.2,color='C1')
    ax.scatter(de_lateralsharpness[1], 0.8,color='C1',label= f'LP 80-20% (left) = {round(de_deltaxL,3)} mm')
    ax.scatter(de_lateralsharpness[2], 0.8,color='C1',label= f'LP 80-20% (right) = {round(de_deltaxR,3)} mm')
    ax.scatter(de_lateralsharpness[3], 0.2,color='C1')
    ax.scatter(de_lateralinflections[0], 0.5,color='C1')
    ax.scatter(de_lateralinflections[1], 0.5,color='C1',label= f'Field Size = {round(de_fieldsize,3)} mm')
    plt.suptitle(f'{algo.nozzle.title}',size='small')
    plt.legend()
    plt.xlabel('x [mm]')
    plt.ylabel('Dose [fraction of prescribed dose]')
    plt.title(f'Dose in x direction through isocentre for patient {ep_model.patient_config["patient_number"]}')
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    #plt.show()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_lateral_profile_at_iso',
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass

    
    out = {'Eyeplan 80-20% distance (left)':ep_deltaxL,'Eyeplan 80-20% distance (right)':ep_deltaxR,'Eyeplan distance between 50% doses':ep_fieldsize,'Dose engine 80-20% distance (left)':de_deltaxL,'Dose engine 80-20% distance (right)':de_deltaxR,'Dose engine distance between 50% doses':de_fieldsize,'x_bin':x_bin,'y_bin':y_bin,'Distance between 20% doses (left)(dose engine - eyeplan)':compare_20L,'Distance between 80% doses (left)(dose engine - eyeplan)':compare_80L,'Distance between 20% doses (right)(dose engine - eyeplan)':compare_20R,'Distance between 80% doses (right)(dose engine - eyeplan)':compare_80R,'Distance between 50% doses (left)(dose engine - eyeplan)':compare_50L,'Distance between 50% doses (right)(dose engine - eyeplan)':compare_50R}        
    
    return out
    

def compare_depth_profile_at_iso(algo, ep_model, path = ""):
  
    p_iso = np.asarray([0,0,0])
    x_bin = np.asarray(np.floor((p_iso[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_iso[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    z_bin = np.asarray(np.floor((p_iso[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
    range_data = compare_data(algo,ep_model)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:],color='C0', label = "Eyeplan model")
    ax.plot(algo.medium.z_voxel_pos,algo.dose[x_bin,y_bin,0:],color='C1', label = "Dose engine")
    ep_R90 = normalization_utils.find_R90_from_left(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    de_R90 = normalization_utils.find_R90_from_left(algo.medium.z_voxel_pos,algo.dose[x_bin,y_bin,0:])

    try:
        ep_R90_R = normalization_utils.find_R90_from_right(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    except:
        ep_R90_R = algo.config['skin_plane_point'][2]
    if ep_R90_R > algo.config['skin_plane_point'][2]:
        ep_R90_R = algo.config['skin_plane_point'][2]
    
    try:
        de_R90_R = normalization_utils.find_R90_from_right(algo.medium.z_voxel_pos,algo.dose[x_bin,y_bin,0:])
    except:
        de_R90_R = algo.config['skin_plane_point'][2]
        
    ep_R10 = normalization_utils.find_R10_from_left(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    de_R10 = normalization_utils.find_R10_from_left(algo.medium.z_voxel_pos,algo.dose[x_bin,y_bin,0:])
    ep_tr = range_data['Eyeplan target range']
    ep_mr = range_data['Eyeplan modulation range']
    de_tr = range_data['Dose engine target range']
    de_mr = range_data['Dose engine modulation range']
    ep_distal_fall_off = ep_R90 - ep_R10
    de_distal_fall_off = float(de_R90) - float(de_R10)
    ep_sobp_range = ep_R90_R - ep_R90
    de_sobp_range = de_R90_R - de_R90
    ep_R100L = ep_R90 + (ep_sobp_range/10) 
    ep_R100R = ep_R90_R - (ep_sobp_range/10)
    de_R100L = de_R90 + (de_sobp_range/10)
    de_R100R = de_R90_R - (de_sobp_range/10)
    compare_R90 = de_R90 - ep_R90
    
    ep_pos_list = ep_model.doseplane_h.z_voxel_pos.tolist()
    ep_list = ep_model.doseplane_h.image[x_bin,0:].tolist()
    ep_R100L_index = min(ep_pos_list, key = lambda x : abs(ep_R100L - x))
    ep_R100R_index = min(ep_pos_list, key = lambda x : abs(ep_R100R - x))
    
    ep_flat_range = ep_list[ep_pos_list.index(ep_R100L_index):ep_pos_list.index(ep_R100R_index)]
    ep_max = max(ep_flat_range)
    ep_min = min(ep_flat_range)
    ep_homogeneity = ep_max/ep_min
    
    de_pos_list = algo.medium.z_voxel_pos.tolist()
    de_list = algo.dose[x_bin,y_bin,0:].tolist()
    de_R100L_index = min(de_pos_list, key = lambda x : abs(de_R100L - x))
    de_R100R_index = min(de_pos_list, key = lambda x : abs(de_R100R - x))
    
    de_flat_range = de_list[de_pos_list.index(de_R100L_index):de_pos_list.index(de_R100R_index)]
    de_max = max(de_flat_range)
    de_min = min(de_flat_range)
    de_homogeneity = de_max/de_min
    
    struct = ep_model.structure_set_clips_registered['target']
    struct.resample_contour_points_in_grid(1)
    image = struct.binary_mask[x_bin,y_bin,0:]
    image_list = image.tolist()
    res = []
    i = 0
    while True:
        try:
            i = image_list.index(True, i)
            res.append(i)
            i += 1
        except ValueError:
            break
    width = (np.max(res)) - (np.min(res))
    width = width*algo.medium.resolution[2]
   # image_list = image.tolist()
   # index = image_list.index(image==True)
    x = algo.medium.z_voxel_pos[np.min(res)]
    xy = (x, 0.8)
    rec = plt.Rectangle(xy, width,0.2,alpha=0.2,label='Target')
    ax.add_patch(rec)
    
    ax.scatter(ep_R90,0.9,color='C0',label=f'Target Range = {round(ep_tr,3)} mm Modulation Range = {round(ep_mr,3)} mm')
    ax.scatter(de_R90,0.9,color='C1',label=f'Target Range = {round(de_tr,3)} mm Modulation Range = {round(de_mr,3)} mm')
    ax.scatter(ep_R10,0.1,color='C0',label=f'Distal Fall-off = {round(ep_distal_fall_off,3)} mm')
    ax.scatter(de_R10,0.1,color='C1',label=f'Distal Fall-off = {round(de_distal_fall_off,3)} mm')
    ax.scatter(ep_R90_R,0.9,color='C0',label=f'Distance between proximal and distal R90 = {round(ep_sobp_range,3)} mm')
    ax.scatter(de_R90_R,0.9,color='C1',label=f'Distance between proximal and distal R90 = {round(de_sobp_range,3)} mm')
    ax.axvline(x=algo.config['skin_plane_point'][2], color = "g", linewidth = 4)
    plt.suptitle(f'{algo.nozzle.title}',size='small')
    plt.legend()
    plt.xlabel('z [mm]')
    plt.ylabel('Dose [fraction of prescribed dose]')
    plt.title(f'Dose in z direction through isocentre for patient {ep_model.patient_config["patient_number"]}')
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    #plt.show()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_depth_profile_at_iso',
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
    
    out = {'Eyeplan distal fall-off':ep_distal_fall_off, 'Dose engine distal fall-off':de_distal_fall_off, 'Eyeplan field size':ep_sobp_range, 'Dose engine field size':de_sobp_range, 'Eyeplan field homogeneity':ep_homogeneity, 'Dose engine field homogeneity':de_homogeneity, 'Dose engine max dose in field':de_max, 'Dose engine min dose in field':de_min, 'Difference between 90% doses':compare_R90}

    return out
  
def compare_horizontal_plane(algo, ep_model, path = ""):
    
    p_iso = np.asarray([0,0,0])
    x_bin = np.asarray(np.floor((p_iso[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_iso[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    z_bin = np.asarray(np.floor((p_iso[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
    
    fig = plt.figure()
    ax = fig.add_subplot(221)
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, algo.dose[0:,y_bin,0:])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    ax2 = ax.twinx()
    
    ax2.plot(algo.medium.z_voxel_pos,algo.dose[x_bin,y_bin,0:],color='C1', label = "Dose engine")
    ax2.axvline(x=algo.config['skin_plane_point'][2], color = "g", linewidth = 4)
    

    for name in ['eyeglobe', 'target']: #ep_model.structure_set.name_list:
        if name == "Clips":
            continue
        struct = ep_model.structure_set_clips_registered[name]
        
        points = struct.structure_points

        
        alpha = 0.2
        color = "k"
        if name == 'target':
            alpha = 1
            color = "C2"
        
    
        ax.scatter(points[0:,2], points[0:,0],  label = name, alpha = alpha, color = color)

    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax2.set_ylabel('Dose [fraction of prescribed dose]')
    ax.set_title('Dose engine')
    
    ax3 = fig.add_subplot(222)
    eyeplan_image = ep_model.doseplane_h.interpolate_from(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos)
    plt.pcolor(algo.medium.Z, algo.medium.X, eyeplan_image[0:,0:])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    ax4 = ax3.twinx()
    
    ax4.plot(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:],color='C1', label = "Dose engine")
    ax4.axvline(x=algo.config['skin_plane_point'][2], color = "g", linewidth = 4)


    for name in ['eyeglobe', 'target']: #ep_model.structure_set.name_list:
        if name == "Clips":
            continue
        struct = ep_model.structure_set_clips_registered[name]
        
        points = struct.structure_points

        
        alpha = 0.2
        color = "k"
        if name == 'target':
            alpha = 1
            color = "C2"
        
    
        ax3.scatter(points[0:,2], points[0:,0],  label = name, alpha = alpha, color = color)

    
    ax3.set_xlabel('z [mm]')
    ax3.set_ylabel('x [mm]')
    ax4.set_ylabel('Dose [fraction of prescribed dose]')
    ax3.set_title('Eyeplan')
    
    ax5 = fig.add_subplot(223)
    diff_image = algo.dose[0:,y_bin,0:] - eyeplan_image[0:,0:]
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])
    plt.pcolor(algo.medium.Z, algo.medium.X, diff_image, cmap = "seismic", vmin = -limit, vmax = limit )
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose Difference', rotation=270, size='large')
    ax6 = ax5.twinx()
    
    ax6.plot(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:],color='C0', label = "Eyeplan model")
    ax6.plot(algo.medium.z_voxel_pos,algo.dose[x_bin,y_bin,0:],color='C1', label = "Dose engine")
    ax6.axvline(x=algo.config['skin_plane_point'][2], color = "g", linewidth = 4)
    
    for name in ['eyeglobe', 'target']: #ep_model.structure_set.name_list:
        if name == "Clips":
            continue
        struct = ep_model.structure_set_clips_registered[name]
        
        points = struct.structure_points

        
        alpha = 0.2
        color = "k"
        if name == 'target':
            alpha = 1
            color = "C2"
        
    
        ax5.scatter(points[0:,2], points[0:,0],  label = name, alpha = alpha, color = color)
    
    plt.suptitle(f'{algo.nozzle.title}',size='small')
    plt.suptitle(f'Horizontal plane through isocentre for patient {ep_model.patient_config["patient_number"]}')
    ax5.set_xlabel('z [mm]')
    ax5.set_ylabel('x [mm]')
    ax6.set_ylabel('Dose [fraction of prescribed dose]')
    ax5.set_title('Dose difference')
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    plt.tight_layout()
   # plt.show()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_horizontal_plane',
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
 
    
def doseplane_structures(algo, ep_model, path = ""):
    p_iso = np.asarray([0,0,0])
    x_bin = np.asarray(np.floor((p_iso[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_iso[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    z_bin = np.asarray(np.floor((p_iso[2] - algo.medium.mesh_origin[2])/algo.medium.resolution[2]), dtype = 'int')
    
    fig = plt.figure()

    ax = fig.add_subplot(221)  
    ax.axis('equal')
    
    
    legend_handles = []
    legend_labels = []
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        image = struct.binary_mask[0:,y_bin,0:]
            
        plt.contour(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, image , [4], linewidths = 3, colors=color, label=structure_name.capitalize())
        legend_handles.append(plt.Line2D([0], [0], color=color))
        legend_labels.append(structure_name.capitalize())
        
    plt.legend(legend_handles, legend_labels,loc='lower left',fontsize='x-small')
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, algo.dose[0:,y_bin,0:])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose [fraction of prescribed dose]', rotation=270, size='large', verticalalignment='baseline')
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax.set_title('Horizontal Doseplane')    
    
    
    plt.ylim(-20,20)
    
    ax = fig.add_subplot(222)
    ax.axis('equal')
   
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        image = struct.binary_mask[0:,0:,z_bin]
            
        plt.contour(algo.medium.x_voxel_pos, algo.medium.y_voxel_pos, image , [4], linewidths = 3, colors=color, label=structure_name.capitalize())

    plt.legend(legend_handles, legend_labels,loc='lower left',fontsize='x-small')
    plt.pcolor(algo.medium.x_voxel_pos, algo.medium.y_voxel_pos, algo.dose[0:,0:,z_bin])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose [fraction of prescribed dose]', rotation=270, size='large', verticalalignment='baseline')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_title('Transversal Doseplane') 
    
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    
    ax = fig.add_subplot(223)
    ax.axis('equal')
   
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        image = struct.binary_mask[x_bin, 0:,0:]
            
        plt.contour(algo.medium.z_voxel_pos, algo.medium.y_voxel_pos, image, [4], linewidths = 3, colors=color, label=structure_name.capitalize())

    plt.legend(legend_handles, legend_labels,loc='lower left',fontsize='x-small')
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.y_voxel_pos, algo.dose[x_bin,0:,0:])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose [fraction of prescribed dose]', rotation=270, size='large', verticalalignment='baseline')
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_title('Vertical Doseplane')
    
    plt.ylim(-20,20)
    
    
    ax = fig.add_subplot(224)
    
    dose_fractions = np.linspace(0, 1, 1000)
    
    # legend_handles = [plt.Line2D([0], [0], color='k'),plt.Line2D([0], [0], color='k', linestyle='--')]
    # legend_labels = ['Dose Engine','Eyeplan']
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
            
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        dvh = DVH(algo.dose, struct.binary_mask)
        volumes = dvh.V(dose_fractions)    
        
        plt.plot(dose_fractions*100, volumes*100,
                  label=structure_name.capitalize(), color=color)
        
        # legend_handles.append(plt.Line2D([0], [0], color=color))
        # legend_labels.append(structure_name.capitalize())


      
    plt.legend( loc='lower left', fontsize='x-small')
    ax.set_xlabel('Dose [% of prescribed dose]')
    ax.set_ylabel('Volume [% of structure volume]')
    ax.set_title('DVH', size='x-large')
    
    
    plt.suptitle(f'Dose for P{ep_model.patient_config["patient_number"]}')
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    plt.tight_layout()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_doseplane_structures',
                        bbox_inches='tight', pad_inches=0.1)

        except:
            pass


def compare_dvh(algo,ep_model,path=''):
    p_iso = np.asarray([0,0,0])
    y_bin = np.asarray(np.floor((p_iso[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    
    fig = plt.figure()
    ax = fig.add_subplot(221)
    legend_handles = [plt.Line2D([0], [0], color='k'),plt.Line2D([0], [0], color='k', linestyle='--')]
    legend_labels = ['Dose Engine','Eyeplan']
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
            
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        dose_2d = np.expand_dims(algo.dose[0:,y_bin,0:], axis=1)
        struct_2d = struct.binary_mask[:,y_bin,:]
        struct_2d = np.expand_dims(struct_2d, axis=1)
        dvh = DVH(dose_2d, struct_2d)

        dose_fractions = np.linspace(0, 1, 1000)
            
        volumes = dvh.V(dose_fractions)    
        
        plt.plot(dose_fractions*100, volumes*100,
                  label=structure_name.capitalize(), color=color)

        eyeplan_image = ep_model.doseplane_h.interpolate_from(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos)
        eyeplan_im = np.expand_dims(eyeplan_image, axis=1)
        ep_struct = struct.binary_mask[:,y_bin,:]
        ep_struct = np.expand_dims(ep_struct, axis=1)
        dvh_e = DVH(eyeplan_im,ep_struct)
        volumes_e = dvh_e.V(dose_fractions)
        plt.plot(dose_fractions*100, volumes_e*100,
                      color=color, linestyle='--', label='Eyeplan')
        
        legend_handles.append(plt.Line2D([0], [0], color=color))
        legend_labels.append(structure_name.capitalize())


      
    plt.legend(legend_handles, legend_labels, loc='lower left', fontsize='x-small')
    ax.set_xlabel('Dose [% of prescribed dose]')
    ax.set_ylabel('Volume [% of structure volume]')
    ax.set_title('DVH', size='x-large')
    
    ax = fig.add_subplot(222)
    
    legend_handles = []
    legend_labels = []
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        image = struct.binary_mask[0:,y_bin,0:]
            
        plt.contour(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, image , [4], linewidths = 3, colors=color, label=structure_name.capitalize())
        legend_handles.append(plt.Line2D([0], [0], color=color))
        legend_labels.append(structure_name.capitalize())
    plt.legend(legend_handles, legend_labels,loc='lower left',fontsize='x-small')
    
    my_rays = list(filter(lambda ray: abs(ray.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--", alpha=0)
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--",alpha=0)
        
    eyeplan_image = ep_model.doseplane_h.interpolate_from(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos)
    diff_image = algo.dose[0:,y_bin,0:] - eyeplan_image[0:,0:]
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])

    plt.pcolor(algo.medium.Z, algo.medium.X, diff_image, cmap = "seismic", vmin = -limit, vmax = limit)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose Difference [fraction of prescribed dose]', rotation=270, size='large', verticalalignment='baseline')
    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    ax.set_xlim(-20,20)
    ax.set_ylim(-15,15)
    ax.set_xlabel('z [mm]', size='x-large')
    ax.set_ylabel('x [mm]', size='x-large')
    ax.set_title('Horizontal Dose Difference')
    
    ax = fig.add_subplot(223)
    
    diff_volume = diff_image
    diff_vector = diff_volume.ravel()

    N_bins = 50

    counts, bins = np.histogram(diff_vector, bins=N_bins)
    counts = counts / np.max(counts)
    centers = bins[0:-1] + (bins[1] - bins[0])

    plt.step(centers, counts, label = "All voxels", color = "k")

    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], color_utils.my_colors): #"target", 
        struct = ep_model.structure_set_clips_registered[structure_name]
    
        if structure_name != "optdisk":
            struct.resample_contour_points_in_grid(1)
    
        vec = diff_volume[struct.binary_mask[0:,y_bin,0:]].ravel()
    
    
        counts, _ = np.histogram(vec, bins=bins)
        counts = counts / np.max(counts)
        plt.step(centers, counts, label = structure_name.capitalize(), color = color, linewidth = 3)

    ax.set_yscale('log')
    ax.set_xlabel("Dose difference [fraction of prescribed dose]")
    ax.set_ylabel("Probability")
    ax.set_title("Voxel-by-Voxel Comparison")

    plt.grid(linewidth = 0.3)
    plt.legend()
    plt.suptitle(f'2D dosimetric comparison for P{ep_model.patient_config["patient_number"]} (dose engine - Eyeplan)')
    
    ax = fig.add_subplot(224)
    x = algo.medium.x_voxel_pos
    z = algo.medium.z_voxel_pos
    coords = (x,z)
    gamma3 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],3,3)
    valid_gamma3 = gamma3[~np.isnan(gamma3)]
    passrate3 = len(valid_gamma3[valid_gamma3 <= 1]) / len(valid_gamma3) * 100
    gamma2 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],2,2)
    valid_gamma2 = gamma2[~np.isnan(gamma2)]
    passrate2 = len(valid_gamma2[valid_gamma2 <= 1]) / len(valid_gamma2) * 100
    gamma1 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],1,1)
    valid_gamma1 = gamma1[~np.isnan(gamma1)]
    passrate1 = len(valid_gamma1[valid_gamma1 <= 1]) / len(valid_gamma1) * 100
    gamma3[np.isnan(gamma3)] = 0
    print(np.min(gamma3), np.max(gamma3) )
    try:
        cmap = color_utils.gamma_index_colormap(gamma3)
    except:
        cmap
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos,gamma3,cmap=cmap)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gamma Index (3%/3mm)', rotation=270, size='large', verticalalignment='baseline')
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], color_utils.my_colors):
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
            
        image = struct.binary_mask[0:,y_bin,0:]
            
        plt.contour(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos, image , [4], linewidths = 3, alpha=0)
       
    
    my_rays = list(filter(lambda ray: abs(ray.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    legend_handle = [plt.Line2D([0], [0], color='k'),plt.Line2D([0], [0], color='k'),plt.Line2D([0], [0], color='k')]
    legend_label = [f"Pass rate (3%/3mm): {round(passrate3,3)}%", f"Pass rate (2%/2mm): {round(passrate2,3)}%",f"Pass rate (1%/1mm): {round(passrate1,3)}%"]
    plt.legend(legend_handle, legend_label, loc='upper left',fontsize='large',handlelength=0, handletextpad=0, fancybox=True) 
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax.set_xlim(-20,20)
    ax.set_ylim(-15,15)
    plt.tight_layout()
    ax.set_title('Gamma Index')
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    plt.tight_layout()
    plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_dvh',
                bbox_inches='tight', pad_inches=0.1)
    
    out = {'Gamma Index (3%/3mm)':passrate3,'Gamma Index (2%/2mm)':passrate2,'Gamma Index (1%/1mm)':passrate1}
    return out

def find_first_valley(curve):
    for i in np.flip(range(1,len(curve)-1)):
        if curve[i] < 0.5:
            continue
        if curve[i - 1] > curve[i] < curve[i + 1]:
            return curve[i]
    return 1

def compare_depth_profile_with_measurements(algo2, ep_model, measurements, path = ""):
    
    p_iso = np.asarray([0,0,0])
    x_bin = np.asarray(np.floor((p_iso[0] - algo2.medium.mesh_origin[0])/algo2.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_iso[1] - algo2.medium.mesh_origin[1])/algo2.medium.resolution[1]), dtype = 'int')
    z_bin = np.asarray(np.floor((p_iso[2] - algo2.medium.mesh_origin[2])/algo2.medium.resolution[2]), dtype = 'int')
    
 
    fig = plt.figure()
    ref_R90 = normalization_utils.find_R90_from_right(measurements.iloc[:,0],measurements.iloc[:,1]*0.01)
    ep_R90 = normalization_utils.find_R90_from_left(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    de_R90 = normalization_utils.find_R90_from_left(algo2.SOBP.xes,algo2.SOBP.yes[::-1])
    ref_R10 = normalization_utils.find_R10_from_right(measurements.iloc[:,0],measurements.iloc[:,1]*0.01)
    ep_R10 = normalization_utils.find_R10_from_left(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    de_R10 = normalization_utils.find_R10_from_left(algo2.SOBP.xes,algo2.SOBP.yes[::-1])
    ref_distal_fall_off = ref_R10 - ref_R90
    ep_distal_fall_off = ep_R90 - ep_R10
    de_distal_fall_off = de_R90 - de_R10
    
    ref_xpos = normalization_utils.find_xpos_and_norm_factor(measurements.iloc[:,0],measurements.iloc[:,1]*0.01)
    ep_xpos = normalization_utils.find_xpos_and_norm_factor(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    de_xpos = normalization_utils.find_xpos_and_norm_factor(algo2.SOBP.xes[::-1],algo2.SOBP.yes)
    sobp_length = ref_xpos[0] - min(measurements.iloc[:,0])
    sobp_centre = sobp_length/2
    ref_R50 = normalization_utils.find_inflection_R(measurements.iloc[:,0],measurements.iloc[:,1]*0.01)
    ep_R50 = normalization_utils.find_inflection_L(ep_model.doseplane_h.z_voxel_pos,ep_model.doseplane_h.image[x_bin,0:])
    de_R50 = normalization_utils.find_inflection_L(algo2.SOBP.xes[::-1],algo2.SOBP.yes)
    
    if ep_model.patient_config["patient_number"] == 17301:
        ref_xpos = [21.303720,1.007504]
 
    shift = 0 - ref_R50
    depth1 = measurements.iloc[:,0] 
    depth3 = algo2.SOBP.xes[::-1]
    shift2 = 0 - de_R50
    shift3 = 0 - ep_R50
    depth4 = ep_model.doseplane_h.z_voxel_pos

    depth2 = []
    for i in range(len(depth1)):
        if depth1[i]>ref_R50:
            depth2.append(ref_R50 - (depth1[i]-ref_R50))
        else:
            depth2.append(ref_R50 + (ref_R50-depth1[i]))
    depth = []    
    for i in range(len(depth2)):
        dep = depth2[i] + shift
        depth.append(dep)

    dose1 = measurements.iloc[:,1]
    dose11 = dose1.tolist()
    dose11 = [0 if str(x).lower() == 'nan' else x for x in dose11]
    #np.isnan(dose11) == 0
    #dose11 = dose11[~np.isnan(dose11)]
# =============================================================================
#     dose = []
# 
#     if ref_xpos[1] > 90:
#         for i in range(len(dose11)):
#             dos = dose11[i]/ref_xpos[1]
#             dose.append(dos)
#     else:
#         for i in range(len(dose11)):
#             dos = dose11[i]/np.max(dose11)
#             dose.append(dos)
# =============================================================================

    depth_value = min(depth, key = lambda x : abs(x-sobp_centre))
    depth_index = depth.index(depth_value)

    dose = []   
    for i in range(len(dose11)):
         dos = dose11[i]/dose11[depth_index]
         dose.append(dos)
         
    de_depth = []    
    for i in range(len(depth3)):
        de_dep = depth3[i] + shift2
        de_depth.append(de_dep)
        
# =============================================================================
#     de_dose = []
#     if de_xpos[1] > 0.9:
#         for i in range(len(algo2.SOBP.yes)):
#             do = algo2.SOBP.yes[i]/de_xpos[1]
#             de_dose.append(do)
#     else:
#         for i in range(len(algo2.SOBP.yes)):
#             do = algo2.SOBP.yes[i]/np.max(algo2.SOBP.yes)
#             de_dose.append(do)
# =============================================================================
    
    de_depth_value = min(de_depth, key = lambda x : abs(x-sobp_centre))
    de_depth_index = de_depth.index(de_depth_value)

    de_dose = []
    for i in range(len(algo2.SOBP.yes)):
         do = algo2.SOBP.yes[i]/algo2.SOBP.yes[de_depth_index]
         de_dose.append(do)

    ep_depth = []    
    for i in range(len(depth4)):
        ep_dep = depth4[i] + shift3
        ep_depth.append(ep_dep)
   
    ax = fig.add_subplot(311)
    ax.plot(ep_depth,ep_model.doseplane_h.image[x_bin,0:],color='C0', label = "Eyeplan model")
    ax.plot(de_depth,de_dose,color='C1', label = "Dose engine")
    ax.plot(depth,dose,color='C2',label='Measurements')
     
    plt.legend(fontsize='large')
    ax.set_xlim(-10,30)
    ax.set_xlabel('Distance in z from distal 50% dose [mm]')
    ax.set_ylabel('Dose [fraction of prescribed dose]')
    ax.set_title(f'Dose in z direction through isocentre for patient {ep_model.patient_config["patient_number"]} with R50 aligned')
    
    ax = fig.add_subplot(312)
    ax.plot(ep_depth,ep_model.doseplane_h.image[x_bin,0:],color='C0', label = "Eyeplan model")
    ax.plot(de_depth,de_dose,color='C1', label = "Dose engine")
    ax.plot(depth,dose,color='C2',label='Measurements')
    
    plt.legend(fontsize='large')
    ax.set_xlim(-10,30)
    ax.set_ylim(0.9,1.1)
    ax.set_xlabel('Distance in z from distal 50% dose [mm]')
    ax.set_ylabel('Dose [fraction of prescribed dose]')
    ax.set_title(f'Dose in z direction through isocentre for patient {ep_model.patient_config["patient_number"]} with R50 aligned')
    
    ax = fig.add_subplot(313)
    ep_points = []
    de_points = []
    m_points = []
    sample_points = np.linspace(min(depth),max(depth),200)
    for point in sample_points:
        point1 = min(de_depth, key=lambda x: abs(point - x))
        index1 = de_depth.index(point1)
        de_points.append(de_dose[index1])
        
        point2 = min(ep_depth, key=lambda x: abs(point - x))
        index2 = ep_depth.index(point2)
        ep_points.append(ep_model.doseplane_h.image[x_bin,0:][index2])
        
        point3 = min(depth, key=lambda x: abs(point - x))
        index3 = depth.index(point3)
        m_points.append(dose[index3])
        
    de_ratio = []
    ep_ratio = []
    for point in range(len(m_points)):
        if m_points[point] == 0:
            if de_points[point]==0:
                de_ratio.append(1)
            else:
                de_ratio.append(de_points[point]+1)
            if ep_points[point]==0:
                ep_ratio.append(1)
            else:
                ep_ratio.append(ep_points[point]+1)
        else:
            de_ratio.append(de_points[point]/m_points[point])
            ep_ratio.append(ep_points[point]/m_points[point])
    ax.plot(sample_points,ep_ratio,label="Eyeplan model/measurement")
    ax.plot(sample_points,de_ratio,label="Dose engine/measurement")
    ax.axhline(y=1,alpha=0.25,color='C2')
    ax.set_xlabel('Distance in z from distal 50% dose [mm]')
    ax.set_ylabel('Dose Ratio')
    ax.set_title(f'Dose Ratio for P{ep_model.patient_config["patient_number"]}')
    ax.set_xlim(-10,30)
    ax.set_ylim(0.9,1.1)
    plt.legend(fontsize='large')
    
    de_diff = []
    ep_diff = []
    for point in range(len(m_points)):
        if sample_points[point] > 0:
            de_diff.append(abs(de_points[point]-m_points[point]))
            ep_diff.append(abs(ep_points[point]-m_points[point]))
        else:
            continue
    de_difference = np.mean(de_diff)
    ep_difference = np.mean(ep_diff)
    
    de_xpos_point  = normalization_utils.find_xpos_and_norm_factor(de_depth,de_dose)
    de_valley_point = find_first_valley(de_dose)
    ep_xpos_point  = normalization_utils.find_xpos_and_norm_factor(ep_depth,ep_model.doseplane_h.image[x_bin,0:])
    ep_valley_point = find_first_valley(ep_model.doseplane_h.image[x_bin,0:])
    m_xpos_point  = normalization_utils.find_xpos_and_norm_factor(depth,dose)
    m_valley_point = find_first_valley(dose)
    de_valley = de_xpos_point[1] - de_valley_point
    ep_valley = ep_xpos_point[1] - ep_valley_point
    m_valley = m_xpos_point[1] - m_valley_point
    de_valley_difference = abs(de_valley - m_valley)
    ep_valley_difference = abs(ep_valley - m_valley)
    
# =============================================================================
#     de_valley_diff = []
#     ep_valley_diff = []
#     for point in range(len(m_points)):
#         if point > 0 and point < 2:
#             de_valley_diff.append(abs(de_points[point]-m_points[point]))
#             ep_valley_diff.append(abs(ep_points[point]-m_points[point]))
#         else:
#             continue
#     de_valley_difference = np.mean(de_valley_diff)
#     ep_valley_difference = np.mean(ep_valley_diff)
# =============================================================================
# =============================================================================
#     de_diff = []
#     ep_diff = []
#     for point in sample_points:
#         ind = sample_points.tolist().index(point)
#         de_diff.append(abs(de_ratio[ind]-1))
#         ep_diff.append(abs(ep_ratio[ind]-1))
#         
#     de_difference = np.mean(de_diff)
#     ep_difference = np.mean(ep_diff)
#     
#     de_valley_diff = []
#     ep_valley_diff = []
#     for point in sample_points:
#         if point > 0 and point < 2:
#             ind = sample_points.tolist().index(point)
#             de_valley_diff.append(abs(de_ratio[ind]-1))
#             ep_valley_diff.append(abs(ep_ratio[ind]-1))
#         else:
#             continue
#         
#     de_valley_difference = np.mean(de_valley_diff)
#     ep_valley_difference = np.mean(ep_valley_diff)
# =============================================================================
    
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    plt.tight_layout()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_depth_profile_measurements',
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
    
    out = {'Dose engine 90% dose':de_R90, 'Eyeplan 90% dose':ep_R90, 'Measurements 90% dose':ref_R90, 'Dose engine distal fall-off':de_distal_fall_off, 'Eyeplan distal fall-off':ep_distal_fall_off, 'Measurements distal fall-off':ref_distal_fall_off,'Dose engine difference':de_difference,'Eyeplan difference':ep_difference, 'Dose engine valley difference':de_valley_difference, 'Eyeplan valley difference':ep_valley_difference}
    return out


        
    

def gamma_index_comparison(algo, ep_model, path=''):
    p_iso = np.asarray([0,0,0])
    x_bin = np.asarray(np.floor((p_iso[0] - algo.medium.mesh_origin[0])/algo.medium.resolution[0]), dtype = 'int')
    y_bin = np.asarray(np.floor((p_iso[1] - algo.medium.mesh_origin[1])/algo.medium.resolution[1]), dtype = 'int')
    
    eyeplan_image = ep_model.doseplane_h.interpolate_from(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos)
    x = algo.medium.x_voxel_pos
    z = algo.medium.z_voxel_pos
    coords = (x,z)
    gamma4 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],3,1)
    valid_gamma4 = gamma4[~np.isnan(gamma4)]
    passrate4 = len(valid_gamma4[valid_gamma4 <= 1]) / len(valid_gamma4) * 100
    gamma3 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],3,0.3)
    valid_gamma3 = gamma3[~np.isnan(gamma3)]
    passrate3 = len(valid_gamma3[valid_gamma3 <= 1]) / len(valid_gamma3) * 100
    gamma2 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],1,2)
    valid_gamma2 = gamma2[~np.isnan(gamma2)]
    passrate2 = len(valid_gamma2[valid_gamma2 <= 1]) / len(valid_gamma2) * 100
    gamma1 = pymedphys.gamma(coords, eyeplan_image, coords, algo.dose[0:,y_bin,0:],1,0.3)
    valid_gamma1 = gamma1[~np.isnan(gamma1)]
    passrate1 = len(valid_gamma1[valid_gamma1 <= 1]) / len(valid_gamma1) * 100
    gamma4[np.isnan(gamma4)] = 0
    gamma3[np.isnan(gamma3)] = 0
    gamma2[np.isnan(gamma2)] = 0
    gamma1[np.isnan(gamma1)] = 0
    
    fig = plt.figure()
    ax = fig.add_subplot(221)
    try:
        cmap4 = color_utils.gamma_index_colormap(gamma4)
    except:
        cmap4
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos,gamma4,cmap=cmap4)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gamma Index', rotation=270, size='large', verticalalignment='baseline')
    
    my_rays = list(filter(lambda ray: abs(ray.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")
        
    legend_handle = [plt.Line2D([0], [0], color='k')]
    legend_label = [f"Pass rate: {round(passrate4,3)}%"]
    plt.legend(legend_handle, legend_label, loc='upper left',fontsize='large',handlelength=0, handletextpad=0, fancybox=True) 
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax.set_xlim(-30,30)
    ax.set_ylim(-30,30)
    ax.set_title('3%/1mm')
        
    ax = fig.add_subplot(222)
    try:
        cmap3 = color_utils.gamma_index_colormap(gamma3)
    except:
        cmap3
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos,gamma3,cmap=cmap3)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gamma Index', rotation=270, size='large', verticalalignment='baseline')
    
    my_rays = list(filter(lambda ray: abs(ray.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    legend_handle = [plt.Line2D([0], [0], color='k')]
    legend_label = [f"Pass rate: {round(passrate3,3)}%"]
    plt.legend(legend_handle, legend_label, loc='upper left',fontsize='large',handlelength=0, handletextpad=0, fancybox=True) 
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax.set_xlim(-30,30)
    ax.set_ylim(-30,30)
    ax.set_title('3%/0.3mm')
    
    ax = fig.add_subplot(223)

    try:
        cmap2 = color_utils.gamma_index_colormap(gamma2)
    except:
        cmap2

    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos,gamma2,cmap=cmap2)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gamma Index', rotation=270, size='large', verticalalignment='baseline')
    
    my_rays = list(filter(lambda ray: abs(ray.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    legend_handle = [plt.Line2D([0], [0], color='k')]
    legend_label = [f"Pass rate: {round(passrate2,3)}%"]
    plt.legend(legend_handle, legend_label, loc='upper left',fontsize='large',handlelength=0, handletextpad=0, fancybox=True) 
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax.set_xlim(-30,30)
    ax.set_ylim(-30,30)
    ax.set_title('1%/2mm')
    
    ax = fig.add_subplot(224)
    try:
        cmap1 = color_utils.gamma_index_colormap(gamma1)
    except:
        cmap1
    plt.pcolor(algo.medium.z_voxel_pos, algo.medium.x_voxel_pos,gamma1,cmap=cmap1)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gamma Index', rotation=270, size='large', verticalalignment='baseline')
    
    my_rays = list(filter(lambda ray: abs(ray.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    legend_handle = [plt.Line2D([0], [0], color='k')]
    legend_label = [f"Pass rate: {round(passrate1,3)}%"]
    plt.legend(legend_handle, legend_label, loc='upper left',fontsize='large',handlelength=0, handletextpad=0, fancybox=True) 
    ax.set_xlabel('z [mm]')
    ax.set_ylabel('x [mm]')
    ax.set_xlim(-30,30)
    ax.set_ylim(-30,30)
    ax.set_title('1%/0.3mm')
    
    plt.suptitle(f'Gamma Index for P{ep_model.patient_config["patient_number"]}')
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    plt.tight_layout()
    plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_gamma_index',
                bbox_inches='tight', pad_inches=0.1)
    
    out = {'Gamma Pass Rate (1%/0.3mm)':passrate1,'Gamma Pass Rate (1%/2mm)':passrate2,'Gamma Pass Rate (3%/0.3mm)':passrate3,'Gamma Pass Rate (3%/1mm)':passrate4}
    return out


def overview_plot(algo, ep_model, path=''):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(f'Overview for P{ep_model.patient_config["patient_number"]}')
    ax.text(0.05,0.9,f'Dose engine nozzle: {algo.nozzle.title}',size='large')
    ax.text(0.05,0.8,f'Eyeplan nozzle: Target range: {ep_model.target_range_from_ptd} mm, Modulation range: {ep_model.modulation_range_from_ptd} mm', size='large')
    if ep_model.symmetrical_eye == True:
        ax.text(0.05,0.7,"Symmetrical",size='large')
    else:
        ax.text(0.05,0.7,"Asymmetrical",size='large')
    if ep_model.use_of_wedge == 1:
        ax.text(0.05,0.6,"Wedge",size='large')
    else:
        ax.text(0.05,0.6,"No wedge",size='large')
    if ep_model.treatment_eye == "R":
        ax.text(0.05,0.5,"Right eye",size='large')
    else:
        ax.text(0.05,0.5,"Left eye",size='large')
    read = pd.read_excel(r'G:\ASM_Groups_Users\User\Learmonth\comparisons\by_patient\Patient Overview.xlsx')
    patient_number = read['Patient Number'].tolist()
    comments = read['Comments'].tolist()
    index = patient_number.index(float(ep_model.patient_config["patient_number"]))
    ax.text(0.05,0.4,f'{comments[index]}',size='large')
    if "lid_point" in algo.config.keys() and "lid_thickness" in algo.config.keys():
        ax.text(0.05,0.3,f'Eyelid thickness: {ep_model.eyelid_thickness}',size='large')
    else:
        ax.text(0.05,0.3,'No eyelid involved',size='large')
    ax.text(0.05,0.2,f'Proximal margin: {ep_model.proximal_margin} mm, Distal margin: {ep_model.distal_margin} mm',size='large')
    ax.axis('off')
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)
    ax.set_xticks([])
    ax.set_yticks([])
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.61)
    plt.tight_layout()
    if type(path) == pathlib.WindowsPath:
        try:
            plt.savefig(path / f'P{ep_model.patient_config["patient_number"]}_a_overview',
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
        