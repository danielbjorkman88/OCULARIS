# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import matplotlib.pyplot as plt
import numpy as np
from pyproton.metrics.dvh import DVH
# from pyproton.volume.grid import Grid

from utils.color_utils import my_colors
import dose_engine
from matplotlib import cm
import math
from scipy.interpolate import interp1d


def plot_lateral_degradation(algo, ep_model, fraction1, tot_acc_dose, slice_idx,  list_sharpness, lateral_fall_offs, path):
    
    image1 = algo.dose[0:, slice_idx[1], 0:]
    image2 = tot_acc_dose[0:, slice_idx[1], 0:]
    x_pos = algo.medium.resolution[0] * \
        (slice_idx[0] + 0.5) + algo.medium.mesh_origin[0]
    y_pos = algo.medium.resolution[1] * \
        (slice_idx[1] + 0.5) + algo.medium.mesh_origin[1]
    z_pos = algo.medium.resolution[2] * \
        (slice_idx[2] + 0.5) + algo.medium.mesh_origin[2]
        
        
    N = fraction1.N_resampled
    
    rainbow_colors = cm.rainbow(np.linspace(0,1,N))
    
    

    
    # ----------- PTCOG 2 ---------------------------------------------
    fig = plt.figure()
    
    
    
    ax = plt.subplot(221)
    
    ax.axis('equal')
    
    
    diff_image = image2 - image1
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])
    plt.pcolor(algo.medium.Z, algo.medium.X, diff_image, cmap = "seismic", vmin = -limit, vmax = limit )
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Motion - no motion [fraction prescribed dose]", rotation=270, size='large', labelpad=10)
    
    x = algo.medium.z_voxel_pos
    y = algo.medium.x_voxel_pos
    X, Y = np.meshgrid(x, y)
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors):
    
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
    
        ax.contour(X, Y, struct.binary_mask[0:, slice_idx[1], 0:], [
            1], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )
    
    my_rays = list(filter(lambda x: abs(x.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")
    
    plt.plot(algo.central_axis.zes, algo.central_axis.xes, color = "k" )
    
    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    
    plt.axvline(x = z_pos, label = "Sampling plane", color = "grey", linestyle = ":" )
    
    
    
    plt.title("Dose difference from intra-fraction motion", fontsize = 12)
    
    
    plt.legend(loc='lower left', prop={'size': 6})
    # plt.xlim(-12.75, skin_plane)
    plt.xlim(algo.medium.minZ, algo.medium.maxZ)
    
    plt.xlabel("Z [mm]", fontsize = 10)
    plt.ylabel("X [mm]", fontsize = 10)
    
    
    
    ax = plt.subplot(223)
    plt.plot(range(len(list_sharpness)), list_sharpness, linewidth = 3, label = f"Fraction {fraction1._fraction_number}")
    
    plt.xlabel("Treatment time [s]", fontsize = 12)
    plt.ylabel("Lateral penumbra (R80 - R20) [mm]", fontsize = 12)
    
    plt.legend()
    plt.grid(linewidth = 0.3)
    
    
    plt.title(f"Lateral sharpness degradation z = {round(z_pos, 3)} mm", fontsize = 10)
    
    # plt.title("Lateral sharpness degradation in iso plane due to intra-fractional motion", fontsize = 10)
    
    
    
    plt.subplot(122)
    
        
    for i in range(N):
        xes, yes = lateral_fall_offs[i]
        plt.plot(xes, yes , label = f"idx = {i}. R80 - R20 = {list_sharpness[i]} mm", color = rainbow_colors[i])
            
        # print( norm_ref , acc_norm , norm_ref + acc_norm)
        
        # print(i, norm_ref)
    
    plt.xlabel("x [mm]", fontsize = 12)
    plt.ylabel("Dose intensity", fontsize = 12)
    
    plt.legend(prop={'size': 6})
    plt.grid(linewidth = 0.3)
    
    plt.xlim(6,13.5)
    
    
    plt.suptitle(f"P{ep_model.tps_config['patient_number']}, fraction {fraction1._fraction_number}", fontsize = 14)
    
    xlength = 14
    fig.set_size_inches(xlength, xlength/1.61)
    plt.show()
    try:
        plt.savefig(path / f"P{ep_model.tps_config['patient_number']}_plot_penumbra_degradation",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass



def plot_target_com_movement(fraction, tps_config, path):
    df = fraction.df
    df_resampled = fraction.df_resampled

    
    fig = plt.figure()
    
    # "ax = plt.subplot(311)
    
    plt.plot([], [] , label = "Moving average filter only", color = "k", linewidth = 0.5) 
    plt.plot([], [] , label = "Derived from downsampled signal", color = "k", linewidth = 2)
    
    xes, yes = df.t, df["target.com.x"]
    plt.plot([], [] , color = my_colors[0]) 
    plt.plot(xes, yes , color = my_colors[0], linewidth = 0.5)
    
    
    xes, yes = df_resampled.t, df_resampled["target.com.x"]
    plt.plot([], [] , label = "x", color = my_colors[0]) 
    plt.plot(xes, yes , color = my_colors[0], linewidth = 2)
    
    
    xes, yes = df.t, df["target.com.y"]
    plt.plot([], [] , color = my_colors[1]) 
    plt.plot(xes, yes , color = my_colors[1], linewidth = 0.5)
    
    
    xes, yes = df_resampled.t, df_resampled["target.com.y"]
    plt.plot([], [] , label = "y", color = my_colors[1]) 
    plt.plot(xes, yes , color = my_colors[1], linewidth = 2)
    
    
    xes, yes = df.t, df["target.com.z"]
    plt.plot([], [] , color = my_colors[2]) 
    plt.plot(xes, yes , color = my_colors[2], linewidth = 0.5)
    
    
    xes, yes = df_resampled.t, df_resampled["target.com.z"]
    plt.plot([], [] , label = "z", color = my_colors[2]) 
    plt.plot(xes, yes , color = my_colors[2], linewidth = 2)
    
    plt.legend()
    
    plt.ylabel("x, y or z [mm]", fontsize= 16)
    plt.xlabel("time [ms]", fontsize= 14)
    
    plt.grid(linewidth = 0.3)
    
    plt.title("Target COM", fontsize = 14)
    
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.6)
    plt.show()
    try:
        plt.savefig(path / f"P{tps_config['patient_number']}_plot_taget_com",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass


def plot_translation_only(fraction, path):

    df = fraction.df
    df_resampled = fraction.df_resampled
    
    fig = plt.figure()
    
    ax = plt.subplot(111)
    
    xes, yes = df.t, df.dist_origin
    xes = np.asarray(xes)
    yes = np.asarray(yes)
    plt.plot(xes, yes , label = "Translation only")
    
    
    xes, yes = df_resampled.t, df_resampled.dist_origin
    xes = np.asarray(xes)
    yes = np.asarray(yes)
    plt.plot(xes, yes , label = "Downsampled translation ")
    
    plt.legend()
    
    plt.xlabel("time [ms]", fontsize= 14)
    plt.ylabel("Distance to structure origin [mm]", fontsize= 14 ) 
    
    # plt.suptitle(f" Showing {algo.config['Slice']} ")
    
    plt.grid(linewidth = 0.3)
    
    # plt.title()
    
    xlength = 12
    fig.set_size_inches(xlength, xlength/1.6)
    plt.show()
    try:
        plt.savefig(path / "P{eyeplan_model.tps_config['patient_number']}_plot_translation",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass


def plot_motion_subfraction_frame(alg, alg_moved_model, eyeplan_model, idx_of_interest, acc_dose_subfraction, slice_idx,  my_path):
    
    x = alg.medium.z_voxel_pos
    y = alg.medium.x_voxel_pos
    X, Y = np.meshgrid(x, y)
    
    fig = plt.figure()
    
    ax = plt.subplot(221)
    
    image1 = alg_moved_model.dose[0:, slice_idx[1], 0:]
    plt.pcolor(alg.medium.Z, alg.medium.X, image1, label = "", cmap = "jet", vmin = 0, vmax = 1)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    
    # plt.title("")
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors):

        struct = eyeplan_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)

        ax.contour(X, Y, struct.binary_mask[0:,slice_idx[1],0:], [
            4], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )



    my_rays = list(filter(lambda x: abs(x.T[1]) < alg.medium.resolution[1]/2 , alg.collimator.aperture_rays_expanded))

    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    plt.plot(alg.central_axis.zes, alg.central_axis.xes, color = "k" )

    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    plt.xlim(alg.medium.minZ, alg.config["skin_plane"] )
    
    plt.title("Motion simulation", fontsize = 12)
    
    
    plt.xlabel("z [mm]", fontsize = 12)
    plt.ylabel("x [mm]", fontsize = 12)
    
    ax = plt.subplot(222)
    image2 = acc_dose_subfraction[0:, slice_idx[1], 0:]
    plt.pcolor(alg.medium.Z, alg.medium.X, image2, label = "", cmap = "jet", vmin = 0, vmax = 1)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors):

        struct = eyeplan_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)

        ax.contour(X, Y, struct.binary_mask[0:,slice_idx[1],0:], [
            4], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )



    my_rays = list(filter(lambda x: abs(x.T[1]) < alg.medium.resolution[1]/2 , alg.collimator.aperture_rays_expanded))

    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    plt.plot(alg.central_axis.zes, alg.central_axis.xes, color = "k" )

    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    plt.xlim(alg.medium.minZ, alg.config["skin_plane"] )
    
    
    plt.title("Subfraction", fontsize = 12)
    
    plt.xlabel("z [mm]", fontsize = 12)
    plt.ylabel("x [mm]", fontsize = 12)
    ax = plt.subplot(223)
    
    
    diff_image = image2 - image1
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])
    plt.pcolor(alg.medium.Z, alg.medium.X, diff_image, cmap = "seismic", vmin = -limit, vmax = limit)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose diff', rotation=270, size='large')
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors):

        struct = eyeplan_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)

        ax.contour(X, Y, struct.binary_mask[0:,slice_idx[1],0:], [
            4], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )



    my_rays = list(filter(lambda x: abs(x.T[1]) < alg.medium.resolution[1]/2 , alg.collimator.aperture_rays_expanded))

    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")

    plt.plot(alg.central_axis.zes, alg.central_axis.xes, color = "k" )

    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    plt.xlim(alg.medium.minZ, alg.config["skin_plane"] )
    
    plt.title("Acc - ref", fontsize = 8)
    
    plt.xlabel("z [mm]", fontsize = 12)
    plt.ylabel("x [mm]", fontsize = 12)
    
    
    
    ax = plt.subplot(224)
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], my_colors):
        
        struct = eyeplan_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
        
        
        dvh = DVH(alg.dose, struct.binary_mask)
        
        dvh2 = DVH(acc_dose_subfraction, struct.binary_mask)
        
        dose_fractions = np.linspace(0, 1, 1000)
        
        volumes = dvh.V(dose_fractions)
        volumes2 = dvh2.V(dose_fractions)
    
    
        plt.plot(dose_fractions*100, volumes*100,
                label=structure_name.capitalize(), color=color)
        plt.plot(dose_fractions*100, volumes2 *
                100, color=color, linestyle="--")
    
    plt.plot([], [], color = "k", linewidth = 2, label = "Reference")
    plt.plot([], [], color = "k", linestyle = "--" , label = "Subfraction")
    
    
    plt.legend()
    
    
    plt.xlabel("% Dose", fontsize = 12)
    plt.ylabel("% Volume or surface", fontsize = 12)
    
    plt.grid(linewidth = 0.3)
    
    
    
    plt.suptitle(
        f"Subfraction idx = {idx_of_interest}", fontsize=14)
    
    xlength = 15
    fig.set_size_inches(xlength, xlength/1.61803398875)
    plt.show()
    try:
        plt.savefig(my_path / f"P{eyeplan_model.tps_config['patient_number']}_plot_idx_{idx_of_interest}_",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass




def plot_dose_field_and_motion(algo, ep_model, fraction1, slice_idx, path):

    #-------------- PTCOG 1 ----------------------------------------------------
    fig = plt.figure()
    
    ax = plt.subplot(121)
    ax.axis('equal')
    
    
    x = algo.medium.z_voxel_pos
    y = algo.medium.x_voxel_pos
    X, Y = np.meshgrid(x, y)
    # ax.contour(X, Y, ep_model.structure_set_clips_registered["eyeglobe"].binary_mask[0:,config2["Slice"][1],0:], [1], colors='black')
    
    # ax.contour(X, Y, ep_model.structure_set_clips_registered["target"].binary_mask[0:,config2["Slice"][1],0:], [1], colors='green')
    
    # plt.pcolor(algo.medium.Z, algo.medium.X, acc_volume[0:,config2["Slice"][1],0:])
    plt.pcolor(algo.medium.Z, algo.medium.X, algo.dose[0:,slice_idx[1],0:], alpha = 1, cmap = "jet")
    
    
    cbar = plt.colorbar()
    cbar.set_alpha(1)
    cbar.draw_all()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    
    
    # x = algo.medium.z_voxel_pos
    # y = algo.medium.x_voxel_pos
    # X, Y = np.meshgrid(x, y)
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], my_colors):
    
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
    
        ax.contour(X, Y, struct.binary_mask[0:,slice_idx[1],0:], [
            4], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )
    
    
    
    my_rays = list(filter(lambda x: abs(x.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")
    
    plt.plot(algo.central_axis.zes, algo.central_axis.xes, color = "k" )
    
    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    
    # plt.axvline(x=z_pos, color = "k", linewidth = 1, linestyle = "--")
    # plt.axhline(y=y_pos, color = "k", linewidth = 1, linestyle = "--")
    # plt.plot([],[], color = "k", label = "Eyeglobe")
    # plt.plot([],[], color = "g", label = "Target")
    
    
    plt.xlim(algo.medium.minZ, algo.medium.maxZ)
    
    
    
    plt.legend(loc=3, fontsize=8)
    
    
    
    
    plt.xlabel("Z [mm]", fontsize = 15)
    plt.ylabel("X [mm]", fontsize = 15)
    
    plt.title("Dose field and patient eye structures", fontsize = 12)
    
    
    ax = plt.subplot(122)
    
    df = fraction1.df
    
    
    # plt.plot([], [] , label = "Moving average filter only", color = "k", linewidth = 0.5) 
    plt.plot([], [] , label = "Target center of mass position", color = "k", linewidth = 1) 
    # plt.plot([], [] , label = "Derived from downsampled signal", color = "k", linewidth = 1)
    
    xes, yes = df.t, df["target.com.x"]
    plt.plot([], [] , color = my_colors[0]) 
    plt.plot(xes, yes , color = my_colors[0], linewidth = 1, label = "x")
    
    
    # xes, yes = df_resampled.t, df_resampled["target.com.x"]
    # plt.plot([], [] , label = "x", color = my_colors[0]) 
    # plt.plot(xes, yes , color = my_colors[0], linewidth = 1)
    
    
    xes, yes = df.t, df["target.com.y"]
    plt.plot([], [] , color = my_colors[1]) 
    plt.plot(xes, yes , color = my_colors[1], linewidth = 1, label = "y")
    
    
    # xes, yes = df_resampled.t, df_resampled["target.com.y"]
    # plt.plot([], [] , label = "y", color = my_colors[1]) 
    # plt.plot(xes, yes , color = my_colors[1], linewidth = 1)
    
    
    xes, yes = df.t, df["target.com.z"]
    plt.plot([], [] , color = my_colors[2]) 
    plt.plot(xes, yes , color = my_colors[2], linewidth = 1, label = "z")
    
    
    # xes, yes = df_resampled.t, df_resampled["target.com.z"]
    # plt.plot([], [] , label = "z", color = my_colors[2]) 
    # plt.plot(xes, yes , color = my_colors[2], linewidth = 1)
    
    plt.legend()
    
    plt.ylabel("x, y or z [mm]", fontsize= 12)
    plt.xlabel("Treatment time [ms]", fontsize= 12)
    
    plt.grid(linewidth = 0.3)
    
    plt.title("Change in target position due to detected intra-fractional motion", fontsize = 12)
    
    
    xlength = 14
    fig.set_size_inches(xlength, xlength/2.7)
    plt.show()
    try:
        plt.savefig(path / "P{eyeplan_model.tps_config['patient_number']}_plot_dose_motion",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass





def plot_gaze_acceptance(algo, ep_model, fractions_dict, my_subfractions, plot_path):
        
    fraction1 = fractions_dict[1]
    df = fraction1.df_resampled
    
    i = 0
    
    gaze_vector_light = ep_model.gaze_vector_light
    centre_of_model = ep_model.centre_of_model
    
    pp = list(zip(df["pps.x"], df["pps.y"] , df["pps.z"]))[i]
    pp = np.asarray(pp)
    
    cp = list(zip(df["cps.x"], df["cps.y"] , df["cps.z"]))[i]
    cp = np.asarray(cp)
    
    length = math.dist(pp, cp)
    
    gaze_point = centre_of_model + gaze_vector_light*length
    
    dose_fractions = np.linspace(0, 1, 1000)
    structure_name = "target"
    struct = ep_model.structure_set_clips_registered[structure_name]
    struct.resample_contour_points_in_grid(1)
    mask = struct.binary_mask
    
    
    
    
    
    p1 = ep_model.centre_of_model + gaze_vector_light*5
    p2 = ep_model.centre_of_model - gaze_vector_light*5
    
    
    ray_ref = dose_engine.Ray(p1,p2)
    ref_point = ray_ref.find_intersect_with_plane(np.asarray([0,0,1]), np.asarray([0,0,algo.config['FID']]) )
    
    
    
    
    fig = plt.figure()
    
    ax = fig.add_subplot(121, projection='3d')
    ax.view_init(-100,90) 
    
    # ref.add_frame_to_plot(ax, ref,20)
    # ax.plot(ep_model.light_point[0], ep_model.light_point[1], ep_model.light_point[2], markerfacecolor='C1', markeredgecolor='C1', marker='o', markersize=10, alpha=1, label = f"Fixation point = {ep_model.light_point}")
    
    plt.quiver(ep_model.centre_of_model[0],ep_model.centre_of_model[1],ep_model.centre_of_model[2],gaze_vector_light[0] , gaze_vector_light[1], gaze_vector_light[2], color = "k", length=2)
    ax.scatter(centre_of_model[0], centre_of_model[1], centre_of_model[2], color = "C0", s = 30 )
    ax.scatter(gaze_point[0], gaze_point[1], gaze_point[2], color = "C0", s = 30)
    
    p1 = ep_model.centre_of_model + gaze_vector_light*5
    p2 = ep_model.centre_of_model - gaze_vector_light*5
    xes_light = [p1[0], p2[0]]
    yes_light = [p1[1], p2[1]]
    zes_light = [p1[2], p2[2]]
    
    
    # ray_ref = dose_engine.Ray(p1,p2)
    
    
    plt.plot(xes_light , yes_light, zes_light, color = "k", linestyle= "--", label = "Eyeplan model axis")
    
    
    
    
    counter = 0
    for fraction_number in fractions_dict.keys():
        
        fraction = fractions_dict[fraction_number]
        
        for idx in range(len(fraction.df_resampled)):
        
            # dvh = DVH(algo.dose, struct.binary_mask)
            dvh2 = DVH(my_subfractions[counter], mask)
        
        
            # volumes = dvh.V(dose_fractions)
            volumes2 = dvh2.V(dose_fractions)
        
        
            f = interp1d(dose_fractions, volumes2 )
            # f_rev = interp1d(volumes2, dose_fractions)
        
            color = "C0"
            
            if f(0.95) < 0.95:
                # print(idx)
                color = "r"
            
            
            new_center = fraction.transform_points_by_idx(centre_of_model, idx)[0]
            new_gaze_point = fraction.transform_points_by_idx(gaze_point, idx)[0]
        
            ax.scatter(new_center[0], new_center[1], new_center[2], color = color, s = 30, alpha = 0.2 )
            ax.scatter(new_gaze_point[0], new_gaze_point[1], new_gaze_point[2], color = color, s = 30, alpha = 0.2 )
        
            vec = new_gaze_point - new_center
        
            xes = [new_center[0], new_gaze_point[0]]
            yes = [new_center[1], new_gaze_point[1]]
            zes = [new_center[2], new_gaze_point[2]]
        
            plt.plot(xes, yes, zes, color = color, alpha = 0.2)
            counter +=1
    
        
            # plt.scatter(0.95*100, f(0.95)*100, color = color)
    
    plt.plot([],[], color = "C0", label = "Within threshold")
    plt.plot([],[], color = "r", label = "Exceeds threshold")
    
    
    
    
    plt.title("Gaze vector by frame")
    
    
    
    
    plt.legend()
    
    ax.set_xlabel('X', fontsize = 12)
    ax.set_ylabel('Y', fontsize = 12)
    ax.set_zlabel('Z', fontsize = 12)
    
    
    
    ax = plt.subplot(122)
    ax.axis('equal')
    
    
    
    circle2 = plt.Circle((0.0, 0.0), 60, color='k', fill = None)
    ax.add_patch(circle2)
    
    counter = 0
    for fraction_number in fractions_dict.keys():
        
        fraction = fractions_dict[fraction_number]
        
        for idx in range(len(fraction.df_resampled)):
            
            
    # for idx in range(N):
            # dvh = DVH(algo.dose, struct.binary_mask)
            dvh2 = DVH(my_subfractions[counter], mask)
        
            # volumes = dvh.V(dose_fractions)
            volumes2 = dvh2.V(dose_fractions)
        
            f = interp1d(dose_fractions, volumes2 )
            # f_rev = interp1d(volumes2, dose_fractions)
        
            color = "C0"
            
            if f(0.95) < 0.95:
                # print(idx)
                color = "r"
            
            
            new_center = fraction.transform_points_by_idx(centre_of_model, idx)[0]
            new_gaze_point = fraction.transform_points_by_idx(gaze_point, idx)[0]
        
        
            ray = dose_engine.Ray(new_center, new_gaze_point)
            point = ray.find_intersect_with_plane(np.asarray([0,0,1]), np.asarray([0,0,algo.config['FID']]) )
            plt.scatter(point[0], point[1], color = color)
            counter += 1
    
    
    ax.scatter(ref_point[0], ref_point[1], color = "k", s = 30, label = "Reference fixation point")
    
    plt.plot([],[], color = "C0", label = " > 95% of dose delivered to 95% of target volume")
    plt.plot([],[], color = "r", label = " < 95% of dose delivered to 95% of target volume")
    
    
    plt.legend()
    
    plt.grid(linewidth = 0.3)
    
    limit = 120
    
    plt.xlim(-limit,limit)
    plt.ylim(-limit,limit)
    
    
    plt.xlabel("x [mm]", fontsize = 12)
    plt.ylabel("y [mm]", fontsize = 12)
    
    plt.title("Fixation light plane", fontsize = 12)
    
    plt.suptitle(f"P{ep_model.tps_config['patient_number']}", fontsize = 14)
    
    
    xlength = 12
    fig.set_size_inches(xlength, xlength/2.2)
    plt.show()
    
    try:
        plt.savefig(plot_path / f"P{ep_model.tps_config['patient_number']}_plot_gaze_acceptance",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass




def plot_dose_comp_3v1(algo, ep_model, tot_acc_dose, tps_config, plot_path):

    slice_idx = algo.config["Slice"]
    
    image1 = algo.dose[0:, slice_idx[1], 0:]
    image2 = tot_acc_dose[0:, slice_idx[1], 0:]
    
    
    x = algo.medium.z_voxel_pos
    y = algo.medium.x_voxel_pos
    X, Y = np.meshgrid(x, y)
    
    
    fig = plt.figure()
    
    ax = plt.subplot(131)
    ax.axis('equal')
    
    
    
    diff_image = image2 - image1
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])
    plt.pcolor(algo.medium.Z, algo.medium.X, diff_image, cmap = "seismic", vmin = -limit, vmax = limit)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose diff', rotation=270, size='large')
    
    plt.title("Accumulated dose - ref", fontsize = 12)
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors):
    
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
    
        ax.contour(X, Y, struct.binary_mask[0:,slice_idx[1],0:], [
            4], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )
    
    
    
    my_rays = list(filter(lambda x: abs(x.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")
    
    plt.plot(algo.central_axis.zes, algo.central_axis.xes, color = "k" )
    
    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    
    # plt.axvline(x = z_pos, label = "Sampling plane", color = "grey", linestyle = ":" )
        
    
    
    plt.xlim(algo.medium.minZ, algo.config["skin_plane_point"][2])
    plt.xlabel("z [mm]", fontsize = 12)
    plt.ylabel("x [mm]", fontsize = 12)


    ax = plt.subplot(132)
    

    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target", "macula"], my_colors):
        
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
        
        
        dvh = DVH(algo.dose, struct.binary_mask)
        
        dvh2 = DVH(tot_acc_dose, struct.binary_mask)
        
        dose_fractions = np.linspace(0, 1, 1000)
        
        volumes = dvh.V(dose_fractions)
        volumes2 = dvh2.V(dose_fractions)
    
    
        plt.plot(dose_fractions*100, volumes*100,
                label=structure_name.capitalize(), color=color)
        plt.plot(dose_fractions*100, volumes2 *
                100, color=color, linestyle="--")
    
    plt.plot([], [], color = "k", linewidth = 2, label = "Reference")
    plt.plot([], [], color = "k", linestyle = "--" , label = "Accumulated dose")
    
    
    plt.legend()
    
    
    plt.xlabel("% Dose", fontsize = 12)
    plt.ylabel("% Volume or surface", fontsize = 12)
    
    plt.grid(linewidth = 0.3)
    


    ax = plt.subplot(133)
    
    diff_volume = tot_acc_dose - algo.dose
    diff_vector = diff_volume.ravel()
    
    N_bins = 50
    
    

    
    counts, bins = np.histogram(diff_vector, bins=N_bins)
    counts = counts / np.max(counts)
    
    centers = bins[0:-1] + (bins[1] - bins[0])
    
    
    plt.step(centers, counts, label = "All voxels", color = "k")
    
    
    
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors): #"target", 
        struct = ep_model.structure_set_clips_registered[structure_name]
        
        if structure_name != "optdisk":
            struct.resample_contour_points_in_grid(1)
        
        vec = diff_volume[struct.binary_mask].ravel()
        
        
        counts, _ = np.histogram(vec, bins=bins)
        counts = counts / np.max(counts)
        plt.step(centers, counts, label = structure_name.capitalize(), color = color, linewidth = 3)
        
    
    
    
    plt.yscale('log')
    plt.xlabel("Dose difference (acc - ref)", fontsize = 12)
    plt.ylabel("Probability", fontsize = 12)
    
    plt.grid(linewidth = 0.3)
    plt.legend()
    
    plt.suptitle(f"Patient {tps_config['patient_number']}", fontsize = 16)
    
    # plt.title(
    #     f"Intra-fractional motion impact on dose difference histograms. Patient {tps_config['patient_number']}", fontsize=16)
    



    xlength = 15
    fig.set_size_inches(xlength, xlength/3.4)
    plt.show()
    try:
        plt.savefig(plot_path / f"Patient_{tps_config['patient_number']}_if_motion",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass
    
    


def plot_penumbra_increase(algo, ep_model, tot_acc_dose, area_ratio_dict, areas_dict, tps_config, plot_path):
    
    slice_idx = algo.config["Slice"]
    
    tot_volume_acc = 0
    tot_volume_ref = 0
    for z in areas_dict.keys():
        tot_volume_acc += areas_dict[z][0]*algo.medium.resolution[2]
        tot_volume_ref += areas_dict[z][1]*algo.medium.resolution[2]
        
    
    volume_increase = tot_volume_acc/tot_volume_ref
    
    
    x = algo.medium.z_voxel_pos
    y = algo.medium.x_voxel_pos
    X, Y = np.meshgrid(x, y)
    
    fig = plt.figure()
    
    ax = plt.subplot(111)
    ax.axis('equal')
    
    image1 = algo.dose[0:, slice_idx[1], 0:]
    image2 = tot_acc_dose[0:, slice_idx[1], 0:]
    diff_image = image2 - image1
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])
    plt.pcolor(algo.medium.Z,algo.medium.X, diff_image, cmap = "seismic", vmin = -limit, vmax = limit)
    cbar = plt.colorbar( pad=0.2)
    cbar.ax.set_ylabel('Dose diff', rotation=270, size='large')
    
    
    for structure_name, color in zip(["eyeglobe", "lens", "cornea", "target" , "macula"], my_colors):
    
        struct = ep_model.structure_set_clips_registered[structure_name]
        struct.resample_contour_points_in_grid(1)
    
        ax.contour(X, Y, struct.binary_mask[0:,slice_idx[1],0:], [
            4], colors=color, linewidths = 3)
        
        plt.plot([],[], label = structure_name.capitalize() )
        
    
    
    my_rays = list(filter(lambda x: abs(x.T[1]) < algo.medium.resolution[1]/2 , algo.collimator.aperture_rays_expanded))
    
    plt.plot([],[], color = "k", label = "Projected aperture", linestyle = "--")
    for ray in my_rays:
        new_T = ray.T + ray.vec_norm*100
        new_ray = dose_engine.Ray(ray.S, new_T)
        plt.plot(new_ray.zes, new_ray.xes, color = "k", linestyle = "--")
    
    plt.plot(algo.central_axis.zes, algo.central_axis.xes, color = "k" )
    
    plt.scatter(0,0, label = "Iso", s = 60, color = "r")
    plt.ylabel("x [mm]", fontsize = 12)
    plt.xlabel("z [mm]", fontsize = 12)
    
    ax2 = ax.twinx()
    
    xes = []
    yes = []
    for z_bin in area_ratio_dict.keys():
        xes.append(z_bin)
        yes.append(area_ratio_dict[z_bin])
        
        
    plt.plot(xes,yes, color = "C1", label = "Acc/ref penumbra area", linestyle = "--")
    
    plt.ylabel("Ratio penumbra area", fontsize = 12, rotation=270)
    plt.xlabel("z [mm]", fontsize = 12)
    
    plt.title(f"Penumbra volume acc/ref = {round(volume_increase, 4)}", fontsize = 12)
    
    ax2.yaxis.set_label_coords(1.2, 0.5)
    
    # ax.spines['bottom'].set_color('red')
    # ax2.spines['top'].set_color('C1')
    # ax.xaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='C1')
    
    plt.legend()
    
    plt.xlim(algo.medium.minZ, algo.medium.maxZ)
    
    
    plt.suptitle(f"Patient_{tps_config['patient_number']}", fontsize = 12)
    
    
    xlength = 8
    fig.set_size_inches(xlength, xlength/1.61)
    plt.show()
    try:
        plt.savefig(plot_path / f"Patient_{tps_config['patient_number']}_if_penumbra_area_ratio",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass


    
def plot_contours_cross_section(algo, ep_model, tot_acc_dose, tps_config, plot_path):
    
    slice_idx = algo.config["Slice"]
    
    
    x_pos = algo.medium.resolution[0] * \
        (slice_idx[0] + 0.5) + algo.medium.mesh_origin[0]
    y_pos = algo.medium.resolution[1] * \
        (slice_idx[1] + 0.5) + algo.medium.mesh_origin[1]
    z_pos = algo.medium.resolution[2] * \
        (slice_idx[2] + 0.5) + algo.medium.mesh_origin[2]
    
    p = np.asarray([x_pos, y_pos, z_pos])
    
    
    
    # N_subfractions = 0
    # for fraction_number in fractions:
    #     filenames_f = filenames_at_path(Path(basepath  + str(fraction_number)))
    #     N_subfractions += len(filenames_f)
    
    
    # tot_acc_dose, my_subfractions = accumulate_dose(basepath, algo, tps_config, N_subfractions)
    
    
    
    
    
    x = algo.medium.x_voxel_pos
    y = algo.medium.y_voxel_pos
    X, Y = np.meshgrid(x, y)
    
    
    
    z_bin = slice_idx[2]
    
    
    fig = plt.figure()
    
    ax = plt.subplot(221)
    ax.axis('equal')
    
    image1 = algo.dose[0:, 0:, z_bin].T
    plt.pcolor(algo.medium.trans_X, algo.medium.trans_Y, image1, label = "", cmap = "jet", vmin = 0, vmax = 1)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    
    
    
    ax.contour(X, Y, image1, [0.2, 0.8], colors= "k", linestyles = [":", "--"])
    plt.plot([],[], color = "k", linestyle = "--", label = "80 % iso")
    plt.plot([],[], color = "k", linestyle = ":", label = "20 % iso")
    
    plt.legend()
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    
    plt.xlabel("x [mm]", fontsize = 12)
    plt.ylabel("y [mm]", fontsize = 12)
    
    plt.title("Reference", fontsize = 14)
    
    
    ax = plt.subplot(222)
    ax.axis('equal')
    
    image2 = tot_acc_dose[0:, 0:, z_bin].T
    plt.pcolor(algo.medium.trans_X, algo.medium.trans_Y, image2, label = "", cmap = "jet", vmin = 0, vmax = 1)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose', rotation=270, size='large')
    
    
    
    ax.contour(X, Y, image2, [0.2, 0.8], colors= "k", linestyles = [":", "--"])
    plt.plot([],[], color = "k", linestyle = "--", label = "80 % iso")
    plt.plot([],[], color = "k", linestyle = ":", label = "20 % iso")
    
    plt.legend()
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    
    plt.xlabel("x [mm]", fontsize = 12)
    plt.ylabel("y [mm]", fontsize = 12)
    
    plt.title("Acc", fontsize = 14)
    
    ax = plt.subplot(223)
    ax.axis('equal')
    
    
    # contour1 = ax.contour(X, Y, image1, [0.2, 0.8], colors= "C0", linestyles = [":", "--"])
    # contour2 = ax.contour(X, Y, image2, [0.2, 0.8], colors= "C1", linestyles = [":", "--"])
    
    
    # plt.plot([],[], color = "C0", label = "Reference")
    # plt.plot([],[], color = "C1", label = "Acc")
    
    diff_image = image2 - image1
    limit = np.max([abs(np.min(diff_image)) , np.max(diff_image) ])
    plt.pcolor(algo.medium.trans_X, algo.medium.trans_Y, diff_image, cmap = "seismic", vmin = -limit, vmax = limit)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dose diff (acc - ref)', rotation=270, size='large')
    
    
    plt.legend()
    
    plt.grid(linewidth = 0.3)
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    
    plt.xlabel("x [mm]", fontsize = 12)
    plt.ylabel("y [mm]", fontsize = 12)
    
    
    # ax = plt.subplot(224)
    # ax.axis('equal')
    
    # contour1 = ax.contour(X, Y, image1, [0.2, 0.8], colors= "C0", linestyles = [":", "--"])
    # contour2 = ax.contour(X, Y, image2, [0.2, 0.8], colors= "C1", linestyles = [":", "--"])
    
    # for angle in range(0,360,5):
        
    # # angle = 58
                
    #     sampling_line = LineString([(-amp*np.cos(math.radians(angle)), -amp*np.sin(math.radians(angle))), (amp*np.cos(math.radians(angle)), amp*np.sin(math.radians(angle)))])
        
    #     inner1, outer1, d1_1, d1_2 = find_intersections(contour1, sampling_line)
        
    #     inner2, outer2, d2_1, d2_2  = find_intersections(contour2, sampling_line)
        
    
    #     for p1, p2 in zip(inner1, outer1):
    
    #         xes = [p1[0], p2[0]]
    #         yes = [p1[1], p2[1]]
    #         plt.plot(xes, yes, color = "C0", linewidth = 2)
        
        
    #     for p1, p2 in zip(inner2, outer2):
    #         xes = [p1[0], p2[0]]
    #         yes = [p1[1], p2[1]]
    #         plt.plot(xes, yes, color = "C1", linewidth = 0.5)
    
    
    
    
    # plt.plot([],[], color = "C0", label = "Reference")
    # plt.plot([],[], color = "C1", label = "Acc")
    
    # plt.plot([],[], color = "k", linestyle = "--", label = "80 % iso")
    # plt.plot([],[], color = "k", linestyle = ":", label = "20 % iso")
    
    # plt.legend()
    
    # plt.grid(linewidth = 0.3)
    # plt.xlim(-10,10)
    # plt.ylim(-10,10)
    # plt.xlabel("x [mm]", fontsize = 12)
    # plt.ylabel("y [mm]", fontsize = 12)
    
    plt.suptitle(f"z index = {z_bin}")
    
    xlength = 15
    fig.set_size_inches(xlength, xlength/1.61)
    plt.show()
    try:
        plt.savefig(plot_path / f"Patient_{tps_config['patient_number']}_if_contours_cross_section",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass
    return ax






def plot_position_and_dose(ep_model, fractions_dict, my_subfractions, plot_path):
    
    
    df = fractions_dict[1].df_resampled
    
    gaze_vector_light = ep_model.gaze_vector_light
    centre_of_model = ep_model.centre_of_model
    
    i = 0
    pp = list(zip(df["pps.x"], df["pps.y"] , df["pps.z"]))[i]
    pp = np.asarray(pp)

    cp = list(zip(df["cps.x"], df["cps.y"] , df["cps.z"]))[i]
    cp = np.asarray(cp)

    length = math.dist(pp, cp)

    gaze_point = centre_of_model + gaze_vector_light*length
    
    dose_fractions = np.linspace(0, 1, 1000)
    structure_name = "target"
    struct = ep_model.structure_set_clips_registered[structure_name]
    struct.resample_contour_points_in_grid(1)
    mask = struct.binary_mask
    
    xes = []
    yes = []
    zes = []
    vec_xes = []
    vec_yes = []
    vec_zes = []
    
    target_com_x = []
    target_com_y = []
    target_com_z = []
    
    threshold_color = []
    
    
    counter = 0
    for fraction_number in fractions_dict.keys():
        
        fraction = fractions_dict[fraction_number]
        fraction.fill_table_resampled()
        
        df = fraction.df_resampled
        
        for idx in range(len(df)):
        
            # dvh = DVH(algo.dose, struct.binary_mask)
            dvh2 = DVH(my_subfractions[counter], mask)
        
        
            # volumes = dvh.V(dose_fractions)
            volumes2 = dvh2.V(dose_fractions)
        
        
            f = interp1d(dose_fractions, volumes2 )
            # f_rev = interp1d(volumes2, dose_fractions)
        
            color = "C0"
            
            if f(0.95) < 0.95:
                # print(idx)
                color = "r"
            
            
            new_center = fraction.transform_points_by_idx(centre_of_model, idx)[0]
    
            new_gaze_point = fraction.transform_points_by_idx(gaze_point, idx)[0]
            
            target_com_x.append( df['target.com.x'][idx])
            target_com_y.append( df['target.com.y'][idx])
            target_com_z.append( df['target.com.z'][idx])
            
            # ax.scatter(new_center[0], new_center[1], new_center[2], color = color, s = 30, alpha = 0.2 )
            # ax.scatter(new_gaze_point[0], new_gaze_point[1], new_gaze_point[2], color = color, s = 30, alpha = 0.2 )
        
            vec = new_gaze_point - new_center
            vec = vec / np.linalg.norm(vec)
                    
            xes.append(new_center[0] - centre_of_model[0])
            yes.append(new_center[1] - centre_of_model[1])
            zes.append(new_center[2] - centre_of_model[2])
            vec_xes.append(vec[0])
            vec_yes.append(vec[1])
            vec_zes.append(vec[2])
            threshold_color.append(color)
            
            
        
            # xes = [new_center[0], new_gaze_point[0]]
            # yes = [new_center[1], new_gaze_point[1]]
            # zes = [new_center[2], new_gaze_point[2]]
        
            # plt.plot(xes, yes, zes, color = color, alpha = 0.2)
            counter +=1
            
            
    
    
    fig = plt.figure()
    
    ax = plt.subplot(231)
    
    plt.scatter(np.ones(len(xes)), xes, color = threshold_color)
    plt.scatter(2*np.ones(len(yes)), yes, color = threshold_color)
    plt.scatter(3*np.ones(len(zes)), zes, color = threshold_color)
    
    plt.boxplot([xes, yes, zes])
    
    # plt.legend()
    # plt.xlabel("% Dose", fontsize = 12)
    plt.ylabel("Relative position [mm]", fontsize = 12)
    
    plt.grid(linewidth = 0.3)
    ax.set_xticklabels(["x", "y", "z"] , fontsize = 20)
    
    plt.title("Translation only", fontsize = 12)
    
    ax = plt.subplot(232)
    
    
    plt.scatter(1*np.ones(len(target_com_x)), target_com_x, color = threshold_color)
    plt.scatter(2*np.ones(len(target_com_y)), target_com_y, color = threshold_color)
    plt.scatter(3*np.ones(len(target_com_z)), target_com_z, color = threshold_color)
    
    plt.boxplot([target_com_x, target_com_y, target_com_z])
    
    # plt.legend()
    # plt.xlabel("% Dose", fontsize = 12)
    plt.ylabel("Target position [mm]", fontsize = 12)
    
    plt.title("Target COM", fontsize = 12)
    
    plt.grid(linewidth = 0.3)
    ax.set_xticklabels(["x", "y", "z"] , fontsize = 20)
    
    ax = plt.subplot(233)
    
    plt.scatter(np.ones(len(vec_xes)), vec_xes, color = threshold_color)
    plt.scatter(2*np.ones(len(vec_yes)), vec_yes, color = threshold_color)
    plt.scatter(3*np.ones(len(vec_zes)), vec_zes, color = threshold_color)
    
    plt.boxplot([ vec_xes, vec_yes, vec_zes])
    
    plt.title("Gaze direction", fontsize = 12)
    
    
    plt.ylabel("Vector component", fontsize = 12)
    
    plt.grid(linewidth = 0.3)
    
    ax.set_xticklabels(["x'", "y'", "z'"] , fontsize = 20)
    
    
    plt.scatter([],[], color = "C0", label = "Within threshold")
    plt.scatter([],[], color = "r", label = "Exceeds threshold")
    
    
    plt.legend()
    
    ax = plt.subplot(2,1,2)
    
    
    
    
    
    # fig = plt.figure()
    
    # ax = plt.subplot(111)
    
    counter = 0
    elapsed_time = 0
    elapsed_idx = 0
    for fraction_number in fractions_dict.keys():
    
            
        # fraction_number = 1
        fraction = fractions_dict[fraction_number]
        fraction.fill_table_resampled()
        
        df = fraction.df_resampled
        local_N = len(df)
        
        local_colors = threshold_color[elapsed_idx:elapsed_idx + local_N]
        local_time = df["t"]/1000 + elapsed_time
        
        ax.scatter(local_time, df['target.com.x'], color = local_colors, s = 80)
        ax.scatter(local_time, df['target.com.y'], color = local_colors, s = 80)
        ax.scatter(local_time, df['target.com.z'], color = local_colors, s = 80)
        
        plt.plot(local_time, df['target.com.x'], color = "C2", linewidth = 3)
        plt.plot(local_time, df['target.com.y'], color = "c", linewidth = 3)
        plt.plot(local_time, df['target.com.z'], color = "grey", linewidth = 3)
        
        for t, color  in zip(local_time, local_colors):
            if color == "r":
                ax.axvline(x = t, color = "r", linewidth = 0.3)
        
        elapsed_time += df.t[local_N - 1]/1000
        elapsed_idx += local_N
        ax.axvline(x = elapsed_time, color = "k")
    # plt.scatter(np.ones(len(vec_xes)), vec_xes, color = threshold_color)
    # plt.scatter(2*np.ones(len(vec_yes)), vec_yes, color = threshold_color)
    # plt.scatter(3*np.ones(len(vec_zes)), vec_zes, color = threshold_color)
    
    # plt.boxplot([ vec_xes, vec_yes, vec_zes])
    # plt.legend()
    # plt.ylabel("Gaze direction", fontsize = 12)
    # plt.ylabel("% Volume or surface", fontsize = 12)
    
    plt.plot([],[], color = "C2", label = "x", linewidth = 3)
    plt.plot([],[], color = "c", label = "y", linewidth = 3)
    plt.plot([],[], color = "grey", label = "z", linewidth = 3)
    
    plt.plot([], [], color = "k",label = "Fraction")
    plt.plot([], [], color = "r", label = "Exceeds threshold")
    
    plt.legend()
    
    plt.grid(linewidth = 0.3)
    
    # ax.set_xticklabels(["x'", "y'", "z'"] , fontsize = 20)
    
    plt.ylabel("Target COM", fontsize = 12)
    plt.xlabel("Treatment time [s]", fontsize = 12)
    
    plt.suptitle(f"P{ep_model.tps_config['patient_number']}", fontsize = 16)
    
    xlength = 15
    fig.set_size_inches(xlength, xlength/1.61)
    plt.show()
    try:
        plt.savefig(plot_path / f"P{ep_model.tps_config['patient_number']}_plot_position_and_dose",
                    bbox_inches='tight', pad_inches=0.1)
    except:
        pass

