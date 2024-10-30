# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


from plot_utils.Plotter import Plotter

from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path
from configuration import config, data
from matplotlib.pyplot import cm
import numpy as np

import matplotlib.pyplot as plt
import reference_frames
import geometry

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M )
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

class PlotProjectedAperture(Plotter):
    


    def plot(self):

        
        ref = reference_frames.Reference_Frame(self.algo.config, 'TreatmentRoom')
        geo = geometry.OptisGeometry(self.algo.config)


        
        p_on_axis = geo.aperture_point + self.algo.central_axis.vec_norm*70
                
        
        iso_projection = self.algo.collimator.project_aperture(p_on_axis)
        
        
        Nplanes = 10
        
        colors = cm.rainbow(np.linspace(0, 1, Nplanes))
        
        dist_sifts = np.linspace(-20, 90, Nplanes)
        
        
        my_plane_hulls = []
        
        
        for val in dist_sifts:
            p_on_axis = geo.aperture_point + self.algo.central_axis.vec_norm*(70 + val)
            projection = self.algo.collimator.project_aperture(p_on_axis)
            my_plane_hulls.append(projection)
        


        
        
        
        fig = plt.figure()
        
        #------------------------------------------
        ax1 = fig.add_subplot(121, projection='3d')
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
        
        
        
        
        plt.plot([],[],[], label = "Expanded aperture", color = "C1", linestyle = "--")
        for edge in self.algo.collimator.expanded_hull.boundary_edges:
            xes = edge[0]
            yes = edge[1]
            zes = edge[2]
            plt.plot(xes, yes, zes, color = "C1", linestyle = "--")
            
        
        
        
        
        # for i in range(len(colors)):
        #     p = self.algo.collimator.point_ray_pair_expanded_aperture[i][0]
        #     plt.scatter(p[0], p[1], 70, color = colors[i])
        
            
        # for i in self.algo.collimator.aperture_hull.boundary_vertices:
            
        #     if i >= len(colors):
        #         break
            
        #     ray = algo.collimator.aperture_rays[i]
        #     # ray = self.algo.collimator.point_ray_pair_expanded_aperture[i][1]
        
        #     a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color=colors[i])
        #     ax1.add_artist(a)
        
        
        
        
        
        
            
        plt.plot([],[],[], label = "Collimator aperture", color = "C0", linestyle = "--")
        for edge in self.algo.collimator.aperture_hull.boundary_edges:
            xes = edge[0]
            yes = edge[1]
            zes = edge[2]
            plt.plot(xes, yes, zes, color = "C0", linestyle = "--")    
            
        # plt.plot(self.algo.central_axis.xes, self.algo.central_axis.yes, self.algo.central_axis.zes,  color = "k", label = "Central axis")
        
        
        plt.plot([],[],[], color = "k", label = "Central axis")
        a = Arrow3D(self.algo.central_axis.xes, self.algo.central_axis.yes, self.algo.central_axis.zes , **arrow_prop_dict, color="k")
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
        
        

        
        
        

        
        ax = plt.subplot(122)
        
        

        
        
        
        plt.plot(self.algo.collimator.aperture_hull.boundary_in_plane_xes, self.algo.collimator.aperture_hull.boundary_in_plane_yes, label = "Collimator aperture", color = "C0")
        
        
        plt.plot(self.algo.collimator.expanded_hull.boundary_in_plane_xes, self.algo.collimator.expanded_hull.boundary_in_plane_yes, label = "Expanded aperture", color = "C1")
        
        
        plt.plot(iso_projection.boundary_in_plane_xes, iso_projection.boundary_in_plane_yes, label = "Projection iso plane", color = "k")
        
        
        for p in self.algo.collimator.points_expanded:
            plt.scatter(p[0], p[1], color = "C1")
        
        
        for i in range(len(self.algo.collimator.aperture_hull.boundary_in_plane_xes)):
            x = self.algo.collimator.aperture_hull.boundary_in_plane_xes[i]
            y = self.algo.collimator.aperture_hull.boundary_in_plane_yes[i]
            plt.scatter(x,y, color = "C0")

        
        
        for idx, projection, color in zip(np.linspace(0,len(dist_sifts)-1, len(dist_sifts)), my_plane_hulls, colors):
            # plt.plot([],[],[], label = f"{dist_sifts[int(idx)]} mm")
        
            plt.plot(projection.boundary_in_plane_xes, projection.boundary_in_plane_yes, color = color, linestyle = "-", linewidth = 1)
        
                
                
        plt.xlabel("x [mm]", fontsize = 16)
        plt.ylabel("y [mm]", fontsize = 16)
        
        plt.xlim(-6.5,6.5)
        plt.ylim(-6.5,6.5)
        
        
        plt.title("In plane" , fontsize = 12)
        
        plt.grid(linewidth = 0.3)
        # xlength = 8
        # fig.set_size_inches(xlength, xlength)
        # plt.show()
        
        plt.legend()
        
        
        if self.algo.wedge.traced:
            plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees, Wedge cover = {self.algo.config["wedge_cover"]}', fontsize = 16)
        else:
            plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees', fontsize = 16)
        
        
        
        
        xlength = 12
        fig.set_size_inches(xlength, xlength/2)
        plt.show()
        

        try:
            plt.savefig(self.algo.plot_path / "Plot_projected_aperture.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
            

