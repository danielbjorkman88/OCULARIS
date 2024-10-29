# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:39:00 2022

@author: bjoerk_c
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

class PlotTarget(Plotter):
    


    def plot(self):

        ref = reference_frames.Reference_Frame(self.algo.config, 'TreatmentRoom')

        fig = plt.figure()
                    
        ax = plt.subplot(211)
        
        
        # ax.add_patch(self.algo.collimator_top)
        # ax.add_patch(self.algo.collimator_bot)
        # ax.add_patch(self.algo.phantom)
        # ax.add_patch(self.algo.wedge_patch)
        
        
        image = self.algo.target_matrix[0:, int(self.algo.dose.shape[1]/2), 0:] * 0.6
        
        
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image, cmap = "jet", vmin=0, vmax=1)

        
        
        
        plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes , color = "k", label = "Central axis")
        
        

        
        plt.plot([],[], color = "C1", label = "Target Rays")
        for ray in self.algo.target_rays:
            plt.plot(ray.zes, ray.xes , color = "C1")
            
        
        
        plt.scatter([], [] , color = "k", label = "Ray intersects")    
        plt.scatter(self.algo.collimator.target_rays_intersects_zes, self.algo.collimator.target_rays_intersects_xes , color = "k")
        # plt.scatter(algo.wedge.intersects_zes, algo.wedge.intersects_xes , color = "k")
        
        
        plt.scatter(0, 0, marker = ".", label = "Isocentre", color = "r", s = 200)
        plt.scatter(self.algo.central_axis.S[2], self.algo.central_axis.S[0], marker = ".", label = "VPS", color = "b", s = 200)    
        
        
        
        

        plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2a[2] ], [self.algo.aperture_point1[0], self.algo.aperture_point2a[0] ], linewidth = 3, color = "C0", label = "Projected aperture from VPS at {} mm".format(config["VPS"]))
        plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2b[2] ], [self.algo.aperture_point1[0], self.algo.aperture_point2b[0] ], color = "C0", linewidth = 3)
        
        
        
        plt.plot([self.algo.target_plane_p1x[2],self.algo.target_plane_p2x[2]], [self.algo.target_plane_p1x[0],self.algo.target_plane_p2x[0]], color = "r", label = "Target plane")
        plt.plot([self.algo.aperture_plane_p1x[2],self.algo.aperture_plane_p2x[2]], [self.algo.aperture_plane_p1x[0],self.algo.aperture_plane_p2x[0]], color = "y", label = "Aperture plane")
        

        
        
        plt.grid(linewidth = 0.3)
        plt.legend(loc='upper right', prop={'size': 8})
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm from isocentre]", fontsize = 14)
        
        
        
        # plt.title(f'Theta = {self.algo.config["Gantry rot theta"]} degrees, Phi = {self.algo.config["Gantry rot phi"]} ', fontsize = 14)
        
        

        
        plt.xlim(-150,110)
        plt.ylim(-100,100)
        
        


                    
        ax = plt.subplot(223)
        
        
 
        plt.scatter(self.algo.collimator.target_rays_intersects_in_plane_xes, self.algo.collimator.target_rays_intersects_in_plane_yes , color = "b", label = "Target rays intersects with aperture plane")
        
        points = np.array(self.algo.collimator.target_rays_aperture_intersect_points_in_plane)
        plt.plot([], [], color = "r", label = "Concave hull")
        for i, j in self.algo.collimator.edges_in_plane:
            plt.plot(points[[i, j], 0], points[[i, j], 1], color = "r")
        # for xes, yes, zes in self.algo.collimator.edges_global:
        #     ax.scatter(xes, yes, zes, c="C0",s=40)
        #     ax.plot(xes, yes, zes, color="C0")
        
        #plt.plot(self.algo.collimator.hull_in_plane_xes , self.algo.collimator.hull_in_plane_yes , 'r--', lw=4, label = "Convex hull")
        
        # plt.plot(collimator.hull2_xes , collimator.hull2_yes , 'b--', lw=2, label = "Convex hull 2")
        
        
        # plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
        
        # plt.scatter(mid_point[0], mid_point[1], label = "" )
        
        plt.grid(linewidth = 0.3)
        plt.legend(loc='upper right', prop={'size': 10})
        
        plt.ylabel("y [mm]", fontsize = 14)
        plt.xlabel("x [mm]", fontsize = 14)
        

        
        plt.title("Aperture plane", fontsize = 20)

        

        
        
        ax1 = fig.add_subplot(224, projection='3d')
        ax1.view_init(-70,90) 
        # plt.title("Treatment Room", fontsize = 18, y=1.15)
        
        # ref.add_frame_to_plot(ax1, ref)
        
        
        arrow_prop_dict = dict(mutation_scale=20, arrowstyle='->', shrinkA=0, shrinkB=0)
        
        step = 10
        
        for ray in self.algo.rays[::step]:
            a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color="C0")
            ax1.add_artist(a)   
        
        
        for ray in self.algo.target_rays[::280]:
            a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color="C1")
            ax1.add_artist(a)
            
        
        ax1.text(ref.origin_ref[0], ref.origin_ref[1], ref.origin_ref[2] -0.1, r'$0_{iso}$', fontsize = 12)
        
        
        ax1.set_xlabel('X')
        ax1.set_xlim(-30, 30)
        ax1.set_ylabel('Y')
        ax1.set_ylim(-30, 30)
        ax1.set_zlabel('Z')
        ax1.set_zlim(-30, 30)
        
        
        plt.suptitle(f'XZ plane in Treatment Room, Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees, Wedge cover = {self.algo.config["wedge_cover"]}', fontsize = 16)
        


            
        xlength = 12
        fig.set_size_inches(xlength, xlength/1.61803398875)
        #fig.set_size_inches(xlength, xlength)
        plt.show()
        try:
            plt.savefig(self.plot_path / "Plot_target.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
            

