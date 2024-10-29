# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:27:51 2022

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

def format_func(value, tick_number):
  
    if value == 0:
        return 0
    else:
        return -value




class PlotCross(Plotter):
    


    def plot(self):
        
        fig = plt.figure()
            
        ax = plt.subplot(131)
        
        
        # ax.add_patch(self.algo.collimator_top)
        # ax.add_patch(self.algo.collimator_bot)
        # ax.add_patch(self.algo.phantom)
        # ax.add_patch(self.algo.wedge_patch)
        
        
        image = self.algo.dose[0:, int(self.algo.dose.shape[1]/2), 0:]
    
    
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image, cmap = "jet", vmin=0, vmax=1)
        cbar = plt.colorbar()
        
        cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')
        
        # print(self.algo.config["Gantry rot phi"])
        
        if self.algo.config["Gantry rot phi"] != 0:
            pass
        else:
            plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes , color = "k", label = "Central axis")
            plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2a[2] ], [self.algo.aperture_point1[0], self.algo.aperture_point2a[0] ], linewidth = 3, color = "C0", label = "Aperture limit")
            plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2b[2] ], [self.algo.aperture_point1[0], self.algo.aperture_point2b[0] ], color = "C0", linewidth = 3)
            
        
        
        # colors = cm.rainbow(np.linspace(0,1,len(selfs)))
        
        # for i,alg in enumerate(selfs):
        #     plt.plot(alg.inflection_xes, alg.inflection_yes, color = colors[i]) # label = alg.nozzle.modulator_wheel
        
        
        # plt.plot([],[], color = "C1", label = "Rays")
        # for ray in self.algo.raytracer.rays:
        #     plt.plot(ray.zes, ray.xes , color = "C1")
            
        
        
        # plt.scatter([], [] , color = "k", label = "Ray intersects")    
        # plt.scatter(self.algo.raytracer.intersects_zes, self.algo.raytracer.intersects_xes , color = "k")
        # plt.scatter(self.algo.wedge.intersects_zes, self.algo.wedge.intersects_xes , color = "k")
        
        
        plt.scatter(0, 0, marker = ".", label = "Isocentre", color = "r", s = 400)
        #plt.scatter(self.algo.central_axis.S[2], self.algo.central_axis.S[0], marker = ".", label = "VPS", color = "b", s = 200)    
        
        
        
        
        #plt.plot(self.algo.beta_xes, self.algo.beta_yes, label = "beta")
        
        # plt.plot(self.algo.inflection_xes2, self.algo.inflection_yes2, label = "Half point for MW901", color = "C1")    
        
        
        # #Projected Aperture
        # plt.plot(self.algo.aperture_line_xes, self.algo.aperture_line_yes, linewidth = 3, color = "C0", label = "Aperture limit")
        # plt.plot(self.algo.aperture_line_xes, [ -x for x in self.algo.aperture_line_yes], color = "C0", linewidth = 3)
        

        
        # plt.plot([self.algo.target_plane_p1x[2],self.algo.target_plane_p2x[2]], [self.algo.target_plane_p1x[0],self.algo.target_plane_p2x[0]], color = "r", label = "Target plane")
        # plt.plot([self.algo.aperture_plane_p1x[2],self.algo.aperture_plane_p2x[2]], [self.algo.aperture_plane_p1x[0],self.algo.aperture_plane_p2x[0]], color = "y", label = "Aperture plane")
        #plt.plot([self.algo.wedge_plane_p1x[2],self.algo.wedge_plane_p2x[2],self.algo.wedge_plane_p3x[2]], [self.algo.wedge_plane_p1x[0],self.algo.wedge_plane_p2x[0],self.algo.wedge_plane_p3x[0]], color = "g", label = f"Wedge plane, {self.algo.wedge_angle} degrees")
        
        
        # plt.scatter(self.algo.aperture_point[2], self.algo.aperture_point[0], label = "Aperture point", color = "y", s = 200, marker = ".")
        # plt.scatter(self.algo.wedge_apex_point[2], self.algo.wedge_apex_point[0], label = "Wedge apex point", color = "k", s = 200, marker = ".")
        # plt.scatter(self.algo.wedge_point[2], self.algo.wedge_point[0], label = "Wedge point", color = "g", marker = ".", s = 200)
        # plt.scatter(self.algo.wedge_top_point[2], self.algo.wedge_top_point[0], label = "Wedge top point", color = "c", marker = ".", s = 200)
        
        
        
        plt.grid(linewidth = 0.3)
        plt.legend(loc='upper right', prop={'size': 8})
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        
        plt.title('XZ slice', fontsize = 16)
        
        #plt.title(f'Optis Geometry, azimuth = {self.algo.config["Gantry rot theta"] } ', fontsize = 16)
        # plt.title(f'XZ plane at Gantry rot, theta = {self.algo.config["Gantry rot theta"]} degrees ', fontsize = 16)
        # plt.title(f'XZ plane with Wedge cover = {self.algo.config["wedge_cover"]} mm, Wedge insert angle = {self.algo.config["wedge_insert_angle"]} ', fontsize = 16)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        # ["wedge_cover"]
        
        plt.xlim(self.algo.medium.mesh_origin[2],self.algo.medium.mesh_apex[2])
        plt.ylim(self.algo.medium.mesh_origin[0],self.algo.medium.mesh_apex[0])
        #plt.ylim(-100,100)
        
        
        
        
        #--------------------------------------------------------
        
        
        
        ax = plt.subplot(132)
        
        
        # ax.add_patch(self.algo.collimator_top)
        # ax.add_patch(self.algo.collimator_bot)
        # ax.add_patch(self.algo.phantom)
        # ax.add_patch(self.algo.wedge_patch)
        
        
        image = self.algo.dose[int(self.algo.dose.shape[1]/2), 0:, 0:]
    
    
        plt.pcolor(self.algo.medium.Z, self.algo.medium.Y, image, cmap = "jet", vmin=0, vmax=1)
        cbar = plt.colorbar()
        
        cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')

        # print(self.algo.config["Gantry rot theta"] )

        if self.algo.config["Gantry rot theta"] != 0:
            pass
        else:
            plt.plot(self.algo.central_axis.zes, self.algo.central_axis.yes , color = "k", label = "Central axis")
            plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2a[2] ], [self.algo.aperture_point1[1], self.algo.aperture_point2a[1] ], linewidth = 3, color = "C0", label = "Aperture limit")
            plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2b[2] ], [self.algo.aperture_point1[1], self.algo.aperture_point2b[1] ], color = "C0", linewidth = 3)
            
            
            

        
        
        
        # colors = cm.rainbow(np.linspace(0,1,len(selfs)))
        
        # for i,alg in enumerate(selfs):
        #     plt.plot(alg.inflection_xes, alg.inflection_yes, color = colors[i]) # label = alg.nozzle.modulator_wheel
        
        
        # plt.plot([],[], color = "C1", label = "Rays")
        # for ray in self.algo.raytracer.rays:
        #     plt.plot(ray.zes, ray.yes , color = "C1")
            
        
        
        # plt.scatter([], [] , color = "k", label = "Ray intersects")    
        # plt.scatter(self.algo.raytracer.intersects_zes, self.algo.raytracer.intersects_xes , color = "k")
        # plt.scatter(self.algo.wedge.intersects_zes, self.algo.wedge.intersects_xes , color = "k")
        
        
        plt.scatter(0, 0, marker = ".", label = "Isocentre", color = "r", s = 400)
        #plt.scatter(self.algo.central_axis.S[2], self.algo.central_axis.S[0], marker = ".", label = "VPS", color = "b", s = 200)    
        
        
        
        
        #plt.plot(self.algo.beta_xes, self.algo.beta_yes, label = "beta")
        
        # plt.plot(self.algo.inflection_xes2, self.algo.inflection_yes2, label = "Half point for MW901", color = "C1")    
        
        
        # #Projected Aperture
        # plt.plot(self.algo.aperture_line_xes, self.algo.aperture_line_yes, linewidth = 3, color = "C0", label = "Aperture limit")
        # plt.plot(self.algo.aperture_line_xes, [ -x for x in self.algo.aperture_line_yes], color = "C0", linewidth = 3)
        

        
        # plt.plot([self.algo.target_plane_p1x[2],self.algo.target_plane_p2x[2]], [self.algo.target_plane_p1x[0],self.algo.target_plane_p2x[0]], color = "r", label = "Target plane")
        # plt.plot([self.algo.aperture_plane_p1x[2],self.algo.aperture_plane_p2x[2]], [self.algo.aperture_plane_p1x[0],self.algo.aperture_plane_p2x[0]], color = "y", label = "Aperture plane")
        #plt.plot([self.algo.wedge_plane_p1x[2],self.algo.wedge_plane_p2x[2],self.algo.wedge_plane_p3x[2]], [self.algo.wedge_plane_p1x[0],self.algo.wedge_plane_p2x[0],self.algo.wedge_plane_p3x[0]], color = "g", label = f"Wedge plane, {self.algo.wedge_angle} degrees")
        
        
        # plt.scatter(self.algo.aperture_point[2], self.algo.aperture_point[0], label = "Aperture point", color = "y", s = 200, marker = ".")
        # plt.scatter(self.algo.wedge_apex_point[2], self.algo.wedge_apex_point[0], label = "Wedge apex point", color = "k", s = 200, marker = ".")
        # plt.scatter(self.algo.wedge_point[2], self.algo.wedge_point[0], label = "Wedge point", color = "g", marker = ".", s = 200)
        # plt.scatter(self.algo.wedge_top_point[2], self.algo.wedge_top_point[0], label = "Wedge top point", color = "c", marker = ".", s = 200)
        
        
        
        plt.grid(linewidth = 0.3)
        plt.legend(loc='upper right', prop={'size': 8})
        
        plt.ylabel("y [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        
        
        plt.title('YZ slice', fontsize = 16)
        
        #plt.title(f'Optis Geometry, azimuth = {self.algo.config["Gantry rot theta"] } ', fontsize = 16)
        # plt.title(f'XZ plane at Gantry rot, theta = {self.algo.config["Gantry rot theta"]} degrees ', fontsize = 16)
        # plt.title(f'XZ plane with Wedge cover = {self.algo.config["wedge_cover"]} mm, Wedge insert angle = {self.algo.config["wedge_insert_angle"]} ', fontsize = 16)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        # ["wedge_cover"]
        
        # plt.xlim(-150,110)
        # plt.ylim(-40,40)
        # plt.ylim(-100,100)
        
        plt.xlim(self.algo.medium.mesh_origin[2],self.algo.medium.mesh_apex[2])
        plt.ylim(self.algo.medium.mesh_origin[1],self.algo.medium.mesh_apex[1])
        
        
        
        
        #--------------------------------------------------------
        
        
        
        ax = plt.subplot(133)
        
        
        # ax.add_patch(self.algo.collimator_top)
        # ax.add_patch(self.algo.collimator_bot)
        # ax.add_patch(self.algo.phantom)
        # ax.add_patch(self.algo.wedge_patch)
        
        
        image = self.algo.dose[0:, 0:, int(self.algo.dose.shape[2]/2)]
        image = np.swapaxes(image, 0, 1)
        
        #print(image.shape)
        #print(len(self.algo.medium.trans_X) , len(self.algo.medium.trans_Y))
    
        plt.pcolor(self.algo.medium.trans_X , self.algo.medium.trans_Y, image, cmap = "jet", vmin=0, vmax=1)
        cbar = plt.colorbar()
        
        cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')


        
        #plt.plot(self.algo.central_axis.zes, self.algo.central_axis.yes , color = "k", label = "Central axis")
        
        
        # colors = cm.rainbow(np.linspace(0,1,len(selfs)))
        
        # for i,alg in enumerate(selfs):
        #     plt.plot(alg.inflection_xes, alg.inflection_yes, color = colors[i]) # label = alg.nozzle.modulator_wheel
        
        
        # plt.plot([],[], color = "C1", label = "Rays")
        # for ray in self.algo.raytracer.rays:
        #     plt.plot(ray.zes, ray.xes , color = "C1")
            
        
        
        # plt.scatter([], [] , color = "k", label = "Ray intersects")    
        # plt.scatter(self.algo.raytracer.intersects_zes, self.algo.raytracer.intersects_xes , color = "k")
        # plt.scatter(self.algo.wedge.intersects_zes, self.algo.wedge.intersects_xes , color = "k")
        
        
        #plt.scatter(0, 0, marker = ".", label = "Isocentre", color = "r", s = 400)
        #plt.scatter(self.algo.central_axis.S[2], self.algo.central_axis.S[0], marker = ".", label = "VPS", color = "b", s = 200)    
        
        
        
        
        #plt.plot(self.algo.beta_xes, self.algo.beta_yes, label = "beta")
        
        # plt.plot(self.algo.inflection_xes2, self.algo.inflection_yes2, label = "Half point for MW901", color = "C1")    
        
        
        # #Projected Aperture
        # plt.plot(self.algo.aperture_line_xes, self.algo.aperture_line_yes, linewidth = 3, color = "C0", label = "Aperture limit")
        # plt.plot(self.algo.aperture_line_xes, [ -x for x in self.algo.aperture_line_yes], color = "C0", linewidth = 3)
        
        # plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2a[2] ], [self.algo.aperture_point1[1], self.algo.aperture_point2a[1] ], linewidth = 3, color = "C0", label = "Aperture limit")
        # plt.plot([self.algo.aperture_point1[2], self.algo.aperture_point2b[2] ], [self.algo.aperture_point1[1], self.algo.aperture_point2b[1] ], color = "C0", linewidth = 3)
        
        
        
        # plt.plot([self.algo.target_plane_p1x[2],self.algo.target_plane_p2x[2]], [self.algo.target_plane_p1x[0],self.algo.target_plane_p2x[0]], color = "r", label = "Target plane")
        # plt.plot([self.algo.aperture_plane_p1x[2],self.algo.aperture_plane_p2x[2]], [self.algo.aperture_plane_p1x[0],self.algo.aperture_plane_p2x[0]], color = "y", label = "Aperture plane")
        #plt.plot([self.algo.wedge_plane_p1x[2],self.algo.wedge_plane_p2x[2],self.algo.wedge_plane_p3x[2]], [self.algo.wedge_plane_p1x[0],self.algo.wedge_plane_p2x[0],self.algo.wedge_plane_p3x[0]], color = "g", label = f"Wedge plane, {self.algo.wedge_angle} degrees")
        
        
        # plt.scatter(self.algo.aperture_point[2], self.algo.aperture_point[0], label = "Aperture point", color = "y", s = 200, marker = ".")
        # plt.scatter(self.algo.wedge_apex_point[2], self.algo.wedge_apex_point[0], label = "Wedge apex point", color = "k", s = 200, marker = ".")
        # plt.scatter(self.algo.wedge_point[2], self.algo.wedge_point[0], label = "Wedge point", color = "g", marker = ".", s = 200)
        # plt.scatter(self.algo.wedge_top_point[2], self.algo.wedge_top_point[0], label = "Wedge top point", color = "c", marker = ".", s = 200)
        
        
        
        plt.grid(linewidth = 0.3)
        # plt.legend(loc='upper right', prop={'size': 8})
        
        plt.ylabel("y [mm]", fontsize = 14)
        plt.xlabel("x [mm]", fontsize = 14)
        
        
        plt.title('YX slice', fontsize = 16)
        
        #plt.title(f'Optis Geometry, azimuth = {self.algo.config["Gantry rot theta"] } ', fontsize = 16)
        # plt.title(f'XZ plane at Gantry rot, theta = {self.algo.config["Gantry rot theta"]} degrees ', fontsize = 16)
        # plt.title(f'XZ plane with Wedge cover = {self.algo.config["wedge_cover"]} mm, Wedge insert angle = {self.algo.config["wedge_insert_angle"]} ', fontsize = 16)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        # ["wedge_cover"]
        
        # plt.xlim(-150,110)
        # plt.ylim(-40,40)
        # plt.ylim(-100,100)
        
        plt.xlim(self.algo.medium.mesh_origin[0],self.algo.medium.mesh_apex[0])
        plt.ylim(self.algo.medium.mesh_origin[1],self.algo.medium.mesh_apex[1])
        
        if self.algo.wedge.traced:
            plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees, Wedge cover = {self.algo.config["wedge_cover"]}', fontsize = 16)
        else:
            plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees', fontsize = 16)
        
        
        xlength = 18
        fig.set_size_inches(xlength, xlength/3)
        #fig.set_size_inches(xlength, xlength)
        plt.show()
        try:
            plt.savefig(self.plot_path / "Plot_Cross.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
        
        
        
   
