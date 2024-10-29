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




class PlotCrossCubes(Plotter):
    

    def __init__(self, algo, patient):
        self.algo = algo
        self.patient = patient


    def plot(self,x_bin, y_bin):
        
        
        # z_bin = 40
        
        # z = algo.medium.minZ + z_bin*algo.medium.resolution[2]
        # points = self.patient.patient_model.structure_set["Sclera"].contour_coordinates
        # points = list(filter(lambda p: abs(p[2] - z) < 1, points))
        
        
        #y_bin = 130
        
        # y_bin = int(self.algo.config["Slice"][1])
        
        
        fig = plt.figure()
        
        ax = plt.subplot(221)
        
        image_structures = np.zeros(self.patient.patient_model.structure_set["Sclera"].binary_mask[0:,y_bin, 0:].shape)
        
        image_structures[self.algo.medium.medium[0:,y_bin, 0:] == 1] = 1
        image_structures[self.patient.patient_model.structure_set["Sclera"].binary_mask[0:,y_bin, 0:]] = 2
        image_structures[self.patient.patient_model.structure_set["GTV"].binary_mask[0:,y_bin, 0:]] = 3
        
        # image_structures = self.algo.target_matrix[0:,y_bin, 0:]
        
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image_structures) #,  image_structures self.algo.medium.trans_X , self.algo.medium.trans_Y,
        
        
        plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes , color = "k", label = "Central axis")
        
        # plt.scatter(self.algo.max_target_point[2], self.algo.max_target_point[0], label = "Modulation depth max", color = "C0")
        # plt.scatter(self.algo.min_target_point[2], self.algo.min_target_point[0], label = "Modulation depth min", color = "C1")
        
        
        
        plt.title("Medium", fontsize = 12) # Sclera + GTV
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)

        plt.xlim(self.algo.medium.minZ,self.algo.medium.maxZ)
        #plt.ylim(self.algo.medium.minX,self.algo.medium.maxX)
        plt.ylim(-20,20)
        plt.legend()
        
        ax = plt.subplot(222)
        
        #image = self.patient.patient_model.structure_set["Lens"].binary_mask[0:,y_bin, 0:]
        
        

        
        # y_bin = 141
        # x_bin = 180
        
        # fig = plt.figure()
        
        # ax = plt.subplot(121)
        
        image_structures = np.zeros(self.patient.patient_model.structure_set["Sclera"].binary_mask[0:,y_bin, 0:].shape)
        image_structures[self.patient.patient_model.structure_set["Sclera"].binary_mask[0:,y_bin, 0:]] = 1
        image_structures[self.patient.patient_model.structure_set["GTV"].binary_mask[0:,y_bin, 0:]] = 2
        
        ax.pcolor( self.algo.medium.Z, self.algo.medium.X , image_structures)
        
        image = self.algo.raytracer.depth[0:,y_bin,0:]
        
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        plt.ylim(-20,20)
        
        # ax2 = ax.twinx()
        
        # ax.pcolor(image)
        ax.pcolor( self.algo.medium.Z, self.algo.medium.X , image, alpha = 0.4, cmap = "jet")
        # cbar = plt.colorbar()
        # cbar.ax.set_ylabel('Radiological depth', rotation=270, size='xx-large')
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        plt.ylim(-20,20)
        
        
        plt.title("WEQ Depth", fontsize = 12) # Sclera + GTV
        # plt.xlim(-10, 80)
        
        ax3 = ax.twinx()
        
        
        
        # plt.plot(self.algo.DD.xes, self.algo.DD.yes)
        
        
        plt.plot(self.algo.medium.z_voxel_pos  ,self.algo.raytracer.depth[x_bin,y_bin,0:], color = "C1")
        
        
        
        # plt.xlim(-10, 80)
        # plt.ylabel("Radiological Depth [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        

        
        # ax = plt.subplot(121)
        
        # image_structures = np.zeros(image.shape)
        
        
        # ax.pcolor( algo.medium.Z, algo.medium.X , image_structures)
        
        
        # xlength = 16
        # fig.set_size_inches(xlength, xlength/1.6)
        # plt.show()
        

        
        # plt.xlim(self.algo.medium.minZ,self.algo.medium.maxZ)
        # plt.ylim(self.algo.medium.minX,self.algo.medium.maxX)
        
        ax = plt.subplot(224)
        
        image_dose = self.algo.dose[0:,y_bin, 0:]
        
        #image = self.patient.patient_model.structure_set["Lens"].binary_mask[0:,y_bin, 0:]
        
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X,  image_structures) #self.algo.medium.trans_X , self.algo.medium.trans_Y,

        masked_data = np.ma.masked_where(image_dose < 0.001, image_dose)
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X,  masked_data, alpha = 0.6, cmap = "jet", vmin = 0, vmax = 1) #self.algo.medium.trans_X , self.algo.medium.trans_Y,
        # plt.pcolor(self.algo.medium.Z, self.algo.medium.X,  [], alpha = 1, cmap = "jet", vmin = 0, vmax = 1)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')

        plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes , color = "k", label = "Central axis")
        
        

        # plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image_dose, cmap = "jet") 
        
        # cbar = plt.colorbar()
        
        # cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')
        
        
        plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes , color = "k", label = "Central axis")
        
        plt.title("Dose", fontsize = 12) 
        # for p in points:
        #     plt.scatter(p[0], p[1], color = "C1")
        
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        
        plt.xlim(self.algo.medium.minZ,self.algo.medium.maxZ)
        plt.ylim(-20,20)
                
        ax2 = ax.twinx()
        
        plt.plot(self.algo.medium.z_voxel_pos  ,self.algo.dose[x_bin,y_bin,0:], color = "w")
        
        
        
        # plt.xlim(-10, 80)
        # plt.ylabel("Depth [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        
        
        plt.xlim(self.algo.medium.minZ,self.algo.medium.maxZ)

        # plt.ylim(self.algo.medium.minX,self.algo.medium.maxX)
        
        
        
        ax = plt.subplot(223)
        
        
        my_cmap = "seismic" # vmax = 1/vmin
        
        #image = self.patient.patient_model.structure_set["Lens"].binary_mask[0:,y_bin, 0:]
        
        
        image_dose = self.algo.dist_aperture[0:,y_bin, 0:]
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image_dose, cmap = my_cmap, vmin = -11, vmax = 11) 
        
        cbar = plt.colorbar()
        
        cbar.ax.set_ylabel('Distance to projected aperture', rotation=270, size='xx-large')
        
        
        plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes , color = "k", label = "Central axis")
        

        plt.title("Distance to projected aperture", fontsize = 12) 
        # for p in points:
        #     plt.scatter(p[0], p[1], color = "C1")
        
        
        plt.ylabel("x [mm]", fontsize = 14)
        plt.xlabel("z [mm]", fontsize = 14)
        
        plt.xlim(self.algo.medium.minZ,self.algo.medium.maxZ)
        plt.ylim(-20,20)
        
        
                
        
        
        
        plt.suptitle( f"{self.algo.nozzle.title},  x, y slice {x_bin}, {y_bin}" )
        
        
        xlength = 20
        fig.set_size_inches(xlength, xlength)
        plt.show()


        # if self.algo.wedge.traced:
        #     plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees, Wedge cover = {self.algo.config["wedge_cover"]}', fontsize = 16)
        # else:
        #     plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees', fontsize = 16)
        
        
        # xlength = 18
        # fig.set_size_inches(xlength, xlength/3)
        # #fig.set_size_inches(xlength, xlength)
        # plt.show()
        try:
            plt.savefig(self.plot_path / "Plot_CrossCubes.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
        
        
        
   
