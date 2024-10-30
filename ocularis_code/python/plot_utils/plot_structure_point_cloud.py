# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:39:00 2022

@author: bjoerk_c
"""



from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path
from configuration import config, data
from matplotlib.pyplot import cm
import numpy as np
from scipy.ndimage import zoom

import matplotlib.pyplot as plt
import reference_frames
import geometry
from plot_utils.Plotter import Plotter




class PlotStructurePointCloud(Plotter):
    
    def __init__(self, config, patient):
        self.config = config
        self.patient = patient

    def plot(self):

        
        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        
        fig = plt.figure()
        
        
        
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(-70,90) 
        
        ref.add_frame_to_plot(ax, ref, 10)
        
        
        for name in ['GTV', 'Lens', 'Sclera']:
            struct = self.patient.patient_model.structure_set[name]
            plt.plot([],[],[], color = struct.color , label = name)
            for i,j,k in struct.contour_coordinates[::2]:
                
                alpha = 1
                
                if struct.organ.name == 'Sclera':
                    alpha = 0.08
                elif struct.organ.name == 'Vitreous':
                    alpha = 0.01
                
                plt.plot(i,j,k, marker = "*",color = self.patient.patient_model.structure_set[name].color , alpha = alpha)
        
        
        gaze_vector = self.patient.patient_model.structure_set.gaze_vector
        
        
        x, y, z = self.patient.patient_model.structure_set["Sclera"].com
        
        plt.quiver(x,y,z, gaze_vector[0], gaze_vector[1],gaze_vector[2], color = "k", length=15, label = "Gaze vector")
        
        
        #plt.title("com " + str(self.patient.patient_model.structure_set["GTV"].com))
        
        plt.title(f"Gaze vector = {self.patient.patient_model.structure_set.gaze_vector}", fontsize = 14)
        
        ax.set_xlabel('X', fontsize = 12)
        # ax.set_xlim(0, 120)
        ax.set_ylabel('Y', fontsize = 12)
        # ax.set_ylim(150, 200)
        ax.set_zlabel('Z', fontsize = 12)
        
        plt.legend()
        

        xlength = 16
        fig.set_size_inches(xlength, xlength/1.61803398875)
        plt.show()
        

        
        

        try:
            plt.savefig(self.algo.plot_path / "Plot_structure_point_cloud.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
            
