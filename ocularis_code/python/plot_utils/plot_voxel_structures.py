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
from scipy import ndimage

import matplotlib.pyplot as plt

import reference_frames

def format_func(value, tick_number):
  
    if value == 0:
        return 0
    else:
        return -value




class PlotVoxelStructures(Plotter):
    
    def __init__(self, config, patient):
        self.config = config
        self.patient = patient
        

    def plot(self):
        
        OG = self.patient.patient_model.structure_set.grid


        x0, y0, z0 = np.indices((OG.size[0]+1, OG.size[1]+1, OG.size[2]+1))
        xOG = x0*OG.spacing[0] + OG.origin[0]
        yOG = y0*OG.spacing[1] + OG.origin[1]
        zOG = z0*OG.spacing[2] + OG.origin[2]

        
        
        
        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        
        
        fig = plt.figure()
        
        
        
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(-70,90) 
        
        ref.add_frame_to_plot(ax, ref, 10)
        
        for i, name in enumerate(['GTV', 'Lens', 'Sclera']): #set1.name_list
            struct = self.patient.patient_model.structure_set[name]
            print(struct)
        
            alpha = 1
        
            if struct.organ.name == 'Sclera':
                alpha = 0.03
            elif struct.organ.name == 'Vitreous':
                alpha = 0.01
            
            # volume_sub = zoom(structure1.binary_mask, (downsampling[0], downsampling[1], downsampling[2]))
            # volume_sub = volume_sub/np.max(volume_sub)
        
            # print(volume_sub.shape)
        
            # if np.nan in volume_sub:
            #     print("Warning")
        
            # struct.set_grid(dcm_volume_patient1.grid)
            ax.voxels(xOG, yOG, zOG ,  struct.binary_mask, facecolors=struct._color, alpha = alpha, label = struct.organ.name)
        
        
        
        plt.title("DICOM", fontsize = 12, y=1.15)
        
        # plt.plot([],[],[], color = "C1", label = "GTV")
        # for p in structure1.contour_coordinates:
        #     plt.plot(p[0], p[1], p[2], linestyle = "", marker = "*", color = "C0")
        
        
        
        ax.set_xlabel('X', fontsize = 12)
        # ax.set_xlim(0, 120)
        ax.set_ylabel('Y', fontsize = 12)
        # ax.set_ylim(150, 200)
        ax.set_zlabel('Z', fontsize = 12)
        
        
        xlength = 16
        fig.set_size_inches(xlength, xlength/1.61803398875)
        plt.show()
        
        

        xlength = 16
        fig.set_size_inches(xlength, xlength/3)
        plt.show()
        

        
        
   
