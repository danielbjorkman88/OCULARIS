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

def format_func(value, tick_number):
  
    if value == 0:
        return 0
    else:
        return -value




class PlotCrossPoints(Plotter):
    
    def __init__(self, patient):
        self.patient = patient
        

    def plot(self, x_bin, y_bin, z_bin, names):
        
        
        struct = self.patient.patient_model.structure_set[names[0]]
        
        x_min =self.patient.patient_model.structure_set.grid.origin[0] + self.patient.patient_model.structure_set.grid.spacing[0]*x_bin
        x_max = self.patient.patient_model.structure_set.grid.origin[0] + self.patient.patient_model.structure_set.grid.spacing[0]*(x_bin + 1)
        
        y_min = self.patient.patient_model.structure_set.grid.origin[1] + self.patient.patient_model.structure_set.grid.spacing[1]*y_bin
        y_max = self.patient.patient_model.structure_set.grid.origin[1] + self.patient.patient_model.structure_set.grid.spacing[1]*(y_bin + 1)
        
        z_min = self.patient.patient_model.structure_set.grid.origin[2] + self.patient.patient_model.structure_set.grid.spacing[2]*z_bin
        z_max = self.patient.patient_model.structure_set.grid.origin[2] + self.patient.patient_model.structure_set.grid.spacing[2]*(z_bin + 1)
        



        mask = struct.binary_mask
        # com_bins = np.asarray(ndimage.measurements.center_of_mass(mask))
        # target_com = com_bins*struct.grid.spacing + struct.grid.origin    

        
        fig = plt.figure()
        
        
        ax = fig.add_subplot(131)
        
        
        
        image = mask[x_bin, 0:, 0:]
        plt.pcolor(struct.grid.meshgrids.Z, struct.grid.meshgrids.Y, image, cmap=plt.cm.bone)
        
        
        
        for row in struct.contour_coordinates: # [::3]
            if x_min < row[0] < x_max:
                plt.scatter(row[2], row[1], color = "C1")
        
        
        ax.set_xlabel('Z', fontsize = 12)
        ax.set_ylabel('Y', fontsize = 12)
        
        plt.title(f"x slice = {x_bin}") #. Between {y_min} and {y_max}
        
        
        
        ax = fig.add_subplot(132)
        
        
        
        image = mask[0:, y_bin, 0:]
        plt.pcolor(struct.grid.meshgrids.Z, struct.grid.meshgrids.X, image, cmap=plt.cm.bone)
        
        
        
        for row in struct.contour_coordinates: # [::3]
            if y_min < row[1] < y_max:
                plt.scatter(row[2], row[0], color = "C1")
        
        
        ax.set_xlabel('Z', fontsize = 12)
        ax.set_ylabel('X', fontsize = 12)
        
        plt.title(f"y slice = {y_bin}") #. Between {y_min} and {y_max}
        
        
        
        
        ax = fig.add_subplot(133)
        
        
        
        
        image = mask[0:, 0:, z_bin]
        image = np.swapaxes(image, 0, 1)
        plt.pcolor(struct.grid.meshgrids.trans_X , struct.grid.meshgrids.trans_Y, image, cmap=plt.cm.bone)
        
        
        
        for row in struct.contour_coordinates: # [::3]
            if z_min < row[2] < z_max:
                plt.scatter(row[0], row[1], color = "C1")
        
        
        ax.set_xlabel('X', fontsize = 12)
        ax.set_ylabel('Y', fontsize = 12)
        
        plt.title(f"z slice = {z_bin}") #. Between {y_min} and {y_max}
        
        
        xlength = 16
        fig.set_size_inches(xlength, xlength/3)
        plt.show()
        

        
        
   
