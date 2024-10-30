# -*- coding: utf-8 -*-




from plot_utils.Plotter import Plotter

from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path
from configuration import config, data
from matplotlib.pyplot import cm
import numpy as np
import reference_frames

import matplotlib.pyplot as plt


def format_func(value, tick_number):

    if value == 0:
        return 0
    else:
        return -value


class PlotEyeplanModelTreatmentRoom(Plotter):

   def __init__(self, config, patient):
       self.config = config
       self.patient = patient

   def plot(self):
       
        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
           
        fig = plt.figure()
        
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(-100,90) 
        
        ref.add_frame_to_plot(ax, ref,20)
        
        for name in ['eyeglobe', 'target']: #patient.patient_model.eyeplan_model.structure_set.name_list:
            if name == "Clips":
                continue
            struct = self.patient.patient_model.eyeplan_model.structure_set[name]
            
            points = struct.contour_coordinates
            # points = self.patient.patient_model.eyeplan_model.rotation_matrix_z_to_gaze_light.dot(points.T).T
            
            # save.append(points)
            
            alpha = 0.4
            color = "k"
            if name == 'target':
                alpha = 1
                color = "C2"
            
        
            ax.scatter(points[0:,0], points[0:,1], points[0:,2],  label = name, alpha = alpha, color = color)
        
        
        # plt.plot([],[],[], color = "g" , alpha = 0.1, label = "Skin plane")
        # ax.plot_surface(xx, yy, z, alpha = 0.1, color = "g")
        
        
        ax.plot(self.patient.patient_model.eyeplan_model.light_point[0], self.patient.patient_model.eyeplan_model.light_point[1], self.patient.patient_model.eyeplan_model.light_point[2], markerfacecolor='C1', markeredgecolor='C1', marker='o', markersize=10, alpha=1, label = f"Fixation point = {self.patient.patient_model.eyeplan_model.light_point}")
        
        gaze_vector_light = self.patient.patient_model.eyeplan_model.gaze_vector_light
        plt.quiver(self.patient.patient_model.eyeplan_model.centre_of_model[0],self.patient.patient_model.eyeplan_model.centre_of_model[1],self.patient.patient_model.eyeplan_model.centre_of_model[2],gaze_vector_light[0] , gaze_vector_light[1], gaze_vector_light[2], color = "r", label = f"Fixation light gaze {gaze_vector_light}", length=25)
        
        # plt.quiver([0],[0],[0],gaze_vector_planes[0] , gaze_vector_planes[1], gaze_vector_planes[2], color = "grey", label = "Plane intersect gaze", length=25)
        
        
        
        # xes_light = [0, self.patient.patient_model.eyeplan_model.light_point[0]]
        # yes_light = [0, self.patient.patient_model.eyeplan_model.light_point[1]]
        # zes_light = [0, self.patient.patient_model.eyeplan_model.light_point[2]]
        
        # p = np.asarray([0,0,0]) + gaze_vector_planes*150
        # xes_planes = [0, p[0]]
        # yes_planes = [0, p[1]]
        # zes_planes = [0, p[2]]
        
        p1 = self.patient.patient_model.eyeplan_model.centre_of_model + gaze_vector_light*150
        p2 = self.patient.patient_model.eyeplan_model.centre_of_model - gaze_vector_light*150
        xes_light = [p1[0], p2[0]]
        yes_light = [p1[1], p2[1]]
        zes_light = [p1[2], p2[2]]
        
        
        plt.plot(xes_light , yes_light, zes_light, color = "k", linestyle= "--", label = "Eyeplan model axis")
        
        # plt.plot(xes_planes , yes_planes, zes_planes, color = "grey", linestyle= "--")
        
        # ax.scatter([],[],[], s= 20, c = "C0", label = "Feature")
        # for p in self.patient.patient_model.eyeplan_model.feature_points:
        #     ax.scatter(p[0], p[1], p[2], color = "C0", s= 20)
        
        
        plt.legend()
        
        zoom= 1.3
        
        ax.set_xlabel('X', fontsize = 12)
        ax.set_xlim(-20*zoom, 20*zoom)
        ax.set_ylabel('Y', fontsize = 12)
        ax.set_ylim(-20*zoom, 20*zoom)
        ax.set_zlabel('Z', fontsize = 12)
        ax.set_zlim(-20*zoom, 20*zoom)


        xlength = 12
        fig.set_size_inches(xlength, xlength/1.61803398875)
        plt.show()
        # try:
        #      plt.savefig(self.plot_path / "plot_dose_on_dicom.pdf",
        #                  bbox_inches='tight', pad_inches=0.1)
        # except:
        #      pass


